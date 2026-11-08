"""
Step 04 - substrate docking with AutoDock Vina.

Docks every prepared ligand into one or more receptors (wild type plus any
mutants built in step 06) and records the best affinity per pair.

The number this pipeline actually cares about is not the raw affinity of any
one ligand - Vina scores are far too crude for that. It is the *difference*
between L-Tyr and L-Phe in the same pocket, and how that difference moves
when a residue is substituted. Systematic errors in the scoring function
largely cancel in that comparison, which is why selectivity is a more
defensible readout than absolute binding.

Usage:
    python src/s04_dock.py --receptor data/structures/2O6Y_receptor.pdbqt \\
                           --center 12.3 4.5 -7.8
    python src/s04_dock.py --all-mutants --center-from-residues 89
"""

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

from config import (BOX_PADDING, DOCK_DIR, LIG_DIR, STRUCT_DIR, VINA_BIN,
                    VINA_EXHAUSTIVENESS, VINA_NUM_MODES, VINA_SEED)
from pocket import box_from_center, center_from_ligand, center_from_residues

AFFINITY_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")


def run_vina(receptor, ligand, box, out_path, exhaustiveness, seed):
    """Run one Vina docking. Returns list of (mode, affinity, rmsd_lb, rmsd_ub)."""
    cmd = [
        VINA_BIN,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--out", str(out_path),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(VINA_NUM_MODES),
        "--seed", str(seed),
        "--cpu", "2",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if res.returncode != 0:
        raise RuntimeError(f"vina failed ({res.returncode}): {res.stderr[-600:]}")

    modes = []
    for line in res.stdout.splitlines():
        m = AFFINITY_RE.match(line)
        if m:
            modes.append((int(m.group(1)), float(m.group(2)),
                          float(m.group(3)), float(m.group(4))))
    if not modes:
        raise RuntimeError(f"vina produced no poses:\n{res.stdout[-600:]}")
    return modes


def resolve_center(args):
    """Work out the box centre from whichever anchor was supplied."""
    if args.center:
        import numpy as np
        return np.array(args.center, dtype=float), "explicit coordinates"
    ref = Path(args.center_source) if args.center_source else None
    if ref is None:
        raise SystemExit("need --center, or --center-from-ligand/--center-from-residues "
                         "together with --center-source")
    if args.center_from_ligand:
        c, n = center_from_ligand(ref, args.center_from_ligand)
        return c, f"centroid of {args.center_from_ligand} ({n} atoms) in {ref.name}"
    if args.center_from_residues:
        c, n = center_from_residues(ref, args.center_from_residues, args.chain)
        return c, (f"centroid of residues {args.center_from_residues} "
                   f"({n} atoms) in {ref.name}")
    raise SystemExit("no pocket anchor given")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--receptor", action="append", default=None,
                    help="receptor PDBQT (repeatable)")
    ap.add_argument("--all-mutants", action="store_true",
                    help="dock every *_receptor.pdbqt in data/structures")
    ap.add_argument("--ligands", nargs="*", default=None,
                    help="ligand PDBQT stems (default: all in data/ligands)")
    ap.add_argument("--center", nargs=3, type=float, default=None)
    ap.add_argument("--center-source", default=None,
                    help="PDB file the pocket anchor is read from")
    ap.add_argument("--center-from-ligand", default=None, metavar="HET")
    ap.add_argument("--center-from-residues", nargs="*", type=int, default=None)
    ap.add_argument("--chain", default=None)
    ap.add_argument("--padding", type=float, default=BOX_PADDING)
    ap.add_argument("--exhaustiveness", type=int, default=VINA_EXHAUSTIVENESS)
    ap.add_argument("--replicates", type=int, default=3,
                    help="independent seeds per pair (Vina is stochastic)")
    ap.add_argument("--out", default="docking_results.csv")
    args = ap.parse_args()

    # --- receptors --------------------------------------------------------
    if args.all_mutants:
        receptors = sorted(STRUCT_DIR.glob("*_receptor.pdbqt"))
    elif args.receptor:
        receptors = [Path(r) for r in args.receptor]
    else:
        raise SystemExit("give --receptor or --all-mutants")
    receptors = [r for r in receptors if r.exists()]
    if not receptors:
        raise SystemExit("no receptor PDBQT files found")

    # --- ligands ----------------------------------------------------------
    if args.ligands:
        ligands = [LIG_DIR / f"{s}.pdbqt" for s in args.ligands]
    else:
        ligands = sorted(LIG_DIR.glob("*.pdbqt"))
    ligands = [l for l in ligands if l.exists()]
    if not ligands:
        raise SystemExit("no ligand PDBQT files - run s03_ligands.py first")

    center, how = resolve_center(args)
    box = box_from_center(center, args.padding)
    print(f"Pocket centre: ({box['center_x']}, {box['center_y']}, {box['center_z']})")
    print(f"  anchor: {how}")
    print(f"  box: {box['size_x']} x {box['size_y']} x {box['size_z']} A")
    print(f"\n{len(receptors)} receptor(s) x {len(ligands)} ligand(s) "
          f"x {args.replicates} replicate(s) = "
          f"{len(receptors)*len(ligands)*args.replicates} runs\n")

    rows = []
    for rec in receptors:
        rec_name = rec.stem.replace("_receptor", "")
        for lig in ligands:
            lig_name = lig.stem
            best_per_rep = []
            for rep in range(args.replicates):
                seed = VINA_SEED + rep * 1000
                out_pose = DOCK_DIR / f"{rec_name}__{lig_name}__rep{rep}.pdbqt"
                try:
                    modes = run_vina(rec, lig, box, out_pose,
                                     args.exhaustiveness, seed)
                except Exception as exc:
                    print(f"  {rec_name:<14} {lig_name:<18} rep{rep} FAILED: {exc}",
                          file=sys.stderr)
                    continue
                best_per_rep.append(modes[0][1])
                rows.append({
                    "receptor": rec_name, "ligand": lig_name, "replicate": rep,
                    "seed": seed, "best_affinity": modes[0][1],
                    "n_modes": len(modes),
                    "pose_file": out_pose.name,
                })
            if best_per_rep:
                mean = sum(best_per_rep) / len(best_per_rep)
                spread = max(best_per_rep) - min(best_per_rep)
                print(f"  {rec_name:<14} {lig_name:<18} "
                      f"best {min(best_per_rep):>6.2f}  "
                      f"mean {mean:>6.2f}  spread {spread:>4.2f} kcal/mol")

    if not rows:
        print("\nNo docking run succeeded. Check that the `vina` binary is on "
              "PATH and that the box overlaps the receptor.", file=sys.stderr)
        return 1

    out_csv = DOCK_DIR / args.out
    fields = ["receptor", "ligand", "replicate", "seed",
              "best_affinity", "n_modes", "pose_file"]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    (DOCK_DIR / "box.json").write_text(json.dumps(
        {"box": box, "anchor": how, "padding": args.padding}, indent=2))

    print(f"\n{len(rows)} docking runs -> {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
