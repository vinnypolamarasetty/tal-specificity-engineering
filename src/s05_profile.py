"""
Step 05 - active-site profiling.

Given a pocket anchor, enumerate the residues that line it and classify them.
This is what turns "the active site" into an explicit, citable list of
positions, and it is how the mutagenesis targets in step 06 are chosen rather
than guessed.

Output:
  * the shell of residues within ACTIVE_SITE_RADIUS of the pocket centre,
  * each residue's distance, chain and physicochemical class,
  * a shortlist of mutable positions, excluding catalytic machinery.

Usage:
    python src/s05_profile.py --pdb 2O6Y --center-from-residues 89 --chain A
    python src/s05_profile.py --pdb 2O7D --center-from-ligand HCI
"""

import argparse
import csv
import json
import sys

from config import ACTIVE_SITE_RADIUS, PROFILE_DIR, PRIMARY_PDB, STRUCT_DIR
from pocket import (center_from_ligand, center_from_residues, find_mio,
                    residues_within)

# Physicochemical classes, used to reason about what a substitution would do
# to the pocket rather than just listing names.
CLASSES = {
    "ALA": "small hydrophobic", "GLY": "flexible", "VAL": "hydrophobic",
    "LEU": "hydrophobic", "ILE": "hydrophobic", "MET": "hydrophobic",
    "PRO": "rigid", "PHE": "aromatic", "TRP": "aromatic", "TYR": "aromatic/polar",
    "SER": "polar", "THR": "polar", "CYS": "polar/reactive", "ASN": "polar",
    "GLN": "polar", "HIS": "ionisable/aromatic", "LYS": "basic", "ARG": "basic",
    "ASP": "acidic", "GLU": "acidic",
}

# Residues that should not be mutated casually: the catalytic machinery.
PROTECTED_HINT = {"ASG(unmodified)", "MDO", "MIO"}


def classify(resname):
    return CLASSES.get(resname, "non-standard")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", default=PRIMARY_PDB)
    ap.add_argument("--structure", default=None,
                    help="explicit PDB path (default: data/structures/<pdb>.pdb)")
    ap.add_argument("--center-from-ligand", default=None, metavar="HET")
    ap.add_argument("--center-from-residues", nargs="*", type=int, default=None)
    ap.add_argument("--chain", default=None)
    ap.add_argument("--radius", type=float, default=ACTIVE_SITE_RADIUS)
    args = ap.parse_args()

    pdb_id = args.pdb.upper()
    path = args.structure or (STRUCT_DIR / f"{pdb_id}.pdb")
    if not str(path):
        return 1

    try:
        if args.center_from_ligand:
            centre, n = center_from_ligand(path, args.center_from_ligand)
            anchor = f"ligand {args.center_from_ligand} ({n} atoms)"
        elif args.center_from_residues:
            centre, n = center_from_residues(path, args.center_from_residues,
                                             args.chain)
            anchor = f"residues {args.center_from_residues} ({n} atoms)"
        else:
            print("ERROR: give --center-from-ligand or --center-from-residues",
                  file=sys.stderr)
            return 1
    except (ValueError, FileNotFoundError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"=== active-site profile: {pdb_id} ===")
    print(f"anchor: {anchor}")
    print(f"centre: ({centre[0]:.2f}, {centre[1]:.2f}, {centre[2]:.2f})")

    mio = find_mio(path, args.chain)
    if mio:
        print(f"MIO group: {mio[0]} at chain {mio[1]} residue {mio[2]}")
    else:
        print("MIO group: not located (check deposition style)")

    shell = residues_within(path, centre, args.radius, args.chain)
    print(f"\n{len(shell)} residues within {args.radius} A:\n")
    print(f"  {'chain':<6}{'resi':<7}{'name':<6}{'dist':<7}class")
    print("  " + "-" * 46)

    rows = []
    for chain, resseq, resname, dist in shell:
        cls = classify(resname)
        print(f"  {chain:<6}{resseq:<7}{resname:<6}{dist:<7.2f}{cls}")
        rows.append({"chain": chain, "resseq": resseq, "resname": resname,
                     "distance": dist, "class": cls})

    # Mutable shortlist: side chains that line the pocket, excluding glycine
    # (no side chain to swap), proline (backbone constraint) and catalytic parts.
    mutable = [r for r in rows
               if r["resname"] not in ("GLY", "PRO")
               and r["resname"] in CLASSES]

    print(f"\n{len(mutable)} candidate positions for mutagenesis "
          "(glycine/proline and non-standard residues excluded)")

    out_csv = PROFILE_DIR / f"{pdb_id}_active_site.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["chain", "resseq", "resname",
                                           "distance", "class"])
        w.writeheader()
        w.writerows(rows)

    out_json = PROFILE_DIR / f"{pdb_id}_active_site.json"
    out_json.write_text(json.dumps({
        "pdb_id": pdb_id,
        "anchor": anchor,
        "centre": [round(float(c), 3) for c in centre],
        "radius": args.radius,
        "mio": mio,
        "shell": rows,
        "mutable_positions": [f"{r['resname']}{r['resseq']}" for r in mutable],
    }, indent=2))

    print(f"\n-> {out_csv}")
    print(f"-> {out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
