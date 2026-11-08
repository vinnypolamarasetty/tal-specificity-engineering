"""
Step 02 - inspect and prepare the receptor.

Crystal structures are not docking inputs. This step reports what is actually
in the file, then produces a clean, protonated receptor and a Vina-ready PDBQT.

Two things matter for an MIO-family ammonia-lyase and are handled explicitly:

  1. The biological unit is a homotetramer and the active site sits at a
     subunit interface, so docking a single isolated chain removes part of
     the pocket wall. Chain selection is therefore explicit, not assumed.

  2. The catalytic electrophile is the MIO group, formed autocatalytically
     from an internal Ala-Ser-Gly tripeptide. Depending on the entry it is
     deposited either as three standard residues or as one modified residue.
     The inspector reports which, because it changes how the pocket centre
     is defined.

Usage:
    python src/s02_receptor.py --inspect
    python src/s02_receptor.py --chains A B
"""

import argparse
import json
import sys
from collections import Counter, defaultdict

from config import ACTIVE_SITE_RADIUS, PRIMARY_PDB, PROFILE_DIR, STRUCT_DIR

# Residues that are crystallisation artefacts rather than biology.
CRYO_ADDITIVES = {
    "HOH", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "ACT", "FMT",
    "CL", "NA", "K", "MG", "CA", "ZN", "TRS", "IMD", "DMS", "NO3", "IOD",
    "BME", "MES", "EPE", "CIT", "TLA", "ACY", "PG4", "1PE", "P6G", "DTT",
}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",
}


def parse_pdb(path):
    """Minimal PDB parser - returns per-chain residue inventory."""
    chains = defaultdict(list)
    het = Counter()
    het_atoms = defaultdict(list)
    seen = set()

    for line in open(path):
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        resname = line[17:20].strip()
        chain = line[21].strip() or "_"
        resseq = line[22:27].strip()
        key = (chain, resseq, resname)

        if rec == "HETATM" and resname not in STANDARD_AA:
            if key not in seen:
                het[resname] += 1
            het_atoms[resname].append((
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            ))
        if key not in seen:
            seen.add(key)
            if resname in STANDARD_AA:
                chains[chain].append((resseq, resname))
    return chains, het, het_atoms


def inspect(pdb_id):
    path = STRUCT_DIR / f"{pdb_id}.pdb"
    if not path.exists():
        print(f"ERROR: {path} not found. Run s01_fetch.py first.", file=sys.stderr)
        return None

    chains, het, het_atoms = parse_pdb(path)

    print(f"\n=== {pdb_id} composition ===")
    print(f"\nProtein chains: {len(chains)}")
    for ch, residues in sorted(chains.items()):
        nums = [int("".join(c for c in r[0] if c.isdigit() or c == "-"))
                for r in residues if any(c.isdigit() for c in r[0])]
        rng = f"{min(nums)}-{max(nums)}" if nums else "?"
        print(f"  chain {ch}: {len(residues):>4d} residues  (numbered {rng})")

    ligands = {k: v for k, v in het.items() if k not in CRYO_ADDITIVES}
    additives = {k: v for k, v in het.items() if k in CRYO_ADDITIVES}

    print(f"\nNon-standard / ligand residues: {len(ligands)}")
    for name, count in sorted(ligands.items(), key=lambda x: -x[1]):
        print(f"  {name:<6} x{count:<4} <- candidate ligand or modified residue")
    if not ligands:
        print("  (none - this is an apo structure)")

    print(f"\nCrystallisation additives / solvent: {len(additives)} types")
    if additives:
        print("  " + ", ".join(f"{k}x{v}" for k, v in sorted(additives.items())))

    report = {
        "pdb_id": pdb_id,
        "n_chains": len(chains),
        "chains": {ch: len(r) for ch, r in chains.items()},
        "ligands": ligands,
        "additives": additives,
    }
    out = PROFILE_DIR / f"{pdb_id}_composition.json"
    out.write_text(json.dumps(report, indent=2))
    print(f"\nReport -> {out}")
    return report


def prepare(pdb_id, chains, keep_ligand=None):
    """Clean, protonate and convert the receptor to PDBQT."""
    import subprocess

    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    src = STRUCT_DIR / f"{pdb_id}.pdb"
    if not src.exists():
        raise FileNotFoundError(f"{src} not found - run s01_fetch.py")

    # --- chain selection -------------------------------------------------
    sel = STRUCT_DIR / f"{pdb_id}_sel.pdb"
    wanted = set(chains) if chains else None
    kept = 0
    with open(sel, "w") as out:
        for line in open(src):
            if line[:6].strip() in ("ATOM", "HETATM"):
                ch = line[21].strip() or "_"
                resname = line[17:20].strip()
                if wanted and ch not in wanted:
                    continue
                if resname in CRYO_ADDITIVES:
                    continue
                if keep_ligand is None and resname not in STANDARD_AA:
                    continue
                if keep_ligand and resname not in STANDARD_AA and resname != keep_ligand:
                    continue
                out.write(line)
                kept += 1
            elif line.startswith("TER"):
                out.write(line)
        out.write("END\n")
    print(f"  chain selection {sorted(wanted) if wanted else 'all'}: {kept} atom records")

    # --- repair and protonate -------------------------------------------
    fixer = PDBFixer(filename=str(sel))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    # Do not model in unresolved loops: they are guesses and would add
    # phantom atoms near the surface.
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    nonstd = list(fixer.nonstandardResidues)
    if nonstd:
        print(f"  replacing {len(nonstd)} non-standard residue(s): "
              + ", ".join(f"{r[0].name}->{r[1]}" for r in nonstd[:6]))
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    n_missing = sum(len(v) for v in fixer.missingAtoms.values())
    if n_missing:
        print(f"  rebuilding {n_missing} missing heavy atoms")
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    clean = STRUCT_DIR / f"{pdb_id}_clean.pdb"
    with open(clean, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    n_atoms = fixer.topology.getNumAtoms()
    print(f"  protonated at pH 7.4 -> {clean.name} ({n_atoms} atoms)")

    # --- PDBQT ------------------------------------------------------------
    base = STRUCT_DIR / f"{pdb_id}_receptor"
    cmd = ["mk_prepare_receptor.py", "--read_pdb", str(clean),
           "-o", str(base), "-p", "--box_size", "20", "20", "20",
           "--box_center", "0", "0", "0"]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    pdbqt = base.with_suffix(".pdbqt")
    if not pdbqt.exists():
        print(res.stdout[-1500:], file=sys.stderr)
        print(res.stderr[-1500:], file=sys.stderr)
        raise RuntimeError("mk_prepare_receptor failed to write a PDBQT")
    print(f"  receptor PDBQT -> {pdbqt.name}")
    return str(pdbqt)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", default=PRIMARY_PDB)
    ap.add_argument("--inspect", action="store_true",
                    help="report composition only, do not prepare")
    ap.add_argument("--chains", nargs="*", default=None,
                    help="chain IDs to keep (default: all)")
    ap.add_argument("--keep-ligand", default=None,
                    help="HET code of a ligand to retain (for pocket definition)")
    args = ap.parse_args()

    pdb_id = args.pdb.upper()
    report = inspect(pdb_id)
    if report is None:
        return 1
    if args.inspect:
        return 0

    print(f"\n=== preparing receptor ({pdb_id}) ===")
    prepare(pdb_id, args.chains, args.keep_ligand)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
