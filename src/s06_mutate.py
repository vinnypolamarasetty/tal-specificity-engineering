"""
Step 06 - build mutant receptors.

Point mutants are generated with PDBFixer, which deletes the old side chain,
grafts the new one and rebuilds missing atoms, then the result is reprotonated
and converted to PDBQT exactly like the wild type so that scores stay
comparable.

A caveat worth stating plainly: this is side-chain grafting, not repacking.
PDBFixer places the new side chain in a default rotamer and does not relax the
surrounding shell, so a substitution that introduces a large residue can leave
a strained pocket. That is acceptable for ranking candidates, and it is the
reason step 08 exists - the shortlisted variants get an MD relaxation before
anything is believed.

Usage:
    python src/s06_mutate.py --pdb 2O6Y --chains A B --mutations H89F
    python src/s06_mutate.py --pdb 2O6Y --chains A B --scan 89 --to A V L F Y W S T
"""

import argparse
import json
import subprocess
import sys

from config import PRIMARY_PDB, STRUCT_DIR

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}
THREE_TO_ONE = {v: k for k, v in ONE_TO_THREE.items()}
ALL_AA = "ACDEFGHIKLMNPQRSTVWY"


def parse_mutation(mut):
    """'H89F' -> ('HIS', 89, 'PHE'). Raises on malformed input."""
    mut = mut.strip().upper()
    if len(mut) < 3 or not mut[0].isalpha() or not mut[-1].isalpha():
        raise ValueError(f"malformed mutation {mut!r}, expected e.g. H89F")
    wt, new, num = mut[0], mut[-1], mut[1:-1]
    if not num.isdigit():
        raise ValueError(f"malformed mutation {mut!r}, residue number not numeric")
    if wt not in ONE_TO_THREE or new not in ONE_TO_THREE:
        raise ValueError(f"unknown amino acid code in {mut!r}")
    return ONE_TO_THREE[wt], int(num), ONE_TO_THREE[new]


def residue_at(pdb_path, resnum, chain):
    """What residue is actually at this position? Guards against silent mismatches."""
    for line in open(pdb_path):
        if line[:6].strip() not in ("ATOM", "HETATM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        if line[22:26].strip() != str(resnum):
            continue
        if chain and (line[21].strip() or "_") != chain:
            continue
        return line[17:20].strip()
    return None


def build_mutant(pdb_id, mutation, chains, verbose=True, relax=True):
    """Produce a PDBQT for one point mutant. Returns its path, or None."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    wt3, num, new3 = parse_mutation(mutation)
    clean = STRUCT_DIR / f"{pdb_id}_clean.pdb"
    if not clean.exists():
        raise FileNotFoundError(f"{clean} missing - run s02_receptor.py first")

    target_chain = chains[0] if chains else None
    actual = residue_at(clean, num, target_chain)
    if actual is None:
        raise ValueError(f"{mutation}: no residue numbered {num} in chain "
                         f"{target_chain or 'any'}")
    if actual != wt3:
        raise ValueError(
            f"{mutation}: position {num} is {THREE_TO_ONE.get(actual, actual)}"
            f" ({actual}) in the structure, not {THREE_TO_ONE[wt3]} ({wt3}). "
            "Check the numbering offset between the construct and UniProt.")

    if wt3 == new3:
        if verbose:
            print(f"  {mutation}: identity substitution, skipped")
        return None

    fixer = PDBFixer(filename=str(clean))
    # Apply to every requested chain so a symmetric oligomer stays symmetric.
    for ch in (chains or [target_chain]):
        try:
            fixer.applyMutations([f"{wt3}-{num}-{new3}"], ch)
        except Exception as exc:
            raise RuntimeError(f"{mutation}: PDBFixer rejected the mutation "
                               f"on chain {ch} - {exc}") from exc

    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    topology, positions = fixer.topology, fixer.positions

    # Relax the grafted side chain and its immediate shell. Without this a
    # bulky substitution is left in an idealised rotamer that often clashes,
    # and badly distorted rings break bond perception downstream.
    if relax:
        from relax import relax_site
        try:
            topology, positions = relax_site(
                topology, positions, num, chains=chains, verbose=verbose)
        except Exception as exc:
            print(f"    warning: relaxation skipped for {mutation} ({exc})",
                  file=sys.stderr)

    tag = f"{pdb_id}_{mutation.upper()}"
    out_pdb = STRUCT_DIR / f"{tag}_clean.pdb"
    with open(out_pdb, "w") as fh:
        PDBFile.writeFile(topology, positions, fh, keepIds=True)

    base = STRUCT_DIR / f"{tag}_receptor"
    res = subprocess.run(
        ["mk_prepare_receptor.py", "--read_pdb", str(out_pdb), "-o", str(base),
         "-p", "--box_size", "20", "20", "20", "--box_center", "0", "0", "0"],
        capture_output=True, text=True, timeout=900)
    pdbqt = base.with_suffix(".pdbqt")
    if not pdbqt.exists():
        raise RuntimeError(f"{mutation}: receptor conversion failed - "
                           f"{res.stderr[-400:]}")
    if verbose:
        print(f"  {mutation}: {THREE_TO_ONE[wt3]}{num}{THREE_TO_ONE[new3]} "
              f"-> {pdbqt.name}")
    return str(pdbqt)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", default=PRIMARY_PDB)
    ap.add_argument("--chains", nargs="*", default=None,
                    help="chains to mutate (all copies of a homo-oligomer)")
    ap.add_argument("--mutations", nargs="*", default=None,
                    help="explicit mutations, e.g. H89F H89L")
    ap.add_argument("--no-relax", action="store_true",
                    help="skip restrained minimisation after grafting")
    ap.add_argument("--scan", type=int, default=None,
                    help="saturate this residue position")
    ap.add_argument("--to", nargs="*", default=None,
                    help="residues to scan to (default: all 19)")
    args = ap.parse_args()

    pdb_id = args.pdb.upper()
    clean = STRUCT_DIR / f"{pdb_id}_clean.pdb"

    mutations = list(args.mutations or [])
    if args.scan is not None:
        chain = args.chains[0] if args.chains else None
        wt3 = residue_at(clean, args.scan, chain)
        if wt3 is None:
            print(f"ERROR: no residue {args.scan} found", file=sys.stderr)
            return 1
        wt1 = THREE_TO_ONE.get(wt3)
        targets = [a.upper() for a in (args.to or list(ALL_AA))]
        mutations += [f"{wt1}{args.scan}{a}" for a in targets if a != wt1]
        print(f"Saturation scan at position {args.scan} "
              f"(wild type {wt1}{args.scan}): {len(mutations)} mutants\n")

    if not mutations:
        print("Nothing to do - give --mutations or --scan", file=sys.stderr)
        return 1

    built, failed = [], []
    for mut in mutations:
        try:
            path = build_mutant(pdb_id, mut, args.chains, relax=not args.no_relax)
            if path:
                built.append({"mutation": mut.upper(), "pdbqt": path})
        except Exception as exc:
            print(f"  {mut}: FAILED - {exc}", file=sys.stderr)
            failed.append(mut)

    manifest = STRUCT_DIR / f"{pdb_id}_mutants.json"
    manifest.write_text(json.dumps(built, indent=2))
    print(f"\n{len(built)} mutant receptors built, {len(failed)} failed")
    print(f"Manifest -> {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
