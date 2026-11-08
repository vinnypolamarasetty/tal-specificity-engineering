"""
Pocket definition helpers.

A docking box is only meaningful if it is anchored to something real. Two
anchors are supported:

  * a ligand already present in a holo crystal structure (best - it marks the
    pocket experimentally), or
  * a set of residue numbers, used when working from an apo or predicted
    structure.

For MIO-family ammonia-lyases the natural residue anchor is the catalytic
MIO group plus the position-89 specificity residue, because together they
define the two ends of the substrate channel.
"""

import numpy as np

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE",
}


def load_atoms(pdb_path):
    """Read a PDB into a list of atom dicts. Deliberately dependency-free."""
    atoms = []
    for line in open(pdb_path):
        if line[:6].strip() not in ("ATOM", "HETATM"):
            continue
        try:
            atoms.append({
                "record": line[:6].strip(),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21].strip() or "_",
                "resseq": line[22:26].strip(),
                "icode": line[26].strip(),
                "xyz": np.array([float(line[30:38]),
                                 float(line[38:46]),
                                 float(line[46:54])]),
                "element": (line[76:78].strip() or line[12:16].strip()[:1]).upper(),
            })
        except ValueError:
            continue
    return atoms


def center_from_ligand(pdb_path, het_code):
    """Centroid of a named HET group."""
    atoms = load_atoms(pdb_path)
    sel = [a for a in atoms if a["resname"] == het_code.upper()
           and a["element"] != "H"]
    if not sel:
        raise ValueError(f"no atoms found for HET code {het_code!r} in {pdb_path}")
    coords = np.array([a["xyz"] for a in sel])
    return coords.mean(axis=0), len(sel)


def center_from_residues(pdb_path, resnums, chain=None):
    """Centroid of the side chains of the given residue numbers."""
    atoms = load_atoms(pdb_path)
    want = {str(r) for r in resnums}
    backbone = {"N", "CA", "C", "O", "OXT"}
    sel = [a for a in atoms
           if a["resseq"] in want
           and (chain is None or a["chain"] == chain)
           and a["element"] != "H"
           and a["name"] not in backbone]
    if not sel:
        # Fall back to all atoms of those residues (glycine has no side chain).
        sel = [a for a in atoms if a["resseq"] in want
               and (chain is None or a["chain"] == chain)
               and a["element"] != "H"]
    if not sel:
        raise ValueError(f"no atoms found for residues {sorted(want)} in {pdb_path}")
    coords = np.array([a["xyz"] for a in sel])
    return coords.mean(axis=0), len(sel)


def find_mio(pdb_path, chain=None):
    """Locate the MIO group.

    Returns (resname, chain, resseq). Handles both deposition styles:
    a modified residue (MDO and relatives) or a bare Ala-Ser-Gly run.
    """
    atoms = load_atoms(pdb_path)
    mio_codes = {"MDO", "MIO", "TOQ", "CRO"}
    for a in atoms:
        if a["resname"] in mio_codes and (chain is None or a["chain"] == chain):
            return a["resname"], a["chain"], a["resseq"]

    # Fall back: scan for a contiguous ALA-SER-GLY run.
    by_chain = {}
    for a in atoms:
        if a["name"] != "CA" or a["resname"] not in STANDARD_AA:
            continue
        if chain is not None and a["chain"] != chain:
            continue
        by_chain.setdefault(a["chain"], []).append(a)

    for ch, cas in by_chain.items():
        def num(x):
            digits = "".join(c for c in x["resseq"] if c.isdigit() or c == "-")
            return int(digits) if digits else None
        cas = [c for c in cas if num(c) is not None]
        cas.sort(key=num)
        for i in range(len(cas) - 2):
            trio = cas[i:i + 3]
            names = [t["resname"] for t in trio]
            nums = [num(t) for t in trio]
            if names == ["ALA", "SER", "GLY"] and nums[2] - nums[0] == 2:
                return "ASG(unmodified)", ch, str(nums[1])
    return None


def residues_within(pdb_path, center, radius, chain=None):
    """All residues with a heavy atom inside `radius` of `center`."""
    atoms = load_atoms(pdb_path)
    hits = {}
    for a in atoms:
        if a["element"] == "H":
            continue
        if chain is not None and a["chain"] != chain:
            continue
        d = float(np.linalg.norm(a["xyz"] - center))
        if d <= radius:
            key = (a["chain"], a["resseq"], a["resname"])
            if key not in hits or d < hits[key]:
                hits[key] = d
    return sorted(((k[0], k[1], k[2], round(v, 2)) for k, v in hits.items()),
                  key=lambda x: x[3])


def box_from_center(center, padding):
    """Vina box parameters from a centre and half-width."""
    return {
        "center_x": round(float(center[0]), 3),
        "center_y": round(float(center[1]), 3),
        "center_z": round(float(center[2]), 3),
        "size_x": round(2 * padding, 1),
        "size_y": round(2 * padding, 1),
        "size_z": round(2 * padding, 1),
    }
