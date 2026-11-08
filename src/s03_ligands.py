"""
Step 03 - build 3D ligands and convert them to PDBQT for AutoDock Vina.

Substrates and products are built from SMILES with RDKit, protonated for
pH ~7.5 (alpha-carboxylate deprotonated, alpha-amine protonated), embedded
with ETKDGv3 and minimised with MMFF94s, then written to PDBQT with Meeko.

The L-Tyr / L-Phe pair is the one that matters: the difference in their
docking scores is the selectivity readout the whole pipeline is built around.

Usage:
    python src/s03_ligands.py
"""

import json
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

from config import LIG_DIR, LIGANDS, RANDOM_SEED

RDLogger.DisableLog("rdApp.warning")


def protonate_for_ph74(mol: "Chem.Mol") -> "Chem.Mol":
    """Apply the dominant charge state at physiological pH.

    Amino acids are zwitterions at pH 7.4: COOH -> COO-, NH2 -> NH3+.
    Free cinnamic acids are simply deprotonated to the carboxylate.
    Phenolic OH (pKa ~10) stays neutral.
    """
    rxn_acid = AllChem.ReactionFromSmarts(
        "[CX3:1](=[OX1:2])[OX2H1:3]>>[C:1](=[O:2])[O-:3]"
    )
    # Only the alpha-amine of an amino acid, not an aromatic amine.
    rxn_amine = AllChem.ReactionFromSmarts(
        "[NX3;H2;!$(N-[!#6]);$(N[CX4][CX3](=O)[OX1,OX2]):1]>>[N+;H3:1]"
    )

    for rxn in (rxn_acid, rxn_amine):
        changed = True
        while changed:
            changed = False
            products = rxn.RunReactants((mol,))
            if products:
                cand = products[0][0]
                try:
                    Chem.SanitizeMol(cand)
                except Exception:
                    break
                if Chem.MolToSmiles(cand) != Chem.MolToSmiles(mol):
                    mol = cand
                    changed = True
    return mol


def build_3d(name: str, smiles: str) -> "tuple[Chem.Mol, dict]":
    """SMILES -> protonated, embedded, minimised 3D molecule."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"{name}: RDKit could not parse SMILES {smiles!r}")

    mol = protonate_for_ph74(mol)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    params.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise RuntimeError(f"{name}: 3D embedding failed")

    # MMFF94s is the variant parameterised for planar/conjugated systems,
    # which matters for the cinnamate double bond and the phenol rings.
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    if props is not None:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)
        ff.Minimize(maxIts=2000)
        energy = ff.CalcEnergy()
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
        energy = float("nan")

    meta = {
        "name": name,
        "input_smiles": smiles,
        "protonated_smiles": Chem.MolToSmiles(Chem.RemoveHs(mol)),
        "formal_charge": Chem.GetFormalCharge(mol),
        "mw": round(Descriptors.MolWt(mol), 2),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "mmff_energy_kcal": round(energy, 3),
    }
    return mol, meta


def to_pdbqt(mol: "Chem.Mol", name: str) -> str:
    """Write a Vina-ready PDBQT via Meeko."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    setup = setups[0] if isinstance(setups, list) else setups

    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise RuntimeError(f"{name}: Meeko failed to write PDBQT - {err}")

    slug = name.replace(" ", "_").replace("-", "_").lower()
    dest = LIG_DIR / f"{slug}.pdbqt"
    dest.write_text(pdbqt)

    sdf = LIG_DIR / f"{slug}.sdf"
    w = Chem.SDWriter(str(sdf))
    w.write(mol)
    w.close()
    return str(dest)


def main() -> int:
    print(f"Building {len(LIGANDS)} ligands (RDKit ETKDGv3 + MMFF94s -> Meeko PDBQT)\n")
    records, failures = [], []

    for name, smiles in LIGANDS.items():
        try:
            mol, meta = build_3d(name, smiles)
            path = to_pdbqt(mol, name)
            meta["pdbqt"] = path.split("/")[-1]
            records.append(meta)
            print(f"  {name:<18} charge {meta['formal_charge']:+d}  "
                  f"MW {meta['mw']:>6.2f}  heavy {meta['heavy_atoms']:>2d}  "
                  f"rot {meta['rotatable_bonds']}  -> {meta['pdbqt']}")
        except Exception as exc:
            print(f"  {name:<18} FAILED - {type(exc).__name__}: {exc}", file=sys.stderr)
            failures.append(name)

    out = LIG_DIR / "ligand_manifest.json"
    out.write_text(json.dumps(records, indent=2))
    print(f"\nManifest -> {out}")
    print(f"{len(records)}/{len(LIGANDS)} ligands prepared")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
