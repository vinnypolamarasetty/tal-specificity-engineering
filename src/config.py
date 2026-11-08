"""
Central configuration for the TAL specificity-engineering pipeline.

Target: Tyrosine ammonia-lyase from Rhodobacter sphaeroides (RsTAL).

Why this enzyme:
  TAL catalyses the non-oxidative deamination of L-tyrosine to p-coumarate.
  In E. coli, p-coumarate is the entry metabolite for the plant polyphenol
  pathways (naringenin, resveratrol, curcuminoids), and the ammonia-lyase
  step is a well-documented flux bottleneck: the enzyme is slow, and most
  bacterial/fungal ammonia-lyases carry significant side activity on
  L-phenylalanine, draining flux into cinnamate.

Why it is a good benchmark:
  Position 89 is a published specificity determinant. His89 hydrogen bonds to
  the substrate p-hydroxyl; the H89F substitution (the residue found in
  phenylalanine ammonia-lyases) inverts the kinetic preference from L-Tyr to
  L-Phe. That gives this pipeline a ground-truth answer to recover before any
  novel variant it proposes should be believed.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
STRUCT_DIR = DATA / "structures"
SEQ_DIR = DATA / "sequences"
LIG_DIR = DATA / "ligands"
RESULTS = ROOT / "results"
DOCK_DIR = RESULTS / "docking"
PROFILE_DIR = RESULTS / "profiles"
MD_DIR = RESULTS / "md"
LIB_DIR = RESULTS / "library"
FIG_DIR = ROOT / "figures"

for _d in (STRUCT_DIR, SEQ_DIR, LIG_DIR, DOCK_DIR, PROFILE_DIR, MD_DIR, LIB_DIR, FIG_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Structures
# --------------------------------------------------------------------------
# RsTAL structures from Louie et al. (2006) Chem Biol, "Structural
# determinants and modulation of substrate specificity in
# phenylalanine-tyrosine ammonia-lyases".
PDB_TARGETS = {
    "2O6Y": "RsTAL, unliganded (reference apo scaffold)",
    "2O7B": "RsTAL complex",
    "2O7D": "RsTAL complex",
    "2O7E": "RsTAL complex",
    "2O7F": "RsTAL H89F mutant complex",
}

# The scaffold every downstream step is built on.
PRIMARY_PDB = "2O6Y"

# UniProt accession for RsTAL (TAL_CERS4 / Cereibacter sphaeroides).
UNIPROT_TAL = "Q3IWB0"

# --------------------------------------------------------------------------
# Catalytic / specificity residues
# --------------------------------------------------------------------------
# The MIO prosthetic group is formed autocatalytically from the internal
# Ala-Ser-Gly tripeptide. It is the electrophile that activates the substrate
# amine, so it defines the centre of the catalytic pocket.
MIO_TRIPEPTIDE = ("ALA", "SER", "GLY")

# Position 89 is the published specificity switch (His in TALs, Phe in PALs).
SPECIFICITY_POSITION = 89
GROUND_TRUTH_MUTATION = "H89F"  # should shift preference L-Tyr -> L-Phe

# Radius used to define the active-site shell around the bound ligand.
ACTIVE_SITE_RADIUS = 5.0  # angstrom

# --------------------------------------------------------------------------
# Ligands - SMILES at pH ~7.5 (carboxylates deprotonated, amines protonated)
# --------------------------------------------------------------------------
LIGANDS = {
    "L-tyrosine":    "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
    "L-phenylalanine": "N[C@@H](Cc1ccccc1)C(=O)O",
    "p-coumarate":   "OC(=O)/C=C/c1ccc(O)cc1",
    "cinnamate":     "OC(=O)/C=C/c1ccccc1",
    "caffeate":      "OC(=O)/C=C/c1ccc(O)c(O)c1",
}

# The pair whose score difference defines selectivity.
SUBSTRATE = "L-tyrosine"
ANTI_SUBSTRATE = "L-phenylalanine"

# --------------------------------------------------------------------------
# Docking
# --------------------------------------------------------------------------
VINA_BIN = "vina"
VINA_EXHAUSTIVENESS = 16
VINA_NUM_MODES = 9
VINA_SEED = 42          # fixed so every run in the repo is reproducible
BOX_PADDING = 11.0      # angstrom, box half-width around pocket centre

# --------------------------------------------------------------------------
# MD refinement
# --------------------------------------------------------------------------
MD_FORCEFIELD = ("amber14-all.xml", "amber14/tip3pfb.xml")
MD_TEMPERATURE = 300      # kelvin
MD_TIMESTEP_FS = 2.0
MD_EQUIL_STEPS = 25_000   # 50 ps
MD_PROD_STEPS = 250_000   # 500 ps
MD_REPORT_EVERY = 2_500

RANDOM_SEED = 42
