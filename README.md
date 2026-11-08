# Computational specificity engineering of tyrosine ammonia-lyase

A reproducible pipeline for **rational enzyme redesign**: structure preparation,
substrate docking, active-site profiling, saturation mutagenesis, MD stability
triage, and degenerate-codon library design — applied to a documented flux
bottleneck in microbial production of plant natural products.

The pipeline is **calibrated against a published result before it is used to
propose anything new.**

---

## Status

**The machinery is built and tested. It has not yet been run on RsTAL, and this
repository therefore contains no variant predictions.** That distinction is
deliberate and is kept visible rather than buried.

What has been verified, by executing it:

| Step | Verified by |
|------|-------------|
| `s02_receptor` | Ran on a crystal structure: additive stripping, missing-atom rebuild, protonation, PDBQT conversion |
| `s03_ligands`  | Built 5/5 substrates; zwitterion states and **(S) stereochemistry confirmed from the generated 3D coordinates** |
| `s05_profile`  | Ran: active-site shell enumeration and mutable-position shortlist |
| `s06_mutate`   | Ran: point mutants built; restrained relaxation resolves grafting clashes (+12,663 → −429 kcal/mol on a His→Trp test graft) |
| `s04_dock`     | Ran: AutoDock Vina 1.2.5, replicate seeds, affinity parsing |
| `s07_selectivity` | Ran: ΔΔ ranking, replicate-noise propagation, automatic control check |
| `s08_md`       | Ran: 20 ps implicit-solvent trajectory, backbone RMSD 1.06 Å, active-site RMSF and pocket radius |
| `s09_library`  | Verified against the genetic code: NNK = 32 codons / 20 aa / 1 stop; NDT = 12 / 12 / 0 |

Those runs used small control proteins to exercise the code paths, not RsTAL.

**What remains:** fetch 2O6Y, run `run_all.sh`, and see whether step 07 recovers
the published H89F specificity shift. `run_all.sh` stops with a non-zero exit
code if it does not, rather than proceeding to a saturation scan on a method
that has failed its own control.

---

## The problem

Plant polyphenols — naringenin, resveratrol, curcuminoids — are increasingly
produced in engineered *Escherichia coli* rather than extracted from plants.
Nearly all of those routes pass through a single entry metabolite,
**p-coumarate**, made by an aromatic amino acid ammonia-lyase.

That step is a bottleneck twice over:

1. **It is slow.** The ammonia-lyase is typically the lowest-flux step in the
   pathway, so strain yield is capped by its turnover.
2. **It is leaky.** Tyrosine ammonia-lyases (TAL) retain significant activity
   on L-phenylalanine, converting it to cinnamate. That drains carbon into an
   off-pathway product and competes with the cell's own aromatic amino acid
   demand.

Improving the ratio of L-Tyr to L-Phe activity is therefore a direct lever on
pathway flux, and it is a *protein engineering* problem rather than a
pathway-balancing one.

## The target

**Tyrosine ammonia-lyase from *Rhodobacter sphaeroides* (RsTAL)** — a
MIO-dependent ammonia-lyase, crystallised in apo and ligand-bound forms.

The catalytic electrophile is the **MIO group** (4-methylidene-imidazole-5-one),
formed autocatalytically from an internal Ala-Ser-Gly tripeptide. The enzyme is
a homotetramer, and the active site sits at a subunit interface — so docking an
isolated chain would delete part of the pocket wall. The pipeline handles both
facts explicitly rather than assuming a simple monomeric globular protein.

## Why this target is a benchmark, not just a case study

Position 89 is a **published specificity determinant**. His89 hydrogen bonds to
the p-hydroxyl of the tyrosine substrate; it is conserved in TALs and replaced
by phenylalanine in phenylalanine ammonia-lyases (PALs). The **H89F**
substitution shifts the kinetic preference from L-Tyr towards L-Phe.

That gives a ground-truth answer. Before any novel variant from this pipeline
is taken seriously, the pipeline must recover H89F as a Tyr→Phe shift. Step 07
performs this check automatically and labels the result `RECOVERED` or
`NOT recovered`. **A pipeline that fails its own control is reporting noise**,
and the code says so rather than quietly ranking variants anyway.

---

## The selectivity metric

Raw docking affinities are not trustworthy to better than a couple of kcal/mol.
This pipeline therefore never uses them directly. For a receptor `R`:

```
S(R)  =  ΔG_dock(L-Tyr, R)  −  ΔG_dock(L-Phe, R)       # negative ⇒ prefers Tyr
ΔS    =  S(mutant) − S(wild type)                      # positive ⇒ shifted to Phe
```

L-Tyr and L-Phe differ by one hydroxyl group, and a point mutant differs from
wild type by one side chain. Scoring-function errors are strongly correlated
across both comparisons, so most systematic error cancels in `ΔS` even though
it dominates the absolute numbers.

`ΔS` is a **ranking signal for prioritising experiments** — not a predicted
kinetic constant. Every variant it nominates is a hypothesis to be screened.

---

## Pipeline

| Step | Script | What it does |
|------|--------|--------------|
| 01 | `s01_fetch.py` | Fetch structures (RCSB) and reference sequence (UniProt), with caching |
| 02 | `s02_receptor.py` | Inventory the crystal file, strip additives, rebuild missing atoms, protonate at pH 7.4, convert to PDBQT |
| 03 | `s03_ligands.py` | Build substrates from SMILES — correct zwitterions, ETKDGv3 embedding, MMFF94s minimisation, Meeko PDBQT |
| 05 | `s05_profile.py` | Enumerate the active-site shell, classify residues, shortlist mutable positions |
| 06 | `s06_mutate.py` | Graft point mutations, then **restrained minimisation** of the mutated site |
| 04 | `s04_dock.py` | AutoDock Vina with fixed seeds and replicate runs |
| 07 | `s07_selectivity.py` | Compute `ΔS`, rank variants, run the ground-truth check |
| 08 | `s08_md.py` | Implicit-solvent MD: backbone RMSD, active-site RMSF, pocket width |
| 09 | `s09_library.py` | Degenerate-codon library design and required screening depth |

### Design decisions worth noting

**Mutant side chains are relaxed, not just grafted.** PDBFixer places a new
side chain in an idealised rotamer without repacking its surroundings. For
bulky substitutions this leaves a badly clashing pocket — in testing, a
His→Trp graft scored **+12,663 kcal/mol** before relaxation and **−429
kcal/mol** after. Unrelaxed, such structures either fail bond perception
outright or, worse, dock *successfully* and produce a meaningless score.
`relax.py` minimises atoms within 8 Å of the mutation while harmonically
restraining the rest, so the mutant stays comparable to the wild type instead
of drifting away from the experimental coordinates.

**Docking is replicated.** Vina is stochastic. Every pair is run with several
fixed seeds and the replicate spread is propagated into a signal-to-noise
column, so a 0.1 kcal/mol "improvement" sitting inside run-to-run scatter is
visibly not a result.

**Residue identity is verified before mutating.** Construct numbering routinely
differs from UniProt numbering. `s06_mutate.py` refuses to apply `H89F` unless
residue 89 really is a histidine, instead of silently mutating whatever is
there.

**MD is triage, not thermodynamics.** Short implicit-solvent runs ask whether a
variant's pocket stays folded and rigid. They are *not* binding free energies;
that would require explicit solvent and an alchemical method.

---

## Install

```bash
pip install -r requirements.txt
```

AutoDock Vina is a separate binary:

```bash
curl -L -o vina https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina && sudo mv vina /usr/local/bin/
```

## Run

```bash
python src/s01_fetch.py                                   # structures + sequence
python src/s02_receptor.py --pdb 2O6Y --inspect           # what is in the file?
python src/s02_receptor.py --pdb 2O6Y --chains A B        # prepare receptor
python src/s03_ligands.py                                 # substrates
python src/s05_profile.py --pdb 2O6Y --center-from-residues 89 --chain A
python src/s06_mutate.py --pdb 2O6Y --chains A B --mutations H89F
python src/s04_dock.py --all-mutants --center-source data/structures/2O6Y_clean.pdb \
                       --center-from-residues 89 --chain A
python src/s07_selectivity.py                             # ranking + control check
python src/s09_library.py --positions H89 --codon NNK
```

---

## Limitations

Stated plainly, because a method is only useful if its failure modes are known.

- **Rigid-receptor docking.** The receptor does not move during docking.
  Step 08 partly compensates by checking pocket stability afterwards, but
  induced fit is not modelled.
- **Empirical scoring.** Vina affinities are regression estimates, not free
  energies. `ΔS` mitigates but does not remove this.
- **Ground state, not transition state.** Ammonia-lyase catalysis proceeds
  through MIO-amine adduct chemistry. Docking scores substrate *binding*;
  selectivity in a real enzyme also depends on k_cat. A variant that binds
  L-Tyr better may still turn it over worse.
- **No repacking.** Only the mutated site is relaxed; distant compensating
  rearrangements are not sampled.
- **Implicit solvent** in step 08 omits explicit water structure, which matters
  in a polar active site.
- **Predictions are hypotheses.** Nothing here substitutes for expressing the
  variants and measuring kinetics.

## Reproducibility

All stochastic steps take fixed seeds (`config.py`). Structures, ligands and
intermediate files are cached under `data/`, and every step writes a
machine-readable record to `results/`.

## References

- Louie, Bowman, Moore & Noel (2006). Structural determinants and modulation of
  substrate specificity in phenylalanine-tyrosine ammonia-lyases. *Chemistry &
  Biology*. — RsTAL structures and the position-89 specificity switch.
- Eberhardt, Santos-Martins, Tillack & Forli (2021). AutoDock Vina 1.2.0.
  *J. Chem. Inf. Model.*
- Eastman et al. (2017). OpenMM 7. *PLoS Comput. Biol.*
