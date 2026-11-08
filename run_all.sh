#!/usr/bin/env bash
#
# End-to-end pipeline driver.
#
# Runs the whole workflow in order and stops at the first failure, so a broken
# step is never papered over by the next one. The ground-truth check in step 07
# is the gate: if H89F is not recovered, the saturation scan is pointless and
# the script says so.
#
# Usage:
#   ./run_all.sh                 # full run
#   ./run_all.sh --quick         # fewer replicates, lower exhaustiveness
#
set -euo pipefail

cd "$(dirname "$0")"
PY=${PYTHON:-python3}
PDB=${PDB:-2O6Y}
CHAINS=${CHAINS:-"A B"}
SITE=${SITE:-89}

QUICK=0
[[ "${1:-}" == "--quick" ]] && QUICK=1

if [[ $QUICK -eq 1 ]]; then
  EXH=8;  REPS=2; NS=0.1
else
  EXH=16; REPS=3; NS=0.5
fi

banner() { printf '\n=========== %s ===========\n' "$1"; }

command -v vina >/dev/null 2>&1 || {
  echo "ERROR: the 'vina' binary is not on PATH. See README for the download." >&2
  exit 1
}

banner "01  fetch structures and sequence"
$PY src/s01_fetch.py

banner "02  inspect and prepare receptor"
$PY src/s02_receptor.py --pdb "$PDB" --inspect
$PY src/s02_receptor.py --pdb "$PDB" --chains $CHAINS

banner "03  build ligands"
$PY src/s03_ligands.py

banner "05  profile the active site"
$PY src/s05_profile.py --pdb "$PDB" \
    --structure "data/structures/${PDB}_clean.pdb" \
    --center-from-residues "$SITE" --chain "${CHAINS%% *}"

banner "06  build the ground-truth mutant (H89F)"
$PY src/s06_mutate.py --pdb "$PDB" --chains $CHAINS --mutations "H${SITE}F"

banner "04  dock wild type and H89F"
$PY src/s04_dock.py --all-mutants \
    --ligands l_tyrosine l_phenylalanine \
    --center-source "data/structures/${PDB}_clean.pdb" \
    --center-from-residues "$SITE" --chain "${CHAINS%% *}" \
    --exhaustiveness "$EXH" --replicates "$REPS" \
    --out control_check.csv

banner "07  ground-truth check"
if ! $PY src/s07_selectivity.py --input control_check.csv | tee /tmp/control.txt; then
  echo "selectivity analysis failed" >&2; exit 1
fi

if ! grep -q "RECOVERED" /tmp/control.txt; then
  cat >&2 <<'MSG'

The pipeline did not reproduce the published H89F specificity shift.
Stopping here on purpose: a saturation scan built on a method that fails
its own control would produce confident-looking noise.

Things worth checking before going further:
  * is the docking box actually centred on the substrate channel?
  * were the interface chains kept, or was a single subunit used?
  * does residue 89 in this construct correspond to the literature numbering?
MSG
  exit 2
fi

banner "06b saturation scan at the specificity position"
$PY src/s06_mutate.py --pdb "$PDB" --chains $CHAINS --scan "$SITE" \
    --to A V L I M F Y W S T N Q

banner "04b dock the scan"
$PY src/s04_dock.py --all-mutants \
    --ligands l_tyrosine l_phenylalanine \
    --center-source "data/structures/${PDB}_clean.pdb" \
    --center-from-residues "$SITE" --chain "${CHAINS%% *}" \
    --exhaustiveness "$EXH" --replicates "$REPS" \
    --out docking_results.csv

banner "07b rank variants"
$PY src/s07_selectivity.py --input docking_results.csv

banner "09  library design"
$PY src/s09_library.py --positions "H${SITE}" --codon NNK

banner "done"
echo "MD triage is not run automatically - it is the expensive step."
echo "Run it on the shortlist, ideally on a GPU:"
echo "  $PY src/s08_md.py --all-mutants --ns $NS --site $SITE"
