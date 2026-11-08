"""
Step 01 - fetch the target structures and sequences.

Pulls RsTAL crystal structures from the RCSB PDB and the reference sequence
from UniProt. Everything is cached under data/ so the rest of the pipeline is
reproducible offline once this step has run.

Usage:
    python src/s01_fetch.py
    python src/s01_fetch.py --pdb 2O6Y 2O7F
"""

import argparse
import sys
import time
import urllib.error
import urllib.request

from config import PDB_TARGETS, PRIMARY_PDB, SEQ_DIR, STRUCT_DIR, UNIPROT_TAL

RCSB_FILE = "https://files.rcsb.org/download/{pdb}.pdb"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb}"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
USER_AGENT = "tal-specificity-engineering/1.0 (academic pipeline)"


def _get(url: str, retries: int = 3, timeout: int = 30) -> bytes:
    """GET with retries. Raises RuntimeError with an actionable message."""
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise RuntimeError(f"not found (404): {url}") from exc
            last = exc
        except Exception as exc:  # network / TLS / proxy
            last = exc
        if attempt < retries:
            time.sleep(2 * attempt)
    raise RuntimeError(
        f"could not reach {url} after {retries} attempts ({last}). "
        "If you are behind a restrictive proxy, allow files.rcsb.org, "
        "data.rcsb.org and rest.uniprot.org, or copy the files into data/ by hand."
    )


def fetch_pdb(pdb_id: str) -> "tuple[str, int]":
    """Download one PDB entry. Returns (path, n_atom_records)."""
    pdb_id = pdb_id.upper()
    dest = STRUCT_DIR / f"{pdb_id}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        text = dest.read_text()
        n = sum(1 for ln in text.splitlines() if ln.startswith(("ATOM", "HETATM")))
        print(f"  {pdb_id}: cached ({n} atom records)")
        return str(dest), n

    raw = _get(RCSB_FILE.format(pdb=pdb_id))
    text = raw.decode("utf-8", errors="replace")

    # Sanity check: a real PDB file has a HEADER and coordinate records.
    n = sum(1 for ln in text.splitlines() if ln.startswith(("ATOM", "HETATM")))
    if n == 0:
        raise RuntimeError(f"{pdb_id}: downloaded file has no coordinate records")

    dest.write_text(text)
    print(f"  {pdb_id}: downloaded ({n} atom records) -> {dest.name}")
    return str(dest), n


def fetch_sequence(acc: str) -> str:
    """Download the UniProt FASTA for the reference sequence."""
    dest = SEQ_DIR / f"{acc}.fasta"
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  {acc}: cached")
        return str(dest)
    raw = _get(UNIPROT_FASTA.format(acc=acc))
    text = raw.decode("utf-8", errors="replace")
    if not text.startswith(">"):
        raise RuntimeError(f"{acc}: response is not FASTA")
    dest.write_text(text)
    seq = "".join(text.splitlines()[1:])
    print(f"  {acc}: downloaded ({len(seq)} aa) -> {dest.name}")
    return str(dest)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", nargs="*", default=list(PDB_TARGETS),
                    help="PDB IDs to fetch (default: all RsTAL entries)")
    ap.add_argument("--skip-seq", action="store_true")
    args = ap.parse_args()

    print("Fetching structures from RCSB PDB")
    ok, failed = [], []
    for pdb_id in args.pdb:
        try:
            fetch_pdb(pdb_id)
            ok.append(pdb_id)
        except RuntimeError as exc:
            print(f"  {pdb_id}: FAILED - {exc}", file=sys.stderr)
            failed.append(pdb_id)

    if not args.skip_seq:
        print("\nFetching reference sequence from UniProt")
        try:
            fetch_sequence(UNIPROT_TAL)
        except RuntimeError as exc:
            print(f"  {UNIPROT_TAL}: FAILED - {exc}", file=sys.stderr)

    print(f"\n{len(ok)}/{len(args.pdb)} structures available in {STRUCT_DIR}")
    if PRIMARY_PDB.upper() not in [p.upper() for p in ok]:
        print(f"WARNING: primary scaffold {PRIMARY_PDB} is missing; "
              "downstream steps will not run.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
