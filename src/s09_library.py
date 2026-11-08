"""
Step 09 - variant library design.

Converts a shortlist of positions into something a wet lab can actually order
and screen, and - just as importantly - says how many clones must be picked
to have a real chance of seeing every variant.

Degenerate codons
-----------------
NNK (N = A/C/G/T, K = G/T) is the standard saturation codon: 32 codons
covering all 20 amino acids with a single stop (TAG), versus 64 codons and
3 stops for NNN. For a focused library, reduced codons cut the screening
burden sharply:

  NDT (12 codons, 12 aa, no stop)  - balanced polar/nonpolar/charged set
  DBK (18 codons, 18 aa, no stop)  - broad, still stop-free
  NNK (32 codons, 20 aa, 1 stop)   - full saturation

Screening depth
---------------
For a library of V distinct DNA sequences, picking N clones at random gives
expected coverage P = 1 - (1 - 1/V)^N. Solving for N at a target P gives the
number of transformants to screen. This grows fast with the number of
simultaneously randomised positions, which is the practical argument for
keeping combinatorial libraries small and rational.

Usage:
    python src/s09_library.py --positions H89 F37 L108 --codon NNK
    python src/s09_library.py --positions H89 --codon NDT --coverage 0.99
"""

import argparse
import itertools
import json
import math

from config import LIB_DIR

BASES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

CODON_TABLE = {}
_bases = "TCAG"
_aas = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for _i, (_a, _b, _c) in enumerate(itertools.product(_bases, repeat=3)):
    CODON_TABLE[_a + _b + _c] = _aas[_i]

PRESETS = {
    "NNK": "NNK", "NNS": "NNS", "NNN": "NNN",
    "NDT": "NDT", "DBK": "DBK", "VHG": "VHG",
}


def expand(degenerate):
    """All concrete codons encoded by a degenerate codon."""
    degenerate = degenerate.upper()
    if len(degenerate) != 3:
        raise ValueError(f"codon {degenerate!r} must be 3 characters")
    pools = []
    for ch in degenerate:
        if ch not in BASES:
            raise ValueError(f"unknown degenerate base {ch!r}")
        pools.append(BASES[ch])
    return ["".join(c) for c in itertools.product(*pools)]


def profile_codon(degenerate):
    """Amino-acid composition and stop fraction of a degenerate codon."""
    codons = expand(degenerate)
    counts = {}
    stops = 0
    for c in codons:
        aa = CODON_TABLE[c]
        if aa == "*":
            stops += 1
        counts[aa] = counts.get(aa, 0) + 1
    aa_only = {k: v for k, v in counts.items() if k != "*"}
    return {
        "codon": degenerate.upper(),
        "n_codons": len(codons),
        "n_amino_acids": len(aa_only),
        "amino_acids": "".join(sorted(aa_only)),
        "n_stop_codons": stops,
        "stop_fraction": round(stops / len(codons), 4),
        "codon_bias": {k: aa_only[k] for k in sorted(aa_only)},
    }


def clones_for_coverage(diversity, coverage):
    """Transformants needed for the target expected coverage."""
    if not 0 < coverage < 1:
        raise ValueError("coverage must be between 0 and 1")
    return math.ceil(math.log(1 - coverage) / math.log(1 - 1 / diversity))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--positions", nargs="+", required=True,
                    help="positions to randomise, e.g. H89 F37")
    ap.add_argument("--codon", default="NNK",
                    help="degenerate codon (NNK, NDT, DBK, NNS, NNN, or custom)")
    ap.add_argument("--coverage", type=float, default=0.95)
    ap.add_argument("--out", default="library_design.json")
    args = ap.parse_args()

    codon = PRESETS.get(args.codon.upper(), args.codon.upper())
    prof = profile_codon(codon)
    n_pos = len(args.positions)

    dna_diversity = prof["n_codons"] ** n_pos
    protein_diversity = prof["n_amino_acids"] ** n_pos
    functional = (1 - prof["stop_fraction"]) ** n_pos

    clones = clones_for_coverage(dna_diversity, args.coverage)

    print(f"=== library design ===\n")
    print(f"Randomised positions : {', '.join(args.positions)}  (n = {n_pos})")
    print(f"Degenerate codon     : {prof['codon']}")
    print(f"  codons per position : {prof['n_codons']}")
    print(f"  amino acids encoded : {prof['n_amino_acids']} "
          f"({prof['amino_acids']})")
    print(f"  stop codons         : {prof['n_stop_codons']} "
          f"({prof['stop_fraction']*100:.1f}% per position)")

    print(f"\nLibrary size")
    print(f"  DNA diversity       : {dna_diversity:,}")
    print(f"  protein diversity   : {protein_diversity:,}")
    print(f"  stop-free fraction  : {functional*100:.1f}%")

    print(f"\nScreening depth")
    print(f"  clones for {args.coverage*100:.0f}% coverage : {clones:,}")
    for cov in (0.90, 0.95, 0.99):
        print(f"    {cov*100:>4.0f}% -> {clones_for_coverage(dna_diversity, cov):>12,} "
              "transformants")

    if n_pos > 2 and prof["n_codons"] >= 32:
        print(f"\n  Note: {n_pos} NNK positions needs {clones:,} clones for "
              f"{args.coverage*100:.0f}% coverage, which is beyond routine "
              "plate-based screening.")
        ndt = profile_codon("NDT")
        ndt_div = ndt["n_codons"] ** n_pos
        print(f"  Switching to NDT would cut this to "
              f"{clones_for_coverage(ndt_div, args.coverage):,} clones "
              f"({ndt['n_amino_acids']} aa per position), or randomise "
              "positions in separate sub-libraries.")

    out = LIB_DIR / args.out
    out.write_text(json.dumps({
        "positions": args.positions,
        "codon_profile": prof,
        "dna_diversity": dna_diversity,
        "protein_diversity": protein_diversity,
        "stop_free_fraction": round(functional, 4),
        "coverage_target": args.coverage,
        "clones_required": clones,
        "clones_by_coverage": {str(c): clones_for_coverage(dna_diversity, c)
                               for c in (0.90, 0.95, 0.99)},
    }, indent=2))
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
