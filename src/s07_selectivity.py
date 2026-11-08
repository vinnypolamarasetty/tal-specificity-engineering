"""
Step 07 - selectivity analysis and variant ranking.

Turns raw docking scores into the quantity of interest.

Definitions
-----------
For a given receptor R:

    S(R) = dG_dock(L-Tyr, R) - dG_dock(L-Phe, R)

S is negative when the pocket favours L-tyrosine. The effect of a mutation is
then the shift relative to wild type:

    dS = S(mutant) - S(wild type)

dS > 0  means the mutation moves the pocket towards L-phenylalanine
dS < 0  means it sharpens tyrosine preference

Why a difference of differences
-------------------------------
Vina's scoring function is an empirical regression, and its absolute
affinities carry errors of a couple of kcal/mol. Those errors are strongly
correlated between two ligands as similar as L-Tyr and L-Phe, which differ by
a single hydroxyl, and between two receptors differing by one side chain. So
most of the systematic error cancels in dS even though it dominates the raw
numbers. dS is a ranking signal, not a predicted kinetic constant - the
pipeline is calibrated by checking that it recovers the known H89F result
before any novel variant is taken seriously.

Usage:
    python src/s07_selectivity.py
    python src/s07_selectivity.py --input docking_results.csv --top 15
"""

import argparse
import csv
import statistics
import sys
from collections import defaultdict

from config import (ANTI_SUBSTRATE, DOCK_DIR, GROUND_TRUTH_MUTATION, SUBSTRATE)


def slug(name):
    return name.replace(" ", "_").replace("-", "_").lower()


def load(path):
    """Group best affinities by (receptor, ligand)."""
    data = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                data[(row["receptor"], row["ligand"])].append(
                    float(row["best_affinity"]))
            except (KeyError, ValueError):
                continue
    return data


def summarise(values):
    mean = statistics.mean(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    return mean, sd, min(values)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="docking_results.csv")
    ap.add_argument("--substrate", default=SUBSTRATE)
    ap.add_argument("--anti", default=ANTI_SUBSTRATE)
    ap.add_argument("--wt", default=None,
                    help="receptor name treated as wild type (default: auto)")
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    path = DOCK_DIR / args.input
    if not path.exists():
        print(f"ERROR: {path} not found - run s04_dock.py first", file=sys.stderr)
        return 1

    data = load(path)
    if not data:
        print("ERROR: no usable rows in the docking CSV", file=sys.stderr)
        return 1

    sub, anti = slug(args.substrate), slug(args.anti)
    receptors = sorted({r for r, _ in data})

    # Identify wild type: the receptor whose name carries no mutation tag.
    wt = args.wt
    if wt is None:
        plain = [r for r in receptors if "_" not in r]
        wt = plain[0] if plain else receptors[0]
    if (wt, sub) not in data:
        print(f"ERROR: wild type {wt!r} has no {args.substrate} docking data",
              file=sys.stderr)
        return 1

    # --- per-receptor selectivity ----------------------------------------
    results = []
    for rec in receptors:
        if (rec, sub) not in data or (rec, anti) not in data:
            continue
        s_mean, s_sd, s_best = summarise(data[(rec, sub)])
        a_mean, a_sd, a_best = summarise(data[(rec, anti)])
        sel = s_mean - a_mean
        # Propagate replicate scatter into the selectivity uncertainty.
        sel_sd = (s_sd ** 2 + a_sd ** 2) ** 0.5
        results.append({
            "receptor": rec,
            "sub_mean": s_mean, "sub_sd": s_sd,
            "anti_mean": a_mean, "anti_sd": a_sd,
            "selectivity": sel, "selectivity_sd": sel_sd,
            "n": len(data[(rec, sub)]),
        })

    wt_row = next(r for r in results if r["receptor"] == wt)
    wt_sel = wt_row["selectivity"]

    for r in results:
        r["delta_selectivity"] = r["selectivity"] - wt_sel
        # A shift is only interesting if it exceeds the replicate noise.
        noise = max(0.1, (r["selectivity_sd"] ** 2 + wt_row["selectivity_sd"] ** 2) ** 0.5)
        r["signal_to_noise"] = r["delta_selectivity"] / noise

    # --- report -----------------------------------------------------------
    print(f"Substrate      : {args.substrate}")
    print(f"Anti-substrate : {args.anti}")
    print(f"Wild type      : {wt}\n")
    print(f"S = dG({args.substrate}) - dG({args.anti});  "
          "negative = prefers substrate\n")

    hdr = (f"  {'receptor':<22}{'dG_sub':>9}{'dG_anti':>9}"
           f"{'S':>8}{'dS vs WT':>10}{'S/N':>7}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    ordered = sorted(results, key=lambda r: -r["delta_selectivity"])
    for r in ordered[:args.top]:
        tag = "  <- WT" if r["receptor"] == wt else ""
        print(f"  {r['receptor']:<22}{r['sub_mean']:>9.2f}{r['anti_mean']:>9.2f}"
              f"{r['selectivity']:>8.2f}{r['delta_selectivity']:>+10.2f}"
              f"{r['signal_to_noise']:>+7.1f}{tag}")

    # --- ground-truth check ----------------------------------------------
    gt = [r for r in results if GROUND_TRUTH_MUTATION.upper() in r["receptor"].upper()]
    print()
    if gt:
        g = gt[0]
        verdict = ("RECOVERED" if g["delta_selectivity"] > 0 else "NOT recovered")
        print(f"Ground-truth check ({GROUND_TRUTH_MUTATION}): "
              f"dS = {g['delta_selectivity']:+.2f} kcal/mol -> {verdict}")
        print("  The published result is that this substitution shifts preference")
        print(f"  from {args.substrate} towards {args.anti}, i.e. dS > 0.")
        if g["delta_selectivity"] <= 0:
            print("  The pipeline does NOT reproduce it. Treat every other "
                  "prediction below as unvalidated.", file=sys.stderr)
    else:
        print(f"Ground-truth check: {GROUND_TRUTH_MUTATION} not present in the "
              "docking set - build it with s06_mutate.py and re-dock.")

    out = DOCK_DIR / "selectivity_ranking.csv"
    fields = ["receptor", "sub_mean", "sub_sd", "anti_mean", "anti_sd",
              "selectivity", "selectivity_sd", "delta_selectivity",
              "signal_to_noise", "n"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in ordered:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v)
                        for k, v in r.items()})
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
