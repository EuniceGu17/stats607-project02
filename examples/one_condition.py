"""Learn the loop for one condition: python -m examples.one_condition."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.dgps import family_rng, generate_family
from src.methods import estimate_fdr, reject_bh, reject_fixed
from src.metrics import evaluate_rejections


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("results/one-condition"))
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("Use at least two replicates to calculate a Monte Carlo SE")

    m, pi0, mu, gamma, lam = 1000, 0.8, 2.0, 0.01, 0.5
    rows = []
    for replicate in range(args.replicates):
        rng = family_rng(20261006, m, pi0, mu, replicate)
        pvalues, is_null = generate_family(m, pi0, mu, rng)
        estimates = estimate_fdr(pvalues, gamma, lam)
        fixed = reject_fixed(pvalues, gamma)
        bh, threshold = reject_bh(pvalues, estimates["fdr_hat"])
        for method, rejected in [("fixed", fixed), ("bh", bh)]:
            rows.append({
                "replicate_id": replicate, "method": method,
                "m": m, "pi0": pi0, "mu": mu, "gamma": gamma,
                "lambda": lam, "alpha_rep": estimates["fdr_hat"],
                "bh_threshold": threshold, **estimates,
                **evaluate_rejections(rejected, is_null),
            })
        if (replicate + 1) % 200 == 0:
            print(f"Completed {replicate + 1}/{args.replicates} families")

    raw = pd.DataFrame(rows)
    raw_dir = args.output / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw.to_csv(raw_dir / "replicates.csv", index=False)

    # Raw output is saved before any aggregation. Each family has two method rows.
    summaries = []
    for method, group in raw.groupby("method"):
        row = {"method": method, "replicates": len(group)}
        for source, label in [("FDP", "fdr"), ("power", "power")]:
            row[label] = group[source].mean()
            row[label + "_mcse"] = group[source].std(ddof=1) / np.sqrt(len(group))
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    summary_dir = args.output / "summaries"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_dir / "summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Saved {len(raw)} raw rows under {args.output}")


if __name__ == "__main__":
    main()
