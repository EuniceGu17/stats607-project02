"""Run from the project root: python -m examples.one_replicate"""

from pprint import pprint

from src.dgps import family_rng, generate_family
from src.methods import estimate_fdr, reject_bh, reject_fixed
from src.metrics import evaluate_rejections


def main():
    # 1. Generate a family of 1000 independent tests: 800 null, 200 alternative.
    m, pi0, mu, gamma, lam = 1000, 0.8, 2.0, 0.01, 0.5
    rng = family_rng(20261006, m, pi0, mu, replicate=0)
    pvalues, is_null = generate_family(m, pi0, mu, rng)

    # 2. The estimator and both rules see only the p-values.
    estimates = estimate_fdr(pvalues, gamma, lam)
    fixed = reject_fixed(pvalues, gamma)
    bh, threshold = reject_bh(pvalues, estimates["fdr_hat"])

    # 3. Truth labels are used only now, to evaluate the decisions.
    pprint({"first_five_pvalues": pvalues[:5].tolist(), "estimates": estimates,
            "fixed_rule": evaluate_rejections(fixed, is_null),
            "bh_rule": evaluate_rejections(bh, is_null),
            "bh_threshold": threshold})
    print("These are single-replicate FDP/power values, not Monte Carlo FDR estimates.")


if __name__ == "__main__":
    main()
