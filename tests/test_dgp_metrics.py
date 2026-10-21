from itertools import product
import numpy as np
import pytest
from scipy.stats import binom, norm

from src.dgps import family_rng, generate_family
from src.metrics import alternative_power, analytic_targets, mean_mcse


def test_dgp_count_reproducibility_and_stream_identity():
    a = generate_family(1000, .3, 2, family_rng(20261006, 1000, .3, 2, 0))
    b = generate_family(1000, .3, 2, family_rng(20261006, 1000, .3, 2, 0))
    for left, right in zip(a, b):
        np.testing.assert_array_equal(left, right)
    assert a[1].sum() == 300
    other = generate_family(1000, .3, 2, family_rng(20261006, 1000, .3, 2, 1))
    assert not np.array_equal(a[0], other[0])


def test_dgp_distribution_properties():
    p, null = generate_family(100_000, .5, 2, np.random.default_rng(20))
    z = norm.isf(p)
    # Loose, prespecified tolerances (more than six standard errors).
    assert abs(p[null].mean() - .5) < .01
    assert abs(z[~null].mean() - 2) < .03
    assert abs(np.var(z[null]) - 1) < .04


def test_fixed_fdr_against_binomial_enumeration():
    m, pi0, mu, gamma = 8, .375, 1.5, .1
    g = alternative_power(gamma, mu)
    exact = sum(binom.pmf(v, 3, gamma) * binom.pmf(s, 5, g) * v / max(v + s, 1)
                for v, s in product(range(4), range(6)))
    result = analytic_targets(m, pi0, mu, gamma)
    assert result["true_fdr"] == pytest.approx(exact, abs=1e-12)
    assert result["quadrature_error"] < 1e-9


def test_all_null_and_all_alternative_targets():
    all_null = analytic_targets(10, 1, 2, .01)
    assert all_null["true_fdr"] == pytest.approx(1 - .99**10)
    assert np.isnan(all_null["true_power"])
    assert analytic_targets(10, 0, 2, .01)["true_fdr"] == 0


def test_bernoulli_fdr_against_multinomial_enumeration():
    # Independently enumerate every outcome for m=3:
    # false rejection, true rejection, no rejection.
    m, pi0, mu, gamma = 3, .4, 1, .1
    g = alternative_power(gamma, mu)
    probabilities = [pi0 * gamma, (1 - pi0) * g, 1 - pi0 * gamma - (1 - pi0) * g]
    exact = 0
    for outcomes in product(range(3), repeat=m):
        prob = np.prod([probabilities[outcome] for outcome in outcomes])
        v = outcomes.count(0)
        r = v + outcomes.count(1)
        exact += prob * v / max(r, 1)
    assert analytic_targets(m, pi0, mu, gamma, "bernoulli")["true_fdr"] == pytest.approx(exact)


def test_mcse_and_missing_values():
    mean, se, n = mean_mcse([1, 2, 3, np.nan])
    assert mean == 2 and n == 3
    assert se == pytest.approx(1 / np.sqrt(3))
