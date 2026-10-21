import numpy as np
import pytest
from scipy.stats import false_discovery_control

from src.methods import estimate_fdr, reject_bh, reject_fixed
from src.metrics import evaluate_rejections


def test_bh_hand_calculation_and_original_order():
    p = np.array([0.2, 0.03, 0.001, 0.01])
    mask, threshold = reject_bh(p, 0.05)
    np.testing.assert_array_equal(mask, [False, True, True, True])
    assert threshold == 0.03


@pytest.mark.parametrize("alpha", [0, .001, .05, .2, 1])
def test_bh_against_independent_scipy_implementation(alpha):
    p = np.random.default_rng(100).uniform(size=80)
    p[:3] = [0.0001, 0.0002, 0.001]
    actual, _ = reject_bh(p, alpha)
    np.testing.assert_array_equal(actual, false_discovery_control(p, method="bh") <= alpha)


def test_bh_ties_and_no_rejections():
    mask, _ = reject_bh([.01, .01, .9], .05)
    np.testing.assert_array_equal(mask, [True, True, False])
    mask, threshold = reject_bh([.6, .8], .05)
    assert not mask.any() and threshold == 0


def test_estimator_by_hand_and_zero_rejections():
    estimate = estimate_fdr([.001, .02, .7, .8], .01, .5)
    assert estimate["pi0_hat_raw"] == 1
    assert estimate["fdr_hat_raw"] == pytest.approx(.04)
    estimate = estimate_fdr([.7, .8], .1, .5)
    assert estimate["pi0_hat_raw"] == 2
    assert estimate["fdr_hat_raw"] == pytest.approx(.4)
    assert estimate_fdr([.7, .8], .6, .5)["fdr_hat"] == 1
    assert estimate_fdr([.7, .8], .1, .5, cap_pi0=True)["pi0_hat_used"] == 1


def test_metrics_zero_rejections_and_no_alternatives():
    result = evaluate_rejections([False, False], [True, False])
    assert result["FDP"] == result["power"] == 0
    result = evaluate_rejections([True, False], [True, True])
    assert result["FDP"] == 1 and np.isnan(result["power"])


@pytest.mark.parametrize("p", [[], [-.1], [1.1], [np.nan], [[.1, .2]]])
def test_invalid_pvalues_rejected(p):
    with pytest.raises(ValueError):
        reject_bh(p, .05)
