import math
import pytest

from experiments.closure_ledger.historical_foundation_audit_probe import (
    threshold_audit, finite_band_focus_audit, readout_audit,
)


def test_legacy_threshold_checks_miss_rest_energy_mismatch():
    d=threshold_audit()
    assert d['legacy_minimum_test_pass'] and d['legacy_threshold_test_pass']
    assert d['minimum_to_electron_rest_energy']==pytest.approx(d['analytic_ratio'], rel=1e-13)
    assert d['computed_pair_minimum_MeV'] < .01*d['asserted_pair_threshold_MeV']


def test_finite_band_pole_limit_and_exact_refocus():
    d=finite_band_focus_audit()
    a,b,c=d['near_pole_errors']
    assert b < .26*a and c < .26*b
    assert d['antipodal_identity_error'] < 1e-10


def test_readout_countermodel_preserves_marginals_not_tsirelson():
    a,b=readout_audit()['records']
    assert a['chsh']==pytest.approx(2*math.sqrt(2), abs=1e-12)
    assert b['chsh']==pytest.approx(8*math.sqrt(2)/3, abs=1e-12)
    assert max(a['max_marginal_error'],b['max_marginal_error']) < 1e-12
