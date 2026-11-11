"""Temporal symmetry of exchange across the seam. Freeze ``95d7a80``."""

import numpy as np
import pytest

from geometrodynamics.bulk import klein_handle_matching as km
from geometrodynamics.bulk import klein_handle_parity as kp
from geometrodynamics.bulk import klein_handle_response as kh
from experiments.closure_ledger import klein_handle_parity_probe as probe

SMALL = dict(n_ext=7, port_a=0, port_b=3, n_tube=3)


@pytest.fixture(scope="module")
def rows():
    return kp.parity_report(kh.Network(**SMALL),
                            np.random.default_rng(kp.SEED), 4)


def test_sheet_values_sum_to_the_noether_form():
    """Sheet 1 plus sheet 2 is the conserved total, for every generator."""
    net = kh.Network(**SMALL)
    model = km.Multiplet(net, components=2, loop_sign=-1)
    z = 0.3 * np.random.default_rng(1).standard_normal(model.dim)
    swapped = model._matrix(lambda psi, mom, q, pq: (
        psi @ net.swap.T, mom @ net.swap.T, q, pq)) @ z
    for g in kp.GENERATORS:
        total = model.noether(g, z[None])[0]
        both = kp.sheet_noether(model, g, z)[0] \
            + kp.sheet_noether(model, g, swapped)[0]
        assert both == pytest.approx(total, abs=1e-12)


def test_every_predicted_relation_holds(rows):
    assert len(rows) == 16
    for r in rows:
        ok = kp.satisfies(r)
        assert ok["rule"] and ok["live"], (r["internal"], r["s_t"],
                                           r["generator"])
        assert r["symmetry_defect"] < 1e-9 and r["max_drift"] < 1e-9


def test_charge_and_spin_have_opposite_parity_on_the_transport_seam(rows):
    pick = {(r["internal"], r["s_t"], r["generator"]): r for r in rows}
    charge = pick[("J", -1, "i")]
    assert charge["max_even_part"] < 1e-9 < charge["min_odd_part"]
    for g in kp.GENERATORS[1:]:
        spin = pick[("J", -1, g)]
        assert spin["max_odd_part"] < 1e-9 < spin["min_even_part"]


def test_time_preserving_transport_seam_freezes_spin_and_frees_charge(rows):
    pick = {(r["internal"], r["s_t"], r["generator"]): r for r in rows}
    assert pick[("J", 1, "sigma_z")]["max_variation"] < 1e-9
    free = pick[("J", 1, "i")]
    assert free["min_even_part"] > 1e-7 and free["min_odd_part"] > 1e-7


def test_classifier_can_fail():
    """A sheet history of the wrong parity is rejected."""
    t = np.linspace(0, 1, 50)
    odd = kp.classify(np.sin(t), -np.sin(t), 0.0, 1.0)
    row = {"predicted": "even", "max_odd_part": odd["odd_part"],
           "min_even_part": odd["even_part"], "max_abs_total": 0.0}
    assert not kp.satisfies(row)["rule"]


def test_verdict_and_probe():
    report = probe.run_probe(sizes={"small": SMALL}, draws=3)
    v = report["verdict"]
    assert v["numerical_validity"] == "PASS"
    for name, value in {**kp.CONFIRMED, **kp.FIXED}.items():
        assert v[name] == value
    bad = kp.verdict({"a": False}, {"parity_rule": ["a"]})
    assert bad["parity_rule"] == kh.UNRESOLVED_NUMERICALLY
    assert bad["non_vacuous"] == kh.INCOMPLETE
