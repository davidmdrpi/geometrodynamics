"""Two routes to the sheet pairing, on equal terms. Freeze ``fd78693``.

Corrected after review of ``76c2067``.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.bulk import klein_handle_electric as ke
from geometrodynamics.bulk import klein_handle_response as kh
from geometrodynamics.bulk import klein_handle_routes as kr
from experiments.closure_ledger import klein_handle_routes_probe as probe

SMALL = dict(n_ext=7, port_a=0, port_b=3, n_tube=3)
ARCHIVE = (Path(__file__).resolve().parents[1] / "experiments"
           / "closure_ledger" / "runs" / "20261111_klein_handle_routes"
           / "probe.json")
REVIEW_ARCHIVE = ARCHIVE.parents[1] / "20261111_klein_handle_routes_review" \
    / "probe.json"


@pytest.fixture(scope="module")
def model():
    return ke.Electric(kh.Network(**SMALL))


@pytest.fixture(scope="module")
def ds(model):
    return kr.draws(model, np.random.default_rng(kr.SEED), 6)


@pytest.fixture(scope="module")
def small_report():
    """The probe on a small network with a short long window."""
    return probe.run_probe(sizes={"small": SMALL}, long_steps=600,
                           numerics=True)


# ---- preparations ---------------------------------------------------------

@pytest.mark.parametrize("name", ["matter", "matter_at_rest", "flux_route",
                                  "both"])
def test_preparations_are_admissible_at_the_requested_energy(model, ds, name):
    deck = model.deck(-1, -1)
    for d in ds:
        z = getattr(kr, name)(model, d, 1e-3)
        assert np.abs(deck @ z - z).max() < 1e-12
        assert np.abs(model.gauss(z)).max() < 1e-12
        assert abs(model.charge_density(z).sum()) < 1e-12
        assert model.energy(z)[0] == pytest.approx(1e-3, rel=1e-10)


def test_routes_differ_in_their_gauge_invariant_content(model, ds):
    d = ds[0]
    m0 = kr.matter_at_rest(model, d, 1e-3)
    f = kr.flux_route(model, d, 1e-3)
    assert abs(model.handle_flux(m0)[0]) < 1e-15
    assert np.abs(model.link_current(m0)).max() > 1e-6
    assert model.handle_flux(f)[0] > 1e-4
    assert np.abs(model.link_current(f)).max() < 1e-15


def test_flux_energy_bookkeeping(model, ds):
    """The flux carries the requested fraction of the energy."""
    for f in kr.FLUX_FRACTIONS:
        phi = kr.flux_for(model, 1e-3, f)
        assert 0.5 * phi ** 2 * kr.loop_edges(model) == pytest.approx(f * 1e-3)
        z = kr.flux_route(model, ds[0], 1e-3, f)
        assert model.handle_flux(z)[0] == pytest.approx(phi)


def test_matter_draw_has_zero_charge_term_by_term(model, ds):
    for d in ds:
        assert abs(d.x @ d.py) < 1e-12 and abs(d.y @ d.px) < 1e-12


# ---- P1 - P9 on the small network ----------------------------------------

def test_flux_with_field_at_rest_gives_no_sheet_charge(model, ds):
    r = kr.rest_flux_report(model, ds)
    assert r["max_sheet_charge_at_rest"] < 1e-12
    assert r["min_sheet_charge_moving"] > 1e-10
    assert r["min_flux_change_at_rest"] > 1e-12


def test_energy_scaling_and_ranking(model, ds):
    sc = kr.scaling_report(model, ds)
    assert abs(sc["M"]["slope_of_median"] - 1) < 0.1
    assert abs(sc["M0"]["slope_of_median"] - 1) < 0.1
    assert abs(sc["F"]["slope_of_median"] - 1.5) < 0.15
    assert abs(sc["ratio"]["slope_of_median"] - 0.5) < 0.2


def test_optimal_flux_fraction(model, ds):
    r = kr.flux_fraction_report(model, ds)
    for got, want in zip(r["median_ratio"], r["predicted_ratio"]):
        assert got == pytest.approx(want, rel=0.1)
    assert all(f == pytest.approx(1 / 3) for f in r["best_fraction_by_draw"])


def test_coupling_dependence(model, ds):
    r = kr.coupling_report(model.net, ds)
    assert all(abs(s) < 0.1 for s in r["M0_slopes"])
    assert all(abs(s - 1) < 0.1 for s in r["F_slopes"])


def test_electric_sectors(model, ds):
    r = kr.sector_report(model, ds)
    for key in ("handle_flux", "loop_holonomy", "loop_current",
                "swap_even_field"):
        assert r["M0"][key] < 1e-12
    assert r["M0"]["sheet_charge"] > 1e-9
    # Not a property of the route in general: a moving matter current does
    # develop a little loop flux.
    assert r["M"]["handle_flux"] > 1e-9


def test_no_first_order_interference(model, ds):
    r = kr.interference_report(model, ds)
    assert all(abs(s - 2) < 0.3 for s in r["defect_slopes"])
    assert all(abs(s - 1) < 0.1 for s in r["reversal_slopes"])
    assert r["joint_reversal_error"] < 1e-12


def test_screening_sign_rule_and_leading_order(model, ds):
    r = kr.screening_report(model, ds)
    assert r["sign_agreement"] == r["resolved_draws"] > 0
    lo = kr.screening_convergence_report(model, ds)
    assert abs(lo["error_slope"] - 2) < 0.4
    assert lo["bond_rate_over_bond"] < 1e-12
    assert lo["median_error"][0] < lo["median_error"][-1]


def test_reversal_requirements(model, ds):
    r = kr.reversal_report(model, ds, np.random.default_rng(3))
    assert r["rest_data_fixed_by_alice_deck"] < 1e-12
    assert r["alice_fixed_set_invariance"] < 1e-9
    assert r["rest_sheet_charge"] > 1e-9
    assert r["general_matter_not_alice_fixed"] > 1e-6
    assert r["ordinary_sheet_charge"] < 1e-12
    assert r["ordinary_local_charge"] > 1e-7


def test_long_observables_match_stored_histories(model, ds):
    """The streaming loop agrees with the stored-history route."""
    z = kr.matter(model, ds[0], 1e-3)
    obs = kr.long_observables(model, z, 40)
    hist = ke.flow(model, z, 40)[:, 0]
    assert np.abs(obs["sheet"][:, 0] - model.sheet_charge(hist)).max() < 1e-15
    assert np.abs(obs["flux"][:, 0] - model.handle_flux(hist)).max() < 1e-15
    assert obs["overlap"][0, 0] == pytest.approx(1.0)


# ---- verdict --------------------------------------------------------------

def test_contradicted_prediction_is_not_a_numerical_failure():
    groups = {name: [name] for name in kr.CONFIRMED}
    checks = {name: True for name in groups}
    checks["screening_rule"] = False
    plain = kr.verdict(checks, groups)
    assert plain["screening_rule"] == kr.UNRESOLVED_NUMERICALLY
    assert plain["numerical_validity"] == kr.UNRESOLVED_NUMERICALLY
    registered = kr.verdict(checks, groups, {"screening_rule": "one draw"})
    assert registered["screening_rule"] == kr.COUNTEREXAMPLE
    assert registered["numerical_validity"] == "PASS"
    assert registered["contradicted_checks"] == ["screening_rule"]
    assert registered["failed_checks"] == []


def test_missing_group_is_incomplete():
    groups = {name: [name] for name in kr.CONFIRMED}
    del groups["interference"]
    v = kr.verdict({name: True for name in groups}, groups)
    assert v["interference"] == kr.INCOMPLETE
    assert v["numerical_validity"] == kr.INCOMPLETE


# ---- probe ----------------------------------------------------------------

def test_small_probe_is_numerically_valid(small_report):
    v = small_report["verdict"]
    assert v["numerical_validity"] == "PASS"
    assert not v["failed_checks"] and not v["missing_checks"]
    for name, value in kr.FIXED.items():
        assert v[name] == value
    assert set(kr.VERDICT_FIELDS) <= set(v)


def test_small_probe_confirms_the_symmetry_predictions(small_report):
    v = small_report["verdict"]
    for name in ("flux_route_requirement", "energy_scaling", "flux_fraction",
                 "coupling_dependence", "electric_sectors", "interference",
                 "reversal_requirements", "screening_leading_order",
                 "seam_cells", "temporal_symmetry"):
        assert v[name] == kr.CONFIRMED[name], name


def test_archived_probe_records_the_contradicted_prediction():
    """The full run: ten predictions, one of them contradicted as frozen."""
    report = json.loads(ARCHIVE.read_text())
    v = report["verdict"]
    assert v["numerical_validity"] == "PASS"
    assert v["screening_rule"] == kr.COUNTEREXAMPLE
    assert "screening_rule" in report["contradicted_predictions"]
    assert len(v["contradicted_checks"]) == 1
    # The first archive is kept as written: it predates the review, so two
    # labels are the frozen ones and three fields do not exist yet.
    added = {"seam_cells", "temporal_symmetry", "flux_retention"}
    for name, value in kr.CONFIRMED.items():
        if name == "screening_rule" or name in added:
            continue
        assert v[name] == kr.FROZEN_LABELS.get(name, value), name
    assert set(report["measurements"]) == {"primary", "replicate"}
    assert report["long_window"] == pytest.approx(200.0)


def test_archived_probe_never_ranks_the_flux_route_first_at_short_times():
    report = json.loads(ARCHIVE.read_text())
    for size in ("primary", "replicate"):
        ratio = report["measurements"][size]["P2_P3"]["ratio"]
        assert all(r < 1 for r in ratio["median"])
        assert ratio["draws_above_one"] == [0, 0, 0, 0]


# ---- Review of 76c2067 ----------------------------------------------------

@pytest.mark.parametrize("cell", ke.CELLS)
def test_cell_matter_state_is_admissible_with_no_sheet_charge(model, cell):
    z = kr.cell_matter_state(model, *cell, 1e-3, np.random.default_rng(5))
    assert np.abs(model.deck(*cell) @ z - z).max() < 1e-12
    assert np.abs(model.gauss(z)).max() < 1e-12
    assert abs(model.sheet_charge(z)[0]) < 1e-12
    assert abs(model.handle_flux(z)[0]) < 1e-12
    assert abs(model.charge_density(z).sum()) < 1e-12
    assert model.energy(z)[0] == pytest.approx(1e-3, rel=1e-10)


def test_matter_pairing_exists_in_every_cell_but_the_ordinary_handle(model):
    """R1. Fibre reversal is sufficient, not necessary."""
    rows = {(r["s_t"], r["s_f"]): r for r in kr.seam_cell_matter_report(
        model, np.random.default_rng(6), 6)}
    assert rows[(1, 1)]["max_peak_sheet_charge"] < 1e-12
    for cell in ((1, -1), (-1, 1), (-1, -1)):
        assert rows[cell]["min_peak_sheet_charge"] > 1e-7, cell
        assert rows[cell]["dop853_disagreement"] < 1e-2


def test_time_parity_of_sheet_charge_is_set_by_the_seam(model):
    """R1. ``Q_1(s_t t) = -s_t s_f Q_1(t)``."""
    rows = kr.seam_cell_matter_report(model, np.random.default_rng(7), 6)
    assert all(r["time_parity_error"] < 1e-9 for r in rows)
    parity = {(r["s_t"], r["s_f"]): r["predicted_parity"] for r in rows}
    assert parity[(-1, 1)] == "even" and parity[(-1, -1)] == "odd"


def test_even_parity_is_not_odd_parity(model):
    """The time-reversing fibre-preserving cell is even, and visibly not odd."""
    z = kr.cell_matter_state(model, -1, 1, 1e-3, np.random.default_rng(8))
    fwd = model.sheet_charge(ke.flow(model, z, 100)[:, 0])
    bwd = model.sheet_charge(ke.flow(model, z, 100, -kr.STEP)[:, 0])
    assert np.abs(fwd - bwd).max() < 1e-12
    assert np.abs(fwd + bwd).max() > 1e-7


def test_additivity_defect_has_an_odd_part(model, ds):
    """R2. Additive only at leading order in the matter amplitude."""
    r = kr.odd_defect_report(model, ds[:4])
    for row in r["rows"]:
        assert row["min_peak_odd_defect"] > 1e-12
        assert row["worst_change_on_halving_flux"] < 1e-2
    assert abs(r["slope_in_energy"] - 2) < 0.2


def test_flux_retention_is_reported_per_draw(model, ds):
    """R3. Extrema over draws, not medians."""
    r = kr.flux_retention_report(model, ds, steps=400)
    for row in r["rows"]:
        assert len(row["per_draw_min"]) == len(ds)
        assert row["worst_min"] <= row["median_min"]
        assert row["worst_max"] >= row["median_max"]
        assert row["draws_changing_sign"] == 0


def test_review_archive_records_the_corrections():
    report = json.loads(REVIEW_ARCHIVE.read_text())
    v = report["verdict"]
    assert v["numerical_validity"] == "PASS"
    assert v["screening_rule"] == kr.COUNTEREXAMPLE
    for name, value in kr.CONFIRMED.items():
        if name != "screening_rule":
            assert v[name] == value, name
    assert len(report["withdrawn_claims"]) == 3
    assert report["frozen_labels_changed"] == kr.FROZEN_LABELS
    worst = max(r["largest_excursion"] for size in ("primary", "replicate")
                for r in report["measurements"][size]["review"]
                ["flux_retention"]["rows"])
    assert 0.015 < worst < 0.07
