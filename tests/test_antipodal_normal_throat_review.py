"""Four-cell identities and failure controls for the public review amendment."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.bulk import antipodal_normal_throat as ant
from experiments.closure_ledger import antipodal_normal_throat_review_probe as probe


@pytest.mark.parametrize("kind", ant.CELLS)
@pytest.mark.parametrize("eta", [-1, 1])
def test_face_parity_and_actual_transfer_derivative(kind, eta):
    model = ant.Throat(16, 1)
    rng = np.random.default_rng(23)
    states = np.array([probe.generic_start(model, kind, eta, rng) for _ in range(4)])
    c, r = ant.CELLS[kind]
    face = model.angular_momentum(states, per_face=True)
    if r == -1:
        assert np.max(abs(face[:, 1]-c*face[:, 0])) < 1e-14
        assert np.max(abs(states[:, -4]-states[:, -3])) < 1e-14
        assert np.max(abs(states[:, -2]-states[:, -1])) < 1e-14
    elif c == -1:
        assert np.max(abs(face)) < 1e-14
    if c == -1:
        assert np.max(abs(face.sum(axis=1))) < 1e-14
    else:
        assert np.max(abs(face.sum(axis=1))) > 1e-4
    eps = 1e-6
    dz = model.rhs(states)
    derivative = (model.angular_momentum(states+eps*dz, per_face=True)
                  - model.angular_momentum(states-eps*dz, per_face=True))/(2*eps)
    flux = model.normal_angular_flux(states)
    assert np.max(abs(derivative-np.column_stack([flux, -flux]))) < 1e-10
    if c*r == -1:
        assert np.max(abs(flux)) < 1e-14
    else:
        assert np.max(abs(flux)) > 1e-4
    hist = ant.flow(model, states, duration=.2, h=.005, sample=20)
    assert np.max(abs(hist[-1, :, -4:-2]-states[:, -4:-2])) > 1e-4
    assert np.max(abs(model.deck(hist[-1], kind, eta)-hist[-1])) < 1e-12


def test_geometry_probe_observes_coordinate_action(monkeypatch):
    rows = probe.geometry_report()
    assert len(rows) == 4
    assert max(r["jacobian_error"] for r in rows) < 1e-8
    assert max(r["pullback_error"] for r in rows) < 1e-10
    identify = ant.identify

    def broken(points, kind="antipodal"):
        out = identify(points, kind)
        out[..., 2] *= 2
        return out

    monkeypatch.setattr(ant, "identify", broken)
    assert min(r["jacobian_error"] for r in probe.geometry_report()) > .9


def test_one_face_source_measures_both_normal_channels():
    r = probe.one_face_response_report()
    assert r["one_face_sources"] == 32
    assert r["max_two_face_response_error"] < 1e-10


def test_phi_symbol_report_covers_both_harmonics():
    rows = probe.phi_convergence_report()
    assert {(r["m"], r["n"]) for r in rows} == {(m, n) for m in (1, 2, 3) for n in (1, 2)}
    assert all(1.8 < r["order"] < 2.2 for r in rows)


@pytest.fixture
def archived_rows():
    path = Path(__file__).resolve().parents[1]/"experiments/closure_ledger/runs/20261112_antipodal_normal_throat_review/probe.json"
    return json.loads(path.read_text())["generic_runs"]


def test_generic_gates_reject_actual_broken_and_frozen_histories(archived_rows):
    assert all(probe.generic_checks(archived_rows).values())
    negatives = probe.negative_controls(archived_rows)
    assert negatives["broken_seam_rejected"]
    assert negatives["frozen_transfer_rejected"]
    assert negatives["broken_history_metrics"]["seam_error"] > .01
    assert negatives["frozen_history_metrics"]["face_change"] == 0


def test_generic_coverage_rejects_missing_cells_duplicate_draws_and_short_runs(archived_rows):
    assert not probe.generic_checks(archived_rows[:-1])["generic_coverage"]
    rows = copy.deepcopy(archived_rows)
    rows[0]["draws"][0] = rows[0]["draws"][1]
    assert not probe.generic_checks(rows)["generic_coverage"]
    rows = copy.deepcopy(archived_rows)
    rows[0]["duration"] = 4.
    assert not probe.generic_checks(rows)["generic_coverage"]


def test_amended_required_checks_are_not_derived_from_report():
    checks = dict.fromkeys(probe.REQUIRED, True)
    assert ant.classify(checks, probe.REQUIRED)["numerical_validity"] == "PASS"
    del checks["four_cell_transfer"]
    assert ant.classify(checks, probe.REQUIRED)["numerical_validity"] == "INCOMPLETE"
    checks["four_cell_transfer"] = False
    assert ant.classify(checks, probe.REQUIRED)["numerical_validity"] == "UNRESOLVED_NUMERICALLY"
