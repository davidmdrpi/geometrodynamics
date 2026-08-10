"""
Tests for the spin-2 field: the things that make it *not* a scalar.

* the solver reproduces three **exact** modes — ``q = 1`` is ``ℓ = 2``,
  ``q = cos d`` is ``ℓ = 3``, ``q = 7cos²d − 1`` is ``ℓ = 4`` — at
  ``ω = √(ℓ(ℓ+1))``, and conserves its invariant;
* the field **vanishes at both poles** and has no ``ℓ < 2`` content, so its
  smallest source is a ring and its focus is a ring;
* the tensor is symmetric, trace-free, and the deformation it produces is
  **area-preserving** — pure shear, not breathing;
* the strain pattern has **spin weight 2**: it repeats every 180° and inverts
  every 90°;
* and at the caustic it turns by a quarter, not a half — the full inversion
  takes the round trip.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.spin2_tidal import (
    ANTIPODAL_TIME,
    RETURN_TIME,
    Spin2WaveSim,
    TidalField,
    measure_area_preservation,
    measure_caustic_phase,
    measure_exact_modes,
    measure_node_at_the_focus,
    measure_round_trip_inversion,
    measure_scalar_contrast,
    measure_spin_weight,
    mode_frequency,
)

N = 600          # coarse but convergent; the exact-mode tests use their own


@pytest.fixture(scope="module")
def field():
    return TidalField(sim=Spin2WaveSim(n=900))


# ── the solver ──────────────────────────────────────────────────────────────
def test_no_modes_below_quadrupole():
    """The defining fact: a spin-2 field has no monopole and no dipole."""
    for ell in (0, 1):
        with pytest.raises(ValueError):
            mode_frequency(ell)
    assert mode_frequency(2) == pytest.approx(math.sqrt(6.0))
    assert mode_frequency(3) == pytest.approx(math.sqrt(12.0))


def test_solver_validates_its_inputs():
    with pytest.raises(ValueError):
        Spin2WaveSim(n=8)
    with pytest.raises(ValueError):
        Spin2WaveSim(pulse_width=0.0)
    with pytest.raises(ValueError):
        Spin2WaveSim(pulse_width=4.0)


def test_exact_modes_keep_their_shape_and_frequency():
    r = measure_exact_modes(n=600, periods=6.0)
    assert r["worst_shape_error"] < 2e-3
    by_ell = {m["ell"]: m for m in r["modes"]}
    assert by_ell[2]["shape_error"] < 1e-6      # q = 1 is exact on any grid
    for m in r["modes"]:
        assert m["omega"] == pytest.approx(math.sqrt(m["ell"] * (m["ell"] + 1)))
        assert m["energy_drift"] < 1e-9


def test_the_quadrupole_oscillates_at_root_six():
    """The one mode whose whole time dependence is checkable by hand."""
    sim = Spin2WaveSim(n=600)
    sim.start_from(np.ones_like(sim.d))
    w = math.sqrt(6.0)
    for frac in (0.25, 0.5, 1.0):
        sim.advance_to(frac * 2.0 * math.pi / w)
        expected = math.cos(w * sim.t)
        assert float(np.max(sim.q)) == pytest.approx(expected, abs=2e-4)


def test_energy_is_conserved_through_the_focus(field):
    field.reset()
    field.advance_to(1.05 * ANTIPODAL_TIME)
    assert field.sim.energy_drift() < 1e-9


def test_start_from_validates_shape():
    sim = Spin2WaveSim(n=200)
    with pytest.raises(ValueError):
        sim.start_from(np.ones(7))


# ── a spin-2 field cannot live at a pole ────────────────────────────────────
def test_the_field_vanishes_at_both_poles(field):
    """``h = sin²d·q`` — no ``q`` can put amplitude where the frame degenerates."""
    field.reset()
    field.advance_to(1.0)
    h = field.sim.h
    assert abs(h[0]) < 1e-3 * float(np.max(np.abs(h)))
    assert abs(h[-1]) < 1e-3 * float(np.max(np.abs(h)))


def test_the_focus_is_a_ring_not_a_peak(field):
    r = measure_node_at_the_focus(field, frames=90)
    assert r["is_a_ring_not_a_peak"] is True
    assert r["ring_radius_at_focus"] > 0.05
    assert r["amplitude_at_the_antipode"] < 0.02 * r["peak_amplitude"]


def test_the_scalar_does_the_opposite(field):
    """The contrast that motivates the whole module."""
    r = measure_scalar_contrast(field, frames=110)
    assert r["scalar_sits_on_the_antipode"] is True
    assert r["tensor_is_a_ring"] is True
    assert r["tensor_peak_distance"] < r["scalar_peak_distance"]


# ── it is a tensor, not a scalar in disguise ────────────────────────────────
def test_the_tensor_is_symmetric_and_trace_free(field):
    field.reset()
    field.advance_to(0.8)
    m = field.matrix(1.0).reshape(2, 2)
    assert float(np.trace(m)) == pytest.approx(0.0, abs=1e-15)
    assert m[0, 1] == pytest.approx(m[1, 0], abs=1e-15)


def test_the_deformation_preserves_area(field):
    """Trace-free means pure shear: no first-order area change at all."""
    r = measure_area_preservation(field, eps=0.05)
    assert r["trace"] == pytest.approx(0.0, abs=1e-15)
    assert r["area_preserved_to_first_order"] is True
    assert r["first_order_area_change"] == pytest.approx(
        r["second_order_prediction"], rel=1e-3)


def test_strain_has_spin_weight_two(field):
    r = measure_spin_weight(field)
    assert r["consistent"] is True
    assert r["repeats_after_180_deg"] < 1e-12     # identical after a half turn
    assert r["differs_after_90_deg"] > 1.5        # inverted after a quarter


def test_only_one_polarisation_for_an_axisymmetric_source():
    even = TidalField(sim=Spin2WaveSim(n=400), parity="even")
    odd = TidalField(sim=Spin2WaveSim(n=400), parity="odd")
    for f in (even, odd):
        f.reset()
        f.advance_to(0.7)
    hp_e, hx_e = even.components(1.0)
    hp_o, hx_o = odd.components(1.0)
    assert float(hx_e) == 0.0 and float(hp_e) != 0.0
    assert float(hp_o) == 0.0 and float(hx_o) != 0.0
    # the odd field is the even one rotated by 45°
    assert float(even.strain(1.0, 0.0)) == pytest.approx(
        float(odd.strain(1.0, 0.25 * math.pi)), rel=1e-12)


def test_parity_is_validated():
    with pytest.raises(ValueError):
        TidalField(parity="sideways")


# ── the ellipses really lie on the sphere ───────────────────────────────────
def test_the_frame_is_orthonormal_and_tangent(field):
    p, e_d, e_psi = field.frame(1.1, 0.8)
    assert np.linalg.norm(p) == pytest.approx(1.0)
    assert np.linalg.norm(e_d) == pytest.approx(1.0)
    assert np.linalg.norm(e_psi) == pytest.approx(1.0)
    assert float(p @ e_d) == pytest.approx(0.0, abs=1e-14)
    assert float(p @ e_psi) == pytest.approx(0.0, abs=1e-14)
    assert float(e_d @ e_psi) == pytest.approx(0.0, abs=1e-14)


def test_geodesic_distance_matches_the_frame(field):
    d, psi = 1.3, 2.0
    p = field.point(d, psi)
    cos_d = float(p @ field.source_direction)
    assert math.acos(max(-1.0, min(1.0, cos_d))) == pytest.approx(d)


def test_ellipse_closes_and_is_centred_on_its_point(field):
    field.reset()
    field.advance_to(0.9)
    ring = field.ellipse(1.0, 0.4, eps=10.0, n=64, size=0.1)
    assert np.allclose(ring[0], ring[-1])
    p, _, _ = field.frame(1.0, 0.4)
    centre = ring.mean(axis=0)
    assert np.linalg.norm(centre - p) < 0.02


# ── the caustic ─────────────────────────────────────────────────────────────
def test_one_caustic_passage_is_a_quarter_turn_not_a_flip(field):
    """The correction to the obvious guess, and it is measured.

    A single focal passage shifts the phase by 90° — the outbound waveform is
    the Hilbert transform of the inbound one, not its negative.
    """
    r = measure_caustic_phase(field, frames=1200)
    assert r["is_a_quarter_turn_not_a_flip"] is True
    assert r["best_match"].startswith("phase_")
    assert r["best_correlation"] > 0.7
    assert abs(r["correlations"]["inverted"]) < 0.5


def test_the_round_trip_does_invert_it(field):
    """Two passages — the antipode, then home — and the axes have swapped."""
    r = measure_round_trip_inversion(field)
    assert r["inverts"] is True
    assert r["corr_after_one_round_trip_with_minus_start"] > 0.95
    assert r["corr_after_two_round_trips_with_start"] > 0.95
    assert r["inversion_residual"] < 0.05
