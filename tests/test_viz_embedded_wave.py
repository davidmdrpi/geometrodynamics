"""
Tests for the continuous ℝ³ realisation of the spin-2 field.

The load-bearing claim is a theorem about the drawn surface, so it is tested
as one: the induced metric perturbation of the embedding **is** the solved
``h_ab`` — trace-free, off-diagonal free, and equal component by component.

Around it:

* a purely radial deformation is shown to be *conformal*, which is why the
  tangential slide is needed at all;
* ``ℓ = 2`` comes out as ``P₂(cos d)`` — the textbook prolate–oblate shape,
  derived;
* the free constant is exactly a rigid translation, and removing it leaves no
  dipole and no monopole;
* and the area is unchanged at first order, because trace-free means that.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geometrodynamics.viz.embedded_wave import (
    EmbeddedTidalSurface,
    measure_area_and_multipoles,
    measure_bulk_reach,
    measure_induced_metric,
    measure_quadrupole_shape,
)
from geometrodynamics.viz.spin2_tidal import Spin2WaveSim


@pytest.fixture(scope="module")
def surface():
    return EmbeddedTidalSurface(sim=Spin2WaveSim(n=900), n_theta=61, n_phi=91)


# ── the obstruction that motivates the construction ─────────────────────────
def test_a_radial_deformation_alone_is_conformal():
    """Why a height field cannot represent a trace-free tensor.

    ``g_ab = r²ĝ_ab + ∂_a r ∂_b r`` — the gradient term is second order, so at
    first order the perturbation is pure trace and the trace-free part is zero.
    """
    eps, d0, a0, step = 1e-4, 1.1, 0.4, 1e-5

    def X(d, a):                       # a radial deformation and nothing else
        r = 1.0 + eps * math.cos(2.0 * d)
        return r * np.array([math.sin(d) * math.cos(a),
                             math.sin(d) * math.sin(a), math.cos(d)])

    Xd = (X(d0 + step, a0) - X(d0 - step, a0)) / (2.0 * step)
    Xa = (X(d0, a0 + step) - X(d0, a0 - step)) / (2.0 * step)
    dg_dd = (float(Xd @ Xd) - 1.0) / eps
    dg_aa = (float(Xa @ Xa) - math.sin(d0) ** 2) / (eps * math.sin(d0) ** 2)
    assert dg_dd == pytest.approx(dg_aa, abs=1e-3)          # equal ⇒ conformal
    assert 0.5 * (dg_dd - dg_aa) == pytest.approx(0.0, abs=1e-3)


# ── the theorem ─────────────────────────────────────────────────────────────
def test_the_drawn_surface_has_the_solved_metric(surface):
    """The whole point: ``δg_ab`` of the embedding equals ``h_ab``."""
    r = measure_induced_metric(surface, t=1.2, gain=1e-4)
    assert r["reproduces_h"] is True
    assert r["worst_relative_error"] < 5e-3
    assert r["worst_relative_trace"] < 1e-2
    assert r["worst_off_diagonal"] < 1e-6


def test_component_by_component(surface):
    surface.reset()
    surface.advance_to(1.2)
    for d0 in (0.7, 1.4, 2.1):
        m = surface.induced_metric_perturbation(d0, gain=1e-4)
        peak = float(np.max(np.abs(surface.profiles()["shear"])))
        assert abs(m["h_plus_measured"] - m["h_plus_solved"]) < 5e-3 * peak
        assert abs(m["trace"]) < 1e-2 * peak


# ── what falls out ──────────────────────────────────────────────────────────
def test_the_quadrupole_is_the_legendre_shape():
    """``ℓ = 2`` ⇒ ``ρ ∝ P₂(cos d)``: the prolate–oblate picture, derived."""
    r = measure_quadrupole_shape(n=1200)
    assert r["is_legendre_p2"] is True
    assert r["shape_error"] < 1e-3
    assert r["amplitude_ratio"] == pytest.approx(1.0, abs=1e-3)
    assert abs(r["residual_dipole"]) < 1e-12


def test_the_free_constant_is_a_rigid_translation():
    """The kernel of the construction is ``ℓ = 1`` — i.e. moving the sphere.

    Adding a dipole to the potential displaces every point by one constant
    vector, which changes no geometry at all; that is why it is removed.
    """
    s = EmbeddedTidalSurface(sim=Spin2WaveSim(n=600), n_theta=41, n_phi=61)
    s.advance_to(1.0)
    p = s.profiles()
    d = np.linspace(0.2, math.pi - 0.2, 25)
    a = np.zeros_like(d)
    base = s.positions(d[None, :], a[None, :], gain=1.0)[0]
    # a dipole in W means slide += C sin d, radial += -C cos d
    c = 0.37
    n_hat, e_d = s._frames(d[None, :], a[None, :])
    shifted = base + c * (np.sin(d)[:, None] * e_d[0] - np.cos(d)[:, None] * n_hat[0])
    delta = shifted - base
    assert np.ptp(delta, axis=0).max() < 1e-12          # one constant vector
    assert np.linalg.norm(delta[0]) == pytest.approx(c, rel=1e-9)


def test_no_monopole_and_no_dipole(surface):
    r = measure_area_and_multipoles(surface, t=1.2, gain=1e-3)
    assert abs(r["monopole"]) < 1e-6
    assert abs(r["dipole"]) < 1e-9


def test_area_is_unchanged_at_first_order(surface):
    """Trace-free means the deformed surface keeps its area to ``O(ε²)``."""
    r = measure_area_and_multipoles(surface, t=1.2, gain=1e-3)
    assert r["area_is_first_order_unchanged"] is True
    assert r["relative_area_change"] < 1e-4


# ── the bulk ────────────────────────────────────────────────────────────────
def test_the_wave_reaches_into_the_bulk_without_touching_it(surface):
    r = measure_bulk_reach(surface, frames=80)
    assert r["stays_between_the_dolls"] is True
    assert r["max_outward_fraction"] > 0.2       # it does reach, visibly
    assert r["max_outward_fraction"] < 1.0


def test_profiles_are_finite_and_share_one_grid(surface):
    surface.reset()
    surface.advance_to(0.8)
    p = surface.profiles()
    assert p["d"].shape == p["shear"].shape == p["radial"].shape == p["slide"].shape
    for key in ("shear", "radial", "slide"):
        assert np.all(np.isfinite(p[key]))


def test_mesh_and_material_curves_are_finite(surface):
    surface.reset()
    surface.advance_to(1.1)
    X, Y, Z = surface.mesh()
    assert np.all(np.isfinite(X)) and np.all(np.isfinite(Y)) and np.all(np.isfinite(Z))
    rings = surface.material_latitudes(n_rings=5, n_points=40)
    assert len(rings) == 5
    for r in rings:
        assert r.shape == (40, 3)
        assert np.all(np.isfinite(r))
    assert len(surface.material_meridians(n_lines=4, n_points=30)) == 4


def test_the_deformation_scales_with_the_gain(surface):
    surface.reset()
    surface.advance_to(1.0)
    d = np.array([[1.0]])
    a = np.array([[0.3]])
    base = surface.positions(d, a, gain=0.0)[0, 0]
    one = surface.positions(d, a, gain=1e-3)[0, 0] - base
    two = surface.positions(d, a, gain=2e-3)[0, 0] - base
    assert np.allclose(two, 2.0 * one, rtol=1e-9)
    assert np.linalg.norm(base) == pytest.approx(surface.shells.r_mid)


# ── principal axes ──────────────────────────────────────────────────────────
def test_the_two_eigenvalues_are_equal_and_opposite(surface):
    """Trace-free in two dimensions means ``λ± = ±|h|`` — always.

    So the two bars of a principal cross have the *same* length; an asymmetric
    cross would mean the tensor had a trace.
    """
    surface.reset()
    surface.advance_to(1.2)
    a = surface.principal_axes(np.array([0.6, 1.0, 1.5, 2.0]))
    assert np.allclose(a["lambda_plus"], -a["lambda_minus"])
    assert np.all(a["lambda_plus"] >= 0.0)
    assert np.allclose(a["lambda_plus"], np.abs(a["h_plus"]))


def test_the_stretch_axis_follows_the_sign_of_h(surface):
    """It swaps between ``ê_d`` and ``ê_ψ`` rather than staying put."""
    surface.reset()
    surface.advance_to(1.2)
    p = surface.profiles()
    d_pos = float(p["d"][int(np.argmax(p["shear"]))])
    d_neg = float(p["d"][int(np.argmin(p["shear"]))])
    assert float(np.interp(d_pos, p["d"], p["shear"])) > 0.0
    assert float(np.interp(d_neg, p["d"], p["shear"])) < 0.0
    assert surface.principal_axes(np.array([d_pos]))["angle"][0] == pytest.approx(0.0)
    assert surface.principal_axes(np.array([d_neg]))["angle"][0] == pytest.approx(
        0.5 * math.pi)
    # ...and the drawn bar swaps with it
    e1, _ = surface.tangent_basis(d_pos, 0.3)
    c_pos = surface.principal_cross(d_pos, 0.3)
    c_neg = surface.principal_cross(d_neg, 0.3)
    assert abs(float(c_pos["stretch_direction"] @ e1)) == pytest.approx(1.0, abs=1e-9)
    e1n, _ = surface.tangent_basis(d_neg, 0.3)
    assert abs(float(c_neg["stretch_direction"] @ e1n)) == pytest.approx(0.0, abs=1e-9)


def test_the_cross_is_orthogonal_tangent_and_centred(surface):
    surface.reset()
    surface.advance_to(1.2)
    c = surface.principal_cross(1.3, 0.8, size=0.09)
    u, v = c["stretch_direction"], c["squeeze_direction"]
    assert float(u @ v) == pytest.approx(0.0, abs=1e-12)
    assert np.linalg.norm(u) == pytest.approx(1.0)
    assert np.linalg.norm(v) == pytest.approx(1.0)
    e1, e2 = surface.tangent_basis(1.3, 0.8)
    normal = np.cross(e1, e2)
    assert float(u @ normal) == pytest.approx(0.0, abs=1e-9)   # lies in the surface
    assert float(v @ normal) == pytest.approx(0.0, abs=1e-9)
    for key in ("stretch", "squeeze"):
        assert np.allclose(c[key].mean(axis=0), c["centre"])
    assert np.linalg.norm(c["stretch"][1] - c["stretch"][0]) == pytest.approx(
        np.linalg.norm(c["squeeze"][1] - c["squeeze"][0]))


def test_bar_length_is_proportional_to_the_eigenvalue(surface):
    surface.reset()
    surface.advance_to(1.2)
    p = surface.profiles()
    ref = float(np.max(np.abs(p["shear"])))
    rows = []
    for d0 in (0.8, 1.2, 1.6):
        c = surface.principal_cross(d0, 0.2, size=0.1, reference=ref)
        rows.append((abs(c["lambda_plus"]), c["half_length"]))
    for lam, half in rows:
        assert half == pytest.approx(0.1 * lam / ref, rel=1e-9)
    # a vanishing eigenvalue draws no bar at all
    assert rows[0][1] < rows[-1][1]
