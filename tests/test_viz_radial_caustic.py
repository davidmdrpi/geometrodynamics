"""
Front-topology tests for the pulse-vs-ring argument.

These pin what the module actually claims:

* in a **flat** bulk a point's front never folds — a statement about the bulk,
  not about point sources, since the same source folds on a closed manifold;
* a circle's front folds at ``t = ρ``, **detected from the front's own area
  element** rather than from a stored radius of curvature, and then *stays*
  singular on the symmetry axis;
* the first caustic is infinitely degenerate: the whole ring at once;
* the ring whose caustic lands on the inner sphere is exactly the ring whose
  rays graze it — the core result;
* acceptance across the bulk is asymmetric in **solid angle** while individual
  rays remain reversible;
* and a curved bulk is refused unless its effective radius is monotone.
"""

from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")  # must precede pyplot import

import numpy as np
import matplotlib.pyplot as plt
import pytest

from geometrodynamics.viz.radial_caustic import (
    PointSource,
    RingSource,
    ShellGeometry,
    detect_fold,
    draw_front,
    draw_shell,
    measure_acceptance_asymmetry,
    measure_critical_ring,
    measure_front_topology,
    plot_acceptance_cone,
    plot_pulse_vs_ring,
    signed_area_element,
)

R_IN, R_OUT = 0.74, 1.26


@pytest.fixture(autouse=True)
def _cleanup_figs():
    yield
    plt.close("all")


@pytest.fixture
def shell():
    return ShellGeometry(R_IN, R_OUT)


# ── the vacuole ─────────────────────────────────────────────────────────────
def test_shell_validates_its_radii():
    with pytest.raises(ValueError):
        ShellGeometry(1.2, 1.0)
    with pytest.raises(ValueError):
        ShellGeometry(0.0, 1.0)


def test_gap_and_grazing_angle(shell):
    assert shell.gap == pytest.approx(R_OUT - R_IN)
    assert shell.critical_sin == pytest.approx(R_IN / R_OUT)
    assert shell.critical_angle == pytest.approx(math.asin(R_IN / R_OUT))


@pytest.mark.parametrize("b", [0.0, 0.3, 0.74, 1.1])
def test_flat_bulk_turning_radius_is_the_impact_parameter(shell, b):
    assert shell.turning_radius(b) == pytest.approx(b)


def test_curved_bulk_uses_the_effective_radius():
    """A monotone metric factor rescales r; the structure is unchanged."""
    flat = ShellGeometry(R_IN, R_OUT)
    curved = ShellGeometry(R_IN, R_OUT, f=lambda r: np.full_like(r, 0.25))
    # f = 1/4 doubles every effective radius, so the *ratio* is untouched
    assert curved.effective_radius(R_IN) == pytest.approx(2.0 * R_IN)
    assert curved.critical_sin == pytest.approx(flat.critical_sin)
    assert curved.turning_radius(2.0 * 0.9) == pytest.approx(0.9, abs=1e-6)


def test_non_monotone_effective_radius_is_refused():
    """A photon sphere breaks every closed form here, so reject it loudly."""
    # Schwarzschild-like: R_eff = r/sqrt(1 - rs/r) has a minimum at 1.5 rs,
    # which for rs = 0.5 sits at 0.75 — inside this shell.
    with pytest.raises(ValueError, match="monotone|increase"):
        ShellGeometry(R_IN, R_OUT, f=lambda r: 1.0 - 0.5 / r)


def test_monotone_curved_bulk_is_accepted():
    sh = ShellGeometry(R_IN, R_OUT, f=lambda r: np.full_like(r, 0.64))
    assert sh.critical_sin == pytest.approx(R_IN / R_OUT)


# ── a point never folds ─────────────────────────────────────────────────────
def test_point_front_does_not_fold_in_a_flat_bulk(shell):
    """Scoped deliberately: this is a fact about the flat bulk.

    The same point source folds on a closed manifold — ``S²``/``S³`` focus a
    point's front on the antipode at ``t = πR``, which
    ``test_viz_throat_wavefront`` measures.  Nothing here contradicts that.
    """
    pulse = shell.point_source()
    assert pulse.first_caustic_point() is None
    assert pulse.fold_time is None
    assert pulse.singular_points(5.0).shape == (0, 3)
    top = measure_front_topology(pulse, t_hi=2.2)
    assert top["detected_folds"] is False
    assert top["detected_fold_time"] is None


def test_point_front_area_element_keeps_one_sign(shell):
    """t² sin θ up to orientation — it never passes through zero."""
    pulse = shell.point_source()
    for t in (0.2, 0.8, 1.6):
        j = signed_area_element(pulse, t)
        assert np.min(j) * np.max(j) > 0.0        # no sign change anywhere


def test_point_front_is_the_metric_sphere(shell):
    pulse = shell.point_source()
    f = pulse.front(0.4)
    d = np.linalg.norm(f - pulse.position, axis=1)
    assert np.allclose(d, 0.4)
    assert pulse.arrival_multiplicity(pulse.position + np.array([0.4, 0, 0]), 0.4) == 1
    assert pulse.arrival_multiplicity(pulse.position + np.array([0.9, 0, 0]), 0.4) == 0


def test_pulse_crosses_the_bulk_at_the_gap(shell):
    assert shell.point_source().crossing_time == pytest.approx(shell.gap)


# ── a ring must fold ────────────────────────────────────────────────────────
def test_ring_first_caustic_is_its_centre(shell):
    ring = shell.critical_ring()
    assert np.allclose(ring.first_caustic_point(), ring.centre)
    assert ring.fold_time == pytest.approx(ring.radius)


def test_ring_stays_singular_after_the_first_caustic(shell):
    """The fold is not an isolated event — two axis points separate as √(t²−ρ²)."""
    ring = shell.critical_ring()
    assert ring.singular_points(0.9 * ring.radius).shape == (0, 3)
    at = ring.singular_points(ring.radius)
    assert at.shape == (2, 3)
    assert at[0, 2] == pytest.approx(at[1, 2])            # coincide at the fold
    t = 1.4 * ring.radius
    after = ring.singular_points(t)
    sep = abs(after[0, 2] - after[1, 2])
    assert sep == pytest.approx(2.0 * math.sqrt(t ** 2 - ring.radius ** 2))
    assert np.allclose(after[:, :2], 0.0)                 # on the axis


@pytest.mark.parametrize("theta0", [0.3, 0.7, 1.1, 1.4])
def test_fold_detected_from_the_area_element_alone(shell, theta0):
    """The topology check must not be seeded by the answer it is testing."""
    ring = RingSource(shell=shell, polar_angle=theta0)
    det = detect_fold(ring, t_hi=1.6)
    assert det["folds"] is True
    assert det["fold_time"] == pytest.approx(ring.radius, abs=1e-4)


def test_detector_is_relative_and_orientation_referenced(shell):
    """Two ways the naive detector gets it wrong, both guarded.

    An absolute ``min J <= 0`` test fires on the direction sphere's polar
    degeneracy; and ``(X_u × X_v)·N`` is negative everywhere for the point
    source purely from parameter ordering.
    """
    pulse = shell.point_source()
    j = signed_area_element(pulse, 0.5)
    # the element nearly vanishes at the trimmed poles (sin θ → 0 there), so an
    # absolute "min <= 0" test would misfire; and it is negative everywhere,
    # purely from the (u, v) ordering, so an unreferenced sign test would too
    assert np.min(np.abs(j)) < 1e-2 * np.max(np.abs(j))
    assert np.max(j) < 0.0
    assert detect_fold(pulse, t_hi=2.0)["folds"] is False  # still no fold


def test_the_whole_ring_reaches_its_centre_at_once(shell):
    """Equidistance is what makes the caustic degenerate rather than smooth."""
    ring = shell.critical_ring()
    d = np.linalg.norm(ring.points(4096) - ring.centre, axis=1)
    assert np.ptp(d) < 1e-12
    assert d[0] == pytest.approx(ring.radius)
    # reported as -1: "degenerate — every ray at once"
    assert ring.arrival_multiplicity(ring.centre, ring.radius) == -1


def test_ring_multiplicity_is_two_just_off_the_focus(shell):
    ring = shell.critical_ring()
    off = ring.centre + np.array([0.05, 0.0, 0.0])
    assert ring.arrival_multiplicity(off, ring.radius) == 2
    far = ring.centre + np.array([3.0, 0.0, 0.0])
    assert ring.arrival_multiplicity(far, ring.radius) == 0


def test_front_topology_compares_detection_to_the_closed_form(shell):
    ring = shell.critical_ring()
    top = measure_front_topology(ring, t_hi=1.6)
    assert top["detected_folds"] is True
    assert top["detection_error"] < 1e-4
    assert top["singular_points_after"] == 2


# ── the coincidence ─────────────────────────────────────────────────────────
def test_critical_ring_puts_its_defect_on_the_inner_sphere(shell):
    r = measure_critical_ring(shell)
    assert r["caustic_radius"] == pytest.approx(shell.r_inner, abs=1e-12)
    assert r["fold_time"] == pytest.approx(
        math.sqrt(R_OUT ** 2 - R_IN ** 2), abs=1e-12)


def test_the_focusing_ring_is_the_grazing_ring(shell):
    """One condition, not two — the load-bearing coincidence."""
    r = measure_critical_ring(shell)
    assert r["launch_sin"] == pytest.approx(shell.critical_sin, abs=1e-12)
    assert r["ray_turning_radius"] == pytest.approx(shell.r_inner, abs=1e-12)


@pytest.mark.parametrize("r_in", [0.2, 0.45, 0.74, 0.95, 1.15])
def test_the_coincidence_holds_at_every_ratio(r_in):
    sh = ShellGeometry(r_in, R_OUT)
    ring = sh.critical_ring()
    assert ring.centre_radius == pytest.approx(r_in, abs=1e-12)
    assert ring.launch_sin == pytest.approx(sh.critical_sin, abs=1e-12)
    # and the ring always folds later than the pulse crosses
    assert ring.fold_time > sh.gap


def test_critical_ring_needs_a_flat_bulk():
    curved = ShellGeometry(R_IN, R_OUT, f=lambda r: np.full_like(r, 0.25))
    with pytest.raises(ValueError):
        curved.critical_ring()


# ── the asymmetry ───────────────────────────────────────────────────────────
def test_acceptance_closed_form_matches_monte_carlo(shell):
    r = measure_acceptance_asymmetry(shell, n_samples=200000, seed=3)
    assert r["inward_closed_form"] == pytest.approx(r["inward_monte_carlo"],
                                                    abs=0.01)
    assert r["outward_closed_form"] == 1.0
    assert r["outward_monte_carlo"] == 1.0
    assert r["solid_angle_ratio"] > 2.0
    # the asymmetry is in the measure of directions, not in propagation
    assert r["rays_reversible"] is True


def test_inward_acceptance_tightens_as_the_inner_sphere_shrinks():
    fracs = [ShellGeometry(r, R_OUT).acceptance_fraction("inward")
             for r in (0.2, 0.5, 0.9, 1.2)]
    assert all(fracs[i] < fracs[i + 1] for i in range(len(fracs) - 1))
    assert fracs[0] < 0.05


def test_acceptance_direction_is_validated(shell):
    with pytest.raises(ValueError):
        shell.acceptance_fraction("sideways")


def test_outward_rays_always_escape(shell):
    """b <= R_eff(r_inner) < R_eff(r_outer), so nothing turns back."""
    for alpha in np.linspace(0.0, 0.5 * math.pi - 1e-9, 30):
        b = shell.impact_parameter(shell.r_inner, alpha)
        assert b <= shell.r_inner + 1e-12
        assert b < shell.r_outer


# ── rendering ───────────────────────────────────────────────────────────────
def test_draw_helpers_render(shell):
    fig, ax = plt.subplots()
    draw_shell(ax, shell)
    draw_front(ax, shell.point_source(), 0.4)
    draw_front(ax, shell.critical_ring(), 0.7)
    assert len(ax.lines) >= 4


def test_summary_figures_build(shell):
    f1 = plot_pulse_vs_ring(shell, figsize=(8.0, 4.0))
    assert len(f1.axes) == 2
    f2 = plot_acceptance_cone(shell, figsize=(8.0, 4.0))
    assert len(f2.axes) == 2
