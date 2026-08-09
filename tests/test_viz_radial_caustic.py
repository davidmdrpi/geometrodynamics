"""
Front-topology tests for the ring-vs-pulse argument.

These pin the claim that a throat needs a *ring*:

* the focal set of a point is empty, so its front is an embedded sphere at
  every time and never touches itself;
* the focal set of a circle is a single point — its centre — which the whole
  ring reaches simultaneously, so the front stops being embedded there;
* the ring whose defect lands on the inner sphere is exactly the ring whose
  rays graze it;
* and the bulk accepts far less going in than coming out, from the ordering
  of the two radii alone.
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
    draw_front,
    draw_shell,
    measure_acceptance_asymmetry,
    measure_critical_ring,
    measure_front_topology,
    plot_acceptance_cone,
    plot_pulse_vs_ring,
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
    """A metric factor rescales r; the structure is unchanged."""
    flat = ShellGeometry(R_IN, R_OUT)
    curved = ShellGeometry(R_IN, R_OUT, f=lambda r: np.full_like(r, 0.25))
    # f = 1/4 doubles every effective radius, so the *ratio* is untouched
    assert curved.effective_radius(R_IN) == pytest.approx(2.0 * R_IN)
    assert curved.critical_sin == pytest.approx(flat.critical_sin)
    assert curved.turning_radius(2.0 * 0.9) == pytest.approx(0.9, abs=1e-6)


# ── a point never folds ─────────────────────────────────────────────────────
def test_point_focal_set_is_empty_and_front_never_self_intersects(shell):
    pulse = shell.point_source()
    assert pulse.focal_set().shape == (0, 3)
    assert pulse.self_intersection_time is None
    top = measure_front_topology(pulse, np.linspace(0.05, 2.2, 40))
    assert top["ever_self_intersects"] is False
    assert top["degenerate_times"] == []
    assert max(r["off_axis"] for r in top["rows"]) <= 1


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
def test_ring_focal_set_is_its_centre(shell):
    ring = shell.critical_ring()
    fs = ring.focal_set()
    assert fs.shape == (1, 3)
    assert np.allclose(fs[0], ring.centre)
    assert ring.self_intersection_time == pytest.approx(ring.radius)


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


def test_front_topology_samples_the_fold_even_off_grid(shell):
    """The fold is measure-zero in t, so the sampler must insert it."""
    ring = shell.critical_ring()
    top = measure_front_topology(ring, np.array([0.1, 0.5, 1.4]))
    assert len(top["degenerate_times"]) == 1
    assert top["degenerate_times"][0] == pytest.approx(ring.radius)


# ── the coincidence ─────────────────────────────────────────────────────────
def test_critical_ring_puts_its_defect_on_the_inner_sphere(shell):
    r = measure_critical_ring(shell)
    assert r["defect_radius"] == pytest.approx(shell.r_inner, abs=1e-12)
    assert r["defect_time"] == pytest.approx(
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
    assert ring.self_intersection_time > sh.gap


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
    assert r["asymmetry_ratio"] > 2.0


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
