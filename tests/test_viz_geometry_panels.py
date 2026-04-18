"""
Smoke tests for the four-panel geometry dashboard.

These verify the plotters wire up cleanly to the underlying package
primitives and render under the headless ``Agg`` backend.  They also
check a few numerical invariants that the panels visualise (A at the
equator, holonomy at the pole, Green-function divergences).
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # must precede pyplot import

import numpy as np
import matplotlib.pyplot as plt
import pytest

from geometrodynamics.hopf import hopf_connection, hopf_curvature, hopf_holonomy
from geometrodynamics.transaction import s3_green_field_kernel, s3_green_potential
from geometrodynamics.viz import (
    draw_green,
    draw_handshake,
    draw_hopf,
    draw_throat,
    plot_green_panel,
    plot_handshake_panel,
    plot_hopf_panel,
    plot_throat_panel,
)


@pytest.fixture(autouse=True)
def _cleanup_figs():
    yield
    plt.close("all")


def test_hopf_connection_vanishes_at_equator():
    assert abs(float(hopf_connection(np.pi / 2))) < 1e-12
    assert abs(float(hopf_curvature(np.pi / 2)) - 0.5) < 1e-12


def test_hopf_holonomy_at_poles_is_plus_minus_pi():
    assert float(hopf_holonomy(0.0)) == pytest.approx(np.pi)
    assert float(hopf_holonomy(np.pi)) == pytest.approx(-np.pi)


def test_green_field_peaks_near_source():
    # The analytic dG/dψ = −[cot ψ + (π−ψ) csc² ψ] / (4π²R) is
    # Coulomb-like near ψ=0 and cancels to a finite value near ψ=π
    # (the two terms flip sign and partially cancel).
    near = s3_green_field_kernel(0.08 + 1e-6)
    mid = s3_green_field_kernel(np.pi / 2)
    far = s3_green_field_kernel(np.pi - 0.08 - 1e-6)
    assert near > 10 * mid  # Coulomb divergence at the source
    assert far < mid        # no divergence at the antipode


def test_hopf_panel_draws_three_axes():
    fig = plot_hopf_panel(chi=0.42)
    assert len(fig.axes) == 3
    # 3D fibre axes + two 1D param plots
    assert any(ax.name == "3d" for ax in fig.axes)


def test_throat_panel_uses_solver():
    fig = plot_throat_panel(l=3)
    assert len(fig.axes) == 2
    # Mode curve should land on the lower axis via the solver — check
    # there is at least one line above a nonzero amplitude threshold.
    ax_plot = fig.axes[1]
    has_nontrivial_line = any(
        np.nanmax(np.abs(line.get_ydata())) > 0.05 for line in ax_plot.get_lines()
    )
    assert has_nontrivial_line


def test_green_panel_y_scale_log_on_field():
    fig = plot_green_panel(src_angle=0.7)
    # Third axis is |dG/dψ| plotted on a log scale
    assert fig.axes[2].get_yscale() == "log"


def test_handshake_panel_wraps_step():
    fig_a = plot_handshake_panel(step=0.5)
    fig_b = plot_handshake_panel(step=4.5)   # 4.5 mod 4 == 0.5
    # Both drew something valid
    assert len(fig_a.axes) == 1
    assert len(fig_b.axes) == 1


def test_draw_functions_rebind_same_figure():
    fig = plt.figure()
    draw_hopf(fig, chi=0.6)
    n1 = len(fig.axes)
    draw_throat(fig, l=1)
    n2 = len(fig.axes)
    assert n1 == 3
    assert n2 == 2  # throat uses 2 axes; old hopf axes should be cleared


def test_numeric_potentials_are_finite():
    for psi in np.linspace(0.1, np.pi - 0.1, 5):
        g = s3_green_potential(psi)
        dg = s3_green_field_kernel(psi)
        assert np.isfinite(g)
        assert np.isfinite(dg)


def test_handshake_all_phases_renderable():
    # One figure for each of the four narrative phases.
    for step in (0.2, 1.5, 2.8, 3.9):
        fig = plot_handshake_panel(step=step)
        assert len(fig.axes) == 1
        plt.close(fig)
