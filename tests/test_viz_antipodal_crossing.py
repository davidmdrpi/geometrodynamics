"""Smoke tests for the S³ antipodal worldline-crossing visualisation."""

import matplotlib

matplotlib.use("Agg")

import numpy as np

from geometrodynamics.transaction.s3_geometry import antipode4, geo4
from geometrodynamics.viz import AntipodalCrossingSim, run_animation
from geometrodynamics.viz.antipodal_crossing import (
    EPS_CROSS,
    EPS_PHASE,
    w_cross,
    w_phase,
)


def test_weights_peak_at_zero():
    assert w_cross(0.0) == 1.0
    assert w_phase(0.0) == 1.0
    # Peak at π (other closure branch)
    assert w_phase(np.pi) == 1.0
    assert w_cross(EPS_CROSS) < 1.0
    assert w_cross(6 * EPS_CROSS) < 1e-7


def test_sim_steps_without_error():
    sim = AntipodalCrossingSim(seed=0)
    sim.run(n_steps=300)
    assert sim.t > 0.0
    # Particles stay on S³ (unit 4-vectors)
    for p in sim.particles:
        assert abs(np.linalg.norm(p.p4) - 1.0) < 1e-9
    # Energy history is recorded each step
    assert len(sim.energy_history) == 300


def test_sim_fires_at_least_one_tx():
    sim = AntipodalCrossingSim(seed=42)
    sim.run(n_steps=1500)
    assert sim.tx_count() >= 1
    ab = sim.absorptions[0]
    # A fired TX satisfies the crossing gate
    assert ab.chi < EPS_CROSS
    assert ab.src_id != ab.dst_id


def test_best_chi_monotone_bounded():
    sim = AntipodalCrossingSim(seed=0)
    sim.run(n_steps=200)
    chi = sim.best_chi_now()
    assert 0.0 <= chi <= np.pi


def test_crossings_geometry_consistent():
    """Each recorded Crossing has xᵢ_bar ≈ −xᵢ(t') on a particle trail."""
    sim = AntipodalCrossingSim(seed=7)
    sim.run(n_steps=2000)
    for cr in sim.crossings:
        src = sim.particles[cr.src_id]
        # xi_bar should be close to antipode of some trail sample on the src
        if not src.trail:
            continue
        min_d = min(
            geo4(cr.xi_bar, antipode4(q)) for _t, q in src.trail
        )
        assert min_d < 1e-6


def test_run_animation_headless(tmp_path):
    sim = AntipodalCrossingSim(seed=1)
    sim.run(n_steps=400)
    ani, sim = run_animation(sim=sim, n_steps=1, interval_ms=16, show=False)
    fig = ani._fig
    fig.canvas.draw()  # triggers the full render pipeline
    out = tmp_path / "frame.png"
    fig.savefig(out, dpi=60, facecolor=fig.get_facecolor())
    assert out.exists() and out.stat().st_size > 1000
