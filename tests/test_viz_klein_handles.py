"""Klein-type and torus-type antipodal handles: topology, geometry and waves."""

import math

import numpy as np
import pytest

from geometrodynamics.viz.klein_handles import (
    HandleWaveSim,
    Neck,
    build_surface,
    default_handles,
    euler_characteristic,
    is_orientable,
    orientation_loop,
)


@pytest.mark.parametrize("kinds,chi,orientable", [
    ((), 2, True),
    ((True,), 0, False),        # Klein bottle
    ((False,), 0, True),        # torus
    ((True, True), -2, False),  # connected sum of two Klein bottles
    ((False, False), -2, True), # genus-two surface
    ((True, False), -2, False),
])
def test_topology(kinds, chi, orientable):
    surf = build_surface(0.3, default_handles(len(kinds), list(kinds)) if kinds else ())
    assert euler_characteristic(surf.tris) == chi
    assert is_orientable(surf.tris) is orientable


def test_neck_matches_finite_mouth():
    nk = Neck(0.3)
    assert nk.b == pytest.approx(math.sin(0.3) ** 2)
    assert nk.S == pytest.approx(math.sin(0.3) * math.cos(0.3))
    assert float(nk.r(nk.S)) == pytest.approx(math.sin(0.3))
    assert nk.S / float(nk.r(nk.S)) == pytest.approx(math.cos(0.3))          # slope r'(S) = cos a
    assert float(nk.gaussian_curvature(nk.S)) == pytest.approx(-1.0)         # rim, neck side
    assert 2 * nk.S < 2 * math.cos(0.3)                                       # shorter than the chord


def test_area_and_operators():
    surf = build_surface(0.3, default_handles(2, [True, True]))
    K, mass = surf.operators()
    nk = surf.neck
    exact = 4 * math.pi - 4 * 2 * math.pi * (1 - math.cos(0.3)) + 2 * nk.area
    assert mass.sum() == pytest.approx(exact, rel=2e-3)
    assert abs(K - K.T).max() < 1e-12
    assert np.abs(np.asarray(K.sum(axis=1))).max() < 1e-10
    assert mass.min() > 0


@pytest.mark.parametrize("klein", [True, False])
def test_drawing_lands_on_the_glued_rim(klein):
    surf = build_surface(0.3, default_handles(2, [klein, klein]))
    n = surf.n_ring
    psi = 2 * math.pi * np.arange(n) / n
    perm = (n // 2 - np.arange(n)) % n if klein else (np.arange(n) + n // 2) % n
    for d, (_, rb) in zip(surf.drawn, surf.rims):
        drawn_end = d.point(np.full(n, -surf.neck.S), psi)
        assert np.abs(drawn_end - surf.points[rb[perm]]).max() < 1e-9


@pytest.mark.parametrize("klein", [True, False])
def test_orientation_loop_matches_topology(klein):
    surf = build_surface(0.3, default_handles(2, [klein, klein]))
    for h in range(2):
        _, _, flipped = orientation_loop(surf, h)
        assert flipped is klein


def test_bare_sphere_first_arrival_at_antipode():
    surf = build_surface(0.3, ())
    sim = HandleWaveSim(surf, source=(0, 0, 1), t_max=3.6, n_frames=181, threshold=0.05).run()
    j = int(np.argmin(np.linalg.norm(surf.points - np.array([0, 0, -1.0]), axis=1)))
    assert math.pi - sim.width <= sim.arrival[j] <= math.pi   # inside the cone, before the focus
    assert sim.energy_drift < 1e-2


def _probe(d, th, ph):
    return math.cos(th) * d.c + math.sin(th) * (math.cos(ph) * d.e1 + math.sin(ph) * d.e2)


def test_wave_takes_the_neck_shortcut():
    surf = build_surface(0.3, default_handles(1, [True]))
    d = surf.drawn[0]
    sim = HandleWaveSim(surf, source=_probe(d, 0.65, 0.0), t_max=2.5, n_frames=126).run()
    j = int(np.argmin(np.linalg.norm(surf.points - _probe(d, math.pi - 0.65, math.pi), axis=1)))
    around = math.pi - sim.width                     # first arrival the long way, over the sphere
    assert sim.arrival[j] < 1.5 < around


@pytest.mark.parametrize("klein", [True, False])
def test_first_arrival_respects_the_causal_bound(klein):
    """No sphere point is reached before min(direct, through-the-neck) distance − w.

    Through the neck the shortest path is at least (angle to mouth − a) + 2S + (angle from
    the other mouth − a); over the sphere it is at least the great-circle distance.
    """
    surf = build_surface(0.3, default_handles(1, [klein]))
    d, nk = surf.drawn[0], surf.neck
    sim = HandleWaveSim(surf, source=_probe(d, 0.65, 0.0), t_max=3.0, n_frames=301).run()
    x = surf.points[surf.on_sphere]
    ang = lambda p, q: np.arccos(np.clip(p @ q, -1, 1))
    via = lambda c_in, c_out: (ang(sim.source, c_in) - 0.3) + 2 * nk.S + (ang(x, c_out) - 0.3)
    bound = np.minimum(ang(x, sim.source), np.minimum(via(d.c, -d.c), via(-d.c, d.c))) - sim.width
    arrival = sim.arrival[surf.on_sphere]
    assert not np.isnan(arrival).any()
    assert np.min(arrival - bound) > 0.0


def test_bare_sphere_arrival_is_causal_away_from_the_caustic():
    surf = build_surface(0.3, ())
    sim = HandleWaveSim(surf, source=(0, 0, 1), t_max=3.0, n_frames=301).run()
    dist = np.arccos(np.clip(surf.points[:, 2], -1, 1))
    away = dist < 2.4
    assert np.min(sim.arrival[away] - (dist[away] - sim.width)) > 0.0


@pytest.mark.parametrize("klein", [True, False])
def test_gluing_class_is_visible_in_the_emerging_front(klein):
    """Klein-type: the front leaves mouth B mirrored (pi − phi0); torus-type: rotated (pi + phi0)."""
    phi0 = 1.0
    surf = build_surface(0.3, default_handles(1, [klein]))
    d = surf.drawn[0]
    sim = HandleWaveSim(surf, source=_probe(d, 0.62, phi0), t_max=2.0, n_frames=201).run()
    win = (sim.times > 0.7) & (sim.times < 1.6)
    peak = {}
    for label, ph in (("rot", math.pi + phi0), ("mirror", math.pi - phi0)):
        j = int(np.argmin(np.linalg.norm(surf.points - _probe(d, math.pi - 0.45, ph), axis=1)))
        peak[label] = float(np.abs(sim.frames[win, j]).max())
    assert (peak["mirror"] > 1.1 * peak["rot"]) if klein else (peak["rot"] > 1.1 * peak["mirror"])
