"""Checks of the rim-kick bookkeeping on a coarse glued sphere (fast versions of the probe checks)."""
import numpy as np
import pytest

from geometrodynamics.viz.klein_handles import HandleWaveSim, Neck
from geometrodynamics.waves.rim_kicks import (AXIS, FluxField, GluedWave, RimKicks, displacements, follow,
                                              glued_surface, offset_history, release_source)

COARSE = dict(n_sphere=3000, n_ring=24, ds=.03)


@pytest.fixture(scope='module')
def klein():
    return glued_surface('klein', **COARSE)


@pytest.fixture(scope='module')
def torus():
    return glued_surface('torus', **COARSE)


@pytest.fixture(scope='module')
def release_kick(klein):
    g = GluedWave(klein)
    u, v, _ = g.packet(1., np.pi/2, 'release')
    return g.run(u, v, 3.7)


def test_bare_sphere_conserves_band_angular_momentum():
    g = GluedWave(glued_surface('none', **COARSE))
    u, v, _ = g.packet(1., .7, 'out')
    h = g.run(u, v, 2.5)
    assert np.abs(h.L-h.L[0]).max() < .02*np.linalg.norm(h.L[0])
    assert np.abs(h.energy-h.energy[0]).max() < 1e-3
    # No mouths: the weighted fluxes cancel the weighted contents up to discretisation error.
    assert np.abs(h.J).max() < .03*np.linalg.norm(h.L[0])


def test_axisymmetric_ring_gives_no_sideways_kick():
    g = GluedWave(glued_surface('klein', n_sphere=6000, n_ring=24, ds=.03))
    u, v = g.ring(.8)
    a, b = displacements(g.run(u, v, 3.7).J[-1])
    u, v, _ = g.packet(1., np.pi/2, 'release')
    a0, _ = displacements(g.run(u, v, 3.7).J[-1])
    # Residual from the mesh's lack of rotational symmetry: about 1% at the probe's 24k vertices.
    assert max(np.linalg.norm(a), np.linalg.norm(b)) < .08*np.linalg.norm(a0)


def test_release_moves_both_mouths_the_same_way(release_kick):
    a, b = displacements(release_kick.J[release_kick.at(np.pi+.5)])
    u = np.array([0., -1., 0.])                 # e2 for the handle side (0, 0, 1) at +x
    assert a @ u < 0 and b @ u < 0              # both pushed away from the release side
    assert follow(a, b) > .8


def test_klein_and_torus_agree_on_the_mirror_axis(klein, torus):
    out = []
    for s in (klein, torus):
        g = GluedWave(s)
        u, v, _ = g.packet(1., 0., 'in')
        out.append(g.run(u, v, 2.4).J[-1])
    assert np.allclose(out[0], out[1], atol=.02*np.abs(out[0]).max())


def test_kicks_are_quadratic_and_the_flux_term_bilinear(klein):
    g = GluedWave(klein, m2=0.)
    F = FluxField.unit_energy()
    gs = F.band_gradient(g.kicks)/np.sqrt(F.energy)
    u, v, _ = g.packet(1., np.pi/4, 'release')
    plus = g.run(u, v, 2.0, static=gs)
    minus = g.run(-u, -v, 2.0, static=gs)
    assert np.allclose(plus.J, minus.J, rtol=0, atol=1e-14)
    assert np.allclose(plus.Jx, -minus.Jx, rtol=0, atol=1e-14)
    f1, c1, _ = g.kicks.cross_rates(u, v, gs)
    f2, c2, _ = g.kicks.cross_rates(u, v, -2*gs)
    assert np.allclose(f2, -2*f1) and np.allclose(c2, -2*c1)


def test_flux_field_period_energy_and_flux(klein):
    F = FluxField(c=.7)
    nk = Neck(.3)
    band = 2*np.log(1/np.tan(.15))
    assert np.isclose(2*np.arcsinh(nk.S/nk.b), band)          # the neck carries half the period
    assert np.isclose(F.period, .7*2*band) and np.isclose(F.energy, np.pi*.7**2*2*band)
    k = RimKicks(klein, m2=0.)
    gs = F.band_gradient(k)
    th = np.arccos(np.clip(k.centroid @ AXIS, -1, 1))
    that = -(AXIS[None, :]-(k.centroid @ AXIS)[:, None]*k.centroid)/np.sin(th)[:, None]
    for th0 in (.6, 1.6, np.pi-.6):
        ring = np.abs(th-th0) < .06
        flux = (np.einsum('ti,ti->t', gs, that)[ring] @ k.area[ring])/.12
        assert abs(flux-2*np.pi*F.c) < .04*2*np.pi*F.c
    phi = F.band_potential(klein)
    f, _, _ = k.rates(phi, 0*phi)
    assert np.abs(f[:, 1:]).max() < 1e-3*F.energy             # no sideways self-kick


def test_cross_term_does_not_depend_on_the_rim_weight(klein):
    F = FluxField.unit_energy()
    out = []
    for weight in ((.35, .9), (.5, 1.3)):
        g = GluedWave(klein, m2=0., weight=weight)
        u, v, _ = g.packet(1., np.pi/4, 'release')
        out.append(g.run(u, v, 3.0, static=F.band_gradient(g.kicks)/np.sqrt(F.energy)).Jx[-1])
    assert np.abs(out[0]-out[1]).max() < .08*np.abs(out[0]).max()


def test_offset_route_matches_direct_evolution_at_zero_offset(klein):
    src = release_source(np.pi/4, 1., klein)
    h, info = offset_history(klein, src, (), t_max=3.0)
    g = GluedWave(klein, m2=0.)
    u0 = HandleWaveSim(klein, source=src, width=.25, t_max=1., n_frames=4).initial()
    u0 = u0/np.sqrt(g.energy(u0, 0*u0))
    d = g.run(u0, 0*u0, 3.0)
    assert info['tilt'] == 0
    assert np.abs(h.J[-1]-d.J[-1]).max() < .03*np.abs(d.J[-1]).max()


def test_negative_offset_reemits_through_mouth_b_earlier(klein):
    src = release_source(0., 1., klein)
    info = {off: offset_history(klein, src, (off,) if off else (), t_max=1.8, frame_dt=.005)[1] for off in (0., -1.)}
    assert info[-1.]['tilt'] < .95
    assert abs(info[-1.]['arrival_rim_A']-info[0.]['arrival_rim_A']) < .05     # entry at A is unchanged
    assert abs((info[-1.]['arrival_rim_B']-info[0.]['arrival_rim_B'])+1.) < .1  # exit at B is one unit earlier


def test_offset_release_is_not_differentiated_across_the_padding(klein):
    """With an offset the external-time frames are zero-padded before the release at t = 0.

    The jump at the source must not enter the rim flux: before any wave reaches a rim the
    flux stays small and does not grow as the frames are refined (it grew as 1/frame_dt).
    """
    src = release_source(0., 1., klein)
    a = []
    for frame_dt in (.005, .0025, .00125):
        h, info = offset_history(klein, src, (-1.,), t_max=.12, frame_dt=frame_dt)
        assert info['t_min'] < 0
        a.append(displacements(h.J[h.at(.02)])[0])
    assert np.abs(np.array(a)).max() < 1e-4
    assert np.abs(a[0]-a[2]).max() < 1e-5


def test_offset_flux_converges_in_frame_spacing_and_rim_weight():
    s = glued_surface('klein', n_sphere=6000, n_ring=24, ds=.03)
    src = release_source(0., 1., s)
    weights = [(.35, .9), (.35, .65)]
    coarse, _ = offset_history(s, src, (-1.,), width=.5, t_max=2., frame_dt=.005, weight=weights)
    fine, _ = offset_history(s, src, (-1.,), width=.5, t_max=2., frame_dt=.0025, weight=weights)
    ref = fine[0].J[fine[0].at(2.)]
    scale = np.abs(ref).max()
    assert np.abs(coarse[0].J[coarse[0].at(2.)]-ref).max() < .02*scale
    assert np.abs(fine[1].J[fine[1].at(2.)]-ref).max() < .08*scale
