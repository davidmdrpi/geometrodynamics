"""Checks of the pair-birth model (two Wheeler-flux mouths on S²) on a coarse mode set."""
import numpy as np
import pytest

from geometrodynamics.waves.moving_mouth_pairs import harmonics
from geometrodynamics.waves.pair_birth import (PairBirth, Parameters, bare_timeline, birth_positions,
                                               circle_angle, front_crossings)

COARSE = Parameters(cutoff=16, footprint=.2)
PLANE = ((0., 0., 1.), (1., 0., 0.))


@pytest.fixture(scope='module')
def model():
    return PairBirth(COARSE)


def independent_start(v=.3):
    x, vel, P, T = birth_positions(.3, v, PLANE)
    vel = vel.copy()
    vel[0] *= 1.15
    vel[1] = .85*vel[1] + .08*np.cross(P, T)
    return vel


def test_vector_potential_and_magnetic_force_match_finite_differences(model):
    rng = np.random.default_rng(3)
    x = rng.normal(size=(2, 3))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    vel = np.cross(x, rng.normal(size=(2, 3)))
    c = 3*model.source(x, vel)[0]/model.omega
    y = np.array([.3, .5, .81])
    y /= np.linalg.norm(y)

    def chi(z):
        return (model.smoothing*harmonics(model.ell, model.order, z/np.linalg.norm(z)) @ c).real

    def A(z):
        z = z/np.linalg.norm(z)
        return model._rot(model.footprints(np.array([z, z])), c)[0]

    eps = 1e-6
    grad = np.array([(chi(y+eps*e)-chi(y-eps*e))/(2*eps) for e in np.eye(3)])
    grad -= (grad @ y)*y
    assert np.allclose(A(y), np.cross(y, grad), atol=1e-7)
    v = np.cross(y, [.2, -.7, .4])
    gv = np.array([(v @ A(y+eps*e) - v @ A(y-eps*e))/(2*eps) for e in np.eye(3)])
    gv -= (gv @ y)*y
    vh = v/np.linalg.norm(v)
    dA = (A(np.cos(eps)*y+np.sin(eps)*vh) - A(np.cos(eps)*y-np.sin(eps)*vh))/(2*eps)*np.linalg.norm(v)
    dA -= (dA @ y)*y
    B = -(model.footprints(np.array([y, y]))[0].conj() @ (model.omega*c*model.omega)).real
    assert np.allclose(gv - dA, B*np.cross(v, y), atol=1e-7)          # q[grad(v.A) - (v.grad)A] = q B v x n


def test_energy_and_angular_momentum_are_conserved_for_independent_mouths(model):
    h, info = model.run(model.initial(velocities=independent_start()), 2.5)
    assert info['energy_drift'] < 1e-6
    assert np.abs(h['J']-h['J'][0]).max() < 1e-8
    assert np.abs(h['J'][0]).max() > 1e-2                     # a real test: the start carries angular momentum


def test_flux_mouths_attract_and_scalar_sources_repel(model):
    d = np.linspace(.6, np.pi, 200)
    U = model.coulomb_interaction(d)
    assert np.all(np.diff(U) > 0)                             # energy grows with separation: attraction
    _, force = model.coulomb(np.array([[0., 0, 1], [0, 0, -1]]))
    _, near = model.coulomb(np.array([[0., 0, 1], [np.sin(.6), 0, np.cos(.6)]]))
    assert np.abs(force).max() < 1e-12*np.abs(near).max()      # the pull vanishes at antipodes
    assert np.all(np.diff(model.scalar_interaction(d)) < 0)   # opposite scalar sources push apart


def test_free_mouths_meet_again_at_the_antipode():
    m = PairBirth(Parameters(cutoff=8, charge=0.))
    h, info = m.run(m.initial(.4), 8.)
    assert info['contact_time'] == pytest.approx((np.pi-.6)/.4, abs=1e-6)
    assert np.dot(info['contact_point'], [0, 0, 1]) < 0


def test_instantaneous_coulomb_turning_point_matches_the_energy_budget():
    m = PairBirth(Parameters(cutoff=16, footprint=.2, transverse=False))
    v = .2
    h, info = m.run(m.initial(v), 8.)
    kinetic = 2*(1/np.sqrt(1-v**2)-1)
    d = np.linspace(.6, np.pi, 20001)
    turn = d[np.argmax(m.coulomb_interaction(d)-m.coulomb_interaction(.6) >= kinetic)]
    assert info['min_separation'] <= .6 + 1e-9 and info['contact_time'] is not None
    assert h['d'].max() == pytest.approx(turn, abs=2e-3)
    assert np.dot(info['contact_point'], [0, 0, 1]) > 0         # back at the birth point


def test_birth_flash_is_the_free_wave_of_the_missing_dressing():
    """For weak charges the transverse force difference between the default and the dressed
    start is the post hoc force of the free flash wave along the trajectory."""
    m = PairBirth(Parameters(cutoff=16, footprint=.2, charge=.05))
    full, _ = m.run(m.initial(.3), 4., output_dt=.05)
    dressed, _ = m.run(m.initial(.3, dressed=True), 4., output_dt=.05)
    diff = (full['exchange']+full['own']) - (dressed['exchange']+dressed['own'])
    flash = m.free_wave_force(full, *m.flash(m.initial(.3))).sum(2)
    scale = np.abs(flash).max()
    assert scale > 0
    assert np.abs(diff-flash).max() < .02*scale
    # Re-timing by zero at any time is the identity; removing the flash late changes only its returns.
    same, info = m.run(m.initial(.3), 4., output_dt=.05, retime=(1., 0.))
    assert info['retime_time'] == pytest.approx(1.) and abs(info['retime_energy']) < 1e-15
    assert np.allclose(same['x'], full['x'], atol=1e-12)


def test_release_carries_energy_and_angular_momentum(model):
    h, info = model.run(model.initial(velocities=independent_start(.25)), 14.)
    assert info['contact_time'] is not None
    before, after = np.flatnonzero(np.abs(h['t']-info['contact_time']) < 1e-12)       # recorded on both sides
    e = h['energy_parts']
    assert e['rest'][before] + e['kinetic'][before] + e['coulomb'][before] == pytest.approx(info['released_energy'])
    assert e['deformation'][after] == pytest.approx(info['released_energy'], rel=1e-9)
    assert e['rest'][after] == e['kinetic'][after] == e['coulomb'][after] == 0
    assert e['transverse'][after] == e['transverse'][before]
    assert np.abs(h['energy']-h['energy'][0]).max() < 1e-6*h['energy'][0]
    assert np.abs(h['J']-h['J'][0]).max() < 1e-7
    assert np.abs(info['unreleased_angular']).max() < 1e-9


def test_released_pulse_refocuses_at_the_antipode_of_contact(model):
    h, info = model.run(model.initial(.2), 6.)
    s = info['final_state']
    c = np.asarray(info['contact_point'])/np.linalg.norm(info['contact_point'])
    pts = np.array([c, -c, np.cross(c, [1., 0, 0])/np.linalg.norm(np.cross(c, [1., 0, 0]))])
    u, ut = model.deformation(s, pts, info['contact_time']+np.pi)
    # The pulse starts as pure u_t; pi later it is pure displacement, focused at -c.
    assert np.abs(ut).max() < 1e-12*abs(u[1])
    assert abs(u[1]) > 5*max(abs(u[0]), abs(u[2]))
    _, back = model.deformation(s, pts, info['contact_time']+2*np.pi)
    _, start = model.deformation(s, pts, info['contact_time'])
    assert np.allclose(back, -start, atol=1e-12*abs(start[0]))      # back at c with the sign flipped


def test_front_crossings_for_a_resting_mouth():
    t = np.linspace(0, 7, 7001)
    phi = np.stack([np.full_like(t, 1.), np.full_like(t, -1.)], 1)
    start = np.array([.3, -.3])
    hits = [e for e in front_crossings(phi, t, start=start) if e['mouth'] == 0]
    times = sorted(round(e['time'], 3) for e in hits)
    # own ray .3 + t at .7; partner ray -.3 + t at 1.3; partner ray -.3 - t wraps round at
    # 2 pi - 1.3; own ray .3 - t at 2 pi - .7; own ray .3 + t again after a full turn, 2 pi + .7.
    assert times == [0.7, 1.3, round(2*np.pi-1.3, 3), round(2*np.pi-.7, 3), round(2*np.pi+.7, 3)]
    assert [e['source'] for e in sorted(hits, key=lambda e: e['time'])] == ['own', 'partner', 'partner', 'own', 'own']
    assert circle_angle(np.array([np.sin(1.), 0, np.cos(1.)]), *[np.array(e) for e in PLANE]) == pytest.approx(1.)


def test_bare_field_refocuses_at_the_antipode():
    th, t, e, _, _ = bare_timeline(t=np.array([0, np.pi, 2*np.pi]))
    assert th[np.argmax(e[1])] == pytest.approx(np.pi)
    assert np.allclose(e[2], e[0], atol=1e-9*e[0].max())                # exact return for omega = l + 1/2
    th, t, e, _, _ = bare_timeline(t=np.array([np.pi]), operator='massless')
    assert th[np.argmax(e[0])] == pytest.approx(np.pi)


def test_free_mouth_clocks_run_slow_by_gamma():
    m = PairBirth(Parameters(cutoff=8, charge=0.))
    h, info = m.run_history(m.initial(.4), 3.)
    assert np.allclose(h['tau'][-1], 3.*np.sqrt(1-.4**2), rtol=1e-9)
    assert np.abs(h['clock'][h['alive']]).max() < 1e-12              # equal speeds: no offset


def test_symmetric_cycles_keep_zero_offset_and_conserve_through_surgeries(model):
    h, info = model.run_history(model.initial(.25), 22., cycles=3)
    kinds = [e['kind'] for e in info['events']]
    assert kinds[:4] == ['contact', 'birth', 'contact', 'birth']
    c0, b1 = info['events'][0], info['events'][1]
    assert b1['time'] == pytest.approx(c0['time'] + np.pi)                     # re-birth at the refocus
    assert np.dot(b1['point'], c0['point']) == pytest.approx(-1.)               # at the antipode
    assert abs(b1['energy']) < 1e-12 and np.abs(b1['angular']).max() < 1e-12
    assert b1['correction'] < 1e-10                                            # the mirror needs no correction
    assert np.nanmax(np.abs(h['clock'])) == 0.
    assert info['energy_drift'] < 1e-6 and info['angular_drift'] < 1e-12


def test_asymmetric_history_offset_is_set_by_velocities_not_positions(model):
    h, info = model.run_history(model.initial(velocities=independent_start(.25)), 16., cycles=2)
    alive = h['alive']
    first = alive & (h['cycle'] == 0)
    g = np.sqrt(1 + np.sum(np.cross(h['L'], h['x'])**2, -1))                  # gamma per mouth (mass 1)
    dtau = np.diff(h['t'][first])[:, None]/g[first][:-1]
    integral = np.r_[0, np.cumsum(dtau[:, 0] - dtau[:, 1])]
    assert np.abs(h['clock'][first] - integral).max() < 2e-3                  # offset = integrated dilation
    birth = info['events'][1]
    after = np.flatnonzero(np.abs(h['t'] - birth['time']) < 1e-9)[-1]
    assert h['clock'][after] == 0. and np.allclose(h['tau'][after], 0.)       # clocks restart at re-birth
    assert abs(birth['energy']) < 1e-10 and np.abs(birth['angular']).max() < 1e-10
    assert info['angular_drift'] < 1e-7 and info['energy_drift'] < 1e-6
    ext, fo = h['external'], h['first_order']
    ok = alive & np.isfinite(ext)
    assert np.abs(ext[ok] - fo[ok]).max() < .05*np.abs(fo[ok]).max() + 1e-6  # small-offset approximation
    assert info['min_chronology_margin'] > 0
