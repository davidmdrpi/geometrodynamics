"""Checks of the pair-birth model (two Wheeler-flux mouths on S²) on a coarse mode set."""
import numpy as np
import pytest

from geometrodynamics.waves.moving_mouth_pairs import harmonics
from geometrodynamics.waves.pair_birth import (NeckChannel, PairBirth, Parameters, bare_timeline, birth_positions,
                                               circle_angle, front_crossings, rotate)

COARSE = Parameters(cutoff=16, footprint=.2)
PLANE = ((0., 0., 1.), (1., 0., 0.))


@pytest.fixture(scope='module')
def model():
    return PairBirth(COARSE)


@pytest.fixture(scope='module')
def necked():
    return PairBirth(COARSE, NeckChannel(port=100., nodes=33))


def lopsided_start(fast=.45, slow=.15):
    _, vel, _, _ = birth_positions(.3, 1., PLANE)
    return np.array([fast*vel[0], slow*vel[1]])


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


# ---- neck channel ----------------------------------------------------------------------------------
def first_peak(series, t):
    series = np.abs(series)
    j = np.argmax(series > .5*series.max())
    while j + 1 < len(series) and series[j+1] >= series[j]:
        j += 1
    return t[j]


def test_neck_signals_cross_in_L_plus_offset_one_way_and_L_minus_offset_the_other():
    from scipy import linalg
    m = PairBirth(Parameters(cutoff=4), NeckChannel(nodes=257, port=0.))
    L, sigma, n = m.neck_length, m.neck_nodes, 257
    assert L == pytest.approx(np.sin(.6))                                     # #338: 2S = sin 2a
    for A in (0., .2, -.2):
        k = -A/L                                                              # the history tilt
        step = linalg.expm(m.neck_generator(k, np.ones(2))[:2*n, :2*n]*1e-3)
        for d in (1, -1):                                                     # A -> B, then B -> A
            u = (1 - d*k)*(sigma - (.15 if d > 0 else .85)*L)
            F = np.exp(-u**2/(2*.02**2))
            z = np.r_[F, d*m._n_mass*(1 - d*k)*u/.02**2*F]                    # one-way pulse in neck time
            far = int(round((.85 if d > 0 else .15)*(n-1)))
            series = []
            for _ in range(700):
                z = step @ z
                series.append(z[far])
            arrival = first_peak(series, 1e-3*np.arange(1, 701))
            assert arrival == pytest.approx(.7*(L + d*A), abs=2e-3)


def test_neck_conserves_energy_and_angular_momentum_and_takes_up_off_plane_waves(necked):
    h, info = necked.run_history(necked.initial(velocities=independent_start()), 3.)
    e = h['energy_parts']
    assert info['energy_drift'] < 2e-6 and info['angular_drift'] < 1e-8
    assert e['neck'].max() > 1e-4 and e['port'].max() > 1e-6                  # a real test: the neck fills
    assert np.abs(h['ends'][:, 0]).max() > 1e-4                               # off-plane motion drives m = 0
    assert np.abs(h['spin']).max() > 1e-8 and np.abs(h['mu']).max() > 1e-6
    total = sum(e.values())
    assert np.allclose(total, h['energy'])


def test_coplanar_births_drive_only_the_sideways_part(necked):
    h, info = necked.run_history(necked.initial(velocities=lopsided_start()), 2.)
    assert np.abs(h['ends'][:, 0]).max() < 1e-15                              # B is odd under the plane mirror
    assert np.abs(h['ends'][:, 2]).max() > 1e-3                               # the out-of-plane m = 1 row
    assert np.abs(h['ends'][:, 1]).max() < 1e-15 and np.abs(h['spin']).max() < 1e-15
    assert np.abs(h['clock']).max() > .05 and np.abs(h['mu']).max() > 1e-5   # the offset acts back


def test_rebirth_starts_a_relaxed_neck_and_conserves(necked):
    h, info = necked.run_history(necked.initial(velocities=lopsided_start()), 17.5, cycles=2)
    birth = [e for e in info['events'] if e['kind'] == 'birth'][0]
    assert abs(birth['energy']) < 1e-10 and np.abs(birth['angular']).max() < 1e-10
    after = np.flatnonzero(np.abs(h['t'] - birth['time']) < 1e-9)[-1]
    e = h['energy_parts']
    assert e['port'][after] < 1e-2*e['transverse'][after]                    # ports start relaxed
    assert h['mu'][after] == 0. and h['clock'][after] == 0.
    assert info['energy_drift'] < 3e-6 and info['angular_drift'] < 1e-8


def test_handle_labels_match_the_glued_mesh_orientability():
    """In sphere-oriented rim frames #338's non-orientable (Klein) mesh maps rim to rim by a
    rotation; the neck's 'klein' map must be a rotation and its 'torus' map a reflection."""
    from geometrodynamics.viz.klein_handles import build_surface, default_handles, is_orientable
    surface = build_surface(a=.24, handles=default_handles(1), n_ring=16, n_sphere=180, ds_target=.06)
    assert not is_orientable(surface.tris)
    ends = [surface.points[surface.tube_ids[0][0]], surface.points[surface.tube_ids[0][-1]]]
    centres = [e.mean(0)/np.linalg.norm(e.mean(0)) for e in ends]
    frames = []
    for c in centres:
        e1 = np.cross(c, [0., 1., 0.]) if abs(c[1]) < .9 else np.cross(c, [1., 0., 0.])
        e1 /= np.linalg.norm(e1)
        frames.append((e1, np.cross(c, e1)))                                  # sphere-oriented: e1 x e2 = c
    local = [np.array([[p @ f[0], p @ f[1]] for p in pts]) for pts, f in zip(ends, frames)]
    rim_map, *_ = np.linalg.lstsq(local[0], local[1], rcond=None)
    assert np.linalg.det(rim_map) == pytest.approx(1., abs=1e-6)              # Klein: a rotation
    for handle, det in (('klein', 1.), ('torus', -1.)):
        m = PairBirth(Parameters(cutoff=4), NeckChannel(nodes=5, handle=handle))
        assert np.prod(m._handle[1]) == det


def test_symmetric_births_keep_zero_offset_with_either_handle():
    for handle in ('torus', 'klein'):
        m = PairBirth(COARSE, NeckChannel(port=100., nodes=33, handle=handle))
        h, info = m.run_history(m.initial(.25), 2.)
        assert np.abs(h['clock']).max() == 0. and np.abs(h['mu']).max() < 1e-15
        assert h['energy_parts']['neck'].max() > 1e-4


def test_decoupled_neck_is_the_bare_pair_and_clock_off_freezes_mu():
    start = independent_start()
    bare = PairBirth(COARSE)
    h0, _ = bare.run_history(bare.initial(velocities=start), 1.5)
    loose = PairBirth(COARSE, NeckChannel(port=0., nodes=33))
    h1, _ = loose.run_history(loose.initial(velocities=start), 1.5)
    assert np.abs(h1['x'] - h0['x']).max() < 1e-13 and np.abs(h1['clock'] - h0['clock']).max() < 1e-13
    off = PairBirth(COARSE, NeckChannel(port=100., nodes=33, clock='off'))
    h2, info = off.run_history(off.initial(velocities=start), 1.5)
    assert np.all(h2['mu'] == 0.) and info['energy_drift'] < 2e-6
    assert np.abs(h2['x'] - h0['x']).max() > 1e-7                            # the neck still acts, without clocks


def test_clock_and_port_forces_match_finite_differences(necked):
    s = necked.initial(velocities=independent_start())
    for _ in range(250):
        necked.step(s, .004)
    m, eps = necked, 1e-6

    def energy(state):
        return np.sum(np.sqrt(m.masses(state)**2 + np.sum(state.L**2, 1)))

    # clock: dH/dA at fixed everything else, against the tilt force that moves mu
    up, down = s.copy(), s.copy()
    up.tau[0] += eps
    down.tau[0] -= eps
    analytic = -m._tilt_force(m.tilt(s), m.rates(s), s.eta, s.pi)/m.neck_length
    assert (energy(up) - energy(down))/(2*eps) == pytest.approx(analytic, rel=1e-5)
    assert abs(analytic) > 1e-6
    # ports: force on each mouth, torque on its spin and pull on the field, against -grad V
    ends = s.eta[:, [0, -1]]
    targets = m.rim_targets(s)
    pull0 = m.neck.port*(ends[0] - targets[0])
    pull1 = .5*m.neck.port*(np.einsum('ck,kcj->kj', ends[1:], m.port_frames(s)) - targets[2])
    dQ, dL, ds = m._port_reaction(s, pull0, pull1, targets)
    force = m.forces(s)['port']

    def V(state):
        return m.port_energies(state).sum()

    for j in range(2):
        for axis in np.eye(3):
            axis = axis - (axis @ s.x[j])*s.x[j]
            axis /= np.linalg.norm(axis)
            pair = []
            for sign in (1, -1):
                t = s.copy()
                t.x[j] = rotate(s.x[j], axis, sign*eps)
                t.frames[j] = rotate(s.frames[j], axis, sign*eps)
                pair.append(V(t))
            assert force[j] @ np.cross(axis, s.x[j]) == pytest.approx(-(pair[0] - pair[1])/(2*eps), rel=1e-5, abs=1e-12)
        pair = []
        for sign in (1, -1):
            t = s.copy()
            t.frames[j] = rotate(s.frames[j], s.x[j], sign*eps)
            pair.append(V(t))
        assert ds[j] == pytest.approx(-(pair[0] - pair[1])/(2*eps), rel=1e-5, abs=1e-12)
    direction = harmonics(m.ell, m.order, np.array([[.3, -.5, .81]])/np.linalg.norm([.3, -.5, .81]))[0].conj()
    pair = []
    for sign in (1, -1):
        t = s.copy()
        t.P[0] = s.P[0] + sign*eps*direction
        pair.append(V(t))
    assert np.vdot(dQ, direction).real == pytest.approx((pair[0] - pair[1])/(2*eps), rel=1e-5)
