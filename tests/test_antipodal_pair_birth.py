"""Checks of the pair-birth history under the antipodal-only (BAM) interaction rule."""
import numpy as np
import pytest

from geometrodynamics.waves.antipodal_pair_birth import (AntipodalPairBirth, AntipodalState, Front,
                                                         from_local, local_coordinates)
from geometrodynamics.waves.pair_birth import Parameters, bare_timeline, birth_frames, birth_positions

COARSE = Parameters(cutoff=12, footprint=.2)
NORTH, SOUTH = np.array([0., 0., 1.]), np.array([0., 0., -1.])
PLANE = ((0., 0., 1.), (1., 0., 0.))


def resting_pair(xa, xb):
    """Two mouths at rest (no static field) with frames pointing to a reference point."""
    x = np.array([xa, xb], float)
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    ref = np.array([0., 1., 0.])
    return AntipodalState(0., x, np.zeros((2, 3)), np.zeros(2), birth_frames(x, ref), np.zeros(2))


def kinds(info):
    return [e['kind'] for e in info['events']]


def test_fronts_pass_mouths_untouched_and_reflect_at_empty_foci_with_the_gouy_phase():
    m = AntipodalPairBirth(COARSE)
    # The ring from the north pole sweeps over both mouths long before its focus at the south pole.
    s = resting_pair([np.sin(.6), 0, np.cos(.6)], [np.sin(1.2), 0, np.cos(1.2)])
    s.fronts = [Front(NORTH.copy(), 0., .01, np.zeros(3), 1., 'probe')]
    h, info = m.run(s, 2*np.pi + .1, dt=.05)
    assert kinds(info) == ['reflection', 'reflection']
    assert [e['time'] for e in info['events']] == pytest.approx([np.pi, 2*np.pi])
    front = info['final_state'].fronts[0]
    assert front.phase == pytest.approx(-1.)                                # -i per focus
    np.testing.assert_allclose(front.origin, NORTH, atol=1e-15)
    _, _, _, u, _ = bare_timeline(t=np.array([0., 2*np.pi]))                # the field's own 2 pi sign flip
    assert np.allclose(u[1], -u[0], atol=1e-9*np.abs(u[0]).max())
    assert info['energy_drift'] < 1e-15


def test_capture_at_the_focus_crosses_the_throat_on_the_proper_clocks():
    m = AntipodalPairBirth(COARSE)
    s = resting_pair(NORTH, SOUTH)
    s.fronts = [Front(NORTH.copy(), 0., .02, np.zeros(3), 1., 'probe')]
    h, info = m.run(s, np.pi + m.length + .5, dt=.01)
    cap, exit_ = [e for e in info['events'] if e['kind'] in ('capture', 'throat_exit')][:2]
    assert cap['mouth'] == 1 and cap['time'] == pytest.approx(np.pi) and cap['into_throat'] == pytest.approx(.02)
    assert exit_['mouth'] == 0 and exit_['time'] == pytest.approx(np.pi + m.length)     # sigma_out = sigma_in + L
    np.testing.assert_allclose(exit_['origin'], NORTH, atol=1e-12)
    assert info['energy_drift'] < 1e-14


def test_moving_mouths_exit_on_their_own_clock_histories():
    m = AntipodalPairBirth(COARSE)
    s = resting_pair(NORTH, SOUTH)
    vA, vB = .5, .2                                     # A heads +x from the north, B -x from the south:
    s.L = np.cross(s.x, np.array([[vA/np.sqrt(1-vA**2), 0, 0], [-vB/np.sqrt(1-vB**2), 0, 0]]))   # same sense
    B_pi = np.array([-np.sin(vB*np.pi), 0., -np.cos(vB*np.pi)])               # where B is at t = pi
    s.fronts = [Front(-B_pi, 0., .01, np.zeros(3), 1., 'probe')]            # so the front focuses on it
    h, info = m.run(s, 6., dt=.01)
    cap = next(e for e in info['events'] if e['kind'] == 'capture')
    ex = next(e for e in info['events'] if e['kind'] == 'throat_exit')
    assert cap['mouth'] == 1
    expected = (np.pi*np.sqrt(1-vB**2) + m.length)/np.sqrt(1-vA**2)        # tau_A(t) = tau_B(pi) + L
    assert ex['time'] == pytest.approx(expected, abs=1e-9)


def test_capture_transfers_angular_momentum_exactly():
    m = AntipodalPairBirth(COARSE)
    s = resting_pair(NORTH, SOUTH)
    J = np.array([.003, -.002, .004])                                         # |J| <= E
    s.fronts = [Front(NORTH.copy(), 0., .01, J, 1., 'twisted')]
    h, info = m.run(s, 5., dt=.01)
    assert 'capture' in kinds(info)
    st = info['final_state']
    assert np.abs(st.spin).max() > 0 and np.abs(st.L).max() > 0             # it pushed and spun the mouth
    assert info['angular_drift'] < 1e-13 and info['energy_drift'] < 1e-14


@pytest.mark.parametrize('handle, sign', [('klein', -1.), ('torus', 1.)])
def test_exit_point_follows_the_handle(handle, sign):
    m = AntipodalPairBirth(COARSE, handle=handle)
    s = resting_pair(NORTH, SOUTH)
    offset = np.array([.1, .07])                                            # aperture coordinates at B
    focus = from_local(s.x[1], s.frames[1], offset)
    s.fronts = [Front(-focus, 0., .01, np.zeros(3), 1., 'probe')]
    h, info = m.run(s, np.pi + m.length + .1, dt=.01)
    ex = next(e for e in info['events'] if e['kind'] == 'throat_exit')
    out = local_coordinates(s.x[0], s.frames[0], np.array(ex['origin']))
    np.testing.assert_allclose(out, [-offset[0], sign*offset[1]], atol=1e-12)   # Klein rotates, torus reflects


def test_strict_rule_removes_the_static_attraction_and_the_pair_coasts():
    m = AntipodalPairBirth(COARSE)
    s = m.initial(.3)
    v = np.linalg.norm(m.velocity(s), axis=1)
    assert v[0] < .3                                                         # the pair paid for its flash
    h, info = m.run(s, 12., dt=.01)
    c = next(e for e in info['events'] if e['kind'] == 'contact')
    assert c['time'] == pytest.approx((np.pi - .6)/v[0], abs=1e-9)          # free flight to the antipode
    np.testing.assert_allclose(c['point'], [0., 0., -1.], atol=1e-12)
    assert info['energy_drift'] < 1e-14 and info['angular_drift'] < 1e-15


def test_cycles_conserve_through_release_and_rebirth_in_both_static_modes():
    _, vel, _, _ = birth_positions(.3, 1., PLANE)
    for static in ('none', 'coulomb'):
        m = AntipodalPairBirth(COARSE, static=static)
        h, info = m.run(m.initial(velocities=np.array([.45*vel[0], .15*vel[1]])), 30., cycles=3, dt=.01)
        births = [e for e in info['events'] if e['kind'] == 'birth']
        assert len(births) == 2
        assert info['angular_drift'] < 1e-12
        assert info['energy_drift'] < (1e-14 if static == 'none' else 1e-6)


def test_advanced_exit_is_kept_while_its_focus_is_ahead():
    _, vel, _, _ = birth_positions(.3, 1., PLANE)
    m = AntipodalPairBirth(COARSE)
    h, info = m.run(m.initial(velocities=np.array([.85*vel[0], .15*vel[1]])), 4., dt=.005)
    cap = next(e for e in info['events'] if e['kind'] == 'capture')
    ex = next(e for e in info['events'] if e['kind'] == 'throat_exit')
    assert cap['exit'] == 'advanced' and ex['advanced']
    assert ex['time'] < cap['time'] < ex['time'] + np.pi                    # earlier exit, focus still ahead
    assert abs(cap['offset']) > m.length
    assert info['energy_drift'] < 1e-14


# ---- signed exchange: the Klein throat as a reversal ----------------------------------------------
@pytest.mark.parametrize('handle', ['klein', 'torus'])
def test_klein_throat_reverses_a_fronts_spin_and_sign_and_the_handle_recoils(handle):
    m = AntipodalPairBirth(COARSE, handle=handle, exchange='signed')
    s = resting_pair(NORTH, SOUTH)
    J = np.array([0., 0., .002])                                              # spin about the A-B axis
    # One empty focus behind it (phase -i): in quadrature it neither pushes nor pulls B at capture,
    # so the throat carries all of its angular momentum.
    s.fronts = [Front(NORTH.copy(), 0., .01, J.copy(), -1j, 'probe', 1, 1.)]
    h, info = m.run(s, np.pi + m.length + .05, dt=.01)
    cap = next(e for e in info['events'] if e['kind'] == 'capture')
    ex = next(e for e in info['events'] if e['kind'] == 'throat_exit')
    assert cap['coupling'] == 0. and cap['mouth'] == 1
    out = next(f for f in info['final_state'].fronts if f.kind.startswith('exit'))
    if handle == 'klein':
        np.testing.assert_allclose(out.J, -J, atol=1e-15)                   # angular momentum reversed
        assert out.sign == -1.                                              # compression -> rarefaction
        np.testing.assert_allclose(ex['recoil'], 2*J, atol=1e-15)           # the handle takes 2J
        assert info['final_state'].spin[0] == pytest.approx(.004)
    else:
        np.testing.assert_allclose(out.J, J, atol=1e-15)                    # preserved through a torus
        assert out.sign == 1.
        np.testing.assert_allclose(ex['recoil'], 0., atol=1e-15)
    assert info['angular_drift'] < 1e-15 and info['energy_drift'] < 1e-14


@pytest.mark.parametrize('sign, phase, share', [(-1., 1., 1.), (1., 1., -1.), (1., -1j, 0.), (1., -1., 1.)])
def test_signed_capture_pushes_like_signs_and_pulls_opposite_ones(sign, phase, share):
    m = AntipodalPairBirth(COARSE, exchange='signed')
    s = resting_pair(NORTH, SOUTH)                                            # B (index 1) carries -q
    J = np.array([.002, -.001, 0.])
    s.fronts = [Front(NORTH.copy(), 0., .01, J.copy(), phase, 'probe', 0, sign)]
    h, info = m.run(s, np.pi + .01, dt=.01)
    cap = next(e for e in info['events'] if e['kind'] == 'capture')
    assert cap['coupling'] == pytest.approx(share)                          # same sign pushes, opposite pulls,
    st = info['final_state']                                                # quadrature does neither; two foci flip
    np.testing.assert_allclose(st.L[1] + st.spin[1]*st.x[1], share*J, atol=1e-15)
    assert info['angular_drift'] < 1e-15


def test_signed_exchange_conserves_through_cycles_for_both_handles():
    _, vel, _, _ = birth_positions(.3, 1., PLANE)
    for handle in ('klein', 'torus'):
        m = AntipodalPairBirth(COARSE, handle=handle, exchange='signed')
        h, info = m.run(m.initial(velocities=np.array([.45*vel[0], .15*vel[1]])), 30., cycles=3, dt=.01)
        caps = [e for e in info['events'] if e['kind'] == 'capture']
        assert caps and any(e['coupling'] != 0 for e in caps)
        assert info['energy_drift'] < 1e-14 and info['angular_drift'] < 1e-13


@pytest.mark.parametrize('handle, turned', [('klein', (-1., 1., -1.)), ('torus', (-1., -1., 1.))])
def test_each_handle_turns_a_fronts_angular_momentum_by_pi_and_neither_reverses_it(handle, turned):
    # Mouths at the poles with e1 = y at both: J along x comes from motion along e1, J along y from
    # motion along e2, and J along z is spin about the mouths' axis.
    for axis in range(3):
        m = AntipodalPairBirth(COARSE, handle=handle, exchange='signed')
        s = resting_pair(NORTH, SOUTH)
        J = .002*np.eye(3)[axis]
        s.fronts = [Front(NORTH.copy(), 0., .01, J.copy(), -1j, 'probe', 1, 1.)]   # quadrature: throat holds J
        h, info = m.run(s, np.pi + m.length + .05, dt=.01)
        out = next(f for f in info['final_state'].fronts if f.kind.startswith('exit'))
        np.testing.assert_allclose(out.J, turned[axis]*J, atol=1e-15)
        assert info['angular_drift'] < 1e-15
    assert np.prod(turned) == 1.                                              # a rotation by pi, not -1


def test_a_pull_and_its_exit_are_capped_so_no_front_outruns_its_energy():
    m = AntipodalPairBirth(COARSE, exchange='signed')
    s = resting_pair(NORTH, SOUTH)
    J = np.array([.009, 0., 0.])                                              # |J| = 0.9 E
    s.fronts = [Front(NORTH.copy(), 0., .01, J.copy(), 1., 'probe', 0, 1.)]  # B carries -q: c = -1 would
    h, info = m.run(s, np.pi + m.length + .05, dt=.01)                        # leave the throat 1.8 E
    cap = next(e for e in info['events'] if e['kind'] == 'capture')
    ex = next(e for e in info['events'] if e['kind'] == 'throat_exit')
    assert cap['saturated'] and -1. < cap['coupling'] < 0.                    # still a pull, but capped
    assert (1 - cap['coupling'])*np.linalg.norm(J) == pytest.approx(cap['into_throat'], abs=1e-14)
    out = next(f for f in info['final_state'].fronts if f.kind.startswith('exit'))
    assert ex['saturated'] and np.linalg.norm(out.J) == pytest.approx(out.energy, abs=1e-14)
    assert info['energy_drift'] < 1e-14 and info['angular_drift'] < 1e-15


@pytest.mark.parametrize('v', [.04, 1/3, .8])
def test_first_flight_captures_only_near_birth_or_contact(v):
    # Coasting at u, source k's n-th focus is its birth point (n even) or that point's antipode
    # (n odd). Before contact the partner catches an odd return only if pi - 3a <= n u pi, the
    # last a/u of the flight; the source's own mouth an even return only if n u pi <= a, the first.
    m = AntipodalPairBirth(COARSE)
    s = m.initial(v)
    u, a = np.linalg.norm(m.velocity(s), axis=1)[0], m.par.mouth_angle
    contact = (np.pi - 2*a)/u
    h, info = m.run(s, contact - 1e-6, dt=.02)
    got = sorted({(round(e['time']/np.pi), e['mouth'], e['front']) for e in info['events']
                  if e['kind'] == 'capture' and e['front'].startswith('flash')})
    expected = []
    for k in range(2):
        for n in range(1, int(contact/np.pi) + 1):
            if n % 2 and np.pi - 3*a <= n*u*np.pi < np.pi - 2*a:
                expected.append((n, 1 - k, f'flash {"AB"[k]}'))
                break
            if not n % 2 and n*u*np.pi <= a:
                expected.append((n, k, f'flash {"AB"[k]}'))
                break
    assert got == sorted(expected)
    assert bool(got) == (v != 1/3)                      # own second return; none; partner first return
