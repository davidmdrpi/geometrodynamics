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
