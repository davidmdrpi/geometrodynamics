import numpy as np
import pytest

from geometrodynamics.waves.history_throat import (
    AntipodalTransport, HistoryBoundary, MouthHistory, inertial_mouth,
)

NORTH = np.array([0., 0., 1.])


def pair(**kwargs):
    return [inertial_mouth(0., **kwargs), inertial_mouth(np.pi, **kwargs)]


def test_surface_crossings_cannot_capture_before_antipodal_focus():
    # The expanding ring crosses BOTH stationary mouths before reaching pi.
    mouths = [inertial_mouth(.6), inertial_mouth(1.2)]
    route = AntipodalTransport(mouths).trace(0., NORTH, max_legs=2)
    assert [e['event'] for e in route['edges']] == ['antipodal_reflection']*2
    assert [e['t1'] for e in route['edges']] == pytest.approx([np.pi, 2*np.pi])
    for e in route['edges']:
        np.testing.assert_allclose(e['target'], -np.array(e['origin']), atol=1e-15)


def test_hit_uses_arrival_worldtube_not_emission_position():
    mouths = [inertial_mouth(0.), inertial_mouth(np.pi, .1, radius=.05)]
    assert AntipodalTransport(mouths).trace(0., NORTH, max_legs=1)['edges'][0]['event'] == 'antipodal_reflection'
    mouths[1] = inertial_mouth(np.pi-.1*np.pi, .1, radius=.05)
    assert AntipodalTransport(mouths).trace(0., NORTH, max_legs=1)['edges'][0]['event'] == 'mouth_entry'


def test_antipodal_entry_is_followed_by_throat_then_new_surface_leg():
    r = AntipodalTransport(pair(), length=.4).trace(0., NORTH, max_legs=3)
    a, b, c = r['edges']
    assert [x['channel'] for x in r['edges']] == ['surface', 'throat', 'surface']
    assert a['mouth'] == 1 and b['source'] == 1 and b['destination'] == 0
    assert b['t0'] == pytest.approx(np.pi) and b['t1'] == pytest.approx(np.pi+.4)
    np.testing.assert_allclose(c['origin'], NORTH, atol=1e-14)
    assert c['t1']-c['t0'] == pytest.approx(np.pi)


def test_empty_focus_reflects_front_without_repeated_source_injection():
    r = AntipodalTransport([inertial_mouth(.8), inertial_mouth(1.2)]).trace(0., NORTH, amplitude=.37, max_legs=3)
    assert len(r['edges']) == 3
    assert [complex(*e['amplitude']) for e in r['edges']] == pytest.approx([.37, -.37j, -.37])
    assert r['status'] == 'leg_budget'  # Pending continuation, never absorption.


def test_finite_aperture_and_spatial_reversal_transport_both_rim_coordinates():
    mouths = pair()
    rin = [.04, .06]
    for handle, expected in [('klein', [-.04, -.06]), ('torus', [.04, -.06])]:
        r = AntipodalTransport(mouths, handle=handle).trace(
            1., mouths[0].point(1., rin), channel='throat', local=rin, max_legs=1)
        e = r['edges'][0]
        np.testing.assert_allclose(mouths[1].local(e['t1'], e['target']), expected, atol=1e-14)
        assert complex(*e['amplitude']) == 1.  # A spatial map is not a scalar sign flip.


def test_actual_proper_clock_maps_control_throat_arrival():
    mouths = [inertial_mouth(0., .4), inertial_mouth(np.pi, -.2)]
    m = AntipodalTransport(mouths, length=.2)
    r = m.trace(.5, mouths[0].point(.5, [0., 0.]), channel='throat', max_legs=1)
    expected = (.5*np.sqrt(1-.4**2)+.2)/np.sqrt(1-.2**2)
    assert r['edges'][0]['t1'] == pytest.approx(expected, abs=1e-12)


def test_closed_destination_reflects_inside_throat_not_across_surface():
    mouths = pair(); mouths[1].opened = False
    r = AntipodalTransport(mouths, length=.4).trace(1., NORTH, channel='throat', max_legs=3)
    assert [e['event'] for e in r['edges'][:2]] == ['throat_reflection', 'mouth_exit']
    assert [e['channel'] for e in r['edges']] == ['throat', 'throat', 'surface']
    assert r['edges'][1]['t1'] == pytest.approx(1.8)
    assert complex(*r['edges'][1]['amplitude']) == -1.


def test_closed_surface_mouth_reflects_only_at_antipode():
    mouths = pair(); mouths[1].opened = False
    e = AntipodalTransport(mouths).trace(0., NORTH, max_legs=1)['edges'][0]
    assert e['event'] == 'antipodal_reflection' and e['reason'] == 'closed_mouth'
    assert e['t1'] == pytest.approx(np.pi)


def test_reversed_clock_keeps_advanced_edge_and_flags_cycle_unsolved():
    mouths = [inertial_mouth(0., horizon=18.),
              inertial_mouth(np.pi, horizon=18., orientation=-1, anchor=16.)]
    r = AntipodalTransport(mouths, time_policy='channel_forward').trace(1., NORTH, max_legs=4)
    assert r['status'] == 'cycle_requires_boundary_solution'
    assert any(e['channel'] == 'throat' and e['t1'] < e['t0'] for e in r['edges'])
    for e in r['edges']:
        if e['channel'] == 'throat':
            assert mouths[e['destination']].sigma(e['t1']) == pytest.approx(e['sigma1'])
    assert r['continuation']['edge'] == 0
    assert complex(*r['continuation']['loop_gain']) == -1.
    assert r['continuation']['homogeneous_solution'] == 'zero_only'
    assert r['continuation']['seed_consistency_residual'] == pytest.approx(2.)


def test_history_boundaries_are_unresolved_not_extrapolated_or_absorbed():
    r = AntipodalTransport(pair(horizon=2.)).trace(0., NORTH)
    assert r['status'] == 'history_boundary' and r['edges'][0]['mouth'] is None
    with pytest.raises(HistoryBoundary):
        pair(horizon=2.)[0].external(3.)
    mouths = pair(horizon=2.)
    r = AntipodalTransport(mouths, length=1.).trace(1.5, NORTH, channel='throat')
    assert r['edges'][0]['t1'] is None
    assert r['status'] == 'history_boundary'


def test_overlapping_apertures_fail_instead_of_picking_a_mouth():
    with pytest.raises(ValueError, match='Overlapping'):
        AntipodalTransport([inertial_mouth(np.pi), inertial_mouth(np.pi)]).trace(0., NORTH)


@pytest.mark.parametrize('bad', [np.nan, np.inf, -1.])
def test_invalid_neck_length_is_rejected(bad):
    with pytest.raises(ValueError):
        AntipodalTransport(pair(), length=bad)


@pytest.mark.parametrize('handle', ['klein', 'torus'])
@pytest.mark.parametrize('sides', [(1, 1), (1, -1), (-1, 1), (-1, -1)])
def test_visual_rim_identification_is_independent_of_chart_convention(handle, sides):
    from geometrodynamics.viz.klein_handles import build_surface, default_handles, is_orientable
    surface = build_surface(a=.24, handles=default_handles(1, klein=handle == 'klein'), n_ring=16,
                            n_sphere=180, ds_target=.06)
    assert is_orientable(surface.tris) == (handle == 'torus')
    mouths = [inertial_mouth(0., radius=.24, chart_side=sides[0]),
              inertial_mouth(np.pi, radius=.24, chart_side=sides[1])]
    transport = AntipodalTransport(mouths, handle=handle)
    for ia, ib in zip(surface.tube_ids[0][0], surface.tube_ids[0][-1]):
        a = mouths[0].local(0., surface.points[ia])
        b = mouths[1].local(0., surface.points[ib])
        np.testing.assert_allclose(b, transport.rim_map(0, 1)@a, atol=1e-14)
        np.testing.assert_allclose(a, transport.rim_map(1, 0)@b, atol=1e-14)


def test_caustic_phase_matches_integrated_free_conformal_field():
    from geometrodynamics.waves.legacy_history_throat import HistoryThroat, Parameters
    m = HistoryThroat(Parameters(lmax=3, nodes=3, coupling=0., charge=0.,
        horizon=2*np.pi, neck_horizon=2*np.pi))
    field = m.solve()
    initial = m.initial()[:m.n]
    route = AntipodalTransport([inertial_mouth(.8), inertial_mouth(1.2)]).trace(
        0., NORTH, amplitude=.37, max_legs=2)
    for k, e in enumerate(route['edges'], 1):
        y = field['exterior'].sol(k*np.pi)
        z = y[:m.n]+1j*y[m.n:2*m.n]/m.omega
        expected = complex(*e['outgoing_amplitude'])/.37 * (-1.)**(k*m.ell)*initial
        np.testing.assert_allclose(z, expected, atol=1e-7)
    assert complex(*route['edges'][0]['outgoing_amplitude']) == -.37j
    assert complex(*route['edges'][1]['outgoing_amplitude']) == -.37


@pytest.mark.parametrize('end_condition,phase', [('dirichlet', -1.), ('neumann', 1.)])
def test_outlet_phase_is_independent_of_surface_caustic(end_condition, phase):
    mouths = pair(); mouths[1].opened = False; mouths[1].end_condition = end_condition
    e = AntipodalTransport(mouths).trace(1., NORTH, channel='throat',
                                      amplitude=.3+.4j, max_legs=1)['edges'][0]
    assert complex(*e['outgoing_amplitude']) == pytest.approx(phase*(.3+.4j))
    assert 'caustic_factor' not in e


def test_captured_focus_gets_quadrature_before_throat():
    e = AntipodalTransport(pair()).trace(0., NORTH, max_legs=2)['edges']
    assert e[0]['event'] == 'mouth_entry'
    assert complex(*e[0]['focus_amplitude']) == -1j
    assert complex(*e[1]['amplitude']) == -1j


def test_transported_time_orientation_reverses_exit_surface_direction():
    mouths = [inertial_mouth(0.), inertial_mouth(np.pi, orientation=-1, anchor=16.)]
    r = AntipodalTransport(mouths).trace(1., NORTH)
    a, b, c = r['edges'][:3]
    assert a['time_direction'] == 1
    assert b['time_direction'] == -1 and b['sigma1'] < b['sigma0']
    assert c['time_direction'] == -1 and c['t1'] < c['t0']
    assert complex(*c['caustic_factor']) == 1j
    assert r['status'] != 'cycle_requires_boundary_solution'
    for prev, nxt in zip(r['edges'], r['edges'][1:]):
        if prev['event'] == 'mouth_entry':
            assert nxt['time_direction'] == mouths[prev['mouth']].orientation*prev['time_direction']
        if prev['event'] == 'mouth_exit':
            assert nxt['time_direction'] == mouths[prev['destination']].orientation*prev['time_direction']


def test_default_handle_is_klein_even_with_both_sphere_oriented_charts():
    m = AntipodalTransport(pair())
    assert m.handle == 'klein'
    np.testing.assert_allclose(m.rim_map(0, 1), -np.eye(2))



def test_near_recurrence_reports_residual_and_honors_explicit_tolerance():
    mouths = [inertial_mouth(0., 1e-8),
              inertial_mouth(np.pi, orientation=-1, anchor=16.)]
    seed = mouths[0].point(1., [0., 0.])
    strict = AntipodalTransport(mouths, time_policy='channel_forward', recurrence_tolerance=1e-12)
    loose = AntipodalTransport(mouths, time_policy='channel_forward', recurrence_tolerance=1e-6)
    a = strict.trace(1., seed, max_legs=8)
    b = loose.trace(1., seed, max_legs=8)
    assert a['status'] == 'leg_budget'
    assert 1e-8 < a['closest_return']['position_difference'] < 1e-6
    assert b['status'] == 'cycle_requires_boundary_solution'
    assert b['continuation']['closure']['position_difference'] > 1e-8
    # A tolerance-based recurrence is reported with residual, not exact closure.
    assert b['continuation']['closure']['tolerance'] == 1e-6


def test_manifest_verifies_from_another_checkout_root(tmp_path):
    import hashlib
    import json
    from pathlib import Path
    from experiments.closure_ledger.antipodal_throat_replay import REPO, ROOT
    manifest = json.loads((ROOT/'sha256.json').read_text())
    assert 'experiments/closure_ledger/antipodal_throat_replay.py' in manifest
    for key, expected in manifest.items():
        assert not Path(key).is_absolute() and '..' not in Path(key).parts
        target = tmp_path/key
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((REPO/key).read_bytes())
        assert hashlib.sha256(target.read_bytes()).hexdigest() == expected



def test_real_pulse_quadrature_and_return_match_repository_bare_timeline():
    from geometrodynamics.waves.pair_birth import bare_timeline
    _, _, _, u, ut = bare_timeline(cutoff=12, t=[0., np.pi, 2*np.pi])
    assert np.max(abs(u[1])) < 1e-12
    assert np.max(abs(ut[1])) > 1.
    np.testing.assert_allclose(u[2], -u[0], atol=1e-12)



def test_precaustic_capture_changes_loop_gain_but_not_empty_focus_phase():
    mouths = [inertial_mouth(0.), inertial_mouth(np.pi, orientation=-1, anchor=16.)]
    r = AntipodalTransport(mouths, time_policy='channel_forward', capture_phase='before_caustic').trace(1., NORTH)
    assert r['status'] == 'cycle_requires_boundary_solution'
    assert complex(*r['continuation']['loop_gain']) == 1.
    assert r['continuation']['homogeneous_solution'] == 'any_amplitude'
    assert r['continuation']['seed_consistency_residual'] == 0.
    r = AntipodalTransport([inertial_mouth(.8), inertial_mouth(1.2)], capture_phase='before_caustic').trace(0., NORTH, max_legs=2)
    assert complex(*r['edges'][0]['outgoing_amplitude']) == -1j
    assert complex(*r['edges'][1]['outgoing_amplitude']) == -1.
