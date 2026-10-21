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
    assert [e['amplitude'] for e in r['edges']] == pytest.approx([.37, -.37, .37])
    assert r['status'] == 'leg_budget'  # Pending continuation, never absorption.


def test_finite_aperture_and_spatial_reversal_transport_both_rim_coordinates():
    mouths = pair()
    rin = [.04, .06]
    for parity in [-1, 1]:
        r = AntipodalTransport(mouths, spatial_parity=parity).trace(
            1., mouths[0].point(1., rin), channel='throat', local=rin, max_legs=1)
        e = r['edges'][0]
        np.testing.assert_allclose(mouths[1].local(e['t1'], e['target']), [.04, parity*.06], atol=1e-14)
        assert e['amplitude'] == 1.  # A spatial map is not a scalar sign flip.


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
    assert r['edges'][1]['amplitude'] == -1.


def test_closed_surface_mouth_reflects_only_at_antipode():
    mouths = pair(); mouths[1].opened = False
    e = AntipodalTransport(mouths).trace(0., NORTH, max_legs=1)['edges'][0]
    assert e['event'] == 'antipodal_reflection' and e['reason'] == 'closed_mouth'
    assert e['t1'] == pytest.approx(np.pi)


def test_reversed_clock_keeps_advanced_edge_and_flags_cycle_unsolved():
    mouths = [inertial_mouth(0., horizon=18.),
              inertial_mouth(np.pi, horizon=18., orientation=-1, anchor=16.)]
    r = AntipodalTransport(mouths).trace(1., NORTH)
    assert r['status'] == 'cycle_requires_boundary_solution'
    assert any(e['channel'] == 'throat' and e['t1'] < e['t0'] for e in r['edges'])
    for e in r['edges']:
        if e['channel'] == 'throat':
            assert mouths[e['destination']].sigma(e['t1']) == pytest.approx(e['sigma1'])
    assert r['continuation']['edge'] == 0


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


def test_visual_klein_rim_identification_matches_transport_chart_map():
    from geometrodynamics.viz.klein_handles import build_surface, default_handles, is_orientable
    surface = build_surface(a=.24, handles=default_handles(1), n_ring=16,
                            n_sphere=180, ds_target=.06)
    assert not is_orientable(surface.tris)
    mouths = [inertial_mouth(0., radius=.24),
              inertial_mouth(np.pi, radius=.24, chart_side=-1)]
    # Follow each actual shared mesh column between rims, in the outer/inner
    # chart convention. This catches a painted handle with the wrong gluing.
    for ia, ib in zip(surface.tube_ids[0][0], surface.tube_ids[0][-1]):
        a = mouths[0].local(0., surface.points[ia])
        b = mouths[1].local(0., surface.points[ib])
        np.testing.assert_allclose(b, a*[1., -1.], atol=1e-14)
