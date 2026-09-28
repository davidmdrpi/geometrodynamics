"""Physical identities and independent region measurements, not trajectory tests."""
import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.waves import two_regions as tr


def test_disk_basis_measure_gradient_and_laplacian():
    xy, weights = tr.disk_grid(16, 40)
    basis = tr.DiskBasis(8)
    Y, grad = basis.evaluate(xy, gradient=True)
    np.testing.assert_allclose(Y.T @ (weights[:, None]*Y)/tr.VOLUME, np.eye(Y.shape[1]), atol=2e-14)
    np.testing.assert_allclose(weights.sum(), tr.VOLUME, atol=1e-13)
    x = np.array([[.24, .37], [.63, -.21]])
    h = 2e-5
    val, derivative = basis.evaluate(x, gradient=True)
    first = np.stack([(basis.evaluate(x+h*e)-basis.evaluate(x-h*e))/(2*h) for e in np.eye(2)], axis=1)
    np.testing.assert_allclose(first, derivative, atol=3e-7)
    second = np.stack([(basis.evaluate(x+h*e)-2*val+basis.evaluate(x-h*e))/h**2 for e in np.eye(2)], axis=1)
    cross = (basis.evaluate(x+h)-basis.evaluate(x+h*np.array([1,-1]))
             -basis.evaluate(x+h*np.array([-1,1]))+basis.evaluate(x-h))/(4*h*h)
    lap = ((1-x[:, 0]**2)[:, None]*second[:, 0]+(1-x[:, 1]**2)[:, None]*second[:, 1]
           -2*(x[:, 0]*x[:, 1])[:, None]*cross-3*np.einsum('pi,pij->pj', x, derivative))
    np.testing.assert_allclose(lap, -val*basis.lambdas, atol=4e-5)


def test_round_exact_solution_and_independent_curvature():
    data = tr.solve(8, 12, 32, (0, 0))
    np.testing.assert_allclose(data.coefficients[0], tr.PSI0, atol=1e-15)
    np.testing.assert_allclose(data.coefficients[1:], 0, atol=1e-15)
    assert data.iterations == 0
    errors = [abs(tr.coordinate_curvature(data, [.55, .8, .2], h)-6/tr.F0) for h in (1e-3, 5e-4)]
    assert errors[1] < errors[0]/3.8


def test_tracker_does_not_impose_equal_opposite_momenta():
    xy, weights = tr.disk_grid(32, 96)
    t = np.sqrt(1-np.sum(xy*xy, axis=1))
    points = np.concatenate([np.c_[xy, t*np.cos(p), t*np.sin(p)] for p in np.arange(4)*np.pi/2])
    weights = np.tile(weights/4, 4)
    X = np.c_[-points[:, 1], points[:, 0], np.zeros((len(points), 2))]
    w, _ = tr.window(points, tr.CENTERS[0])
    rows = tr.track_snapshot(points, weights, np.ones(len(points)), np.ones(len(points)), w[:, None]*X)
    assert rows[0]['momentum'][0] > .001
    np.testing.assert_allclose(rows[1]['momentum'], 0, atol=1e-15)
    assert np.arccos(np.dot(rows[0]['centroid'], rows[1]['centroid'])) > 1.19
    with pytest.raises(ValueError, match='overlap'):
        tr.track_snapshot(points, weights, np.ones(len(points)), np.ones(len(points)), w[:, None]*X,
                          centers=np.array([tr.CENTERS[0], -tr.CENTERS[0]]))


def test_source_antipodal_parity_and_independent_amplitudes():
    xy, _ = tr.disk_grid(12, 32)
    a, b = tr.sources(xy, (.02, .03)), tr.sources(-xy, (.02, .03))
    for key in ('q', 'f', 'S', 'U', 'lap'):
        np.testing.assert_allclose(a[key], b[key], atol=1e-14)
    np.testing.assert_allclose(a['dq'], -b['dq'], atol=1e-14)
    assert np.max(abs(tr.sources(xy, (.021, .03))['q']-tr.sources(xy, (.02, .031))['q'])) > .0005


def test_saved_coefficients_replay_constraint_and_independent_balance():
    path = Path(__file__).resolve().parents[1]/'experiments/closure_ledger/runs/20260927_two_regions/initial_data.json'
    record = json.loads(path.read_text())
    row = record['cases']['pair'][-1]
    data = tr.InitialData(row['degree'], (.02, .03), np.array(row['coefficients']),
                          row['iterations'], row['projected_residual'])
    d = tr.diagnostics(data)
    assert d['H_normalized_max'] < 1e-7
    assert max(abs(v['balance_error']) for v in d['ledger']) < 1e-5
    rates = [r['direct_rate_01'] for r in d['ledger']]
    assert rates[0] > 1e-5 and rates[1] < -1e-5
    assert abs(rates[0]+rates[1]) > 1e-6  # no imposed pairwise cancellation
    np.testing.assert_allclose(sum(rates), d['total_direct_rate_01'], atol=2e-15)
    np.testing.assert_allclose(d['total_direct_rate_01'], d['total_metric_work_01'], atol=2e-14)
    np.testing.assert_allclose([r['energy'] for r in d['regions']],
                               [r['energy'] for r in row['diagnostics']['regions']], rtol=1e-10)


def test_nonlinear_solve_cannot_be_replaced_by_adding_single_region_metrics():
    kwargs = dict(degree=12, radial=20, angular=64)
    pair = tr.solve(amplitudes=(.02,.03), **kwargs)
    a = tr.solve(amplitudes=(.02,0), **kwargs)
    b = tr.solve(amplitudes=(0,.03), **kwargs)
    difference = pair.coefficients-a.coefficients-b.coefficients
    difference[0] += tr.PSI0
    assert np.linalg.norm(difference) > 1e-6
