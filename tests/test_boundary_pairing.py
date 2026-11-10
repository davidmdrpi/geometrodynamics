"""Independent oracles and failure controls for boundary-pairing qualification."""
import json

import numpy as np
import pytest
from numpy.testing import assert_allclose

from geometrodynamics.waves.boundary_pairing import (
    conditional_spinor_checks, context_gram_counterexample, direct_waist_load,
    interval_dtn, join_dtn, output_gram, polarization_matrix,
    scattering_from_dtn, third_order_response, transfer_to_dtn,
)
from geometrodynamics.waves.maxwell_throat_response import Bridge, boundary_dtn
from experiments.closure_ledger import boundary_pairing_baseline_probe as probe


def free_transfer(k, length):
    c, s = np.cos(k*length), np.sin(k*length)
    return np.array([[c, s/k], [-k*s, c]])


def test_free_interval_boundary_signs_and_coherent_transmission():
    k, length = 2.3, .7
    d = transfer_to_dtn(free_transfer(k, length))
    expected = k / np.sin(k*length) * np.array([[np.cos(k*length), -1], [-1, np.cos(k*length)]])
    assert_allclose(d, expected, atol=1e-14)
    assert_allclose(scattering_from_dtn(d, k),
                    [[0, np.exp(1j*k*length)], [np.exp(1j*k*length), 0]], atol=1e-14)
    # Actual on-shell action is -u.D.u/2 with these outward-flux conventions.
    u = np.array([.2, -.7])
    direction = np.array([.3, .1])
    h = 1e-5
    action = lambda v: -.5 * (v @ d @ v).real
    derivative = (action(u+h*direction)-action(u-h*direction))/(2*h)
    assert_allclose(derivative, -(d@u)@direction, atol=2e-11)


def test_transfer_retains_nonunit_determinant_as_a_reciprocity_failure():
    d = transfer_to_dtn([[1, 1], [0, 2]])
    assert_allclose(d, [[1, -1], [-2, 2]])
    assert not np.allclose(d, d.T)


def test_three_free_segments_reproduce_analytic_whole_map():
    k = 1.7
    lengths = [.17, .29, .41]
    pieces = [transfer_to_dtn(free_transfer(k, length)) for length in lengths]
    expected = transfer_to_dtn(free_transfer(k, sum(lengths)))
    assert_allclose(join_dtn(pieces), expected, atol=2e-14)
    assert_allclose(join_dtn([join_dtn(pieces[:2]), pieces[2]]), expected, atol=2e-14)
    # Wrong interface orientation must not be an invisible control.
    assert np.linalg.norm(join_dtn([pieces[0], -pieces[1], pieces[2]])-expected) > .1


@pytest.mark.parametrize('bridge', [Bridge(), Bridge(lapse_depth=.4),
    Bridge(lapse_depth=.3, waist=.04, stretch=1.1, charge_fraction=.5)])
def test_noncentral_cut_on_inherited_and_deformed_bridge(bridge):
    x = bridge.limit
    cuts = [-x, -.33*x, .18*x, x]
    pieces = [interval_dtn(bridge, 23., a, b) for a, b in zip(cuts[:-1], cuts[1:])]
    whole = boundary_dtn(bridge, 23.)
    assert_allclose(join_dtn(pieces), whole, rtol=2e-7, atol=2e-7)


def test_singular_maps_fail_without_regularizing_or_clipping():
    with pytest.raises(ValueError, match='pole'):
        transfer_to_dtn(np.eye(2))
    with pytest.raises(ValueError, match='interface'):
        join_dtn([np.array([[1., 1.], [1., 0.]]), np.array([[0., 1.], [1., 1.]])])
    with pytest.raises(ValueError, match='interface'):
        join_dtn([np.array([[1., 1.], [1., 1e-15]]), np.array([[0., 1.], [1., 1.]])])


@pytest.mark.parametrize('omega,left,right,ell', [(0, -.1, .1, 1), (2, .1, -.1, 1),
    (2, -100, .1, 1), (2, -.1, .1, 0), (2, -.1, .1, 1.5), (np.nan, -.1, .1, 1)])
def test_bad_physical_inputs_rejected(omega, left, right, ell):
    with pytest.raises(ValueError):
        interval_dtn(Bridge(), omega, left, right, ell)


def test_polarization_recovers_imaginary_cross_term_with_correct_sign():
    matrix = np.array([[2, .3+.7j], [.3-.7j, 3]])
    assert_allclose(polarization_matrix(lambda v: np.vdot(v, matrix@v).real), matrix, atol=1e-15)
    # Direct load obeys polarization even for incoincident complex inputs.
    bridge, omega = Bridge(lapse_depth=.4), 20.
    load = lambda v: direct_waist_load(bridge, omega, v, order=128)
    w = polarization_matrix(load)
    v = np.array([.4+.3j, -.2+.8j])
    assert_allclose(load(v), np.vdot(v, w@v).real, rtol=1e-11, atol=1e-12)


def test_fixed_background_gram_i3_and_relative_phase_interference():
    s = np.array([[1, 1j], [1j, 1]])/np.sqrt(2)
    inputs = np.array([[1, 0], [0, 1], [.3+.2j, -.1]])
    gram = output_gram(s, inputs)
    assert np.linalg.eigvalsh(gram).min() > -1e-14
    assert abs(third_order_response(gram)) < 1e-14
    constructive = output_gram(s, [[1, 0], [0, 1j]])
    destructive = output_gram(s, [[1, 0], [0, -1j]])
    assert_allclose(constructive.sum(), 2., atol=1e-14)
    assert_allclose(destructive.sum(), 0., atol=1e-14)
    # A quartic detector response does not inherit the fixed quadratic I3 identity.
    q = lambda v: abs(s[1]@v)**4
    a, b, c = inputs
    quartic_i3 = q(a+b+c)-q(a+b)-q(a+c)-q(b+c)+q(a)+q(b)+q(c)
    assert abs(quartic_i3) > .01


def test_positive_contexts_do_not_imply_common_quantum_bound():
    control = context_gram_counterexample()
    assert control['chsh'] == 4
    assert control['common_joint_functional'] == 'NOT_SUPPLIED'
    for row in control['contexts']:
        assert row['min_eigenvalue'] >= 0
        assert row['normalization_error'] == row['marginal_error'] == 0


def test_conditional_spinor_identity_does_not_select_efficiency_or_pin_sign():
    result = conditional_spinor_checks(samples=128)
    assert max(result['residuals'].values()) < 1e-12
    assert result['unequal_efficiency_gap'] > .1
    assert result['scope'] == 'CONDITIONAL_SPINOR_IDENTITY'


@pytest.mark.parametrize('bad', [None, float('nan'), float('inf'), -1., True, 1.])
def test_invalid_gate_metric_cannot_pass(bad):
    values = {key: 0. for key in probe.LIMITS}
    values['dtn_joining'] = bad
    result = probe.evaluate_gates(values)
    assert result['verdict'] == 'UNRESOLVED_NUMERICALLY'
    assert not result['checks']['dtn_joining']['passed']
    json.dumps(result, allow_nan=False)


def test_every_missing_gate_fails():
    for key in probe.LIMITS:
        result = probe.evaluate_gates({k: 0. for k in probe.LIMITS if k != key})
        assert result['verdict'] == 'UNRESOLVED_NUMERICALLY'


def test_failed_cli_writes_failure_and_exits_nonzero(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(probe, 'run_qualification', lambda: probe.evaluate_gates({}))
    output = tmp_path/'nested'/'failed.json'
    assert probe.main(['--output', str(output)]) == 1
    assert json.loads(output.read_text())['verdict'] == 'UNRESOLVED_NUMERICALLY'
    assert 'UNRESOLVED_NUMERICALLY' in capsys.readouterr().out


def test_early_solver_exception_cannot_leave_zero_maxima_passing(monkeypatch):
    def fail(*args, **kwargs):
        raise ValueError('injected interface singularity')
    monkeypatch.setattr(probe, 'interval_dtn', fail)
    result = probe.run_qualification()
    assert result['verdict'] == 'UNRESOLVED_NUMERICALLY'
    assert result['errors'] and not any(c['passed'] for c in result['checks'].values())


def test_probe_qualifies_baseline_only():
    result = probe.run_qualification()
    assert result['verdict'] == 'BASELINE_QUALIFIED_DYNAMIC_EXPERIMENT_PENDING', result
    assert result['dynamic_mouth_response'] == 'NOT_IMPLEMENTED'
    assert result['physical_spinor_matching'] == 'NOT_DERIVED'
    assert result['detector_frequency_law'] == 'NOT_DERIVED'
    assert result['common_all_setting_history_kernel'] == 'NOT_DERIVED'
    assert len(result['cases']) == 12 and len(result['loads']) == 2
