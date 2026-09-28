"""Canonical normalization, nonlinear ledger and explicit receiver limitations."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.waves import action_selection as a
from experiments.closure_ledger.action_selection_probe import decode_states, summarize, validate_schedule

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT/'experiments/closure_ledger/runs/20260928_action_selection'


def test_canonical_one_form_matches_kinetic_legendre_derivative():
    y = a.initial_loop(.08)[8]
    A, Ap, q, qp, M, L = a.dynamics.unpack(y)
    H = A*A-q @ q/6
    Mp = 2*M @ L
    variation = a.U @ M+M @ a.U
    def kinetic(delta):
        rate = np.linalg.solve(M, Mp+delta*variation)/2
        return H*np.trace(rate @ rate)/2
    h = 1e-6
    measured = (kinetic(h)-kinetic(-h))/(2*h)
    Pi = H/2*L @ np.linalg.inv(M)
    np.testing.assert_allclose(measured, np.trace(Pi @ variation), rtol=1e-8, atol=1e-11)
    constraints = a.dynamics.constraints(y)
    assert max(constraints['normalized']) < 1e-13
    assert np.linalg.norm(constraints['matter_momentum']) > 1e-6


def test_exact_initial_action_scaling_orientation_and_phase_parameterization():
    for eps in (.01,.02,.08):
        forward = a.circulation(a.initial_loop(eps,32)[:,None,:])
        np.testing.assert_allclose(forward['total'][0], a.initial_prediction(eps), rtol=1e-12)
        np.testing.assert_allclose(forward['shape'], forward['total'], atol=1e-14)
        shifted = a.circulation(a.initial_loop(eps,32,phase_offset=.31)[:,None,:])
        reverse = a.circulation(a.initial_loop(eps,32)[::-1,None,:])
        np.testing.assert_allclose(shifted['total'], forward['total'], rtol=1e-12)
        np.testing.assert_allclose(reverse['total'], -forward['total'], rtol=1e-12)
    assert a.initial_prediction(.02)/a.initial_prediction(.01) == 4


def test_periodic_derivative_and_invalid_inputs():
    theta = 2*np.pi*np.arange(32)/32
    np.testing.assert_allclose(a.loop_derivative(np.cos(3*theta)), -3*np.sin(3*theta), atol=4e-14)
    for eps,n in [(-.1,32),(.1,3),(.1,31),(float('nan'),32)]:
        with pytest.raises(ValueError):
            a.initial_loop(eps,n)
    with pytest.raises(ValueError):
        a.circulation(np.zeros((16,29)))


def test_archived_states_replay_full_ledger_and_readiness():
    raw = decode_states(RUN/'states.json.gz.b64')
    saved = json.loads((RUN/'action.json').read_text())
    replay = summarize(raw)
    assert replay['numerical_gates'] == saved['numerical_gates']
    assert replay['numerical_verdict'] == saved['numerical_verdict']
    assert replay['receiver_selection']['verdict'] == 'NOT_READY_FOR_RECEIVER_ACTION_SELECTION'
    assert len(replay['receiver_selection']['missing']) == 6
    assert 'selection mechanism' in replay['receiver_selection']['missing'][0]
    assert replay['numerical_verdict'] == 'INTEGRAL_INVARIANT_IMPLEMENTATION_CHECK'
    assert replay['numerical_check_passed'] is True
    for row, original in zip(replay['cases'],saved['cases']):
        for key in a.SECTORS:
            values = row['quadratures'][-1]['action'][key]
            np.testing.assert_allclose(values,original['quadratures'][-1]['action'][key],rtol=1e-9,atol=1e-12)
        fine = row['quadratures'][-1]['action']
        np.testing.assert_allclose(np.array(fine['scale'])+fine['scalar']+fine['shape'],fine['total'],atol=1e-14)
    for path, digest in saved['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    assert hashlib.sha256((RUN/'states.json.gz.b64').read_bytes()).hexdigest() == saved['archive_sha256']


def test_fresh_nonlinear_evolution_and_projection_control():
    states = a.evolve_loop(.04,16)
    action = a.circulation(states)
    prediction = a.initial_prediction(.04)
    assert np.max(abs(action['total']-prediction))/prediction < 1e-6
    assert a.diagnostics(states)['sampled_constraint_max'] < 1e-9
    # Sector-only circulation is not the full conserved canonical action.
    assert np.max(abs(action['shape']-prediction)) > 1e-8*prediction
    assert a.receiver_readiness()['verdict'] != 'ACTION_SELECTION_ESTABLISHED'


def test_replay_rejects_relabelled_preparation_or_duration():
    raw = decode_states(RUN/'states.json.gz.b64')
    raw['times'][-1] = 3.
    with pytest.raises(ValueError, match='schedule'):
        validate_schedule(raw)
    raw['times'][-1] = 2.
    raw['cases'][1]['states'][0][0][0] += .01
    with pytest.raises(ValueError, match='preparation'):
        validate_schedule(raw)


def test_review_relabel_preserves_frozen_results_and_rejects_damaged_flow():
    saved = json.loads((RUN/'action.json').read_text())
    initial = json.loads((RUN/'action_initial.json').read_text())
    assert saved['registered_numerical_verdict'] == initial['numerical_verdict']
    assert saved['numerical_gates'] == initial['numerical_gates']
    assert saved['archive_sha256'] == initial['archive_sha256']
    assert saved['cases'] == initial['cases']
    raw = decode_states(RUN/'states.json.gz.b64')
    raw['cases'][-1]['states'][0][-1][1] *= 1.1
    damaged = summarize(raw)
    assert damaged['numerical_check_passed'] is False
    assert damaged['numerical_verdict'] == 'REGISTERED_NUMERICAL_FAILURE'
    assert damaged['receiver_selection']['verdict'] == 'NOT_READY_FOR_RECEIVER_ACTION_SELECTION'
