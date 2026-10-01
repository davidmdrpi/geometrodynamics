"""A selection criterion must accept a plateau and reject scaling, noise and missing data."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.waves import selection_plateau as p
from experiments.closure_ledger import selection_plateau_probe as probe

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'experiments/closure_ledger/runs/20260928_selection_plateau'


def test_work_identities_from_independent_directional_differences():
    y=p.preparation(.18,.63,.47)
    dy=p.dynamics.conformal_rhs(y)
    action,work,_=p.quantities(y)
    h=1e-6
    derivative=(p.quantities(y+h*dy)[0]-p.quantities(y-h*dy)[0])/(2*h)
    np.testing.assert_allclose(derivative,work,rtol=2e-6,atol=2e-8)
    assert np.max(abs(work))>1e-4  # not a conserved-loop identity
    assert np.all(action>0)


def test_plateau_detector_positive_control_and_conservative_failures():
    constant=np.full((8,12,3,2),.01*p.J_BG)
    assert all(row['passes'] for row in p.classify(constant))
    assert p.selection_verdict({'valid':True},p.classify(constant))=='CANDIDATE_FINITE_TIME_PLATEAU_NOT_QUANTIZATION'
    scaling=constant*(p.AMPLITUDES[1:]/p.AMPLITUDES[1])[:,None,None,None]**2
    assert not any(row['passes'] for row in p.classify(scaling))
    for values in (np.zeros_like(constant),constant*1e-6,-constant):
        assert not any(row['passes'] for row in p.classify(values))
    altered=constant.copy();altered[:,0,:,0]*=2
    assert not any(row['passes'] for row in p.classify(altered))
    assert p.selection_verdict({'valid':False},p.classify(constant))=='INCONCLUSIVE_NUMERICAL_FAILURE'
    assert p.selection_verdict({'valid':True},[])=='INCONCLUSIVE_NUMERICAL_FAILURE'
    with pytest.raises(ValueError):p.classify(constant[:-1])
    constant[0,0,0,0]=np.nan
    with pytest.raises(ValueError):p.classify(constant)


def test_fresh_background_and_coupled_work_balance():
    background=p.summarize_history(p.evolve(0.,0.,.4))
    assert np.max(abs(np.array(background['changes'])))/p.J_BG<1e-9
    coupled=p.summarize_history(p.evolve(.08,.4,.2))
    assert coupled['constraint_max']<1e-8
    assert coupled['work_error_over_Jbg']<1e-8
    assert np.max(abs(np.array(coupled['changes'])))>1e-5


def test_raw_archive_replay_and_provenance():
    raw=probe.load_raw(RUN/'states.json.gz.b64')
    saved=json.loads((RUN/'plateau.json').read_text())
    fresh=probe.summarize(raw)
    assert fresh['selection_verdict']==saved['selection_verdict']
    assert fresh['numerical_gates']==saved['numerical_gates']
    assert fresh['localized_receiver_verdict']=='NOT_TESTED'
    for a,b in zip(fresh['cases'],saved['cases']):
        if 'failure' in a:
            assert a==b
        else:
            np.testing.assert_allclose(a['diagnostics']['window_means'],b['diagnostics']['window_means'],rtol=1e-9,atol=1e-12)
    for path,digest in saved['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    assert hashlib.sha256((RUN/'states.json.gz.b64').read_bytes()).hexdigest()==saved['archive_sha256']


def test_missing_failed_or_relabelled_case_cannot_pass():
    raw=probe.load_raw(RUN/'states.json.gz.b64')
    removed=raw['cases'].pop()
    with pytest.raises(ValueError,match='schedule'):probe.validate_raw(raw)
    raw['cases'].append(dict(parameters=removed['parameters'],failure='synthetic solver failure'))
    report=probe.summarize(raw)
    assert report['selection_verdict']=='INCONCLUSIVE_NUMERICAL_FAILURE'
    raw['cases'][-1]=removed
    raw['cases'][0]['states'][0][0]+=.01
    with pytest.raises(ValueError,match='preparation'):probe.validate_raw(raw)


def test_refinement_preserves_original_failure_and_rejects_damaged_data():
    from experiments.closure_ledger import selection_plateau_refinement as ref
    original=json.loads((RUN/'plateau.json').read_text())
    raw=json.loads((RUN/'refinement_states.json').read_text())
    saved=json.loads((RUN/'refinement.json').read_text())
    replay=ref.analyze(original,raw['rows'])
    assert original['selection_verdict']=='INCONCLUSIVE_NUMERICAL_FAILURE'
    assert original['numerical_gates']['window_quadrature'] is False
    assert replay['selection_verdict']==saved['selection_verdict']=='NO_ROBUST_PLATEAU_IN_REGISTERED_FAMILY'
    assert all(replay['numerical_gates'].values())
    for path,digest in saved['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    assert hashlib.sha256((RUN/'refinement_states.json').read_bytes()).hexdigest()==saved['raw_sha256']
    with pytest.raises(ValueError):ref.analyze(original,raw['rows'][:1])
    raw['rows'][0]['states'][-1][1]*=1.1
    assert ref.analyze(original,raw['rows'])['selection_verdict']=='INCONCLUSIVE_NUMERICAL_FAILURE'
