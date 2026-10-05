import numpy as np
import pytest
from geometrodynamics.waves import r3_nonlinear_budget as nb


def test_horizon_from_twist():
    assert nb.horizon(4.908582988640952,.001)==102
    assert nb.horizon(4.908582988640952,-.0005)==204
    with pytest.raises(ValueError): nb.horizon(0.,.001)
    with pytest.raises(ValueError): nb.horizon(4.9,1e-8)


def test_failure_labels_are_reachable():
    assert nb.decision(2.1,1,0,0,True,1,0,100)=='TUBE_ESCAPE'
    assert nb.decision(.5,1,.021,0,True,1,0,100)=='CONTROL_BUDGET_EXCEEDED'
    assert nb.decision(.5,1,.011,.24,True,1,0,100)=='CONTROL_BUDGET_EXCEEDED'
    assert nb.decision(.5,1,0,0,True,100,.49,100)=='DRIFT_HORIZON_NOT_RESOLVED'
    assert nb.decision(.5,1,.01,.2,True,100,.6,100)=='HORIZON_COMPLETED'


def test_intervention_not_charged_to_unforced_arm():
    assert nb.decision(.5,1,0,0,False,1,0,100)=='CONTINUE'
    assert nb.decision(2.1,1,0,0,False,1,0,100)=='TUBE_ESCAPE'


@pytest.mark.parametrize('x',[float('nan'),float('inf'),-1.])
def test_invalid_diagnostics_fail_closed(x):
    with pytest.raises(ValueError): nb.decision(x,1,0,0,True,1,0,100)

import copy
import json
import shutil
from experiments.closure_ledger import r3_nonlinear_budget_probe as probe
from experiments.closure_ledger import r3_nonlinear_budget_replay as portable

MANIFEST_SHA = '4ea7db2a0f0dc002bb1a4d0f7a12375d9c9f5a38d719d395741efdb4d51e44e1'


@pytest.fixture(scope='module')
def measured_family():
    return probe.family()


@pytest.fixture(scope='module')
def first_case():
    return probe.read(probe.RUN/'case_0.json.gz.b64')


def test_clock_completion_is_included_in_norm(first_case):
    y=np.array(first_case['initial']['initial']);yy=y.copy();yy[7]+=.01
    assert np.linalg.norm(nb.coords(yy)-nb.coords(y))==pytest.approx(.01)


def test_reject_fabricated_pass_or_free_reset(measured_family,first_case):
    r=copy.deepcopy(first_case['runs'][1])
    r['terminal']='HORIZON_COMPLETED'
    with pytest.raises(ValueError): probe.validate_run(measured_family,first_case['initial'],r)
    r=copy.deepcopy(first_case['runs'][1]);r['steps'][0]['proposed_cost']=0.
    with pytest.raises(ValueError): probe.validate_run(measured_family,first_case['initial'],r)


def test_reject_changed_history_and_omitted_return(measured_family,first_case):
    r=copy.deepcopy(first_case['runs'][1]);r['steps'][0]['returns'][0]['states'][0][0]+=.001
    with pytest.raises(ValueError): probe.validate_run(measured_family,first_case['initial'],r)
    r=copy.deepcopy(first_case['runs'][1]);r['steps'][0]['returns'].pop()
    with pytest.raises(ValueError): probe.validate_run(measured_family,first_case['initial'],r)


def test_over_budget_proposal_was_not_applied(measured_family,first_case):
    r=first_case['runs'][1];ini=first_case['initial'];s=r['steps'][0]
    out=probe.validate_run(measured_family,ini,r)
    assert out['terminal']=='CONTROL_BUDGET_EXCEEDED'
    assert s['proposed_cost']>.02*ini['d0']
    assert not s['applied'] and s['cumulative_cost']==0
    assert np.array_equal(s['after'],s['pre'])


def test_authenticated_replay():
    assert portable.MANIFEST_SHA==MANIFEST_SHA
    result=portable.replay()
    assert result['label']=='CONTROLLED_NONLINEAR_BOUND_FAILED'
    assert len(result['cases'])==4
    for case in result['cases']:
        assert case['controlled']=='CONTROL_BUDGET_EXCEEDED'
        assert case['unforced']=='NUMERICALLY_UNRESOLVED'
        assert case['integrator_comparison']['controlled']['agrees']
        for r in (case['runs'][1],case['runs'][3]):
            assert r['steps']==1
            assert r['applied_cumulative_cost_over_d0']==0
            assert r['peak_proposed_cost_over_d0']>.7


@pytest.mark.parametrize('change',['bytes','missing','label','manifest'])
def test_reject_corrupted_evidence(tmp_path,change):
    dst=tmp_path/'archive';shutil.copytree(probe.RUN,dst)
    if change=='bytes':
        with (dst/'case_0.json.gz.b64').open('a') as f: f.write('X')
    elif change=='missing':
        (dst/'case_3.json.gz.b64').unlink()
    elif change=='label':
        p=dst/'result.json';obj=json.loads(p.read_text())
        obj['label']='CONTROLLED_NONLINEAR_BOUND_SUPPORTED_ON_TESTED_HORIZON'
        p.write_text(json.dumps(obj))
    else:
        p=dst/'manifest.json';p.write_text(p.read_text()+'\n')
    with pytest.raises((ValueError,FileNotFoundError)):
        portable.replay(dst)


@pytest.mark.parametrize('change',['missing_case','reorder_cases','reorder_runs'])
def test_reject_incomplete_or_reordered_schedule(measured_family,change):
    raw=json.loads((probe.RUN/'provenance.json').read_text())
    raw['cases']=[probe.read(probe.RUN/f'case_{k}.json.gz.b64') for k in range(4)]
    if change=='missing_case': raw['cases'].pop()
    elif change=='reorder_cases': raw['cases'].reverse()
    else: raw['cases'][0]['runs'].reverse()
    with pytest.raises(ValueError): probe.assess(measured_family,raw)


def test_action_offset_and_initial_linear_tuning(measured_family,first_case):
    from geometrodynamics.waves import r3_family as rf
    f=measured_family;ini=first_case['initial'];loop=np.array(ini['loop'])
    assert rf.loop_action(loop)-f.action==pytest.approx(-.001,abs=1e-10)
    projector=f.projector(0.)
    assert np.linalg.norm(projector@(loop[0]-f.Z[0]))<1e-11
    assert ini['tuning_norm']>0
    assert ini['log10_autonomous_unstable_tolerance']<-395


@pytest.mark.parametrize('change',['phase','state'])
def test_portable_replay_checks_recorded_initial_setup(measured_family,first_case,change):
    case=copy.deepcopy(first_case)
    if change=='phase': case['initial']['initial_phase']+=.001
    else: case['initial']['initial'][0]+=.001
    f=portable.RecordedInitialFamily(measured_family,[case])
    with pytest.raises(ValueError): f.initial(case['detuning'])
