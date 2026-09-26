import copy
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from geometrodynamics.waves import handle_evolution as e
from experiments.closure_ledger import handle_evolution_probe as p

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'experiments/closure_ledger/runs/20260923_handle_evolution'


@pytest.fixture(scope='module', autouse=True)
def restored_evidence():
    from experiments.closure_ledger.restore_handle_evidence import restore
    restore(RUN)


def test_independent_equations_and_nonconstant_norm():
    result=e.validate_equations()
    assert result['target_metric_compatible'] and result['equivariant_closure']
    assert result['round_evolution_error']<1e-10
    assert result['coordinate_geometry_error']<1e-10
    assert result['initial_norm_acceleration']==pytest.approx(-3.95897327444315)


def test_twisted_derivatives_converge_through_identification():
    errors=[]
    for n in (64,128):
        s=np.linspace(-e.L,e.L,n,endpoint=False)
        y=np.ones((8,n))
        w=np.pi/(2*e.L)
        y[4]=np.sin(w*s);y[6]=np.cos(w*s)
        first,second=e.derivatives(y,2*e.L/n)
        errors.append(np.max(abs(first[4]-w*np.cos(w*s))))
        assert np.max(abs(second[4]+w*w*y[4]))<1e-7
    assert 15<errors[0]/errors[1]<17


def test_geodesic_crossing_does_not_insert_a_momentum_sign_flip():
    n=64;y=np.zeros((8,n));y[:2]=1
    dy=np.zeros_like(y)
    g={'f':np.ones(n)}
    particles=np.array([[e.L-1e-8,e.L+1e-8,-e.L-1e-8],[.5,.5,-.5]])
    result=e.particle_rhs(y,dy,g,particles)
    np.testing.assert_allclose(result[0],particles[1]/np.sqrt(1+particles[1]**2))
    assert np.max(abs(result[1]))==0


def test_evolution_domain_stop():
    y=np.zeros((8,16));y[:2]=1;y[4]=3
    with pytest.raises(ArithmeticError,match='positivity'):
        e.geometry(y,np.zeros_like(y),np.zeros_like(y))


def test_homogeneous_controls_keep_growing_mode_and_constraint():
    controls=e.homogeneous_controls()
    assert max(r['constraint_max'] for r in controls)<1e-9
    assert np.max(abs(np.asarray(controls[0]['states'])[0]-1))<1e-10
    for row in controls[1:]:
        scale=np.asarray(row['states'])[0]
        assert abs(scale[-1]-1)>abs(scale[0]-1)


@pytest.fixture(scope='module')
def evidence():
    return json.loads((RUN/'evolution.json').read_text())


def test_full_evolution_evidence_replay(evidence):
    measured=p.score(evidence)
    recorded=json.loads((RUN/'evolution_verdict.json').read_text())
    assert p.agreement(measured,recorded)
    assert measured['evidence']
    assert not measured['verdicts']['DISCRETE_RECIPROCAL_MOMENTUM_EXCHANGE']


@pytest.mark.parametrize('damage',['momentum','metric','source','missing','nan'])
def test_corrupted_evolution_clears_verdicts(evidence,monkeypatch,damage):
    # Full replay above recomputes the baseline; these isolate evidence damage.
    monkeypatch.setattr(p,'run',lambda:evidence)
    bad=copy.deepcopy(evidence)
    if damage=='momentum':bad['groups'][0]['runs'][-1]['trajectories'][50][1][0]+=1
    elif damage=='metric':bad['groups'][0]['runs'][-1]['final_fields'][0][0]+=1
    elif damage=='source':bad['sources'][p.SOURCES[0]]='0'*64
    elif damage=='missing':bad['groups'][0]['runs'][-1]['crossings'].pop()
    else:bad['groups'][0]['runs'][-1]['final_fields'][0][0]=float('nan')
    result=p.score(bad)
    assert not result['evidence']
    assert not any(result['verdicts'].values())


def test_geometry_work_is_not_an_optional_flux_term(evidence):
    for group in evidence['groups']:
        scored=p.group_score(group)
        assert scored['omitted_geometry_defects'][-1]>scored['balance_defects'][-1]


def test_impulse_classifier_can_detect_a_resolution_independent_jump(evidence):
    group=copy.deepcopy(evidence['groups'][0])
    for run in group['runs'][-2:]:
        for window in run['crossings'][0]['windows']:
            window['p_hat']=.01
    # A positive classifier control, not fabricated evidence: archive replay
    # would reject these altered windows against the computed worldlines.
    assert p.group_score(group)['finite_impulse']


def test_failed_cli_withdraws_stale_verdict(tmp_path):
    target=tmp_path/'evolution_verdict.json'
    target.write_text(json.dumps({'verdicts':dict.fromkeys(p.VERDICTS,True)}))
    command=[sys.executable,'-m','experiments.closure_ledger.handle_evolution_probe',
             '--output-dir',str(tmp_path),'--rescore',str(tmp_path/'absent.json')]
    run=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
    assert run.returncode!=0
    assert not any(json.loads(target.read_text())['verdicts'].values())
