"""Independent identities, both asymptotic columns and fail-closed verdicts."""
import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from geometrodynamics.waves import tt_turning_transport as t
from experiments.closure_ledger import tt_turning_transport_probe as probe

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/'experiments/closure_ledger/runs/20260912_tt_turning_transport/probe.json'

@pytest.fixture(scope='module')
def evidence():return json.loads(ARCHIVE.read_text())


def result(r,checks=None):
    return t.verdict(r['checks'] if checks is None else checks,r['certificate'],r['transports'],r,r)


def test_dependency_table_matches_freeze():
    frozen={}
    for line in (ROOT/'docs/tt_turning_transport_prereg.md').read_text().splitlines():
        cells=[s.strip() for s in line.split('|')]
        if len(cells)==4 and cells[1] in t.DEPENDENCIES:frozen[cells[1]]=cells[2].replace(', ','')
    assert frozen==t.DEPENDENCIES


@pytest.mark.parametrize('gate',list(t.DEPENDENCIES))
@pytest.mark.parametrize('missing',[False,True])
def test_each_gate_only_invalidates_dependent_targets(evidence,gate,missing):
    assert all(result(evidence)[name]!='UNRESOLVED' for name in t.TARGETS.values())
    checks=evidence['checks'].copy()
    if missing:checks.pop(gate)
    else:checks[gate]=False
    outcome=result(evidence,checks)
    for target,name in t.TARGETS.items():
        assert (outcome[name]=='UNRESOLVED')==(target in t.DEPENDENCIES[gate])


@pytest.mark.parametrize('key,targets',[('certificate','B'),('transports','T'),('series','FT'),('endpoint','A')])
@pytest.mark.parametrize('damage',['missing','nonfinite'])
def test_malformed_evidence_is_target_specific(evidence,key,targets,damage):
    r=copy.deepcopy(evidence)
    if damage=='missing':r[key]=None
    elif key=='certificate':r[key]['facts']['P_at_1']=float('nan')
    elif key=='transports':r[key][0]['comparison_maps'][0][0][0]=float('nan')
    elif key=='series':r[key][0]['coefficients'][0][0]=float('nan')
    else:r[key][0]['ratios'][0]=float('nan')
    outcome=result(r)
    for target,name in t.TARGETS.items():assert (outcome[name]=='UNRESOLVED')==(target in targets)


def test_certificate_rejects_scan_and_mutated_facts(evidence):
    assert t.verify_certificate(evidence['certificate'])
    assert not t.verify_certificate({'verified':True,'roots':evidence['roots']})
    altered=t.turning_certificate();altered['facts']['P_at_1']='-120'
    assert not t.verify_certificate(altered)
    assert t.verify_certificate(t.turning_certificate())


@pytest.mark.parametrize('alpha',[0.,.47,1.4,2.8])
def test_potential_and_eta_sign_independent(alpha):
    for x in (.2,.5,1.,3.):
        m,mp,mpp,k,_=t.coefficients(x,alpha)
        W,Wp=t.potential(x,alpha)
        assert W==pytest.approx(k/m-mpp/(2*m)+mp**2/(4*m*m),abs=1e-12)
        h=1e-5
        dx=(t.potential(x+h,alpha)[0]-t.potential(x-h,alpha)[0])/(2*h)
        assert Wp==pytest.approx(-dx,rel=2e-8,abs=1e-8)
    root=t.turning_point(alpha)
    assert 2.97<=root['R']<=3.03 and root['Wprime']<-20


def test_bare_root_and_attainable_supported_lower_endpoint():
    # Bare m=R^2, k=8R^2; independent reduction gives W0=9-R^2.
    for R in (1.1,2.,3.,4.):
        mp=math.sqrt(2)*R*(R*R-1);mpp=3*R**4-4*R*R+1
        assert 8-mpp/(2*R*R)+mp*mp/(4*R**4)==pytest.approx(9-R*R)
    x=t.x_from_radius(3.)
    alpha=math.atan2(-4,3*math.sqrt(2))+2*x
    assert abs(t.potential(x,alpha)[0])<1e-13


@pytest.mark.parametrize('alpha',[0.,.37,math.pi/2,2.1])
def test_series_resonance_and_review_correction(alpha):
    b,resonance=t.series_coefficients(alpha,12)
    assert np.array_equal(resonance,[0,0])
    assert np.array_equal(b[1],[0,0])
    assert np.array_equal(b[2],[4,0])
    assert b[4,0]==pytest.approx(math.cos(alpha)**2/16-28/3)
    assert b[4,0]!=pytest.approx(-8)
    data=t.series_basis(.01,alpha,12)
    assert np.linalg.det(data['matrix'])==pytest.approx(1,abs=1e-11)
    assert max(data['normalized'])<1e-11
    # Independent second-order equation, evolved in x; compare beta and p.
    x0=.02;x1=.4;start=t.series_basis(x0,alpha,12)['matrix']
    m=t.coefficients(x0,alpha)[0]
    initial=np.array([start[0],-start[1]/m]).ravel()
    def rhs(x,y):
        b,v=y.reshape(2,2);m,mp,_,k,_=t.coefficients(x,alpha)
        return np.array([v,(mp*v-k*b)/m]).ravel()
    sol=solve_ivp(rhs,(x0,x1),initial,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.001)
    b,v=sol.y[:,-1].reshape(2,2)
    expected=np.array([b,-t.coefficients(x1,alpha)[0]*v])
    assert np.allclose(t.output_basis(alpha,x0=x0,match=x1),expected,rtol=1e-10,atol=1e-10)


def test_complete_map_independent_reproduction_and_rank_one_control(evidence):
    row=evidence['transports'][3];alpha=row['alpha']
    matrix=np.linalg.solve(t.output_basis(alpha),t.input_basis(20,alpha))
    assert np.allclose(matrix,row['comparison_maps'][0],rtol=1e-10,atol=1e-11)
    assert np.linalg.norm(matrix.T@t.J@matrix-t.J)<1e-10
    damaged=copy.deepcopy(evidence['transports'])
    for row in damaged:
        for M in row['comparison_maps']:M[1]=[0,0]
    assert not t.valid_transport(damaged)
    # Keeping C alone cannot distinguish this nonzero initial datum from zero.
    initial=np.array([-matrix[0,1],matrix[0,0]])
    assert abs((matrix@initial)[0])<1e-12
    assert abs((matrix@initial)[1])>1e-3


def test_actions_and_scope_controls(evidence):
    assert probe.action_identity()=='0'
    assert t.valid_actions(evidence) and t.valid_future(evidence)
    assert max(evidence['constant_oscillator_errors'])<1e-12
    assert max(r['pulled_back_invariant_error'] for r in evidence['actions'])<1e-10
    assert max(max(r['J_beta']) for r in evidence['actions'])>100
    for row in evidence['tuned']:
        vals=[r['J_y'] for r in row['samples']]
        assert (vals[-1]<vals[0])==(row['turn_yprime']==0)
    assert result(evidence)['all_classical_invariants_excluded'] is False
    assert result(evidence)['quantization']=='NOT_DERIVED'


@pytest.mark.parametrize('failure',['exception','nonfinite','missing','failed_gate'])
def test_cli_overwrites_stale_success_with_failure(evidence,tmp_path,monkeypatch,failure):
    (tmp_path/'probe.json').write_text('{"checks_passed":true}')
    r=copy.deepcopy(evidence)
    if failure=='exception':
        def broken():raise RuntimeError('injected integration failure')
        monkeypatch.setattr(probe,'run_probe',broken)
    else:
        if failure=='nonfinite':r['transports'][0]['comparison_maps'][0][0][0]=float('nan')
        elif failure=='missing':r.pop('series')
        else:r['certificate']['facts']['P_at_1']='1'
        monkeypatch.setattr(probe,'run_probe',lambda:r)
    assert probe.main(['--output-dir',str(tmp_path)])==1
    output=json.loads((tmp_path/'probe.json').read_text())
    assert not output['checks_passed']
    assert output['verdict']['all_phase_turning_bound']=='UNRESOLVED'


def test_finalize_accepts_reproduced_archive(evidence):
    r=probe.finalize(copy.deepcopy(evidence))
    assert r['checks_passed']
    assert all(v!='UNRESOLVED' for k,v in r['verdict'].items() if k in t.TARGETS.values())
