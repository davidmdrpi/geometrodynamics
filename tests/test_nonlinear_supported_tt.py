"""Nonlinear full-field completion, prospective refinement and scoped failures."""
import copy
import hashlib
import json
import gzip
import math
from pathlib import Path
import numpy as np
import pytest
from scipy.linalg import expm

from geometrodynamics.waves import nonlinear_supported_tt as n
from experiments.closure_ledger import nonlinear_supported_tt_probe as p
from experiments.closure_ledger import nonlinear_supported_tt_refinement_probe as r

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'experiments/closure_ledger/runs/20260913_nonlinear_supported_tt'

@pytest.fixture(scope='module')
def original():return json.loads(gzip.decompress((RUN/'probe.json.gz').read_bytes()))
@pytest.fixture(scope='module')
def refined():return json.loads(gzip.decompress((RUN/'refinement.json.gz').read_bytes()))


def test_original_misses_are_retained_and_refinement_is_separate(original,refined):
    assert sum(original['checks'].values())==11
    assert not original['checks']['linear_recovery'] and not original['checks']['quadratic_response']
    assert original['verdict']['finite_amplitude_response']=='UNRESOLVED'
    assert original['verdict']['future_persistence']=='UNRESOLVED'
    assert p.evidence_gates(original)==original['checks']
    assert p.verdict(p.evidence_gates(original),original)==original['verdict']
    assert r.valid_rows(refined,original)
    assert r.finalize(copy.deepcopy(refined),original)==refined
    assert refined['checks_passed']


@pytest.mark.parametrize('level',[1,2])
@pytest.mark.parametrize('time_index',[0,3,6])
def test_refinement_constraints_detect_invisible_velocity_damage(original,refined,level,time_index):
    damaged=copy.deepcopy(refined)
    # A' is absent from the measured observables: variation checks alone pass.
    damaged['rows'][0]['levels'][level]['plus'][time_index][1]+=.1
    assert r.valid_rows(damaged,original)
    assert not r.constraints_valid(damaged)
    checks=r.evidence_checks(damaged,original)
    assert [key for key,passed in checks.items() if not passed]==['constraint_propagation']
    out=r.verdict(refined['checks'],damaged,original)
    for target,name in p.TARGETS.items():
        assert (out[name]=='UNRESOLVED')==(target in 'NF')


def test_pre_refinement_source_hashes_are_archived():
    freeze=(ROOT/'docs/nonlinear_supported_tt_refinement_prereg.md').read_text()
    for file in ('original_module.py.txt','original_probe.py.txt','probe.json'):
        data=gzip.decompress((RUN/(file+'.gz')).read_bytes()) if file=='probe.json' else (RUN/file).read_bytes()
        assert hashlib.sha256(data).hexdigest() in freeze


def test_dependency_table_matches_unchanged_scientific_freeze():
    frozen={}
    for line in (ROOT/'docs/nonlinear_supported_tt_prereg.md').read_text().splitlines():
        cells=[v.strip() for v in line.split('|')]
        if len(cells)==4 and cells[1] in p.DEPENDENCIES:frozen[cells[1]]=cells[2].replace(', ','')
    assert frozen==p.DEPENDENCIES


@pytest.mark.parametrize('gate',list(p.DEPENDENCIES))
@pytest.mark.parametrize('missing',[False,True])
def test_each_gate_fails_only_its_dependent_targets_from_affirmative_evidence(original,refined,gate,missing):
    checks=refined['checks'].copy()
    if missing:checks.pop(gate)
    else:checks[gate]=False
    result=r.verdict(checks,refined,original)
    for target,name in p.TARGETS.items():assert (result[name]=='UNRESOLVED')==(target in p.DEPENDENCIES[gate])


@pytest.mark.parametrize('damage,targets',[('future','F'),('initial','CNF'),('coarse','NF'),('unknown','DCNF'),('refinement','NF')])
def test_raw_evidence_cannot_be_replaced_by_passing_flags(original,refined,damage,targets):
    a=copy.deepcopy(original);b=copy.deepcopy(refined);checks=b['checks'].copy()
    if damage=='future':a['future']['rational_margins']['matter_margin']='0'
    elif damage=='initial':a['initial'][0]['state'][1]+=.1
    elif damage=='coarse':a['comparisons'][0]['coarse_states'][2][11]+=.1
    elif damage=='unknown':checks['not_frozen']=True
    else:b['rows'][0]['levels'][-1]['plus'][3][11]+=.01
    out=r.verdict(checks,b,a)
    for target,name in p.TARGETS.items():assert (out[name]=='UNRESOLVED')==(target in targets)


@pytest.mark.parametrize('kind',['missing','nonfinite'])
def test_malformed_future_certificate_preserves_independent_targets(original,refined,kind):
    a=copy.deepcopy(original)
    if kind=='missing':a.pop('future')
    else:a['future']['entry']['A_min']=float('nan')
    out=r.verdict(refined['checks'],refined,a)
    assert out['future_persistence']=='UNRESOLVED'
    assert all(out[p.TARGETS[t]]!='UNRESOLVED' for t in 'DCN')


def test_exact_certificates_are_recomputed_and_not_mutable_shared_flags():
    c=n.exact_certificate();assert c['all_zero'] and set(c['identities'].values())=={'0'}
    c['identities']['completion_current']='1'
    assert n.exact_certificate()['identities']['completion_current']=='0'
    f=n.future_certificate();assert f['all_positive']
    from fractions import Fraction
    assert all(Fraction(v)>0 for v in f['rational_margins'].values())
    f['rational_margins']['shape_margin']='0'
    assert n.future_certificate()['rational_margins']['shape_margin']!='0'


@pytest.mark.parametrize('phase',[0.,math.pi/4,math.pi/2,2.3])
@pytest.mark.parametrize('pair',[1,2])
def test_complete_initial_data_and_independent_full_fields(pair,phase):
    U,V=p.pairs()[pair]
    y=n.initial_data(U,V,.071,departure=.12,phase=phase,a=.9,kappa=.7)
    assert max(n.constraints(y,.9,.7)['normalized'])<1e-12
    for angles in ((.73,.92,.57),(1.2,.6,1.4)):
        left=n.full_geometry(y,a=.9,kappa=.7,angles=angles)
        right=n.full_geometry(y,a=.9,kappa=.7,angles=angles,clock='proper')
        assert left['normalized']<1e-11
        assert np.linalg.norm(left['KG'])<1e-11
        assert np.linalg.norm(left['residual']-right['residual'])<1e-11


def test_shift_variation_recovers_nonzero_momentum_constraint():
    U,V=p.pairs()[2];y=n.initial_data(U,V,.1,phase=.4,rigid=True)
    A,Ap,q,qp,M,L=n.unpack(y);H,_,rr,ell,tt,_,_=n.ingredients(y)
    Mp=2*M@L
    def lagrange(shift,lapse=1.):
        D=2*n.rt.cross_matrix(shift)
        Lcov=np.linalg.solve(M,Mp-(D@M-M@D))/2
        velocity=n.B(qp)-n.B(q)@np.einsum('i,ijk->jk',shift,n.S)
        kinetic=-3*Ap*Ap+np.trace(velocity.T@velocity)/8+H*np.trace(Lcov@Lcov)/2
        potential=H*rr/2-(q@q)*tt/2-1.5*A**4
        return kinetic/lapse+lapse*potential
    h=1e-5
    derivative=np.array([(lagrange(h*e)-lagrange(-h*e))/(2*h) for e in np.eye(3)])
    C=n.constraints(y)['residual']
    assert np.linalg.norm(C[1:])>.01
    assert np.allclose(derivative,C[1:],atol=1e-9)
    assert (lagrange(np.zeros(3),1+h)-lagrange(np.zeros(3),1-h))/(2*h)==pytest.approx(-C[0],abs=1e-8)


def test_off_shell_tensor_and_scalar_accelerations_are_not_assumed_on_shell():
    U,V=p.pairs()[4];y=n.initial_data(U,V,.19,phase=.63)
    dy=n.conformal_rhs(y);dy[1]+=.4;dy[7:11]+=[.2,-.1,.3,.1]
    M=n.unpack(y)[4];root=expm(.19*U)
    dy[20:29]+=(np.linalg.solve(root,n.rt.STF_BASIS[3]@root)).ravel()
    g=n.full_geometry(y,dy,angles=(.63,1.21,.88));expected=n.expected_residual(y,dy,angles=(.63,1.21,.88))
    assert np.linalg.norm(g['residual'])>.1
    assert np.linalg.norm(g['residual']-expected['residual'])<1e-10
    assert np.linalg.norm(g['KG']-expected['KG'])<1e-11


def test_pure_rigid_obstruction_does_not_exclude_responsive_fields(original):
    assert original['controls'][0]['rigid_momentum']<1e-12
    assert original['controls'][1]['rigid_momentum']>.03
    assert original['controls'][1]['responsive_momentum']<1e-12
    assert min(v['minimal_stress_error'] for v in original['controls'])>.1
    assert min(v['wrong_clock_error'] for v in original['controls'])>.1


def test_variational_coefficient_reproduced_at_unfrozen_phase():
    U,V=p.pairs()[4];times=[0.,.5,1.,2.];eps=.003
    prediction=n.second_variations(U,V,times,phase=.27)
    vals={e:np.array([n.observables(y) for y in n.evolve(n.initial_data(U,V,e,phase=.27),times)]) for e in (0.,-eps,eps)}
    second=(vals[eps]+vals[-eps]-2*vals[0])/(2*eps*eps)
    assert max(p.relative(v,row['second']) for v,row in zip(second,prediction))<1e-3


def test_future_bootstrap_inequalities_on_independent_regular_data():
    # These estimates do not use the momentum constraint; the expansion bounds
    # use a positive root of the energy constraint and the same bootstrap box.
    rng=np.random.default_rng(19)
    for _ in range(40):
        coeff=rng.normal(size=5);coeff*=.07/np.linalg.norm(coeff)
        root=expm(n.rt.tensor(coeff));M=root@root
        V=n.rt.tensor(rng.normal(size=5));V*=.4/np.linalg.norm(V);L=np.linalg.solve(root,V@root)
        q=rng.normal(size=4);qp=rng.normal(size=4)
        fac=7/(np.linalg.norm(q)+np.linalg.norm(qp));q*=fac;qp*=fac
        A=float(rng.uniform(32,80));y=n.pack(A,0,q,qp,M,L)
        E=n.constraints(y)['residual'][0];y[1]=math.sqrt(E/3)
        H,Hp,rr,ell,tt,force,inv=n.ingredients(y)
        assert H>=.98*A*A and 2/3<=y[1]/A**2<=1
        assert abs(rr)<=15 and tt<=4 and abs(tt+(rr+ell)/6)<7
        assert np.linalg.norm(force)<=14*np.linalg.norm(M-n.I)


def test_frozen_tail_samples_are_diagnostics_not_the_continuation_certificate(original):
    assert len(original['tails'])==30
    for row in original['tails']:
        eligible=[v for v in row['diagnostics'] if v['A']>=32]
        assert eligible
        assert any(v['shape']<=.1 and v['velocity']<=.1 and v['matter']<=4 for v in eligible)
        assert abs(row['diagnostics'][-1]['proper_Hubble']-1/math.sqrt(2))<1e-7
    bad=copy.deepcopy(original);bad['future']={'all_positive':True,'tails':original['tails']}
    assert not p.evidence_gates(bad)['future_continuation_bound']


@pytest.mark.parametrize('bad',[dict(departure=0),dict(a=-1),dict(kappa=0),dict(epsilon=float('nan'))])
def test_invalid_chart_data_rejected(bad):
    args=dict(epsilon=.01);args.update(bad)
    with pytest.raises(ValueError):n.initial_data(*p.pairs()[0],**args)


@pytest.mark.parametrize('kind',['exception','nonfinite','missing'])
def test_primary_cli_overwrites_stale_success(original,tmp_path,monkeypatch,kind):
    a=copy.deepcopy(original)
    if kind=='nonfinite':a['future']['entry']['A_min']=float('nan')
    elif kind=='missing':a.pop('initial')
    def run():
        if kind=='exception':raise ArithmeticError('injected failure')
        return a
    monkeypatch.setattr(p,'run_probe',run)
    (tmp_path/'probe.json').write_text('{"checks_passed":true}')
    assert p.main(['--output-dir',str(tmp_path)])==1
    result=json.loads((tmp_path/'probe.json').read_text())
    assert not result['checks_passed']
    if kind=='nonfinite':
        assert result['verdict']['exact_reduction']=='EXACT_HOMOGENEOUS_REDUCTION_VERIFIED'
        assert result['nonfinite_evidence_paths']


@pytest.mark.parametrize('kind',['exception','nonfinite','tampered'])
def test_refinement_cli_overwrites_stale_success(original,refined,tmp_path,monkeypatch,kind):
    b=copy.deepcopy(refined)
    if kind=='nonfinite':b['rows'][0]['levels'][0]['plus'][0][0]=float('nan')
    elif kind=='tampered':b['rows'][0]['levels'][0]['plus'][2][11]+=.1
    def run(a):
        if kind=='exception':raise ArithmeticError('injected extension failure')
        return b
    monkeypatch.setattr(r,'run_probe',run)
    (tmp_path/'refinement.json').write_text('{"checks_passed":true}')
    assert r.main(['--original',str(RUN/'probe.json.gz'),'--output-dir',str(tmp_path)])==1
    result=json.loads((tmp_path/'refinement.json').read_text())
    assert not result['checks_passed']
    if kind=='nonfinite':
        assert result['verdict']['exact_reduction']=='EXACT_HOMOGENEOUS_REDUCTION_VERIFIED'
        assert result['nonfinite_evidence_paths']
    assert result['verdict']['finite_amplitude_response']=='UNRESOLVED'
