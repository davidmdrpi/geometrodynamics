"""Round 12 independent stress, exact-rank and review failure controls."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.waves import odd_multiplet_support as ms
from experiments.closure_ledger import odd_multiplet_support_probe as probe


@pytest.fixture(scope='module')
def report():
    return probe.finalize(probe.run_probe())


def test_gate_schema_matches_published_review_amendment():
    text=(Path(__file__).parents[1]/'docs/odd_multiplet_support_prereg.md').read_text()
    symbols=dict(zip('BKMDP',ms.PHYSICAL_VERDICTS))
    rows={}
    for line in text.splitlines():
        if line.startswith('| `'):
            parts=line.split('|')
            rows[parts[1].strip().strip('`')]=tuple(symbols[x.strip()] for x in parts[2].split(','))
    assert tuple(rows)==ms.REQUIRED_CHECKS
    assert rows==ms.CHECK_TARGETS
    assert all(ms.CHECK_TARGETS.values())


def test_all_frozen_gates_pass_and_scope_stays_open(report):
    assert report['checks_passed'],report['verdict']
    assert {k:report['verdict'][k] for k in ms.SCOPE}==ms.SCOPE
    assert all(c['all_zero'] for c in report['exact_addition'])
    assert report['exact_background']['all_zero']


@pytest.mark.parametrize('k',[1,3,5])
def test_componentwise_odd_parity_and_actual_even_control(k,report):
    row=next(r for r in report['backgrounds'] if r['degree']==k)
    assert row['parity_absolute']<1e-12
    assert row['even_parity_control']>1
    with pytest.raises(ValueError):ms.MultipletSupport(2)


@pytest.mark.parametrize('k',[1,3,5])
def test_full_stress_is_constant_and_matches_even_background(k,report):
    b=next(r for r in report['backgrounds'] if r['degree']==k)
    assert b['individual_anisotropy_min']>1e-8
    assert b['coherent_failure_min']>.1
    assert b['inherited_stress_relative']<1e-10
    assert b['euclidean_stress_relative']<1e-10
    for row in b['radii']:
        assert row['stress_relative']<1e-10
        assert row['einstein_scaled']<1e-10
        assert row['KG_absolute']<1e-10
        assert row['wrong_lambda']>2.9
        assert row['wrong_amplitude']>.3
        assert row['even_control_density_difference']==row['even_control_Lambda_difference']==0
    for row in b['kinetic']:
        assert row['f_min']==pytest.approx(1-1/(2*(k+1)**2),abs=1e-12)
        assert row['K_min']>0
        assert row['eigenvalue_error']<1e-10


def test_minimality_uses_different_phases_and_not_zero_phase_anisotropy():
    m=ms.MultipletSupport(1)
    points=np.eye(4)
    C=np.diag([.1,.2,.3,.4])
    T0,_=ms.component_stress(m,np.pi/2,points,C)
    T1,_=ms.component_stress(m,0.,points,C)
    def tf(T):
        Q=T.sum(axis=1)[:,1:,1:].copy()
        return Q-np.trace(Q,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
    assert np.linalg.norm(tf(T0))<1e-14
    assert np.linalg.norm(tf(T1))>.01
    Y,dY,_=ms.basis_jets(1,np.array([[.5,.5,.5,.5]]))
    assert np.linalg.matrix_rank(dY[0])==3
    assert np.linalg.norm(Y[0]@dY[0])<1e-14


@pytest.mark.parametrize('k',[1,3,5])
def test_exact_certificate_independent_of_numeric_singular_cutoff(k,report):
    c=next(c for c in report['certificates'] if c['degree']==k)
    op=ms.polynomial_operators(k)
    for name,M in [('A',op['A']),('L',np.vstack([op['Z'],op['A']]))]:
        constrained=np.vstack([op['trace'],M])
        certificate=c['operators'][name]
        assert ms.verify_certificate(constrained,certificate)
        broken=copy.deepcopy(certificate); broken['rank']-=1
        assert not ms.verify_certificate(constrained,broken)
        if certificate['kernel_columns']:
            broken=copy.deepcopy(certificate); broken['kernel_columns'][0][0]+=1
            assert not ms.verify_certificate(constrained,broken)
            broken=copy.deepcopy(certificate); broken['kernel_columns'][0][0]=.5
            assert not ms.verify_certificate(constrained,broken)


def test_modular_rank_alone_does_not_certify_a_kernel():
    # Rank drops modulo the prime, but the rational matrix is invertible.
    A=np.diag([1,1000003])
    assert ms.modular_rank(A)==1
    assert not ms.verify_certificate(A,dict(rank=1,kernel_columns=[[0,1]]))


@pytest.mark.parametrize('k',[1,3,5])
def test_normalized_sensitivity_against_two_spatial_rules_and_component_stress(k,report):
    p=next(p for p in report['preparations'] if p['degree']==k)
    assert p['quadrature_norm_error']<1e-9
    assert p['refined_norm_error']<1e-9
    assert p['reconstruction_relative']<1e-9
    assert p['linearity_relative']<1e-9
    assert p['trace_direction']['full']==pytest.approx(np.sqrt(4/3))
    assert p['trace_direction']['anisotropy']<1e-12
    assert len([r for r in p['directions'] if r['direction'].startswith('random_')])==20
    for row in p['directions']:
        for step in row['steps']:
            assert step['gram_mismatch']==pytest.approx(abs(step['epsilon']),abs=1e-12)
            assert step['psd_min']>0
            assert step['field_reconstruction_relative']<1e-9
        if row['direction'].startswith('kernel_'):
            assert row['full_per_gram']<1e-10
    for d in p['diagonal']:
        assert d['gram_mismatch']/abs(d['epsilon'])==pytest.approx(2,rel=.02)


def test_internal_rotations_preserve_gram_but_coherent_sum_does_not():
    rng=np.random.default_rng(142)
    C=rng.normal(size=(4,4))*.1
    O=np.linalg.qr(rng.normal(size=(4,4)))[0]
    points=rng.normal(size=(7,4));points/=np.linalg.norm(points,axis=1)[:,None]
    model=ms.MultipletSupport(1)
    t,_=ms.component_stress(model,.31,points,C)
    rotated,_=ms.component_stress(model,.31,points,O@C)
    coherent,_=ms.component_stress(model,.31,points,C.sum(axis=0)[None,:])
    assert np.max(abs(t.sum(axis=1)-rotated.sum(axis=1)))<1e-14
    assert np.max(abs(t.sum(axis=1)-coherent.sum(axis=1)))>.01


def test_psd_gate_precedes_factorization_and_zero_field_kinetic_limit():
    with pytest.raises(ValueError): ms.factor_gram(np.diag([1.,-1.]))
    with pytest.raises(ValueError): ms.factor_gram(np.array([[1.,1.],[0.,1.]]))
    model=ms.MultipletSupport(1,kappa=.4)
    F,K=model.kinetic(np.zeros(4))
    assert F==2.5
    assert np.array_equal(K,np.eye(4))
    with pytest.raises(ValueError): model.kinetic(np.ones(4)*10)


@pytest.mark.parametrize('key',ms.REQUIRED_CHECKS)
@pytest.mark.parametrize('mode',['missing','failed'])
def test_every_gate_failure_replaces_stale_reports_and_invalidates_targets(report,tmp_path,monkeypatch,key,mode):
    broken=copy.deepcopy(report)
    if mode=='missing':del broken['checks'][key]
    else:broken['checks'][key]=False
    monkeypatch.setattr(probe,'run_probe',lambda **kwargs:broken)
    (tmp_path/'probe.json').write_text('{"stale":true}')
    (tmp_path/'probe.md').write_text('STALE_SUCCESS')
    assert probe.main(['--output-dir',str(tmp_path)])==1
    out=json.loads((tmp_path/'probe.json').read_text())
    assert key in out['verdict']['failed_checks']
    assert all(out['verdict'][target]=='UNRESOLVED' for target in ms.CHECK_TARGETS[key])
    assert 'STALE_SUCCESS' not in (tmp_path/'probe.md').read_text()


@pytest.mark.parametrize('checks',[{}, {'unrelated':True}])
def test_empty_or_unrelated_gate_sets_fail(checks):
    v=ms.verdict(checks)
    assert all(v[key]=='UNRESOLVED' for key in ms.PHYSICAL_VERDICTS)


def test_missing_evidence_never_becomes_an_affirmative_kernel():
    v=ms.verdict(dict.fromkeys(ms.REQUIRED_CHECKS,True))
    assert v['full_preparation_kernel']=='UNRESOLVED'
    assert 'kernel_evidence' in v['failed_checks']


def test_partial_target_map_cannot_leave_untargeted_success(monkeypatch):
    monkeypatch.setattr(ms,'CHECK_TARGETS',{**ms.CHECK_TARGETS,'kinetic_matrix':()})
    v=ms.verdict(dict.fromkeys(ms.REQUIRED_CHECKS,True),{k:dict(certified=True,nullity=0) for k in (1,3,5)})
    assert all(v[key]=='UNRESOLVED' for key in ms.PHYSICAL_VERDICTS)


def test_computation_exception_also_overwrites_stale_success(tmp_path,monkeypatch):
    def broken(**kwargs):raise ArithmeticError('deliberate certificate failure')
    monkeypatch.setattr(probe,'run_probe',broken)
    (tmp_path/'probe.md').write_text('STALE_SUCCESS')
    assert probe.main(['--output-dir',str(tmp_path)])==1
    r=json.loads((tmp_path/'probe.json').read_text())
    assert all(r['verdict'][key]=='UNRESOLVED' for key in ms.PHYSICAL_VERDICTS)
    assert 'ArithmeticError' in r['error']
    assert 'STALE_SUCCESS' not in (tmp_path/'probe.md').read_text()


def test_post_freeze_divergence_identity_explains_kernel_equality(report):
    for c in report['post_freeze_divergence']:
        assert c['residuals']==[0,0,0]
        assert c['minimum_spectral_gap']>0
        assert c['post_freeze']
    for c in report['certificates']:
        assert c['operators']['A']['nullity']==c['operators']['L']['nullity']
