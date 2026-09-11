"""Physical constraint completion, cover/restricted dipoles, and report gates."""
import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm, null_space

from geometrodynamics.waves import multiplet_scalar_stability as m
from experiments.closure_ledger import multiplet_scalar_stability_probe as probe


@pytest.fixture(scope='module')
def report():return probe.finalize(probe.run_probe())


def test_frozen_gates_and_exact_identities(report):
    freeze=(Path(__file__).parents[1]/'docs/multiplet_scalar_stability_prereg.md').read_text()
    block=freeze.split('Required gate names are fixed:')[1].split('All ten gates')[0]
    names=tuple(x.strip().rstrip('.') for x in block.split(','))
    assert names==m.REQUIRED_CHECKS
    assert report['checks_passed'],report['checks']
    assert set(report['exact']['residuals'])==set(m.EXACT_IDENTITIES)
    assert set(report['exact']['residuals'].values())=={'0'}
    assert report['homogeneous_evidence_passed']


def test_physical_growth_has_an_exact_constraint_satisfying_family(report):
    assert len(report['continuations'])==72
    assert {r['sign'] for r in report['continuations']}=={-1.,1.}
    assert {r['epsilon'] for r in report['continuations']}=={-.01,.01}
    for row in report['continuations']:
        assert max(row['full_residual'],row['KG'],abs(row['constraint']))<1e-9
        assert row['f']>0
    # Exact constraint branch has the predicted growing tangent, including
    # the nonzero velocity required beyond the degenerate linear constraint.
    tau=.4;eps=1e-5
    plus=m.exact_scale(tau,eps);minus=m.exact_scale(tau,-eps)
    tangent=(plus-minus)/(2*eps)
    assert np.allclose(tangent,np.array([1,math.sqrt(2),2])*np.exp(math.sqrt(2)*tau),rtol=1e-8)


def test_off_shell_full_geometry_and_finite_variations(report):
    assert len(report['frw_off_shell'])==10
    assert max(max(r.values()) for r in report['frw_off_shell'])<1e-9
    for name in ('homogeneous','dipole','gauges','cover'):
        assert report[name]
        for r in report[name]:
            assert max(r['steps'][-1]['errors'].values())<2e-4
            assert all(3.5<x<4.5 for x in r['ratios'])


@pytest.mark.parametrize('phase',[0.,.31,math.pi/4])
def test_two_clocks_agree_and_missing_scalar_phase_conversion_is_detected(report,phase):
    rows=[r for r in report['clocks'] if r['phase']==phase]
    assert len(rows)==3
    assert rows[-1]['error']<2e-4
    assert rows[0]['error']/rows[1]['error']==pytest.approx(4,rel=.01)
    assert rows[1]['error']/rows[2]['error']==pytest.approx(4,rel=.01)
    assert rows[-1]['omitted_clock_term']>.05


def test_full_field_period_and_half_period_are_not_confused(report):
    h=np.asarray(report['period_evidence']['homogeneous'][-1])
    half=np.asarray(report['half_period'])
    assert np.allclose(h,half@half,rtol=1e-10,atol=1e-10)
    assert np.linalg.det(h)==pytest.approx(1,abs=1e-7)
    assert np.trace(h)==pytest.approx(2*math.cosh(math.sqrt(2)*math.pi),rel=1e-12)
    assert sorted(np.linalg.eigvals(h))==pytest.approx(sorted(np.exp([-math.sqrt(2)*math.pi,math.sqrt(2)*math.pi])),rel=1e-10)
    assert np.allclose(m.field_jets(.37+math.pi),m.field_jets(.37))
    assert np.allclose(m.field_jets(.37+math.pi/2),-m.field_jets(.37))


@pytest.mark.parametrize('tau',[0.,math.pi/4,.31,math.pi/2])
def test_dipole_constraint_rank_never_divides_by_field_or_velocity_zero(tau):
    C=m.constraints(tau)
    assert np.linalg.matrix_rank(C)==2
    # Check propagation directly by differentiating C, separately from the
    # exact symbolic coefficient certificate.
    step=1e-5;Cp=(m.constraints(tau+step)-m.constraints(tau-step))/(2*step)
    assert np.allclose(Cp+C@m.DIPOLE,np.array([[0.,-3.],[1/3,0.]])@C,atol=1e-8)


def test_dipole_physical_subspace_is_stress_free_and_not_a_jordan_block(report):
    assert len(report['cover'])==96
    for row in report['cover']:
        assert row['constraint']<1e-9
        assert max(abs(x) for x in row['stress_coefficients'].values())<1e-9
    M=expm(m.DIPOLE*math.pi)
    assert np.allclose(M,-np.eye(4),atol=1e-12)
    # At full period the same constraint subspace returns. A multiplier -1
    # alone would not exclude a Jordan block; the full matrix does.
    for phase in (0.,.31,math.pi/4):
        basis=null_space(m.constraints(0.,phase))
        assert np.allclose(M@basis,-basis,atol=1e-12)


def test_parity_restriction_excludes_nonzero_even_dipoles_without_counting_gauge(report):
    assert max(r['norm'] for r in report['parity'])>.1
    assert all(r['even_error']==0 for r in report['parity'])
    assert max(r['metric_even_defect'] for r in report['parity'])>.1
    assert report['controls']['gauge_without_metric']>1e-3
    assert report['verdict']['dipole_cover_block']=='CONSTRAINED_LINEAR_NEUTRAL_SEMISIMPLE'
    assert report['verdict']['dipole_restricted_admissibility']=='EXCLUDED_UNDER_STATED_ANTIPODAL_RESTRICTIONS'


@pytest.mark.parametrize('key',m.REQUIRED_CHECKS)
@pytest.mark.parametrize('missing',[False,True])
def test_every_missing_or_failed_gate_starts_from_an_affirmative_report(report,key,missing):
    checks=report['checks'].copy()
    assert m.verdict(checks,report['period_evidence'])['homogeneous_physical_block']=='HYPERBOLIC_GROWING_MODE'
    if missing:checks.pop(key)
    else:checks[key]=False
    v=m.verdict(checks,report['period_evidence'])
    assert all(v[k]=='UNRESOLVED' for k in ('homogeneous_physical_block','homogeneous_constraint_completion','dipole_cover_block','dipole_restricted_admissibility'))


@pytest.mark.parametrize('kind',['missing','nan','shape','wrong_map','disagree','unknown_gate'])
def test_evidence_cannot_be_replaced_with_true_gate_flags(report,kind):
    checks=report['checks'].copy();evidence=copy.deepcopy(report['period_evidence'])
    if kind=='missing':evidence.pop('dipole')
    if kind=='nan':evidence['homogeneous'][0][0][0]=float('nan')
    if kind=='shape':evidence['homogeneous']=[]
    if kind=='wrong_map':evidence['dipole']=[np.eye(2).tolist()]*2
    if kind=='disagree':evidence['homogeneous'][0][0][0]+=.1
    if kind=='unknown_gate':checks['unregistered']=True
    assert m.verdict(checks,evidence)['homogeneous_physical_block']=='UNRESOLVED'


def test_failed_dipole_does_not_erase_independent_homogeneous_evidence(report):
    bad=copy.deepcopy(report);bad['period_evidence']['dipole']=[np.eye(2).tolist()]*2
    result=probe.finalize(bad)
    assert not result['checks_passed']
    assert result['homogeneous_evidence_passed']
    assert result['verdict']['homogeneous_physical_block']=='UNRESOLVED'


def test_missing_exact_identity_is_not_vacuously_certified(report):
    bad=copy.deepcopy(report);bad['exact']['residuals'].pop('density')
    assert not probe.finalize(bad)['checks']['field_reduction']


@pytest.mark.parametrize('kind',['exception','missing_section','nonfinite','empty_geometry','bad_period'])
def test_cli_failure_does_not_leave_a_success_artifact(report,tmp_path,monkeypatch,kind):
    bad=copy.deepcopy(report)
    if kind=='missing_section':bad.pop('continuations')
    if kind=='nonfinite':bad['period_evidence']['homogeneous'][0][0][0]=float('nan')
    if kind=='empty_geometry':bad['dipole']=[]
    if kind=='bad_period':bad['period_evidence']['dipole']=[np.eye(2).tolist()]*2
    def failed():
        if kind=='exception':raise ArithmeticError('injected independent-route failure')
        return bad
    monkeypatch.setattr(probe,'run_probe',failed)
    assert probe.main(['--output',str(tmp_path)])==1
    result=json.loads((tmp_path/'probe.json').read_text())
    assert not result['checks_passed']
    assert result['verdict']['homogeneous_physical_block']=='UNRESOLVED'
    assert 'UNRESOLVED' in (tmp_path/'probe.md').read_text()
