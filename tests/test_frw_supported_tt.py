"""Independent full-field closure, both clocks, and fail-closed FRW transport."""
import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.waves import frw_supported_tt as f
from geometrodynamics.waves import reciprocal_scalar_tt as rt
from experiments.closure_ledger import frw_supported_tt_probe as probe


@pytest.fixture(scope='module')
def report():return probe.finalize(probe.run_probe())


def test_frozen_gate_schema_and_exact_action(report):
    freeze=(Path(__file__).parents[1]/'docs/frw_supported_tt_prereg.md').read_text()
    block=freeze.split('Required gate names:')[1].split('Machine verdicts:')[0]
    names=tuple(x.strip().rstrip('.') for x in block.split(','))
    assert names==f.REQUIRED_CHECKS
    assert report['checks_passed'],report['checks']
    assert set(report['exact']['residuals'])==set(f.IDENTITIES)
    assert set(report['exact']['residuals'].values())=={'0'}
    assert report['exact']['bad_linear_lambda_identity']=='0'
    assert report['exact']['bad_linear_lambda_value']==2


@pytest.mark.parametrize('basis',range(5))
@pytest.mark.parametrize('jet',range(3))
def test_all_five_components_and_independent_time_coefficients(report,basis,jet):
    rows=[r for r in report['geometry'] if r['case']==f'basis_{basis}_jet_{jet}']
    assert len(rows)==1
    row=rows[0]
    assert max(row['steps'][-1]['errors'].values())<2e-4
    assert all(3.5<x<4.5 for x in row['ratios'])


def test_nonlinear_background_departures_and_complete_constraints(report):
    assert len(report['geometry'])==55
    on=[r for r in report['geometry'] if r['case']=='on_equation']
    assert len(on)==18
    for row in report['geometry']:
        assert row['steps'][-1]['constraint']<2e-4
        assert max(row['steps'][-1]['errors'].values())<2e-4
    for row in on:
        assert row['expected_residual_norm']<1e-10
        assert row['steps'][-1]['errors']['conformal_residual']<2e-4
    assert max(r['commutator'] for r in report['geometry'] if r['case']=='noncommuting')>.01
    assert {r['model']['departure'] for r in report['geometry']}=={.1,.3,-.1}
    for bg in report['backgrounds']:
        assert max(bg['residual'],bg['KG'],bg['density_error'],bg['parity_error'])<1e-9


def test_physical_field_velocity_zero_is_distinct_from_conformal_qprime_zero(report):
    rows=list(zip(report['geometry'],report['backgrounds']))
    physical=[bg for row,bg in rows if row['case']=='physical_velocity_zero']
    assert len(physical)==2
    assert all(abs(bg['field_velocity'])<1e-12 for bg in physical)
    assert all(abs(bg['qp'])>.01 for bg in physical)
    qzero=[bg for row,bg in rows if row['case']=='conformal_velocity_zero'][0]
    assert abs(qzero['qp'])<1e-12
    assert abs(qzero['field_velocity'])>.01


def test_field_zero_decouples_instantaneous_stress_not_the_normal_potential():
    model=f.FRWSupport(departure=.3);eta=math.pi/4
    A,Ap,App,q,qp,qpp=model.jets(eta)
    assert abs(q)<1e-12
    b=.2*rt.STF_BASIS[0];v=.3*rt.STF_BASIS[1]
    bare=-2*Ap*v/A-8*b
    expected=f.expected_response(model,eta,b,v,bare)
    assert np.linalg.norm(expected['stress'])<1e-12
    assert np.linalg.norm(expected['residual'])<1e-12
    # q'^2 still enters m'': sqrt(m) normalization carries a second derivative.
    offset=model.normal_potential(eta)-(8-App/A)
    assert offset==pytest.approx(model.kappa*qp*qp/(6*A*A),abs=1e-12)
    assert offset>.01


def test_both_clocks_and_all_transport_backgrounds(report):
    assert len(report['map_evidence'])==108
    for row in report['geometry']:
        assert max(s['clock_error'] for s in row['steps'])<1e-9
    for row in report['map_evidence']:
        matrices=np.asarray(row['maps'])
        assert row['eta_error']<1e-9 and row['scale_error']<1e-9
        assert row['kinetic_bound']>.2
        assert max(np.linalg.norm(x-matrices[0])/max(1.,np.linalg.norm(matrices[0])) for x in matrices[1:])<1e-8
        assert max(abs(np.linalg.det(x)-1) for x in matrices)<1e-8
        assert max(np.linalg.norm(x.T@f.J@x-f.J) for x in matrices)<1e-8
    assert len(report['composition'])==18
    assert all(max(r['composition_error'],r['inverse_error'])<1e-8 for r in report['composition'])


def test_static_control_recovers_supported_esu_and_limit(report):
    for row in report['static']:
        assert row['map_error']<1e-8 and row['coefficient_error']<1e-12
        assert row['trace']==pytest.approx(-.096306540195,abs=1e-10)
    for row in report['static_limit']:
        assert row['errors'][0]>row['errors'][1]>row['errors'][2]>0


def test_wrong_operators_separate_from_full_geometry_floor(report):
    assert all(v>1e-5 for v in report['controls'].values())
    good=max(r['steps'][-1]['errors']['conformal_residual'] for r in report['geometry'] if r['case']=='on_equation')
    assert min(report['controls'].values())>100*good


def test_regular_domain_is_enforced():
    with pytest.raises(ValueError):f.FRWSupport(departure=-1)
    with pytest.raises(ValueError):f.FRWSupport(departure=.3).jets(3.)
    with pytest.raises(ValueError):f.FRWSupport(departure=-.9).coefficients(0.)
    with pytest.raises(ValueError):f.evolve(f.FRWSupport(departure=-.9))
    with pytest.raises(ValueError):f.proper_evolve(f.FRWSupport(departure=-.9))


@pytest.mark.parametrize('key',f.REQUIRED_CHECKS)
@pytest.mark.parametrize('missing',[False,True])
def test_every_missing_or_failed_gate_rejects_an_initially_passing_report(report,key,missing):
    checks=report['checks'].copy()
    assert f.verdict(checks,report['map_evidence'])['supported_operator']=='FULLY_SUPPORTED_FRW_EQUATION_VERIFIED'
    if missing:checks.pop(key)
    else:checks[key]=False
    verdict=f.verdict(checks,report['map_evidence'])
    for name in ('linear_tensor_sector','supported_operator','clock_agreement','finite_interval_transport'):
        assert verdict[name]=='UNRESOLVED'


@pytest.mark.parametrize('kind',['missing','shape','nonfinite','different','determinant','endpoint','domain','unknown_gate'])
def test_map_evidence_is_required_and_checked(report,kind):
    checks=report['checks'].copy();evidence=copy.deepcopy(report['map_evidence'])
    if kind=='missing':evidence=[]
    if kind=='shape':evidence[0]['maps']=[]
    if kind=='nonfinite':evidence[0]['maps'][0][0][0]=float('nan')
    if kind=='different':evidence[0]['maps'][0][0][0]+=.1
    if kind=='determinant':evidence[0]['maps']=(np.asarray(evidence[0]['maps'])*2).tolist()
    if kind=='endpoint':evidence[0]['eta_error']=.1
    if kind=='domain':evidence[0]['kinetic_bound']=0.
    if kind=='unknown_gate':checks['not_frozen']=True
    assert f.verdict(checks,evidence)['supported_operator']=='UNRESOLVED'


def test_missing_identity_is_not_vacuously_certified(report):
    bad=copy.deepcopy(report);bad['exact']['residuals'].pop('action')
    assert not probe.finalize(bad)['checks']['action_derivation']


@pytest.mark.parametrize('kind',['exception','missing_geometry','empty_geometry','nonfinite','bad_map'])
def test_cli_fails_closed_and_overwrites_stale_success(report,tmp_path,monkeypatch,kind):
    bad=copy.deepcopy(report)
    if kind=='missing_geometry':bad.pop('geometry')
    if kind=='empty_geometry':bad['geometry']=[]
    if kind=='nonfinite':bad['controls']['bare_frw']=float('nan')
    if kind=='bad_map':bad['map_evidence'][0]['maps']=[]
    def run():
        if kind=='exception':raise ArithmeticError('injected full-field failure')
        return bad
    (tmp_path/'probe.json').write_text('{"checks_passed":true}')
    (tmp_path/'probe.md').write_text('stale success')
    monkeypatch.setattr(probe,'run_probe',run)
    assert probe.main(['--output-dir',str(tmp_path)])==1
    result=json.loads((tmp_path/'probe.json').read_text())
    assert result['checks_passed'] is False
    assert result['verdict']['supported_operator']=='UNRESOLVED'
    assert 'UNRESOLVED' in (tmp_path/'probe.md').read_text()
