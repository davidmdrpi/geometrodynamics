"""Independent Einstein/stress response and fail-closed tensor evolution."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.waves import coupled_multiplet_response as cm
from geometrodynamics.waves import reciprocal_scalar_tt as rt
from experiments.closure_ledger import coupled_multiplet_response_probe as probe


@pytest.fixture(scope='module')
def report():return probe.finalize(probe.run_probe())


def test_published_gate_names_and_all_checks(report):
    text=(Path(__file__).parents[1]/'docs/coupled_multiplet_response_prereg.md').read_text()
    names=tuple(line.split('|')[1].strip().strip('`') for line in text.splitlines() if line.startswith('| `'))
    assert names==cm.REQUIRED_CHECKS
    assert report['checks_passed'],report['verdict']
    assert report['exact']['all_zero']
    assert set(report['exact']['residuals'].values())=={'0'}


def test_background_four_fields_remain_odd_at_zeros_and_generic_phases(report):
    assert report['exact']['components']==4
    for row in report['backgrounds']:
        assert row['einstein_scaled']<1e-10
        assert row['KG_absolute']<1e-10
        assert row['odd_absolute']<1e-10
        assert row['density_error']<1e-10


@pytest.mark.parametrize('basis',range(5))
def test_every_tensor_component_matches_independent_full_geometry(report,basis):
    rows=[r for r in report['geometry'] if r['case']==f'basis_{basis}']
    assert len(rows)==4
    for row in rows:
        last=row['steps'][-1]
        assert last['G_absolute']/row['scales']['G']<2e-4
        assert last['T_absolute']/row['scales']['T']<2e-4
        assert last['constraints_absolute']<2e-4
        assert last['KG_absolute']<2e-4
        assert all(3.5<r['ratio']<4.5 for r in row['convergence'])


@pytest.mark.parametrize('case',['pure_beta','pure_velocity','pure_acceleration','generic','on_equation'])
def test_independent_coefficients_noncommuting_data_and_complete_equations(report,case):
    rows=[r for r in report['geometry'] if r['case']==case]
    assert rows
    for r in rows:
        assert r['steps'][-1]['residual_absolute']/r['scales']['residual']<2e-4
        assert r['steps'][-1]['scalar_curvature_absolute']<2e-4
    if case=='on_equation':assert rows[0]['steps'][-1]['complete_Einstein_response']<2e-4


def test_pressure_metric_term_is_removed_by_mixed_components():
    m=cm.CoupledSupport();E=rt.STF_BASIS[2]
    zero=np.zeros((3,3))
    # At field zero all remaining matter stress is isotropic even when the
    # finite metric is anisotropic. A covariant-component STF would be wrong.
    result=cm.full_geometry(m,np.pi/4,E,zero,zero,.01)
    assert np.linalg.norm(cm.stf(result['stress'][1:,1:]))<1e-12


def test_bare_oscillator_failure_is_not_a_damping_or_instability_test(report):
    c=report['controls']
    assert c['bare_turning_residual']==pytest.approx(1.5,abs=1e-12)
    assert c['reversed_fprime_residual']>.1
    assert c['dropped_metric_response_residual']>.1
    assert c['bare_history_max']>1
    v=report['verdict']
    assert v['bare_frequency_transfer']=='BARE_EQUATION_FAILS_FOR_THIS_SUPPORT'
    assert v['full_dynamical_stability']=='NOT_ESTABLISHED'


def test_periodic_hamiltonian_and_normal_form_are_independent_checks(report):
    e=report['evolution']
    assert e['refinement']<1e-8
    assert e['symplectic_error']<1e-8
    assert abs(e['determinant']-1)<1e-8
    assert e['normal_form_error']<1e-8
    assert e['twenty_period_error']<1e-8
    assert e['five_component_error']<1e-8
    assert e['bare_control_error']<1e-8
    assert all(r['error']<1e-8 for r in e['physical_time'])
    assert all(r['trace_difference']<1e-8 for r in e['phase_controls'])
    assert abs(e['trace']-e['bare_trace'])>.01


def test_time_coefficients_and_canonical_generator_have_correct_units():
    for a in (.7,1.,2.):
        m=cm.CoupledSupport(a,.4)
        tau=.29
        f,fp,fpp,g=m.coefficients(tau)
        P,Pd,_=m.field_jets(a*tau)
        assert f==pytest.approx(1-m.kappa*P*P/6)
        assert fp/a==pytest.approx(-m.kappa*P*Pd/3)
        assert g==pytest.approx(8+2*m.kappa*P*P/3)
        A=m.generator(tau);J=np.array([[0.,1.],[-1.,0.]])
        assert np.linalg.norm(A.T@J+J@A)<1e-15
        assert np.linalg.norm(m.generator(tau,eta=0)-np.array([[0.,1.],[-8.,0.]]))==0


@pytest.mark.parametrize('traces,expected',[
    ([0.,0.],'NUMERICALLY_ELLIPTIC'),([3.,3.],'NUMERICALLY_HYPERBOLIC'),
    ([2.,2.],'UNRESOLVED'),([-2.,-2.],'UNRESOLVED'),
    ([2-1e-8,2-1e-8],'UNRESOLVED'),([float('nan'),0.],'UNRESOLVED'),
    ([float('inf'),float('inf')],'UNRESOLVED'),([0.,.001],'UNRESOLVED')])
def test_floquet_classification_has_a_band_edge_and_convergence_gate(traces,expected):
    assert cm.classify_period(traces)==expected


@pytest.mark.parametrize('key',cm.REQUIRED_CHECKS)
@pytest.mark.parametrize('mode',['missing','failed'])
def test_each_frozen_gate_invalidates_all_physical_verdicts_and_overwrites(report,tmp_path,monkeypatch,key,mode):
    broken=copy.deepcopy(report)
    if mode=='missing':del broken['checks'][key]
    else:broken['checks'][key]=False
    monkeypatch.setattr(probe,'run_probe',lambda:broken)
    (tmp_path/'probe.json').write_text('{"stale":true}')
    (tmp_path/'probe.md').write_text('STALE_SUCCESS')
    assert probe.main(['--output-dir',str(tmp_path)])==1
    out=json.loads((tmp_path/'probe.json').read_text())
    assert key in out['verdict']['failed_checks']
    assert all(out['verdict'][k]=='UNRESOLVED' for k in cm.PHYSICAL)
    assert {k:out['verdict'][k] for k in cm.SCOPE}==cm.SCOPE
    assert 'STALE_SUCCESS' not in (tmp_path/'probe.md').read_text()


@pytest.mark.parametrize('damage',['absent','nan','inconsistent_trace','wrong_determinant','wrong_shape','missing_geometry','malformed_evolution'])
def test_period_evidence_is_validated_independently_of_success_flags(report,tmp_path,monkeypatch,damage):
    broken=copy.deepcopy(report)
    if damage=='missing_geometry':broken.pop('geometry')
    elif damage=='malformed_evolution':broken['evolution']={}
    elif damage=='absent':broken.pop('period_evidence')
    elif damage=='nan':broken['period_evidence']['maps'][0][0][0]=float('nan')
    elif damage=='inconsistent_trace':broken['period_evidence']['traces'][0]+=.5
    elif damage=='wrong_determinant':broken['period_evidence']['maps'][0][0][0]+=1
    else:broken['period_evidence']['maps']=[[1.,0.],[0.,1.]]
    monkeypatch.setattr(probe,'run_probe',lambda:broken)
    (tmp_path/'probe.md').write_text('STALE_SUCCESS')
    assert probe.main(['--output-dir',str(tmp_path)])==1
    out=json.loads((tmp_path/'probe.json').read_text())
    assert all(out['verdict'][k]=='UNRESOLVED' for k in cm.PHYSICAL)
    assert 'STALE_SUCCESS' not in (tmp_path/'probe.md').read_text()


def test_unknown_checks_and_exceptions_cannot_reuse_success(tmp_path,monkeypatch):
    v=cm.verdict({'unrelated':True})
    assert all(v[k]=='UNRESOLVED' for k in cm.PHYSICAL)
    def bad():raise ArithmeticError('deliberate geometry failure')
    monkeypatch.setattr(probe,'run_probe',bad)
    (tmp_path/'probe.md').write_text('STALE_SUCCESS')
    assert probe.main(['--output-dir',str(tmp_path)])==1
    assert 'STALE_SUCCESS' not in (tmp_path/'probe.md').read_text()


def test_no_general_stability_preparation_or_causality_promotion(report):
    assert {k:report['verdict'][k] for k in cm.SCOPE}==cm.SCOPE
    assert 'modulo-4' in report['evolution']['quasifrequency']
