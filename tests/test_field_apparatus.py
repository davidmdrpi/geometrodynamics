"""Physical normalization, frozen evidence and scoped failure regressions."""
import copy
import gzip
import json
from pathlib import Path
import numpy as np
import pytest
from geometrodynamics.waves import field_apparatus as f
from experiments.closure_ledger import field_apparatus_probe as p

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/'experiments/closure_ledger/runs/20260914_field_apparatus'


@pytest.fixture(scope='module')
def report():
    return json.loads(gzip.decompress((RUN/'probe.json.gz').read_bytes()))


def test_frozen_result_recomputes_from_raw_evidence(report):
    assert report['prereg']==f.PREREG
    assert p.evidence_gates(report)==report['checks']
    assert p.verdict(report['checks'],report)==report['verdict']
    assert report['verification_passed']
    assert [k for k,v in report['checks'].items() if not v]==[
        'map_and_interface_equations','reciprocal_action_and_full_residuals']
    assert report['verdict']['field_generated_amplitude'] is None
    assert report['verdict']['field_response']=='BLOCKED_BY_UNSPECIFIED_MAP'


def test_kinetic_weight_cancels_nonzero_stiffness_in_the_same_degenerate_pair():
    c=f.circle(1.,.1,2,0.,128)
    K,W=f.decode(c['K1']),f.decode(c['W1'])
    assert abs(K[0,1]-.5)<1e-14
    assert abs(W[0,1]-.5)<1e-14
    assert np.linalg.norm(K-W)<1e-14
    # Missing W predicts a spurious first-order eigenvalue splitting of one.
    assert np.ptp(np.linalg.eigvalsh(K))==pytest.approx(1.)
    assert np.ptp(np.linalg.eigvalsh(K-W))<1e-14


def test_prescribed_graph_and_potential_are_sensitive_controls():
    D=f.graph(0.,0.,2,derivative=True)
    off=f.graph_block(D)[0,1]
    assert off.real==pytest.approx(2-np.sqrt(2)) and abs(off.imag)<1e-14
    assert abs(f.potential(.1,2,0.,128)[0,1]-.05)<1e-14
    for m in (1,3):
        assert abs(f.graph_block(f.graph(0.,0.,m,derivative=True))[0,1])<1e-14
        assert abs(f.potential(.1,m,0.,128)[0,1])<1e-14


def test_selected_component_is_not_the_complete_quartet():
    c=f.quartet(np.array([1.,0.,0.,0.]),1.,[0,1],128)
    assert np.max(abs(f.decode(c['fourier'])))<1e-14
    chi=f.grid(128)
    assert np.mean(np.cos(chi)**2*np.exp(-2j*chi))==pytest.approx(.25)


def test_exact_certificates_cannot_be_changed_through_a_cached_return():
    c=f.certificate();c['circle_blocks'][1]['generalized'][1]='1'
    assert f.certificate()['circle_blocks'][1]['generalized'][1]=='0'


@pytest.mark.parametrize('kind,gate',[
    ('missing_W','kinetic_measure_and_coordinate_control'),
    ('circle_matrix','kinetic_measure_and_coordinate_control'),
    ('circle_nonfinite','kinetic_measure_and_coordinate_control'),
    ('symbolic','quartet_identity'),('graph_matrix','graph_control'),
    ('potential_nonfinite','physical_potential_control'),
    ('inventory','source_inventory'),('forged_map','source_inventory'),
    ('duplicate_quartet','quartet_identity'),('missing_circle','kinetic_measure_and_coordinate_control')])
def test_raw_tampering_fails_its_gate_preserving_other_results(report,kind,gate):
    r=copy.deepcopy(report)
    if kind=='missing_W':r['circle'][0].pop('W1')
    elif kind=='circle_matrix':r['circle'][0]['K1'][0][1][0]+=.01
    elif kind=='circle_nonfinite':r['circle'][0]['lhs'][0][0][0]=float('nan')
    elif kind=='symbolic':r['exact']['intensity']='1'
    elif kind=='graph_matrix':r['graph'][0]['levels'][0]['plus'][0][0]+=.01
    elif kind=='potential_nonfinite':r['potential'][0]['block'][0][0][0]=float('nan')
    elif kind=='inventory':r['source_audit']['rows'].pop()
    elif kind=='forged_map':r['source_audit']['map']='graph_control_passed'
    elif kind=='duplicate_quartet':r['quartet'][1]=copy.deepcopy(r['quartet'][0])
    else:r.pop('circle')
    checks=p.evidence_gates(r)
    assert not checks[gate]
    assert [k for k in p.GATES if checks[k]!=report['checks'][k]]==[gate]
    out=p.verdict(report['checks'],r)
    if gate=='source_inventory':
        assert out['field_response']=='UNRESOLVED'
        assert out['operator_status']==report['verdict']['operator_status']
    else:
        assert out['bulk_mouth_map']=='BULK_MOUTH_MAP_UNSPECIFIED'
        assert out['field_response']=='BLOCKED_BY_UNSPECIFIED_MAP'


@pytest.mark.parametrize('gate',p.GATES)
def test_missing_gate_never_creates_a_physical_response(report,gate):
    checks=report['checks'].copy();checks.pop(gate)
    v=p.verdict(checks,report)
    assert v['physical_mouth_matrix_element'] is None
    assert v['field_response'] in ('UNRESOLVED','BLOCKED_BY_UNSPECIFIED_MAP')
    if gate in ('scope_and_provenance','raw_evidence_and_failure_paths'):
        assert v['quartet_intensity']=='UNRESOLVED'
        assert set(v['operator_status'].values())=={'UNRESOLVED'}


def test_forced_flags_unknown_gates_and_empty_evidence_fail_closed(report):
    checks={k:True for k in p.GATES}
    assert p.verdict(checks,report)['field_response']=='BLOCKED_BY_UNSPECIFIED_MAP'
    checks['extra_gate']=True
    assert p.verdict(checks,report)['quartet_intensity']=='UNRESOLVED'
    assert p.verdict({}, {})['bulk_mouth_map']=='UNRESOLVED'


@pytest.mark.parametrize('kind',['exception','failed_control','nonfinite'])
def test_cli_failure_overwrites_stale_success(report,tmp_path,monkeypatch,kind):
    if kind=='exception':
        def run():raise ArithmeticError('injected')
    else:
        damaged=copy.deepcopy(report)
        if kind=='nonfinite':damaged['circle'][0]['lhs'][0][0][0]=float('nan')
        else:damaged['graph'][0]['levels'][0]['plus'][0][0]+=.1
        def run():return damaged
    monkeypatch.setattr(p,'run_probe',run)
    (tmp_path/'probe.json.gz').write_bytes(gzip.compress(b'{"verification_passed":true}'))
    assert p.main(['--output-dir',str(tmp_path)])==1
    out=json.loads(gzip.decompress((tmp_path/'probe.json.gz').read_bytes()))
    assert not out['verification_passed']
    if kind=='nonfinite':
        assert out['verdict']['quartet_intensity']=='NO_ANGULAR_SOURCE_IN_EXACT_QUARTET'
        assert out['verdict']['operator_status']['intrinsic_circle']=='UNRESOLVED'
        assert out['nonfinite_evidence_paths']
