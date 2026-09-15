"""Regression tests for scientific false positives in the frozen experiment."""
import copy
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from geometrodynamics.waves import mouth_momentum as m
from experiments.closure_ledger import mouth_momentum_probe as p

ARCHIVE=Path(__file__).parents[1]/'experiments/closure_ledger/runs/20260915_mouth_momentum/probe.json.gz'


@pytest.fixture(scope='module')
def archive():
    return json.loads(gzip.decompress(ARCHIVE.read_bytes()))


def test_frozen_archive_passes_without_importing_quantum_claims(archive):
    result=p.score(archive)
    assert result['passed']==8
    assert all(result['verdicts'].values())
    assert all(result['unestablished'].values())
    assert result['metrics']['bad_seam_max']>.1
    assert result['metrics']['physical_sample_max_K_norm']>1


@pytest.mark.parametrize('damage',[
    'missing_solution','duplicate_solution','nan_solution','changed_solution',
    'removed_K','altered_K','altered_R','altered_momentum','changed_rhs',
    'missing_physical_point','duplicate_physical_point','wrong_coordinate_step',
    'missing_gate','unknown_gate','missing_control','forged_seam','changed_section',
    'changed_parameter','changed_iterations','changed_prereg',
])
def test_damaged_evidence_withdraws_both_claims(archive,damage):
    data=copy.deepcopy(archive)
    case=data['physical'][-1]['cases'][0]
    if damage=='missing_solution':data['solutions'].pop()
    elif damage=='duplicate_solution':data['solutions'][-1]=copy.deepcopy(data['solutions'][-2])
    elif damage=='nan_solution':data['solutions'][-1]['psi'][0][0]=float('nan')
    elif damage=='changed_solution':data['solutions'][-1]['psi'][0][0]+=.001
    elif damage=='removed_K':case.pop('extrinsic_curvature')
    elif damage=='altered_K':case['extrinsic_curvature'][0][0]+=.01
    elif damage=='altered_R':case['R']+=.01
    elif damage=='altered_momentum':case['momentum'][0]+=.01
    elif damage=='changed_rhs':data['momentum_fd'][0]['rhs'][0]+=.01
    elif damage=='missing_physical_point':data['physical'][-1]['cases'].pop()
    elif damage=='duplicate_physical_point':data['physical'][-1]['cases'][0]=copy.deepcopy(data['physical'][-1]['cases'][1])
    elif damage=='wrong_coordinate_step':case['h']=.01
    elif damage=='missing_gate':data['gate_schema'].pop()
    elif damage=='unknown_gate':data['gate_schema'].append('invented_success')
    elif damage=='missing_control':data['controls'].pop('untwisted')
    elif damage=='forged_seam':data['bad_seam']['value_residuals']=(np.zeros((30,3,3))).tolist()
    elif damage=='changed_section':data['sections'][-1]['values'][0]['area']-=1
    elif damage=='changed_parameter':data['solutions'][-1]['C']=.4
    elif damage=='changed_iterations':data['solutions'][-1]['iterations']=29
    elif damage=='changed_prereg':data['prereg']='not-the-freeze'
    result=p.score(data)
    assert not any(result['verdicts'].values()),damage


def test_real_omitted_momentum_correction_fails_physical_divergence(archive):
    r=next(x for x in archive['solutions'] if x['epsilon']==.1 and x['ns']==40)
    f=m.Interpolant(r)
    x=(.77,.83,.37)
    correct=m.coordinate_constraints(f,x,2.5e-4)
    missing=m.coordinate_constraints(f,x,2.5e-4,correction=0)
    damaged=m.coordinate_constraints(f,x,2.5e-4,correction=1.1)
    assert correct['momentum_normalized']<1e-7
    assert missing['momentum_normalized']>1e-4
    assert damaged['momentum_normalized']>1e-4


def test_coordinate_checker_on_exact_nonzero_K_product_geometry(archive):
    r=next(x for x in archive['solutions'] if x['epsilon']==0 and x['ns']==24)
    f=m.Interpolant(r)
    c=m.coordinate_constraints(f,(.7,1.2,.37),5e-4)
    assert abs(c['R']-2)<1e-6
    assert abs(c['K2']-1.5)<1e-13
    assert abs(c['trace_K'])<1e-13
    assert c['momentum_normalized']<1e-7


def test_antipodal_tensor_pullback_is_not_scalar_periodicity():
    A=m.conformal_tensor(.7,.8,.1)
    B=m.conformal_tensor(.7+2*np.pi,np.pi-.8,.1)
    J=np.diag([1.,-1.,1.])
    assert np.max(abs(A-J@B@J))<1e-14
    assert np.max(abs(A-B))>1e-3


def test_symbolic_certificate_is_immutable():
    assert isinstance(m.symbolic_checks(),tuple)
    assert len(m.symbolic_checks())==14
    assert set(m.symbolic_checks())=={'0'}


def test_cli_parse_failure_replaces_stale_success(tmp_path):
    (tmp_path/'probe.md').write_text('PASSED ALL GATES')
    (tmp_path/'verdict.json').write_text('{"success":true}')
    bad=tmp_path/'bad.json';bad.write_text('{not json')
    result=subprocess.run([sys.executable,'-m','experiments.closure_ledger.mouth_momentum_probe',
        '--output-dir',str(tmp_path),'--rescore',str(bad)],capture_output=True,text=True,
        env={**os.environ,'OPENBLAS_NUM_THREADS':'1'})
    assert result.returncode!=0
    assert 'PASSED ALL GATES' not in (tmp_path/'probe.md').read_text()
    assert not any(json.loads((tmp_path/'verdict.json').read_text())['verdicts'].values())


def test_cli_damaged_raw_evidence_returns_failure(tmp_path,archive):
    data=copy.deepcopy(archive);data['solutions'].pop()
    source=tmp_path/'damaged.json';source.write_text(json.dumps(data))
    result=subprocess.run([sys.executable,'-m','experiments.closure_ledger.mouth_momentum_probe',
        '--output-dir',str(tmp_path),'--rescore',str(source)],capture_output=True,text=True,
        env={**os.environ,'OPENBLAS_NUM_THREADS':'1'})
    assert result.returncode==1
    assert not any(json.loads((tmp_path/'verdict.json').read_text())['verdicts'].values())
