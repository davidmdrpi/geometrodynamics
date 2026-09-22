"""Protect the localized-mouth experiment's failed gates and physical sources."""
import ast
import copy
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pytest
from geometrodynamics.waves import localized_mouth as m
from experiments.closure_ledger import localized_mouth_probe as p
from experiments.closure_ledger import localized_mouth_refinement_probe as r
from experiments.closure_ledger.evidence_archive import read_bytes

ROOT=Path(__file__).parents[1]
RUN=ROOT/'experiments/closure_ledger/runs/20260916_localized_mouth'


@pytest.mark.parametrize('name',['probe.json.gz','refinement.json.gz','bernstein_refinement.json.gz'])
def test_lossless_parts_preserve_original_archives(name):
    manifest=json.loads((RUN/'source_manifest.json').read_text())['files']
    expected=manifest[str((RUN/name).relative_to(ROOT))]
    raw=read_bytes(RUN/name)
    assert hashlib.sha256(raw).hexdigest()==expected['sha256']
    assert hashlib.sha256(gzip.decompress(raw)).hexdigest()==expected['decompressed_sha256']


@pytest.mark.parametrize('damage',['part','missing','order','path'])
def test_archive_parts_reject_damage(tmp_path,damage):
    archive=tmp_path/'sample.gz'
    parts=[]
    for i,data in enumerate((b'first',b'second')):
        name=f'sample.gz.part{i:03d}';(tmp_path/name).write_bytes(data)
        parts.append(dict(name=name,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()))
    manifest=dict(bytes=11,sha256=hashlib.sha256(b'firstsecond').hexdigest(),parts=parts)
    if damage=='part':(tmp_path/parts[0]['name']).write_bytes(b'wrong')
    elif damage=='missing':(tmp_path/parts[0]['name']).unlink()
    elif damage=='order':parts.reverse()
    elif damage=='path':parts[0]['name']='../outside'
    archive.with_name('sample.gz.parts.json').write_text(json.dumps(manifest))
    with pytest.raises((ValueError,FileNotFoundError)):read_bytes(archive)


@pytest.fixture(scope='module')
def original():return r.read_original(RUN/'probe.json.gz')


@pytest.fixture(scope='module')
def refined():return json.loads(gzip.decompress(read_bytes(RUN/'refinement.json.gz')))


def test_original_failure_remains_six_of_eight(original):
    result=p.score(original)
    assert result['passed']==6
    assert result['gates']==r.ORIGINAL_GATES
    assert not any(result['verdicts'].values())


def test_refinement_absolute_accuracy_is_not_a_convergence_pass(original,refined):
    result=r.score(refined,original)
    assert result['passed']==7
    assert result['gates']['evidence']
    assert not result['gates']['physical']
    assert result['metrics']['H_offgrid_max']<1e-7
    assert result['metrics']['physical_H'][-1]<1e-5
    assert result['metrics']['physical_ratios'][0]<8
    assert not any(result['verdicts'].values())


@pytest.mark.parametrize('damage',[
    'missing_case','duplicate_case','nonfinite','changed_polynomial','wrong_frame_Pi',
    'removed_Pi','removed_K','changed_K','changed_source','changed_gradient',
    'changed_curvature','changed_momentum','changed_profile','missing_seam',
    'missing_gate','unknown_gate','changed_freeze','changed_original_hash',
    'changed_localization','changed_physical_step',
])
def test_corrupted_refinement_evidence_is_rejected(original,refined,damage):
    data=copy.deepcopy(refined);case=data['physical'][-1]['cases'][0]
    if damage=='missing_case':data['solutions'].pop()
    elif damage=='duplicate_case':data['solutions'][-1]=copy.deepcopy(data['solutions'][-2])
    elif damage=='nonfinite':case['Pi_J'][0]=float('nan')
    elif damage=='changed_polynomial':data['solutions'][-1]['solution']['c'][-1][0][0]+=.001
    elif damage=='wrong_frame_Pi':case['Pi_J']=(np.array(case['Pi_J'])/np.sqrt(m.F)).tolist()
    elif damage=='removed_Pi':case.pop('Pi_J')
    elif damage=='removed_K':case['geometry'].pop('extrinsic_curvature')
    elif damage=='changed_K':case['geometry']['extrinsic_curvature'][0][0]+=.001
    elif damage=='changed_source':case['current'][0]+=.001
    elif damage=='changed_gradient':case['gradient'][0][0]+=.001
    elif damage=='changed_curvature':case['geometry']['R']+=.001
    elif damage=='changed_momentum':case['geometry']['momentum'][0]+=.001
    elif damage=='changed_profile':data['profiles']['3.5']['theta']['c'][-1][0]+=.001
    elif damage=='missing_seam':data['seams'].pop()
    elif damage=='missing_gate':data['gate_schema'].pop()
    elif damage=='unknown_gate':data['gate_schema'].append('invented')
    elif damage=='changed_freeze':data['refinement_freeze']='unknown'
    elif damage=='changed_original_hash':data['original_sha256']='wrong'
    elif damage=='changed_localization':data['diagnostics'][-1]['data']['radii'][-1]*=2
    elif damage=='changed_physical_step':data['physical'][-1]['h']=.0002
    result=r.score(data,original)
    assert not result['gates']['evidence'],damage
    assert not any(result['verdicts'].values())


def test_changed_original_cannot_become_the_registered_baseline(original,refined):
    changed=copy.deepcopy(original)
    changed['solutions'][0]['message']='changed evidence'
    assert not r.score(refined,changed)['gates']['evidence']


def test_quintic_preserves_knots_without_bernstein_cancellation(original):
    for before in original['solutions']:
        if before['n_initial']!=513:continue
        profile=original['profiles'][str(before['L'])]
        after=m.reconstruct(before,profile)
        a=m.Data(after,profile);b=m.Data(before,profile)
        assert np.max(abs(a.sol(b.sol.x)-b.sol(b.sol.x)))<1e-12
        assert after['solution']['axis']==1


def test_reconstruction_allows_roundoff_without_allowing_metadata_changes(original):
    before=next(row for row in original['solutions'] if row['n_initial']==513)
    expected=m.reconstruct(before,original['profiles'][str(before['L'])])
    actual=copy.deepcopy(expected)
    value=actual['solution']['c'][-1][0][0]
    actual['solution']['c'][-1][0][0]=float(np.nextafter(value,np.inf))
    assert actual!=expected
    agrees,used=r.reconstruction_agrees(actual,expected)
    assert agrees and 0<used<1
    actual['message']='changed provenance'
    assert not r.reconstruction_agrees(actual,expected)[0]


@pytest.mark.parametrize('damage',['mesh','axis','shape','velocity','curvature'])
def test_reconstruction_rejects_structural_and_derivative_damage(original,damage):
    before=next(row for row in original['solutions'] if row['n_initial']==513)
    expected=m.reconstruct(before,original['profiles'][str(before['L'])])
    actual=copy.deepcopy(expected);sol=actual['solution']
    if damage=='mesh':sol['x'][1]=float(np.nextafter(sol['x'][1],np.inf))
    elif damage=='axis':sol['axis']=0
    elif damage=='shape':sol['c'].pop()
    elif damage=='velocity':sol['c'][-1][0][1]+=1e-10
    elif damage=='curvature':
        # Tiny value perturbation, zero endpoint values/slopes, but curvature
        # well above arithmetic roundoff throughout the interval: a*t^2(1-t)^2.
        h=sol['x'][1]-sol['x'][0];amplitude=1e-8*h*h
        sol['c'][1][0][0]+=amplitude/h**4
        sol['c'][2][0][0]-=2*amplitude/h**3
        sol['c'][3][0][0]+=amplitude/h**2
        for i in range(5):sol['c'][i+1][0][1]=(5-i)*sol['c'][i][0][0]
        a=m.restore(sol);e=m.restore(expected['solution'])
        ends=np.array([sol['x'][0],sol['x'][1]])
        assert np.max(abs(a(ends)-e(ends)))<1e-12
    assert not r.reconstruction_agrees(actual,expected)[0]


def test_archived_refinement_rescores_with_cpu_dispatch_disabled(tmp_path):
    # NumPy's dispatch groups differ between supported versions. Disable the
    # groups exposed by this build in a fresh process, before NumPy imports.
    try:from numpy._core import _multiarray_umath as cpu
    except ImportError:from numpy.core import _multiarray_umath as cpu
    mask=','.join(cpu.__cpu_dispatch__)
    proc=subprocess.run([sys.executable,'-m','experiments.closure_ledger.localized_mouth_refinement_probe',
        '--original',str(RUN/'probe.json.gz'),'--rescore',str(RUN/'refinement.json.gz'),
        '--output-dir',str(tmp_path)],capture_output=True,text=True,
        env={**os.environ,'OPENBLAS_NUM_THREADS':'1','NPY_DISABLE_CPU_FEATURES':mask})
    assert proc.returncode==1,proc.stderr
    result=json.loads((tmp_path/'refinement_verdict.json').read_text())
    assert result['passed']==7,result
    assert result['gates']['evidence']
    assert not result['gates']['physical']
    assert not any(result['verdicts'].values())


def test_original_antipodal_bulk_scalar_profile_is_retained(original):
    for L in m.LENGTHS:
        theta=m.restore(original['profiles'][str(L)]['theta'])
        s=np.linspace(0,L-1,30)
        assert np.max(abs(theta(s)-np.arcsin(np.tanh(s))))<1e-10


def test_field_transport_requires_the_declared_line_bundle(original,refined):
    d=m.Data(refined['solutions'][-1],original['profiles']['5.5'])
    x=(5.5,.91,.37);z=(-5.5,np.pi-.91,.37+np.pi)
    a,pa=d.fields(x);b,pb=d.fields(z)
    assert np.max(abs(a+b))<1e-12
    assert np.max(abs(pa+pb))<1e-12
    assert np.max(abs(a-b))>1


@pytest.mark.parametrize('module',[m,p,r])
def test_python310_grammar(module):
    ast.parse(Path(module.__file__).read_text(),feature_version=(3,10))


def test_cli_refinement_returns_failure_despite_small_absolute_residuals(tmp_path,original):
    proc=subprocess.run([sys.executable,'-m','experiments.closure_ledger.localized_mouth_refinement_probe',
        '--original',str(RUN/'probe.json.gz'),'--rescore',str(RUN/'refinement.json.gz'),
        '--output-dir',str(tmp_path)],capture_output=True,text=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1'})
    assert proc.returncode==1
    result=json.loads((tmp_path/'refinement_verdict.json').read_text())
    assert result['passed']==7
    assert not any(result['verdicts'].values())


def test_cli_bad_original_replaces_stale_success(tmp_path):
    (tmp_path/'refinement.md').write_text('ALL GATES PASS')
    (tmp_path/'refinement_verdict.json').write_text('{"passed":8}')
    bad=tmp_path/'bad.json';bad.write_text('{}')
    proc=subprocess.run([sys.executable,'-m','experiments.closure_ledger.localized_mouth_refinement_probe',
        '--original',str(bad),'--output-dir',str(tmp_path)],capture_output=True,text=True)
    assert proc.returncode!=0
    assert 'ALL GATES PASS' not in (tmp_path/'refinement.md').read_text()
    assert not any(json.loads((tmp_path/'refinement_verdict.json').read_text())['verdicts'].values())
