"""Archive gates must survive source evolution and still detect lost evidence."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

import pytest

SPEC = importlib.util.spec_from_file_location(
    'bam_release_verifier', Path(__file__).resolve().parents[1]/'scripts/verify_bam_foundations_release.py')
verifier = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verifier)


@pytest.fixture
def evidence(tmp_path):
    def write(name, content):
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        return dict(path=name, bytes=len(content), sha256=hashlib.sha256(content).hexdigest(),
                    git_blob=hashlib.sha1(b'blob '+str(len(content)).encode()+b'\0'+content).hexdigest())

    run = 'experiments/closure_ledger/runs/20261023_deformable_boundary/'
    state = write(run+'case.npz', b'replay bytes')
    record = write(run+'case.json', json.dumps({'state_sha256': state['sha256']}).encode())
    source = write('geometrodynamics/waves/pair_birth.py', b'original = True\n')
    test = write('tests/test_pair_birth.py', b'assert True\n')
    report = write('docs/report.md', b'Original interpretation\n')
    image = write('docs/images/probe.png', b'figure bytes')
    files = [state, record, source, test, report, image]
    manifest = write(run+'manifest.json', json.dumps({'files': files}).encode())

    def git(*args):
        return subprocess.check_output(['git', '-C', str(tmp_path), *args], stderr=subprocess.PIPE).decode().strip()

    git('init', '-q')
    git('add', '.')
    git('-c', 'user.name=Verifier Test', '-c', 'user.email=verifier@example.invalid',
        'commit', '-qm', 'Original evidence')
    ref = git('rev-parse', 'HEAD')
    source['origin_commit'] = ref
    index = dict(release='test', source_commit=ref, imported_files=files+[manifest],
                 component_manifests=[manifest['path']])
    index_path = tmp_path/verifier.INDEX_PATH
    index_path.parent.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(index))
    return tmp_path, ref, index, index_path


def test_current_source_tests_and_reports_can_evolve(evidence):
    root, ref, _, _ = evidence
    for name in ('geometrodynamics/waves/pair_birth.py', 'tests/test_pair_birth.py', 'docs/report.md'):
        (root/name).write_text('later upstream correction\n')
    result = verifier.verify(root)
    assert result['unique_files_checked'] == 4
    assert result['source_files_verified'] == 0
    assert verifier.verify(root, ref=ref)['unique_files_checked'] == 7
    assert verifier.verify(root, source_ref=ref)['source_files_verified'] == 7


@pytest.mark.parametrize('name', ['case.json', 'case.npz', 'manifest.json'])
def test_changed_archive_is_rejected(evidence, name):
    root, _, _, _ = evidence
    (root/'experiments/closure_ledger/runs/20261023_deformable_boundary'/name).write_bytes(b'changed')
    with pytest.raises(ValueError, match='Evidence mismatch'):
        verifier.verify(root)


def test_changed_figure_is_rejected(evidence):
    root, _, _, _ = evidence
    (root/'docs/images/probe.png').write_bytes(b'changed')
    with pytest.raises(ValueError, match='Evidence mismatch'):
        verifier.verify(root)


def test_historical_check_uses_ref_even_with_missing_working_files(evidence):
    root, ref, index, _ = evidence
    for entry in index['imported_files']:
        (root/entry['path']).unlink()
    assert verifier.verify(root, ref=ref, source_ref=ref)['unique_files_checked'] == 7
    with pytest.raises(FileNotFoundError):
        verifier.verify(root)


def test_provenance_checks_git_blob_not_just_local_hashes(evidence):
    root, ref, index, path = evidence
    index['imported_files'][2]['git_blob'] = '0'*40
    path.write_text(json.dumps(index))
    verifier.verify(root)  # Internal consistency is explicitly the default scope.
    with pytest.raises(ValueError, match='Source provenance mismatch'):
        verifier.verify(root, source_ref=ref)


def test_missing_ref_fails_explicitly(evidence):
    root, _, _, _ = evidence
    with pytest.raises(ValueError, match='Git ref unavailable'):
        verifier.verify(root, ref='unavailable-release')


def test_wrong_source_commit_rejected(evidence):
    root, ref, index, path = evidence
    index['source_commit'] = '0'*40
    path.write_text(json.dumps(index))
    with pytest.raises(ValueError, match='Source ref does not resolve'):
        verifier.verify(root, source_ref=ref)
