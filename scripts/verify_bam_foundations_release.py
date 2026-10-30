"""Verify saved v1.3.5 evidence without freezing evolving source or tests.

Default: check working-tree run archives and images against the release index.
--ref: check all indexed imports at a historical Git ref (including source).
--source-ref: additionally verify the indexed blobs against the source commit.
These checks establish consistency/provenance, not physical truth. Git history
is required only for the explicit ref checks; nothing is fetched automatically.
"""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
INDEX_PATH = 'docs/releases/bam-foundations-evidence.json'


def archive_path(name):
    return name.startswith(('experiments/closure_ledger/runs/', 'docs/images/'))


def verify(root=ROOT, ref=None, source_ref=None):
    root = Path(root).resolve()
    index = json.loads((root/INDEX_PATH).read_text())

    def git(*args):
        return subprocess.check_output(['git', '-C', str(root), *args], stderr=subprocess.PIPE)

    def commit(value):
        try:
            return git('rev-parse', '--verify', '--end-of-options', value+'^{commit}').decode().strip()
        except subprocess.CalledProcessError as exc:
            raise ValueError('Git ref unavailable; fetch its history first: '+value) from exc

    snapshot = commit(ref) if ref else None
    source = commit(source_ref) if source_ref else None
    if source and source != index['source_commit']:
        raise ValueError('Source ref does not resolve to the indexed source_commit')

    def read(name, revision=snapshot):
        path = PurePosixPath(name)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('Evidence path escapes the repository')
        if revision:
            return git('show', revision+':'+name)
        local = (root/name).resolve()
        if not local.is_relative_to(root):
            raise ValueError('Evidence path escapes the repository')
        return local.read_bytes()

    checked = set()
    snapshot_catalog = {}
    if snapshot and 'capstone' in index:
        snapshot_catalog = json.loads(read(INDEX_PATH)).get('capstone', {})
    replacements = {entry['path']: entry for entry in snapshot_catalog.get('source_imports', [])}

    def check(entry):
        entry = replacements.get(entry['path'], entry)
        name = entry['path']
        if not snapshot and not archive_path(name):
            return
        content = read(name)
        if len(content) != entry['bytes'] or hashlib.sha256(content).hexdigest() != entry['sha256']:
            raise ValueError('Evidence mismatch: '+name)
        checked.add(name)

    provenance_checked = 0
    for entry in index['imported_files']:
        check(entry)
        if source:
            content = read(entry['path'], source)
            blob = hashlib.sha1(b'blob '+str(len(content)).encode()+b'\0'+content).hexdigest()
            if (blob != entry['git_blob'] or len(content) != entry['bytes']
                    or hashlib.sha256(content).hexdigest() != entry['sha256']):
                raise ValueError('Source provenance mismatch: '+entry['path'])
            if 'origin_commit' in entry:
                origin = commit(entry['origin_commit'])
                if read(entry['path'], origin) != content:
                    raise ValueError('Origin provenance mismatch: '+entry['path'])
            provenance_checked += 1

    for entry in snapshot_catalog.get('source_imports', []) + snapshot_catalog.get('artifacts', []):
        check(entry)
    component_entries = 0
    for name in index['component_manifests']:
        entries = json.loads(read(name))['files']
        component_entries += len(entries)
        for entry in entries:
            check(entry)
    # Use indexed records, not a working-tree glob, so historical ref checks
    # also work when replay files are absent or changed in the current checkout.
    for entry in index['imported_files']:
        name = entry['path']
        if name.startswith('experiments/closure_ledger/runs/20261023_deformable_boundary/') and name.endswith('.json'):
            value = json.loads(read(name))
            if 'state_sha256' in value:
                state = str(PurePosixPath(name).with_suffix('.npz'))
                if hashlib.sha256(read(state)).hexdigest() != value['state_sha256']:
                    raise ValueError('Full-state replay mismatch: '+state)
    catalog_checks = {}
    if not snapshot and 'capstone' in index:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from bam_release_catalog import verify_catalog
        catalog_checks = verify_catalog(index, root, provenance=bool(source))
    return dict(**catalog_checks, release=index['release'], mode='historical imports' if snapshot else 'working-tree archives',
                checked_ref=snapshot, imported_files=len(index['imported_files']),
                component_manifests=len(index['component_manifests']),
                component_entries=component_entries, unique_files_checked=len(checked),
                source_commit=index['source_commit'], source_files_verified=provenance_checked,
                archive_status='all verified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ref', help='Check all indexed imports at this release Git ref')
    parser.add_argument('--source-ref', help='Verify source blobs at this ref; must resolve to source_commit')
    parser.add_argument('--tests', action='store_true', help='Also run integration tests in the CURRENT checkout')
    args = parser.parse_args()
    if args.ref and args.tests:
        parser.error('--ref and --tests cannot be combined; check out the ref to test its code')
    print(json.dumps(verify(ref=args.ref, source_ref=args.source_ref), indent=2), flush=True)
    if args.tests:
        index = json.loads((ROOT/INDEX_PATH).read_text())
        subprocess.run([sys.executable, '-m', 'pytest', '-q', *index['validation']['test_files']],
                       cwd=ROOT, check=True)
