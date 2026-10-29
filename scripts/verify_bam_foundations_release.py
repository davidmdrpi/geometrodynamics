"""Verify the v0.44.0 imported evidence and its unchanged component manifests.

Standard library only. This verifies provenance, not the truth of a physical
model. Run the test files listed in the evidence index for numerical checks.
"""
import hashlib
import json
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT/'docs/releases/bam-foundations-evidence.json'


def verify():
    index = json.loads(INDEX.read_text())
    checked = set()
    def check(entry):
        path = (ROOT/entry['path']).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError('Evidence path escapes the repository')
        content = path.read_bytes()
        if len(content) != entry['bytes'] or hashlib.sha256(content).hexdigest() != entry['sha256']:
            raise ValueError('Evidence mismatch: '+entry['path'])
        checked.add(entry['path'])
    for entry in index['imported_files']:
        check(entry)
    component_entries = 0
    for name in index['component_manifests']:
        entries = json.loads((ROOT/name).read_text())['files']
        component_entries += len(entries)
        for entry in entries:
            check(entry)
    for record in ROOT.glob('experiments/closure_ledger/runs/20261023_deformable_boundary/*.json'):
        value = json.loads(record.read_text())
        if 'state_sha256' in value:
            state = record.with_suffix('.npz')
            if hashlib.sha256(state.read_bytes()).hexdigest() != value['state_sha256']:
                raise ValueError('Full-state replay mismatch: '+str(state))
    return dict(release=index['release'], imported_files=len(index['imported_files']),
                component_manifests=len(index['component_manifests']),
                component_entries=component_entries, unique_files_checked=len(checked),
                source_commit=index['source_commit'], archive_status='all verified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tests', action='store_true', help='Also run the release integration tests')
    args = parser.parse_args()
    print(json.dumps(verify(), indent=2), flush=True)
    if args.tests:
        index = json.loads(INDEX.read_text())
        subprocess.run([sys.executable, '-m', 'pytest', '-q', *index['validation']['test_files']],
                       cwd=ROOT, check=True)
