#!/usr/bin/env python
"""Restore the exact large assets from checked pieces; no simulation rerun.

Large individual upload requests did not complete. These pieces contain the
original bytes, not regenerated or rounded numerical data.
"""
import hashlib
import json
from pathlib import Path


def load_asset(path):
    path = Path(path)
    parts = Path(str(path)+'.parts')
    if not parts.exists():
        return path.read_bytes()
    manifest = json.loads((parts/'manifest.json').read_text())
    if path.exists():
        data = path.read_bytes()
    else:
        chunks = []
        for entry in manifest['parts']:
            chunk = (parts/entry['file']).read_bytes()
            if len(chunk) != entry['size'] or hashlib.sha256(chunk).hexdigest() != entry['sha256']:
                raise ValueError(f'Corrupt asset piece: {entry["file"]}')
            chunks.append(chunk)
        data = b''.join(chunks)
    if len(data) != manifest['size'] or hashlib.sha256(data).hexdigest() != manifest['sha256']:
        raise ValueError(f'Asset checksum mismatch: {path}')
    return data


def restore():
    root = Path(__file__).resolve().parents[1]
    for name in ['docs/images/moving_mouth_pairs/moving_mouth_pairs.gif',
                 'experiments/closure_ledger/runs/20261016_moving_mouth_pairs/render_states.npz']:
        path = root/name
        data = load_asset(path)
        path.write_bytes(data)
        print(f'Restored and verified {path} ({len(data)} bytes)')


if __name__ == '__main__':
    restore()
