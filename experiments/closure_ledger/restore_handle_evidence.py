"""Restore the exact registered JSON from small lossless transport parts."""
import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path

RAW_SHA256 = '8a858b6d7ecce0de3047f92a8baa77fd9c262ec9475d6886d539eac26cadf815'
RAW_SIZE = 2204136
DEFAULT = Path(__file__).resolve().parent/'runs/20260923_handle_evolution'


def restore(directory=DEFAULT):
    directory = Path(directory)
    target = directory/'evolution.json'
    if target.exists():
        raw = target.read_bytes()
        if len(raw) != RAW_SIZE or hashlib.sha256(raw).hexdigest() != RAW_SHA256:
            raise ValueError('existing evolution evidence hash mismatch')
        return target
    manifest = json.loads((directory/'evolution.transport.json').read_text())
    if manifest['raw_sha256'] != RAW_SHA256 or manifest['raw_size'] != RAW_SIZE:
        raise ValueError('transport manifest changed registered evidence identity')
    chunks = []
    for i, part in enumerate(manifest['parts']):
        if part['name'] != f'evolution.json.gz.b64.part{i:03d}':
            raise ValueError('transport part sequence mismatch')
        chunk = (directory/part['name']).read_bytes()
        if len(chunk) != part['size'] or hashlib.sha256(chunk).hexdigest() != part['sha256']:
            raise ValueError('transport part hash mismatch')
        chunks.append(chunk)
    compressed = base64.b64decode(b''.join(chunks), validate=True)
    if hashlib.sha256(compressed).hexdigest() != manifest['gzip_sha256']:
        raise ValueError('compressed evidence hash mismatch')
    raw = gzip.decompress(compressed)
    if len(raw) != RAW_SIZE or hashlib.sha256(raw).hexdigest() != RAW_SHA256:
        raise ValueError('restored evolution evidence hash mismatch')
    temporary = target.with_suffix('.json.restoring')
    temporary.write_bytes(raw)
    temporary.replace(target)
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=Path, default=DEFAULT)
    print(restore(parser.parse_args().directory))
