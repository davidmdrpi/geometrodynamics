"""Read an archive directly or reassemble hash-checked, lossless file parts."""
import hashlib
import json
from pathlib import Path


def read_bytes(path):
    path = Path(path)
    if path.is_file():
        return path.read_bytes()
    manifest = json.loads(path.with_name(path.name + '.parts.json').read_text())
    parts = []
    for item in manifest['parts']:
        name = item['name']
        if Path(name).name != name:
            raise ValueError('archive part must be a sibling file')
        data = path.with_name(name).read_bytes()
        if len(data) != item['bytes'] or hashlib.sha256(data).hexdigest() != item['sha256']:
            raise ValueError('archive part hash/size mismatch')
        parts.append(data)
    data = b''.join(parts)
    if len(data) != manifest['bytes'] or hashlib.sha256(data).hexdigest() != manifest['sha256']:
        raise ValueError('reassembled archive hash/size mismatch')
    return data
