"""Data-bound release headlines and snapshot provenance (standard library)."""
import hashlib
import json
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INDEX_PATH = 'docs/releases/bam-foundations-evidence.json'
NOTES_PATH = 'docs/releases/bam-foundations-1.3.5.md'
START, END = '<!-- BEGIN VERIFIED RESULTS -->', '<!-- END VERIFIED RESULTS -->'


def pointer(value, path):
    for part in path.strip('/').split('/'):
        key = part.replace('~1', '/').replace('~0', '~')
        value = value[int(key)] if isinstance(value, list) else value[key]
    return value


def metric_value(data, metric):
    value = pointer(data, metric['pointer'])
    if metric.get('operation') == 'max_samples':
        samples = [float(x) for item in value for x in pointer(item, metric['item_pointer'])]
        if not samples or not all(math.isfinite(x) for x in samples):
            raise ValueError('Invalid sampled bound')
        value = max(samples)
    if metric.get('operation') == 'difference':
        if len(value) != 2:
            raise ValueError('Separation requires exactly two centres')
        value = value[1]-value[0]
    if metric.get('operation') == 'ratio':
        value = float(value) / float(pointer(data, metric['denominator_pointer']))
    if metric['format'] == 's':
        return str(value).lower()
    value = float(value) * metric.get('scale', 1)
    if not math.isfinite(value):
        raise ValueError('Nonfinite release metric')
    if metric['format'] == 'd':
        if not value.is_integer():
            raise ValueError('Nonintegral count')
        return str(int(value))
    return format(value, metric['format'])


def table(index, root=ROOT):
    lines = ['| Result | Verified measurement | Scope |', '| --- | --- | --- |']
    for row in index['capstone']['headlines']:
        data = json.loads((root/row['file']).read_text())
        values = '; '.join(m['label']+': **'+metric_value(data,m)+'**' for m in row['metrics'])
        url = 'https://github.com/davidmdrpi/geometrodynamics/blob/'+row['source_commit']+'/'+row['report']
        lines.append(f"| [{row['title']}]({url}) | {values} | {row['scope']} |")
    return '\n'.join(lines)


def verify_catalog(index, root=ROOT, provenance=False):
    if 'capstone' not in index:
        return {}
    root = Path(root)
    catalog = index['capstone']
    for entry in catalog['artifacts']:
        data = (root/entry['path']).read_bytes()
        if len(data) != entry['bytes'] or hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError('Capstone archive mismatch: '+entry['path'])
    source_checks = 0
    if provenance:
        for entry in catalog['artifacts']+catalog['source_imports']:
            data = subprocess.check_output(['git', '-C', str(root), 'show',
                                            entry['source_commit']+':'+entry['source_path']])
            blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
            if (blob != entry['git_blob'] or len(data) != entry['bytes']
                    or hashlib.sha256(data).hexdigest() != entry['sha256']):
                raise ValueError('Capstone provenance mismatch: '+entry['path'])
            source_checks += 1
    rendered = table(index, root)
    notes = (root/NOTES_PATH).read_text()
    expected = START+'\n'+rendered+'\n'+END
    if notes.count(START) != 1 or notes.count(END) != 1 or expected not in notes:
        raise ValueError('Release headline drift; regenerate the notes table')
    return dict(capstone_archives=len(catalog['artifacts']),
                headlines_checked=len(catalog['headlines']), capstone_source_files_verified=source_checks)


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write',action='store_true',help='Regenerate only the delimited notes table')
    args=p.parse_args()
    index=json.loads((ROOT/INDEX_PATH).read_text())
    if args.write:
        path=ROOT/NOTES_PATH
        notes=path.read_text()
        before,rest=notes.split(START)
        _,after=rest.split(END)
        path.write_text(before+START+'\n'+table(index)+'\n'+END+after)
    print(json.dumps(verify_catalog(index),indent=2))
