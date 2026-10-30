"""Build a wheel and deterministic, self-verifying evidence bundle; never publish."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile

from verify_bam_foundations_release import verify
from bam_release_catalog import INDEX_PATH, NOTES_PATH, ROOT


def evidence_paths(index):
    paths = {e['path'] for e in index['imported_files']
             if e['path'].startswith(('experiments/closure_ledger/runs/', 'docs/images/'))}
    paths.update(e['path'] for e in index['capstone']['artifacts'])
    paths.update([INDEX_PATH, NOTES_PATH, 'docs/releases/bam-foundations-research-record.md',
                  'docs/images/releases/bam_capstone.png', 'docs/images/releases/bam_foundations.png',
                  'scripts/verify_bam_foundations_release.py', 'scripts/bam_release_catalog.py',
                  'CITATION.cff', 'CHANGELOG.md'])
    return sorted(paths)


def write_bundle(target, root, paths):
    # Fix gzip and tar timestamps/owners so identical evidence gives identical bytes.
    with target.open('wb') as output, gzip.GzipFile(filename='', fileobj=output, mode='wb', mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode='w', format=tarfile.PAX_FORMAT) as archive:
            for name in sorted(paths):
                data = (root/name).read_bytes()
                info = tarfile.TarInfo(name)
                info.size, info.mode, info.mtime = len(data), 0o644, 0
                archive.addfile(info, io.BytesIO(data))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'dist/bam-1.3.5')
    args = parser.parse_args()
    index = json.loads((ROOT/INDEX_PATH).read_text())
    verify()
    target = args.output_dir.resolve()
    target.mkdir(parents=True, exist_ok=True)
    version = index['release']
    subprocess.run([sys.executable, '-m', 'pip', 'wheel', str(ROOT), '--no-deps',
                    '--no-build-isolation', '--wheel-dir', str(target)], check=True)
    wheel = target/f'geometrodynamics-{version}-py3-none-any.whl'
    if not wheel.is_file():
        raise ValueError('Wheel version does not match evidence index')
    bundle = target/f'geometrodynamics-{version}-evidence.tar.gz'
    write_bundle(bundle, ROOT, evidence_paths(index))
    checksums = target/'SHA256SUMS'
    checksums.write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n'
                                for p in (wheel,bundle)))
    print(json.dumps({'wheel':str(wheel),'evidence':str(bundle),'checksums':str(checksums)},indent=2))


if __name__ == '__main__':
    main()
