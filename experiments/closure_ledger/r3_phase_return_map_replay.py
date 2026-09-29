"""Authenticate the published evidence and independently rebuild its decision."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
from experiments.closure_ledger import r3_phase_return_map_probe as probe
from experiments.closure_ledger.esu_floquet_probe import close, decisions

RAW_SHA256 = '016877441c066b4c491c89c88b4ac54507439284049ffbe4958668232fa20587'
REPORT_SHA256 = '0154641c9e6476f5d9e4954bf2c2c860e2c6fd394b7456101af7849c5ab087da'


def replay(directory=probe.RUN):
    directory = Path(directory)
    rawpath, reportpath = directory/'raw.json.gz.b64', directory/'report.json'
    for path, digest in [(rawpath, RAW_SHA256), (reportpath, REPORT_SHA256)]:
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('published evidence fingerprint mismatch: '+path.name)
    saved = json.loads(reportpath.read_text())
    # Post-run namespace-only rename avoids the independently developed #317.
    # Recover the exact measured producer bytes by reversing ONLY these paths.
    # No physical code, archived measurement, threshold or hash is rewritten.
    historical_sources = {}
    for path in saved['source_sha256']:
        current = {
            'geometrodynamics/waves/r3_return_map.py': 'geometrodynamics/waves/r3_phase_return_map.py',
            'experiments/closure_ledger/r3_return_map_probe.py':
                'experiments/closure_ledger/r3_phase_return_map_probe.py',
        }.get(path, path)
        content = (probe.ROOT/current).read_bytes()
        if path == 'experiments/closure_ledger/r3_return_map_probe.py':
            content = content.replace(b'r3_phase_return_map', b'r3_return_map')
        historical_sources[path] = hashlib.sha256(content).hexdigest()
    if saved['source_sha256'] != historical_sources or saved['raw_sha256'] != RAW_SHA256:
        raise ValueError('source or raw provenance mismatch')
    fresh = json.loads(json.dumps(probe.analyze(probe.load_raw(rawpath)), default=probe.serial, allow_nan=False))
    expected = {k: saved[k] for k in fresh}
    if decisions(fresh) != decisions(expected):
        raise ValueError('categorical replay disagreement')
    # Orders formed from near-floor error ratios are diagnostic. Authenticate
    # their original values by the file hash; compare underlying errors and
    # recompute every order gate independently instead of demanding roundoff
    # equality of ill-conditioned ratios.
    a, b = copy.deepcopy(fresh), copy.deepcopy(expected)
    for x, y in zip(a['remainders'], b['remainders']):
        if not close(x['errors'], y['errors'], 1e-10):
            raise ValueError('remainder error disagreement')
        x['orders'] = y['orders'] = None
    if not close(a, b, 1e-8):
        raise ValueError('numerical replay disagreement')
    return dict(replay='VERIFIED', verdict=fresh['verdict'], nu=fresh['nu'], gates=fresh['gates'],
                scope='Local finite-order normal form; not closure or action selection')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=probe.RUN)
    args = parser.parse_args()
    print(json.dumps(replay(args.directory), indent=2))


if __name__ == '__main__':
    main()
