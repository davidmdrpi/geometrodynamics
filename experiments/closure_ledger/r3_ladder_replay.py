"""Authenticated replay of the ladder archives (docs/r3_ladder_prereg.md).

Fail-closed: archive bytes must match the pinned SHA-256, the bound source hashes must match the files
on disk, every usable high-precision row must be well formed, and the registered labels are re-derived
with the frozen scoring code in a temporary directory. With --full, the high-precision solution at
every 10th phase of every rung is re-closed with the unchanged 160-bit map.
Usage: python -m experiments.closure_ledger.r3_ladder_replay [--full]
"""
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path
import gmpy2
from gmpy2 import mpfr
import numpy as np
from geometrodynamics.waves import lrs_taylor as lt
from experiments.closure_ledger import r3_ladder_probe as probe

SHA256 = {
    'started.json': 'd770e64407a9e9eda858cc9bd486a57fe9c674740097bcdb3d265b452f5ab65e',
    'resume_r1.json': '32e798160a10f42ad65cd18e7413f75f4b585f862c315a117de2946a94d0e1b9',
    'resume_r2.json': '780bb0fd2f7157212bdbbb1edd9a530d2d261a63783f45be84ec1545dea422b7',
    'result.json': '4cb3e6d028e621ad08a978e5807a490e706625aef5d5eb499d901a68d84641a7',
    'scan_5_11.json': '243dc174d6dd11c691aa7eba194b4c40896511aa810a7e555fdd22abb52d61d3',
    'hp_5_11.json': 'f878a8a2a0d13379bff99876c332a74e0019de81979e3656686f9874c83beb0c',
    'hpnoise_5_11.json': '4ab867d94e7d74d104115084df041f362f53732870194d8cf5ede5ab800892f0',
    'scan_4_9.json': '995062c6d0ed9e0b6dd638f2db6542fe5d1f1b80fac2b08dc4cad606858460df',
    'hp_4_9.json': '4e5009ef8ff88f3cdde7e12a0427b13d60782720bd609481aeb2928f1d89898e',
    'hpnoise_4_9.json': '94026205fbb95b8432e9ea456bfa6d76a990e358a45a3e4bc850c2c32aaa4f2f',
    'scan_3_7.json': 'a91756fa9a55768013e4b2149e88431e5fb24cf9999090bb1872a84688151b9a',
    'hp_3_7.json': 'ea1f54af3181d67c8f6f98879f0faf758791b74e786e3441ef23cdf764b854c1',
    'hpnoise_3_7.json': 'e4864afe34214a68b70b9e2483ed990454441ffeb58af01b788187d1cbda9fe5',
    'scan_5_12.json': '95e216fd9e32d4009492f9f08bf917515340c240531b6d1ae2607230f5b103df',
    'hp_5_12.json': 'aad5e7e7a0514845b525c6c4c14b777e118c84c057c48b71412c44ed8e0b7152',
    'hpnoise_5_12.json': '723c1330e0f56a57206b07b3861a86f827e23b92bb614d0235834608549d5178',
    'scan_2_5.json': 'eec80874c4fefb1ce3d79e0dbb98da83edf598024b7ab0c069c4c8aef1b14640',
    'hp_2_5.json': '2173736c004b1081f7a85f082476385f68935bc9423e55ff0f153079fca18f49',
    'hpnoise_2_5.json': '29d6f01fe98c9d47ae4e987f9c9522bd54beab5c0839be2237b8f535497bb8b4',
    'scan_3_8.json': '89731e7a8d716296cf0c9a190e6d684f3179188b040811b013a8b8b4c6dd8269',
    'hp_3_8.json': 'fb37d92d5ba101de18890bd075fa11599a17baa1970e59fb7b03e97830058ff0',
    'hpnoise_3_8.json': 'db48bb18f87ee83080e9ea4546212824548e64c411874eca54d07491ead8a332',
    'scan_4_11.json': 'b7e377bd35f00a994af74a602e707bd01db786cd95c694048449bccdf3e41948',
    'hp_4_11.json': 'c8f05ee174102ebadb00640934d286a599f9bec27ea0207b960448d49ba97967',
    'hpnoise_4_11.json': 'd25ead2960282b342b3bbe177907f5c75fee28206b76abc4e027f8e3ea5c7d12',
}
FULL_TOL = 1e-33
LABELS = ('primary', 'signal_2_5', 'harmonic_selection', 'exact_integrability')


def names():
    out = ['started.json', 'resume_r1.json', 'resume_r2.json', 'result.json']
    for p, q in probe.RUNGS:
        out += [f'{s}_{probe.tag(p, q)}.json' for s in ('scan', 'hp', 'hpnoise')]
    return out


def _load(directory):
    if set(SHA256) != set(names()):
        raise ValueError('pinned archive list incomplete')
    out = {}
    for name, h in SHA256.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != h:
            raise ValueError(f'fingerprint mismatch: {name}')
        out[name] = json.loads(raw)
    return out


def validate_rows(scan, rows, tol):
    if len(rows) != probe.N_PHASE or [r['j'] for r in rows] != list(range(probe.N_PHASE)):
        raise ValueError('rows incomplete or reordered')
    q = scan['q']
    for r, pt in zip(rows, scan['points']):
        if not r['ok']:
            continue
        if not pt['ok']:
            raise ValueError('high-precision row without a converged scan point')
        if len(r['nodes']) != q or any(len(z) != 4 for z in r['nodes']):
            raise ValueError('malformed node set')
        if not (r['residual'] < tol and r['history'][-1] == r['residual']):
            raise ValueError('ok flag inconsistent with residual')
        if abs(float(mpfr(r['lam']))-r['lam_float']) > 1e-15*max(1., abs(r['lam_float'])):
            raise ValueError('lambda string and float disagree')


def reclose(scan, rows, every=10):
    worst = 0.
    P = probe.hp_P(lt.C1)
    with gmpy2.context(gmpy2.get_context(), precision=lt.C1['bits']):
        for j in range(0, probe.N_PHASE, every):
            r, pt = rows[j], scan['points'][j]
            if not r['ok']:
                continue
            Z = [[mpfr(v) for v in z] for z in r['nodes']]
            g, t0, c0 = ([mpfr(v) for v in pt[k]] for k in ('g', 't0', 'c0'))
            F = probe.hp_residual(P, Z, mpfr(r['lam']), g, t0, c0)
            worst = max(worst, float(max(abs(f) for f in F)))
            if worst > FULL_TOL:
                raise ValueError(f'row {j} does not re-close: {worst:.2e}')
    return worst


def replay(directory=None, full=False):
    directory = Path(directory or probe.RUN_DIR)
    recs = _load(directory)
    src = probe.sources()
    for name, rec in recs.items():
        if rec['sources'] != src:
            raise ValueError(f'source hashes differ: {name}')
    for p, q in probe.RUNGS:
        t = probe.tag(p, q)
        validate_rows(recs[f'scan_{t}.json'], recs[f'hp_{t}.json']['rows'], probe.HP_TOL['C1'])
        validate_rows(recs[f'scan_{t}.json'], recs[f'hpnoise_{t}.json']['rows'], probe.HP_TOL['C2'])
    old = probe.RUN_DIR
    with tempfile.TemporaryDirectory() as tmp:
        for name in recs:
            if name != 'result.json':
                shutil.copy(directory/name, Path(tmp)/name)
        probe.RUN_DIR = Path(tmp)
        try:
            again = json.loads(json.dumps(probe.score()))
        finally:
            probe.RUN_DIR = old
    ref = recs['result.json']['result']
    for k in LABELS:
        if again[k] != ref[k]:
            raise ValueError(f'{k} does not re-derive')
    if again['exponent_fit']['label'] != ref['exponent_fit']['label']:
        raise ValueError('exponent_fit does not re-derive')
    out = dict(replay='VERIFIED', **{k: again[k] for k in LABELS}, exponent_fit=again['exponent_fit']['label'])
    if full:
        out['worst_reclosure'] = max(reclose(recs[f'scan_{probe.tag(p, q)}.json'], recs[f'hp_{probe.tag(p, q)}.json']['rows'])
                                     for p, q in probe.RUNGS)
    return out


if __name__ == '__main__':
    print(json.dumps(replay(full='--full' in sys.argv), indent=1))
