"""Authenticate the diagonal Bianchi IX archive and rebuild every decision."""
import argparse
import hashlib
import json
from pathlib import Path
from . import diagonal_bianchi_probe as p
from .esu_floquet_probe import close

MANIFEST_SHA256 = 'eebc287e12240f2d293f99a7206c10465a089825dd93b115a05c4bd1b6536dda'


def decisions(report):
    return dict(crossing=report['crossing'],closure=report['closure'],accepted=report['accepted'],
                termination=report['termination'],
                circles=[m['accepted'] for m in report['circles']],
                refinement=[(m['accepted'],m['resolved_side']) for m in report.get('refinements',[])],
                shooting=[(m.get('candidate',False),m['verified']) for m in report.get('shooting',[])])


def replay(directory=p.RUN):
    directory=Path(directory)
    data=(directory/'manifest.json').read_bytes()
    if hashlib.sha256(data).hexdigest()!=MANIFEST_SHA256:
        raise ValueError('manifest fingerprint mismatch')
    manifest=json.loads(data)
    if {x.name for x in directory.iterdir() if x.is_file()} != set(manifest)|{'manifest.json'}:
        raise ValueError('unexpected or missing evidence file')
    for name,digest in manifest.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('evidence fingerprint mismatch: '+name)
    saved=json.loads((directory/'report.json').read_text())
    if saved['sources']!=p.source_hashes():raise ValueError('measured source fingerprint mismatch')
    raw=p.load_run(directory)
    fresh=json.loads(json.dumps(p.analyze(raw),default=p.serial,allow_nan=False))
    if decisions(fresh)!=decisions(saved):raise ValueError('categorical replay disagreement')
    expected={k:saved[k] for k in fresh}
    if not close(fresh,expected,1e-8):raise ValueError('numerical replay disagreement')
    return dict(replay='VERIFIED',**decisions(fresh))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=p.RUN)
    print(json.dumps(replay(parser.parse_args().directory),indent=2))
