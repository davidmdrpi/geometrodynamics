"""Serialise the frozen r3_family_probe.score() result (added 2026-10-01).

The frozen probe's `score` stage crashed when printing: numpy bools are not
JSON-serialisable. Editing the probe would change its bound source hash, so
this wrapper calls the unchanged score() and converts numpy scalars. It also
authenticates that both stage archives were produced by the committed sources.
"""
import json
import numpy as np
from experiments.closure_ledger import r3_family_probe as probe


def _plain(o):
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError(type(o).__name__)


def result():
    F = json.loads((probe.RUN_DIR/'stage_F.json').read_text())
    S = json.loads((probe.RUN_DIR/'stage_S.json').read_text())
    src = probe.sources()
    if F['sources'] != src or S['sources'] != src:
        raise RuntimeError('stage archives do not match the committed sources')
    return json.loads(json.dumps(dict(sources=src, **probe.score(F, S)), default=_plain))


def main():
    path = probe.RUN_DIR/'result.json'
    if path.exists():
        raise FileExistsError('append-only: '+str(path))
    rec = result()
    path.write_text(json.dumps(rec, indent=1, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in rec.items() if k not in ('diagonal_pairs', 'off_diagonal_pairs', 'sources')}, indent=1))


if __name__ == '__main__':
    main()
