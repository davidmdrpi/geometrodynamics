"""Authenticated replay using the verified, recorded initial state and phase.

Added after measurement; the producer and its source hashes remain unchanged.
Recomputing the trial-loop root can move the initial state by roundoff, which
can move its fitted phase by several nanoradians across SciPy builds. Verify
the setup independently, then replay from the actual archived starting state.
No budget, trajectory or terminal-decision threshold changes.
"""
import copy
import json
import numpy as np
from experiments.closure_ledger import r3_nonlinear_budget_probe as p

MANIFEST_SHA = '4ea7db2a0f0dc002bb1a4d0f7a12375d9c9f5a38d719d395741efdb4d51e44e1'


class RecordedInitialFamily:
    def __init__(self, family, cases):
        self.family = family
        self.cases = cases

    def __getattr__(self, name):
        return getattr(self.family, name)

    def initial(self, delta):
        fresh = self.family.initial(delta)
        matches = [c for c in self.cases if c['detuning'] == delta]
        if len(matches) != 1:
            raise ValueError('nonunique case')
        saved = matches[0]['initial']
        for key in ('detuning', 'scale', 'achieved_action_offset', 'initial',
                    'd0', 'target', 'loop', 'tuning_vector', 'tuning_norm',
                    'log10_autonomous_unstable_tolerance'):
            a = np.asarray(saved[key], float)
            b = np.asarray(fresh[key], float)
            if a.shape != b.shape or not np.isfinite(a).all():
                raise ValueError('invalid initial setup')
            p.close(a, b)
        # Fit the recorded state, not a roundoff-different regenerated state.
        phase, distance = self.family.nearest(p.nb.coords(np.array(saved['initial'])))
        if not np.isfinite(saved['initial_phase']):
            raise ValueError('invalid initial phase')
        if abs(np.angle(np.exp(1j*(phase-saved['initial_phase'])))) > 1e-8:
            raise ValueError('initial phase differs')
        p.close(distance, saved['d0'])
        return copy.deepcopy(saved)


def replay(directory=p.RUN):
    if p.digest(directory/'manifest.json') != MANIFEST_SHA:
        raise ValueError('manifest fingerprint mismatch')
    manifest = json.loads((directory/'manifest.json').read_text())
    names = [f'case_{k}.json.gz.b64' for k in range(4)]+['provenance.json', 'result.json']
    if set(manifest) != set(names):
        raise ValueError('archive inventory')
    for name in names:
        if p.digest(directory/name) != manifest[name]:
            raise ValueError('evidence fingerprint mismatch')
    raw = json.loads((directory/'provenance.json').read_text())
    raw['cases'] = [p.read(directory/f'case_{k}.json.gz.b64') for k in range(4)]
    fresh = p.assess(RecordedInitialFamily(p.family(), raw['cases']), raw)
    saved = json.loads((directory/'result.json').read_text())
    if fresh['label'] != saved['label'] or [(c['controlled'], c['unforced']) for c in fresh['cases']] != [(c['controlled'], c['unforced']) for c in saved['cases']]:
        raise ValueError('decisions changed')
    return fresh


if __name__ == '__main__':
    print(json.dumps(replay(), indent=2))
