#!/usr/bin/env python
"""Post-review replication, NOT new pre-registered evidence.

Reproduce the off-focus/emitter suggestions in PR #339 comment 6206173894.
The original first-transit protocol and archive are not rewritten.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.antipodal_kick import AntipodalKick, Parameters
from experiments.closure_ledger.antipodal_kick_probe import save_trace


def emitter_state(model, moving, emitting):
    s = model.initial()
    s.a[:] = 0
    s.p[:] = 0
    s.q = .2 if emitting else 0.
    s.angular[2] = .03 if moving else 0.
    return s


def event_index(a, time):
    i = np.argmin(abs(a['t']-time))
    if abs(a['t'][i]-time) > 1e-12:
        raise ValueError('Diagnostic endpoint was not integrated exactly')
    return i


def emitter_summary(model, a):
    def impulse(lo, hi):
        return (a['angular'][event_index(a, hi)]-a['angular'][event_index(a, lo)]).tolist()
    budget = {}
    for name, time in [('initial', 0.), ('early', 2.), ('final', a['t'][-1])]:
        i = event_index(a, time)
        budget[name] = dict(time=float(time),
            wave=model.wave_energy(a['a'][i], a['p'][i]),
            shape=float(a['shape_energy'][i]),
            kinetic=float(a['angular'][i]@a['angular'][i]/(2*model.par.mass)),
            interaction=float(model.par.coupling*(a['q'][i]-a['ubar'][i])**2/2),
            total=float(a['energy'][i]))
    return dict(early_impulse=impulse(0., 2.),
                antipode_window_impulse=impulse(np.pi-.6, np.pi+.6),
                return_window_impulse=impulse(2*np.pi-.6, 2*np.pi+.6), energy_budget=budget)


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    summary = dict(status='POST_REVIEW_REPLICATION',
        review='https://github.com/davidmdrpi/geometrodynamics/pull/339#issuecomment-6206173894',
        off_focus={}, emitter={}, refinement={})
    model = AntipodalKick()
    for packet, eta in [('symmetric', 0.), ('directional', .7)]:
        for label, position in [('pole', [0., 0., 1.]), ('equator', [0., 1., 0.])]:
            name = f'{packet}_{label}'
            s = model.initial(eta)
            s.x = np.array(position)
            a, stats = model.evolve(s)
            save_trace(output/f'{name}.csv', a)
            summary['off_focus'][name] = stats
            print(name, stats['final_impulse'], flush=True)
    events = [2., 2*np.pi-.6, 2*np.pi, 2*np.pi+.6]
    for name, moving, emitting in [('rest_emitter', False, True),
                                   ('moving_emitter', True, True),
                                   ('moving_unexcited', True, False)]:
        s = emitter_state(model, moving, emitting)
        a, stats = model.evolve(s, end=2*np.pi+1.2, event_times=events)
        stats.update(emitter_summary(model, a))
        summary['emitter'][name] = stats
        save_trace(output/f'{name}.csv', a)
        print(name, stats['final_impulse'], stats['return_window_impulse'], flush=True)
    # Check the new late-time result rather than rerun the already reproduced
    # first-transit experiment. Same widths, initial state and event endpoints.
    for label, cutoff, dt in [('L40', 40, .002), ('dt_half', 32, .001)]:
        model = AntipodalKick(replace(Parameters(), cutoff=cutoff))
        a, stats = model.evolve(emitter_state(model, True, True), dt=dt,
                                end=2*np.pi+1.2, event_times=events)
        stats.update(emitter_summary(model, a))
        summary['refinement'][label] = stats
        save_trace(output/f'emitter_{label}.csv', a)
        print(label, stats['final_impulse'], flush=True)
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    files = sorted(output.glob('*.csv'))+[output/'summary.json']
    (output/'sha256.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                 for p in files}, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
        default=Path('experiments/closure_ledger/runs/20261016_antipodal_kick_review'))
    run(parser.parse_args().output)
