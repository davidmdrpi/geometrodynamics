"""Summarize histories, make a compact replay, and chunk verified binary archives."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from geometrodynamics.waves.coupled_mouth_throat import real_basis
from geometrodynamics.waves.legacy_history_throat import Parameters
from experiments.closure_ledger.history_throat_probe import CASES, ROOT


def assemble(root=ROOT):
    index = json.loads((root/'replay_parts.json').read_text())
    for name, entry in index.items():
        if Path(name).name != name or not name.endswith('.npz'):
            raise ValueError('Invalid replay name')
        target = root/name
        if target.exists():
            if hashlib.sha256(target.read_bytes()).hexdigest() != entry['sha256']:
                raise ValueError(f'Existing replay differs: {name}')
            continue
        pieces = []
        for part in entry['parts']:
            data = (root/part['path']).read_bytes()
            if hashlib.sha256(data).hexdigest() != part['sha256']:
                raise ValueError(f'Part checksum failed: {part["path"]}')
            pieces.append(data)
        data = b''.join(pieces)
        if hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise ValueError(f'Replay checksum failed: {name}')
        target.write_bytes(data)


def chunk_replays(root):
    directory = root/'replay_parts'; directory.mkdir(exist_ok=True)
    for old in directory.glob('*.npz.part[0-9][0-9][0-9]'):
        old.unlink()
    index = {}
    for path in sorted(root.glob('*.npz')):
        data = path.read_bytes(); parts = []
        for i, start in enumerate(range(0, len(data), 480000)):
            target = directory/f'{path.name}.part{i:03d}'
            block = data[start:start+480000]; target.write_bytes(block)
            parts.append(dict(path=str(target.relative_to(root)), bytes=len(block),
                              sha256=hashlib.sha256(block).hexdigest()))
        index[path.name] = dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), parts=parts)
    (root/'replay_parts.json').write_text(json.dumps(index, indent=2)+'\n')
    (root/'.gitignore').write_text('*.npz\n')
    (root/'README.md').write_text('''# History-coupled throat replay data

CSV and JSON files are directly readable. The complete NPZ archives are split
into lossless parts of at most 480,000 bytes so interrupted uploads can resume.
`replay_parts.json` records every part hash and the original archive hash.

From the repository root:

```bash
python -m experiments.closure_ledger.history_throat_report --assemble
```

This verifies and reassembles the exact NPZ files; it does not rerun a trajectory.
The solver and test data are classical. Geometric event tables are ray predictions,
not measured tagged pulses. Read the protocol and results before interpreting them.
''')


def write_viewer(results, destination):
    template = Path('experiments/closure_ledger/history_throat_viewer.html').read_text()
    theta = np.linspace(.06, np.pi-.06, 13)
    phi = np.arange(24)*2*np.pi/24
    th, ph = np.meshgrid(theta, phi, indexing='ij')
    xyz = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], axis=-1)
    payload = dict(points=np.round(xyz.reshape(-1, 3), 5).tolist(), rows=len(theta), cols=len(phi), cases={})
    for name in ['baseline', 'asymmetric', 'reversed']:
        a = np.load(ROOT/f'{name}.npz'); cfg = results[name]['parameters']; n = (cfg['lmax']+1)**2
        indices = np.unique(np.linspace(0, len(a['t'])-1, 61).astype(int))
        t = a['t'][indices]; ys = a['states'][indices]; q = ys[:, 2*n:2*n+2]
        s = ys[:, 2*n+4:2*n+6]
        basis = real_basis(cfg['lmax'], th, ph).reshape(-1, n)
        field = ys[:, :n]@basis.T; scale = float(np.max(abs(field)))
        quantized = np.rint(field/max(scale, 1e-20)*127).astype(int)
        channels = []
        for i in range(2):
            sigma = s[:, i] if i == 0 or cfg['time_orientation'] == 1 else cfg['neck_horizon']-s[:, i]
            channels.append(np.array([np.interp(sigma, a['sigma'], a['neck_states'][:, j])
                                     for j in range(cfg['nodes'])]).T)
        mouth_values = []
        for row, angles in zip(ys, q):
            y = real_basis(cfg['lmax'], np.arccos(np.cos(angles)),
                           np.where(np.sin(angles) >= 0, 0., np.pi))
            mouth_values.append((y@row[:n]).tolist())
        payload['cases'][name] = dict(t=np.round(t, 5).tolist(), q=np.round(q, 6).tolist(),
            clocks=np.round(s, 6).tolist(), momenta=np.round(ys[:, 2*n+2:2*n+4], 6).tolist(),
            field=quantized.tolist(), field_scale=scale,
            mouth_field=np.round(mouth_values, 6).tolist(),
            neck_a=np.round(channels[0], 7).tolist(), neck_b=np.round(channels[1], 7).tolist(),
            orientation=cfg['time_orientation'], anchor=cfg['neck_horizon'], length=cfg['length'],
            energy_error=results[name]['energy_relative_residual'],
            clock_work_error=results[name]['neck_clock_work_residual'])
    destination.write_text(template.replace('/*REPLAY_DATA*/{}', json.dumps(payload, separators=(',', ':'))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assemble', action='store_true')
    args = parser.parse_args()
    if args.assemble:
        assemble(); return
    results = {name: json.loads((ROOT/f'{name}.json').read_text()) for name in CASES}
    data = {name: np.genfromtxt(ROOT/f'{name}.csv', delimiter=',', names=True) for name in CASES}
    comparisons = {}
    for name, reference in [('baseline', 'no_throat'), ('baseline', 'blocked'),
            ('reversed', 'baseline'), ('higher_resolution', 'high_no_throat'),
            ('higher_resolution', 'high_blocked'), ('high_reversed', 'higher_resolution'),
            ('tight_tolerance', 'baseline'), ('higher_resolution', 'baseline'),
            ('high_asymmetric', 'asymmetric')]:
        a, b = data[name], data[reference]
        differences = {key: a[key]-np.interp(a['t'], b['t'], b[key])
                       for key in ['qa', 'qb', 'pa', 'pb', 'separation']}
        comparisons[name+'__minus__'+reference] = dict(
            final={k: float(v[-1]) for k, v in differences.items()},
            max_absolute={k: float(np.max(abs(v))) for k, v in differences.items()})
    summary = dict(cases=results, comparisons=comparisons,
        all_converged=all(r['converged'] for r in results.values()),
        all_ledgers_pass=all(r['ledger_pass'] for r in results.values()),
        interpretation='Numerical consistency of a constitutive action; not GR, calibrated QED-like exchange, or topology surgery.')
    (ROOT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    images = Path('docs/images/history_throat'); images.mkdir(parents=True, exist_ok=True)
    plt.style.use('dark_background')
    plt.rcParams['svg.hashsalt'] = 'history-throat-20261020'
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = dict(baseline='#5ad8e6', blocked='#b4b8ce', no_throat='#e0b467',
                  reversed='#c893ff', spatial_even='#a3e9aa', absorbing_load='#f29da5')
    b = data['baseline']
    for name, color in colors.items():
        d = data[name]
        axes[0, 0].plot(d['t'], d['separation'], label=name.replace('_', ' '), color=color)
        axes[0, 1].plot(d['t'], d['qa']-np.interp(d['t'], b['t'], b['qa']), color=color,
                        label=name.replace('_', ' '))
    axes[0, 0].axhline(np.pi, color='white', alpha=.3, ls=':')
    for name in ['baseline', 'asymmetric', 'reversed']:
        d = data[name]
        axes[1, 0].plot(d['t'], d['sb']-d['sa'], label=name)
        axes[1, 1].plot(d['t'], d['ka']+d['kb'], label=name)
    titles = ['Mouth separation (rad)', 'Mouth A position relative to baseline (rad)',
              'Simultaneous proper-clock difference sB − sA', 'Clock energies kA + kB (exterior ledger)']
    for ax, title in zip(axes.flat, titles):
        ax.set(title=title, xlabel='Exterior time t'); ax.grid(alpha=.15); ax.legend(fontsize=8)
    fig.suptitle('Reciprocal histories: a resolved throat changes continuous mouth motion\n'
                 'Constitutive fixed-neck model; temporal reversal uses an explicit future anchor', fontsize=14)
    fig.savefig(images/'comparison.svg', metadata={'Date': None})
    svg = images/'comparison.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    fig.savefig(images/'comparison.png', dpi=120); plt.close(fig)
    write_viewer(results, images/'legacy_scalar_replay.html')
    chunk_replays(ROOT)
    paths = [p for p in ROOT.rglob('*') if p.is_file() and p.suffix != '.npz' and p.name != 'sha256.json']
    paths += list(images.glob('*')) + [Path('geometrodynamics/waves/legacy_history_throat.py'),
        Path('tests/test_history_throat.py'), Path('docs/history_coupled_throat_protocol.md'),
        Path('docs/history_coupled_throat_results.md')]
    paths += list(Path('experiments/closure_ledger').glob('history_throat*'))
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    (ROOT/'sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')
    print(json.dumps({'all_converged': summary['all_converged'],
        'all_ledgers_pass': summary['all_ledgers_pass'], 'comparisons': comparisons}, indent=2))


if __name__ == '__main__':
    main()
