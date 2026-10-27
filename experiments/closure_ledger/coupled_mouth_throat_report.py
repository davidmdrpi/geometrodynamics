"""Render actual replay states and record controls, recovery checks and hashes."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from geometrodynamics.waves.coupled_mouth_throat import CoupledMouthThroat, Parameters, real_basis
from experiments.closure_ledger.coupled_mouth_throat_probe import CASES, DEFAULT_ROOT


def main():
    root = DEFAULT_ROOT
    image_dir = Path('docs/images/coupled_mouth_throat')
    image_dir.mkdir(parents=True, exist_ok=True)
    results = {name: json.loads((root/f'{name}.json').read_text()) for name in CASES}
    data = {name: np.genfromtxt(root/f'{name}.csv', names=True, delimiter=',') for name in CASES}
    original = Path('experiments/closure_ledger/runs/20261019_coupled_mouth_throat')
    recovery = {}
    for name in ['baseline', 'dephased', 'fixed_mouths']:
        before, after = np.load(original/f'{name}.npz'), np.load(root/f'{name}.npz')
        recovery[name] = {key: float(np.max(abs(before[key]-after[key])))
                          for key in ['initial', 'final', 'states', 't']}
    b = data['baseline']; r0 = results['baseline']['parameters']['radius0']
    refinements = {}
    for name in ['tight_tolerance', 'higher_resolution']:
        c = data[name]
        refinements[name] = {key: float(np.max(abs(b[key]-np.interp(b['tau'], c['tau'], c[key]))))
                             for key in ['separation', 'radius']}
    summary = dict(cases=results, recovered_archive_max_absolute_differences=recovery,
                   refinement_max_absolute_differences=refinements)
    (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    plt.style.use('dark_background')
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for name, d in data.items():
        args = dict(label=name.replace('_', ' '), lw=2 if name == 'baseline' else 1,
                    alpha=1 if name == 'baseline' else .7)
        axes[0, 0].plot(d['tau']/np.pi, d['separation'], **args)
        axes[0, 1].plot(d['tau']/np.pi, d['radius']/r0, **args)
        axes[1, 0].plot(d['tau']/np.pi, d['speed'], **args)
        keys = ['sphere', 'throat', 'port', 'mouth_kinetic', 'radial_kinetic', 'neck_potential']
        e = sum(d[k] for k in keys)
        axes[1, 1].plot(d['tau']/np.pi, abs(e-e[0])/max(abs(e[0]), 1e-30), **args)
    axes[0, 1].axhline(.75, color='gold', ls='--', label='25% shrinkage gate')
    axes[0, 1].axhline(.25, color='red', ls=':', label='positive cutoff')
    axes[1, 1].set_yscale('symlog', linthresh=1e-12)
    for ax, title in zip(axes.flat, ['Mouth separation (radians)', 'Radius / initial radius',
                                   'Combined mouth speed', 'Relative energy residual']):
        ax.set(title=title, xlabel=r'Clock time $\tau/\pi$'); ax.grid(alpha=.15)
    axes[0, 0].legend(fontsize=8, ncol=2)
    axes[0, 1].legend(handles=axes[0, 1].get_legend_handles_labels()[0][-2:], fontsize=8)
    fig.suptitle('Scalar port candidate: no sustained approach or collapse in four returns\n'
                 'Finite constitutive coupling on an unholed sphere; no repeated launch', fontsize=15)
    fig.savefig(image_dir/'comparison.png', dpi=130); plt.close(fig)

    model = CoupledMouthThroat()
    archive = np.load(root/'baseline.npz')
    theta, phi = np.meshgrid(np.linspace(.015, np.pi-.015, 45), np.linspace(0, 2*np.pi, 72), indexing='ij')
    xyz = np.stack([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)], axis=-1)
    basis = real_basis(model.par.lmax, theta, phi)
    sampled = np.unique(np.linspace(0, len(archive['t'])-1, 68).astype(int))
    fields = np.einsum('ijk,tk->tij', basis, archive['states'][sampled, :model.ns])
    scale = 3.
    vmax = np.max(abs(fields))
    fig = plt.figure(figsize=(10.4, 5.8))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], left=.02, right=.97, bottom=.12, top=.85, wspace=.12, hspace=.5)
    ax = fig.add_subplot(grid[:, 0], projection='3d', computed_zorder=False)
    neck = fig.add_subplot(grid[0, 1]); ledger = fig.add_subplot(grid[1, 1])
    times = archive['t'][sampled]
    neckmax = max(1e-4, float(np.max(abs(archive['states'][:, model.ns:model.n]))))
    def frame(i):
        ax.clear(); neck.clear(); ledger.clear()
        state = archive['states'][sampled[i]]; t = times[i]
        q = state[2*model.n:2*model.n+3]
        rad = 1+scale*fields[i]
        surf = rad[..., None]*xyz
        ax.plot_surface(*np.moveaxis(surf, -1, 0), facecolors=plt.cm.coolwarm(.5+.5*fields[i]/vmax),
                        rstride=1, cstride=1, linewidth=0, antialiased=False, shade=False)
        for qi, color, label in zip(q[:2], ['cyan', 'gold'], ['A', 'B']):
            fp = real_basis(model.par.lmax, np.arccos(np.cos(qi)), 0. if np.sin(qi) >= 0 else np.pi)
            radius = 1+scale*(fp@state[:model.ns])
            point = radius*np.array([np.sin(qi), 0., np.cos(qi)])
            ax.scatter(*point, s=65, color=color, edgecolor='black', depthshade=False, zorder=10)
            ax.text(*(point*1.1), label, color=color, fontsize=12, zorder=11)
        ax.set(xlim=(-1.4, 1.4), ylim=(-1.4, 1.4), zlim=(-1.4, 1.4))
        ax.set_box_aspect((1, 1, 1)); ax.set_axis_off(); ax.view_init(15, -75)
        ax.set_title('Scalar surface displacement ×3\nMouth angles shown without magnification', fontsize=10)
        neck.plot(np.linspace(0, model.par.length, model.par.nodes), state[model.ns:model.n], color='orange')
        neck.set(xlim=(0, model.par.length), ylim=(-1.1*neckmax, 1.1*neckmax),
                 xlabel='Resolved throat coordinate', ylabel='Scalar field', title='Same field, finite transmission path')
        neck.grid(alpha=.2)
        mask = b['tau'] <= t+1e-10
        ledger.plot(b['tau'][mask]/np.pi, b['separation'][mask]-b['separation'][0], label='Separation change (rad)', color='cyan')
        ledger.plot(b['tau'][mask]/np.pi, b['radius'][mask]/r0-1, label='Fractional radius change', color='orange')
        ledger.set(xlim=(0, 8), ylim=(-.04, .1), xlabel=r'Clock time $\tau/\pi$')
        ledger.axhline(0, lw=.5, color='gray'); ledger.legend(fontsize=8, loc='upper left'); ledger.grid(alpha=.2)
        fig.suptitle(f'One continuous packet • τ/π = {t/np.pi:.2f} • r/r₀ = {q[2]/r0:.4f}', fontsize=15)
        return []
    fig.text(.5, .025, 'Equal-τ frames in a prescribed tilted clock chart. Constitutive ports; no glued geometry or topology surgery.',
             ha='center', fontsize=9)
    frame(len(sampled)-1); fig.savefig(image_dir/'snapshot.png', dpi=130)
    movie = FuncAnimation(fig, frame, frames=len(sampled), interval=100, blit=False)
    movie.save(image_dir/'evolution.gif', writer=PillowWriter(fps=10)); plt.close(fig)
    paths = list(root.glob('*'))+list(original.glob('*'))+list(image_dir.glob('*'))+[
        Path('geometrodynamics/waves/coupled_mouth_throat.py'),
        Path('experiments/closure_ledger/coupled_mouth_throat_probe.py'),
        Path('experiments/closure_ledger/coupled_mouth_throat_report.py'),
        Path('tests/test_coupled_mouth_throat.py')]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths) if p.name != 'sha256.json'}
    (root/'sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')
    print(json.dumps({'recovery': recovery, 'refinements': refinements}, indent=2), flush=True)


if __name__ == '__main__':
    main()
