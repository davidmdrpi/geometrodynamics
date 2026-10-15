"""Regenerate the report figure from authenticated production and validation."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from experiments.closure_ledger import r3_half_clock_probe as probe
from experiments.closure_ledger import r3_half_clock_replay as replay


def main(preview=None):
    replay.replay()
    read = lambda name: json.loads((probe.DIRECTORY/name).read_text())
    primary, validation = read('result.json'), read('validation.json')
    points = read('scan.json')['points']
    phase = np.array([p['phase']/(2*np.pi) for p in points])
    lam = np.array([p['lam'] for p in points])
    alt = np.array([p['lam'] for p in validation['rows'][:84]])
    result = validation['result']
    plt.rcParams.update({'font.size': 10, 'svg.hashsalt': 'r3-half-clock'})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    ax = axes[0]
    ax.plot(np.r_[phase, 1], np.r_[lam, lam[0]], color='#166b9d', lw=1.3, label='Scalar half-clock map squared')
    ax.plot(phase, alt, '.', color='#c85a17', ms=3, label='Full matrix return map')
    ax.axhline(0, color='.7', lw=.6)
    ax.set(xlabel='Phase / 2π', ylabel='Obstruction λ (section units)', title='Held-out 3/7 resonance')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax.legend(fontsize=8, loc='upper center', bbox_to_anchor=(.5, -.20), frameon=False)
    ax = axes[1]
    k = np.arange(1, 22)
    ax.semilogy(k, primary['harmonics'][1:22], 'o-', color='#166b9d', ms=3, lw=.8, label='Primary spectrum')
    ax.semilogy(k, result['alternative_harmonics'][1:22], 'x', color='#c85a17', ms=4, label='Independent spectrum')
    ax.axvline(7, color='#877321', lw=.8, alpha=.7, label='Predicted harmonic 7')
    ax.axhline(10*result['combined_resolution'], color='.3', ls='--', lw=.8, label='Detection gate: 10r')
    ax.set(xlabel='Fourier harmonic k', ylabel='|FFT(λ)ₖ| / 84', title='Spectrum and frozen prediction', xticks=[1, 5, 7, 10, 14, 21])
    ax.legend(fontsize=8, loc='upper right')
    path = probe.ROOT/'docs/figures/r3_half_clock.svg'
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, metadata={'Date': None})
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    if preview is not None:
        fig.savefig(preview, dpi=140)
    plt.close(fig)
    print(path.relative_to(probe.ROOT))


if __name__ == '__main__':
    main()
