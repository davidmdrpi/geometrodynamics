"""PR #343 (cb4b5ec): the density formula lets the walls pass through each other once A > b + r_half.

    PYTHONPATH=scripts python docs/reviews/pr343/crossing_figure.py docs/reviews/pr343/crossing.png
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from neck_wave_contact import ContactNeck, Parameters, LEVEL

E = .025
class Clamped(ContactNeck):
    def radius(self, x, t):
        r = super().radius(x, t); d = E/4
        return (r+np.sqrt(r*r+d*d))/2
x = np.linspace(-.3, .3, 901); rho = np.linspace(-.2, .2, 601)
Xg, Rg = np.meshgrid(x, rho)
fig, axs = plt.subplots(2, 2, figsize=(11, 7.4), constrained_layout=True)
for row, A in enumerate([.18, .30]):
    for col, (cls, name) in enumerate([(ContactNeck, 'PR: signed R'), (Clamped, 'smooth positive part R₊ = (R + √(R² + (ε/4)²))/2')]):
        m = cls(Parameters(amplitude=A)); ax = axs[row, col]
        p = m.density(Xg, Rg, 4.25)
        ax.imshow(p, extent=(x[0], x[-1], rho[0], rho[-1]), origin='lower', cmap='magma', vmin=0, vmax=1, aspect='auto')
        ax.contour(Xg, Rg, p, [LEVEL], colors='w', linewidths=1)
        ax.axhline(0, color='w', lw=.5, ls=':')
        ax.set(title=f'{name}\nA = {A:.2f}, t = 4.25', xlabel='x', ylabel='ρ (signed)')
axs[1, 0].text(0, -.17, 'walls cross the axis and reappear:\nan open gap between two bridges', color='w', ha='center', fontsize=9)
fig.suptitle('Signed R lets the diffuse walls pass through each other once A > b + r½ = 0.1877 (slider reaches 0.30)', fontsize=12)
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else 'crossing.png', dpi=100)
