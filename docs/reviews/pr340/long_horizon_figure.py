"""Review figure for PR #340: long-horizon pair momenta and angles.

Reads long_{directional,symmetric}.npz written by review_controls.py."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BG, FG = '#0f1722', '#d6e1e8'
C1, C2 = '#6edbd5', '#ffcc66'
fig, axes = plt.subplots(2, 2, figsize=(12, 7.2), dpi=110, facecolor=BG)
for ax in axes.flat:
    ax.set_facecolor(BG)
    ax.tick_params(colors=FG, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color('#506070')
    ax.grid(alpha=.12)
    for k in range(5):
        ax.axvspan(k*np.pi-.6, k*np.pi+.6, color='white', alpha=.07, lw=0)
    ax.set_xlim(0, 4*np.pi+1.2)
    ax.set_xticks([k*np.pi for k in range(5)], ['0', 'π', '2π', '3π', '4π'])
for row, tag in enumerate(['directional', 'symmetric']):
    o = np.load(f'long_{tag}.npz')
    t = o['t']
    ax = axes[row, 0]
    if tag == 'directional':
        ax.plot(t, o['L'][:, 0, 0], color=C1, lw=2, label='pair 1, L_x (target of z-packet)')
        ax.plot(t, o['L'][:, 1, 2], color=C2, lw=2, label='pair 2, L_z (target of x-packet)')
        ax.plot(t, o['L'][:, 0, 1], color=C1, lw=1.2, ls='--', label='pair 1, L_y (packet cross term)')
        ax.plot(t, o['L'][:, 1, 1], color=C2, lw=1.2, ls='--', label='pair 2, L_y (packet cross term)')
    else:
        ax.plot(t, o['L'][:, 0, 1], color=C1, lw=2, label='pair 1, L_y')
        ax.plot(t, o['L'][:, 1, 1], color=C2, lw=2, label='pair 2, L_y')
    ax.axhline(0, color='#506070', lw=.8)
    ax.set_title(f'{tag.capitalize()} packets: pair angular momenta (PR model, M=0.1)', color='white', fontsize=10)
    ax.legend(facecolor=BG, labelcolor='white', fontsize=7, framealpha=.6,
              loc='lower right' if tag == 'directional' else 'upper left')
    ax = axes[row, 1]
    x0 = o['x'][0]
    disp = np.degrees(np.arccos(np.clip(np.sum(o['x']*x0[None], -1), -1, 1)))
    sep = np.degrees(np.arccos(np.clip(np.sum(o['x'][:, 0]*o['x'][:, 1], -1), -1, 1)))-90
    ax.plot(t, disp[:, 0], color=C1, lw=2, label='pair 1 axis displacement')
    ax.plot(t, disp[:, 1], color=C2, lw=2, label='pair 2 axis displacement')
    ax.plot(t, sep, color='#ad9cff', lw=2, label='change in pair separation')
    ax.axhline(0, color='#506070', lw=.8)
    ax.axvline(np.pi+1.2, color='white', lw=.8, ls=':')
    ax.text(np.pi+1.3, .03 if row == 0 else .9, 'PR horizon', color=FG, fontsize=7,
            transform=ax.get_xaxis_transform())
    ax.set_title(f'{tag.capitalize()}: true angles (degrees)', color='white', fontsize=10)
    ax.legend(facecolor=BG, labelcolor='white', fontsize=7, framealpha=.6, loc='upper left' if row == 0 else 'lower left')
for ax in axes[1]:
    ax.set_xlabel('time (R/c); shaded = focus windows kπ ± 0.6 (even k: source side, odd k: target side)', color=FG, fontsize=8)
fig.suptitle('PR #340 model unchanged, run past its horizon to 4π + 1.2 (review control, not a frozen result)',
             color='white', fontsize=11)
fig.tight_layout(rect=(0, 0, 1, .96))
fig.savefig('pr340_long_horizon.png', facecolor=BG)
