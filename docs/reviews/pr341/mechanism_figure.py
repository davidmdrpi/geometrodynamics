"""Figure for the PR #341 review: where the leading mixed exchange on pair 1 comes from.

Run from this folder after fixed_axis_mechanism.py, PR #341 checkout on PYTHONPATH.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import fixed_axis_mechanism as mech  # reuses the fixed-axis Strang propagation

x, z = mech.components(.7)
end, n = np.pi+1.2, 2172
h = end/n
t = np.arange(n+1)*h
ez = np.zeros((n+1, 4)); gx = np.zeros((n+1, 4))
for i in range(n+1):
    e_z, _ = mech.fields(z)
    _, g_x = mech.fields(x)
    ez[i], gx[i] = e_z, g_x[:, 1]
    if i < n:
        for st in (x, z):
            mech.kick(st, h/2); mech.drift(st, h); mech.kick(st, h/2)
d = np.load('mechanism_directional.npz')
cum, tc = d['cum'], d['t']

BG, FG = '#ffffff', '#222222'
fig, axes = plt.subplots(3, 1, figsize=(9, 8.2), sharex=True, constrained_layout=True)
c_src, c_tgt = '#166c9e', '#cb6633'
for ax in axes:
    ax.axvspan(np.pi-.6, np.pi+.6, color='gray', alpha=.15, lw=0)
    ax.axvline(np.pi/2, color='#555', ls=':', lw=1)
    ax.grid(alpha=.2)
axes[0].plot(t, ez[:, 1]*1e3, color=c_src, label='−z mouth (source side of the z-packet)')
axes[0].plot(t, ez[:, 0]*1e3, color=c_tgt, label='+z mouth (target side of the z-packet)')
axes[0].set_ylabel('z-packet shape mismatch\n$e^z$ at pair 1 ($10^{-3}$)')
axes[0].legend(frameon=False, fontsize=8, loc='upper left')
axes[0].set_title('Directional preparation, mouth axes fixed (leading order)', fontsize=10)
axes[1].plot(t, gx[:, 1]*1e3, color=c_src, label='at −z')
axes[1].plot(t, gx[:, 0]*1e3, color=c_tgt, ls='--', label='at +z')
axes[1].set_ylabel('x-packet gradient\n$D_y U^x$ at pair 1 ($10^{-3}$)')
axes[1].legend(frameon=False, fontsize=8, loc='upper left')
axes[1].text(np.pi/2+.04, axes[1].get_ylim()[1]*.82, 'x-packet front crosses ±z (t = π/2)', fontsize=8, color='#555')
axes[2].plot(tc, cum[:, 0, 1, 1]*1e4, color=c_src, label='−z mouth: $\\kappa e^z D_yU^x$')
axes[2].plot(tc, cum[:, 0, 0, 1]*1e4, color=c_tgt, label='+z mouth: $\\kappa e^z D_yU^x$')
axes[2].plot(tc, cum[:, :, :2, 1].sum(axis=(1, 2))*1e4, color='k', lw=2, label='pair 1 total (other term ≡ 0)')
axes[2].set_ylabel('cumulative mixed $L_y$\nof pair 1 ($10^{-4}$)')
axes[2].set_xlabel('time (unit sphere); shaded: target focus window π ± 0.6')
axes[2].legend(frameon=False, fontsize=8, loc='upper left')
axes[2].set_xticks([0, 1, np.pi/2, 2, np.pi, 4], ['0', '1', 'π/2', '2', 'π', '4'])
fig.suptitle('PR #341 review: the leading exchange is a mid-transit crossing, not a focus event', fontsize=11)
fig.savefig('mechanism.png', dpi=120)
print('wrote mechanism.png')
