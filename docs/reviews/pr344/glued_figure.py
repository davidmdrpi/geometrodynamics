"""Figure for the PR #344 review from glued_drive.npz (run glued_drive.py first).

    PYTHONPATH=.:scripts python docs/reviews/pr344/glued_figure.py docs/reviews/pr344/glued.png
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from klein_release_model import KleinRelease

m = KleinRelease(); d = np.load('glued_drive.npz'); t = d['t']; pk = float(d['peak']); E0 = d['klein_E'][0]
fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), constrained_layout=True)
ax = axs[0]
ax.plot(t, d['modal']/pk, color='0.6', lw=3, label='S² modes: u(A, t)/peak (PR mouth trace)')
ax.plot(t, d['bare_pole']/pk, color='k', lw=1, ls='--', label='bare-sphere FE at the pole (validation)')
ax.plot(t, m.transmitted(0., t), color='#c0392b', lw=2, label='PR readout at the waist, w(0, t)')
ax.plot(t, d['klein_waist']/pk, color='#1f5fa8', lw=2, label='glued Klein surface: u at the waist / peak')
thr = (m.neck.par.waist-m.neck.level_radius)/m.gain
ax.axhline(thr, color='#c0392b', lw=.8, ls=':'); ax.text(1.62, thr+.03, f'contact needs w > {thr:.3f} at gain {m.gain}', color='#c0392b', fontsize=9)
ax.axvline(m.tc, color='0.5', lw=.8); ax.text(m.tc+.02, -.45, 'PR t$_c$', fontsize=9, color='0.4')
ax.set(xlabel='t', ylabel='normalised field', title='(a) what drives the throat', xlim=(0, 2.6)); ax.legend(fontsize=8.5, loc='upper left'); ax.grid(alpha=.3)
ax = axs[1]
ax.plot(t, d['klein_E_neck']/E0, color='#1f5fa8', lw=2, label='glued: energy inside the neck')
ax.plot(t, d['klein_E_B']/E0, color='#138d75', lw=2, label='glued: energy on the sphere within 1 rad of B')
ax.plot(t, d['bare_E_B']/d['bare_E'][0], color='0.5', lw=1.5, ls='--', label='bare S²: within 1 rad of B (arrives near t = 4.3)')
ax.axvspan(m.tc, m.tc+m.withdrawal, color='#c0392b', alpha=.1); ax.text(m.tc+.02, .9, 'PR contact →\nwithdrawal', fontsize=9, color='#c0392b')
ax.set(xlabel='t', ylabel='fraction of incident energy', title='(b) on the glued surface the packet goes through the throat', xlim=(0, 2.6), ylim=(0, 1)); ax.legend(fontsize=8.5, loc='center left'); ax.grid(alpha=.3)
fig.suptitle("PR #344's incident packet on #338's glued Klein surface (same initial data, conformal operator, repository cotangent FE)", fontsize=11.5)
fig.savefig(sys.argv[1] if len(sys.argv) > 1 else 'glued.png', dpi=105)
