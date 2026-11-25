"""Render the vacuum and matter radius histories from the gravity probe."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    root = Path(__file__).resolve().parents[3]
    archive = root/'experiments/closure_ledger/runs/20261125_nariai_loop/probe.json'
    data = json.loads(archive.read_text())
    vacuum = data['vacuum'][1]
    scalar = data['gaussian'][2]['frames']
    plt.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'nariai',
                         'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(1,2,figsize=(9.4,3.7),sharey=True,layout='constrained')
    for i,(label,color) in enumerate((('Contracting center','#b64632'),('Inverted antipode','#146b73'))):
        t = np.linspace(0,vacuum['centers'][i]['proper_time_to_future_pinch'],450)
        time_on_cover = t+(-1)**i*vacuum['tilt']
        r = np.sqrt(np.maximum(0,1-vacuum['eps']*np.cosh(np.sqrt(2)*time_on_cover)))
        axes[0].plot(t,r,color=color,label=label)
        axes[1].plot([f['t'] for f in scalar],[f['centers'][i]['radius'] for f in scalar],color=color)
    axes[1].plot([f['t'] for f in scalar],[f['minimum']['radius'] for f in scalar],
                 color='#403e45',ls='--',label='Minimum on the slice')
    axes[0].set_title('Vacuum control · exact solution',loc='left')
    axes[1].set_title('Scalar packet · 512 points',loc='left')
    for ax in axes:
        ax.set(xlabel='Proper time since initial slice',ylim=(-.02,1.05),xlim=(0,2.5))
        ax.grid(axis='y',alpha=.25)
        ax.legend(frameon=False,loc='lower left')
    axes[0].set_ylabel('Areal radius / Nariai radius')
    axes[1].annotate('Stopped at r = 0.1',xy=(scalar[-1]['t'],.1),
                     xytext=(1.23,.37),arrowprops={'arrowstyle':'->','color':'#484848'})
    fig.savefig(root/'docs/figures/nariai_loop_radii.svg',metadata={'Date':None})
    plt.close(fig)


if __name__ == '__main__':
    main()
