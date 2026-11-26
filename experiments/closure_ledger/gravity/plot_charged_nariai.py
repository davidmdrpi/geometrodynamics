"""Exact neck turnaround versus the scalar packet's matched-radius events."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .charged_nariai import ROOTS, shell


def main():
    root=Path(__file__).resolve().parents[3]
    archive=root/'experiments/closure_ledger/runs/20261125_charged_nariai/probe.json'
    data=json.loads(archive.read_text())
    turn=data['electrovac']['trajectories'][0]['time_from_outer']
    times=np.linspace(turn-.07,turn+.07,250)
    _,speed,_,_=shell(np.sqrt(ROOTS[1]),0.,turn+.08).sol(times)
    plt.rcParams.update({'svg.fonttype':'none','svg.hashsalt':'charged-nariai',
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(9.7,3.8),layout='constrained')
    ax[0].plot(times-turn,speed,color='#176a70')
    ax[0].scatter([0],[0],s=30,color='#176a70',zorder=3)
    ax[0].set(title='Exact electrovac neck',xlabel='Proper time relative to the inner horizon',
              ylabel='Areal velocity')
    for n,mark,color in ((8192,'o','#555555'),(16384,'+','#d05627'),(32768,'x','#176a70')):
        run=next(r for r in data['characteristics'] if r['n']==n)
        events=[c for c in run['crossings'] if .197<=c['target']<=.23]
        ax[1].plot([c['r'] for c in events],[c['u'] for c in events],marker=mark,
                   fillstyle='none',ls='none',color=color,label=f'{n:,} points')
    ax[1].invert_xaxis()
    ax[1].set(title='Coupled scalar packet',xlabel='Areal radius during contraction',ylabel='Areal velocity',
              ylim=(-1.65,.12))
    ax[1].legend(frameon=False)
    for panel in ax:
        panel.axhline(0,color='#7d7d7d',ls='--',lw=1)
        panel.grid(axis='y',alpha=.2)
    fig.savefig(root/'docs/figures/charged_nariai_turn.svg',metadata={'Date':None})
    plt.close(fig)


if __name__=='__main__':
    main()
