"""Render six research lines directly from the pinned capstone evidence."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]

def archived(n):
    return json.loads((ROOT/f'docs/releases/evidence-v1.3.5/pr{n}.json').read_text())

def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig, axes=plt.subplots(3,2,figsize=(12,12),constrained_layout=True)
    ax=axes[0,0];d=archived(304)['metrics']
    for key,label,color in [('physical_H','Hamiltonian','#226b98'),('physical_M','Momentum','#c15127')]:
        ax.semilogy([1,2,3],d[key],'o-',label=label,color=color)
    ax.set(title='A. Vacuum handle constraints',xlabel='Successive coordinate refinements',ylabel='Normalized residual',xticks=[1,2,3]);ax.legend(frameon=False)
    ax.text(.04,.06,'Compact twisted initial data; eight gates pass',transform=ax.transAxes,fontsize=9)
    ax=axes[0,1];rows=archived(311)['result']['packets'];rows=[r for r in rows if not r['even_only'] and r['phase']==0]
    ax.plot([r['center'] for r in rows],[100*r['final_weyl_cap_fraction'] for r in rows],'o-',lw=2,color='#226b98')
    ax.set(title='B. Finite tensor packet refocusing',xlabel='Spectral window centre',ylabel='Target-cap squared-Weyl fraction (%)',ylim=(65,102))
    ax.text(.04,.06,'Linear S³-cover packets; zero initial phase\nCurvature diagnostic, not absorbed energy',transform=ax.transAxes,fontsize=9)
    ax=axes[1,0];d=archived(325);h=np.asarray(d['harmonics'])
    ax.semilogy(np.arange(1,16),h[1:16],'o-',color='#55834b')
    ax.axhline(100*d['noise'],color='#c15127',ls='--',label='Registered detection gate')
    ax.set(title='C. Held-out half-clock prediction',xlabel='Angular harmonic',ylabel='Obstruction Fourier amplitude',xticks=[1,3,5,7,9,11,13,15]);ax.legend(frameon=False,fontsize=9)
    ax.set_title('C. Held-out half-clock prediction',pad=28)
    ax.text(.5,1.02,'Prediction: dominant seventh harmonic; 14 sign changes',ha='center',transform=ax.transAxes,fontsize=9)
    ax=axes[1,1];d=archived(345);band=d['bands'][0];rows=d['router'][:3]
    ax.axvspan(band['two_port'],band['u_k'],color='#55834b',alpha=.12,label='Derived closure band')
    ax.semilogy([r['speed'] for r in rows],[r['speed_drift'] for r in rows],'o',color='#226b98',ms=7)
    ax.set(title='D. Strict-cycle band confirmation',xlabel='Birth speed',ylabel='Cycle-to-cycle speed drift')
    ax.legend(frameon=False,fontsize=9,loc='upper right')
    ax.text(.04,.45,'Below: drift ≈ 0.030\nInside and upper edge: roundoff closure\nStated packet router, unit lapse, k=1',transform=ax.transAxes,fontsize=9)
    ax=axes[2,0];data=ROOT/'experiments/closure_ledger/runs';rows=json.loads((data/'20261028_maxwell_throat/scattering.json').read_text())
    for depth,color in [(0.,'#226b98'),(.4,'#55834b'),(.8,'#c15127')]:
        sub=[r for r in rows if r['angle']==.3 and r['ell']==1 and r['charge_fraction']==0 and r['lapse_depth']==depth]
        ax.semilogy([r['omega'] for r in sub],[r['transmission'] for r in sub],'o-',color=color,label=f'Neck lapse {1-depth:.1f}')
    ax.set(title='E. Geometry-derived throat filtering',xlabel='Frequency (R=c=1)',ylabel='Finite-port transmission fraction');ax.legend(frameon=False,fontsize=9)
    ax=axes[2,1];rows=json.loads((data/'20261028_einstein_throat/histories.json').read_text())
    for fraction,color in [(0.,'#226b98'),(.9,'#c15127')]:
        row=next(r for r in rows if r['angle']==.3 and r['charge_fraction']==fraction);h=row['history']
        ax.plot(np.array(h['tau'])/row['endpoint_proper_time'],np.array(h['radius'])/row['initial_neck'],lw=2,color=color,label='Neutral' if fraction==0 else 'Charged, q/qmax=0.9')
    ax.set(title='F. Conditional Einstein radius evolution',xlabel='Central proper time / (πM)',ylabel='Neck radius / initial radius');ax.legend(frameon=False,fontsize=9)
    ax.text(.03,.04,'Unsupported electrovacuum completion\nEndpoints: singularity / inner horizon',transform=ax.transAxes,fontsize=9)
    fig.suptitle('BAM v1.3.5 · constrained geometry, transport and response',fontsize=17)
    target=ROOT/'docs/images/releases/bam_capstone.png'
    pending=target.with_suffix('.tmp.png')
    fig.savefig(pending,dpi=155)
    pending.replace(target)
    plt.close(fig)

if __name__=='__main__':main()
