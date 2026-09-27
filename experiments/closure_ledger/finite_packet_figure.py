"""Reproduce the finite-packet report figure from archived raw solutions."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometrodynamics.waves import finite_packet as p
from . import finite_packet_probe as probe


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','svg.hashsalt':'finite-packet-20260927'})
    M=probe.load_modes(probe.RUN,'DOP853')
    norms=p.normalization()
    c=p.coefficients(24,6,norms)
    state=p.packet_state(M,c,0.)
    chi=np.linspace(.00001,np.pi-.00001,1600)
    H=p.radial(chi)/norms[:,None,None]
    fig,axes=plt.subplots(2,2,figsize=(11,7.2),constrained_layout=True)
    colors=('#395f9f','#b6771b','#34836b')
    for eta,color,label in zip((0.,np.pi/2,np.pi),colors,('Initial','Mid-transit','First antipode')):
        j=np.argmin(abs(p.TIMES-eta)); t=p.TIMES[j]
        h=state[j,:,0]; hp=state[j,:,1]*(p.DEGREES+1)
        R,_,f,fp,_=p.fl.background(t)
        k=p.DEGREES*(p.DEGREES+2)
        hpp=-fp*hp/f-(k+2*R*R/f)*h
        E=(k*h-hpp)/(4*f)
        radial=np.einsum('n,ncp->cp',E,H)
        density=np.sum(radial**2,axis=0)*np.sin(chi)**2/(E@E)
        axes[0,0].plot(chi/np.pi,density*np.pi,color=color,label=label)
    axes[0,0].set(xlabel=r'$\chi/\pi$',ylabel='Normalized radial Weyl power',title='Cover packet: center 24, width 6, phase 0')
    axes[0,0].legend(frameon=False)
    with np.load(probe.RUN/'powers.npz') as a:
        for region,color in zip(('north','south','belt'),colors):
            key='cover_24_0.00000000_weyl_'
            axes[0,1].plot(p.TIMES/np.pi,a[key+region]/a[key+'full'],color=color,label=region.title())
        axes[0,1].set(xlabel=r'$\eta/\pi$',ylabel='Fraction of total Weyl power',title='A localized signal arrives on the cover')
        axes[0,1].legend(frameon=False)
        key='paired_24_0.00000000_weyl_'
        for region,color,style in (('north',colors[0],'-'),('south',colors[2],'--')):
            axes[1,0].plot(p.TIMES/np.pi,a[key+region]/a[key+'full'],color=color,linestyle=style,label=region.title())
        axes[1,0].set(xlabel=r'$\eta/\pi$',ylabel='Fraction of total Weyl power',title='Even modes: both lobes exist initially')
        axes[1,0].legend(frameon=False)
        take=abs(p.TIMES-np.pi)<=.1
        for phase,color,label in ((0.,colors[0],'Phase 0'),(np.pi/4,colors[1],r'Phase $\pi/4$')):
            power=a[f'cover_12_{phase:.8f}_weyl_south'][take]
            axes[1,1].plot(p.TIMES[take]-np.pi,power/power.max(),color=color,label=label)
        axes[1,1].axvspan(-.05,.05,color='#34836b',alpha=.10,label='Frozen arrival window')
        axes[1,1].set(xlabel=r'$\eta-\pi$',ylabel='South-cap Weyl power / peak',title='Center 12: phase-dependent peak delay')
        axes[1,1].legend(frameon=False,fontsize=9)
    fig.suptitle('Finite-packet first transit on the supported ESU',fontsize=15)
    output=probe.ROOT/'docs/figures/finite_packet_first_transit.svg'
    output.parent.mkdir(exist_ok=True)
    fig.savefig(output,metadata={'Date':None})
    output.write_text('\n'.join(line.rstrip() for line in output.read_text().splitlines())+'\n')
    fig.savefig('/tmp/finite_packet_first_transit.png',dpi=140)
    print(output)


if __name__=='__main__':
    main()
