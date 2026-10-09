"""Plot archived MTY packet histories; does not rerun a trajectory."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from experiments.closure_ledger import mty_packet_probe as p
from geometrodynamics.transaction import mty_packet as m


def plot():
    connected=p.read(p.RUN/'D3_A1_on_fine.npz.b64')
    disconnected=p.read(p.RUN/'D3_A1_off_fine.npz.b64')
    c=m.Config(**connected['config']);t=m.grid(c)
    fig, axes=plt.subplots(2,1,figsize=(9,6.8),sharex=True)
    fig.subplots_adjust(top=.92,bottom=.18,hspace=.1,left=.11,right=.98)
    ax=axes[0]
    ax.plot(t,m.packet(t,c.amplitude),color='0.6',lw=1.3,label='Incoming source')
    ax.plot(t,connected['outgoing'][:,0,2],color='#1769aa',lw=1.6,label='A output: shifted handle')
    ax.plot(t,disconnected['outgoing'][:,0,2],color='#dc7b22',lw=1.2,ls='--',label='A output: handle disconnected')
    ax.set_ylabel('Characteristic derivative')
    ax.set_title('Finite packet with a time-shifted handle: reduced scalar model')
    ax.legend(loc='upper right',fontsize=8)
    tb=c.start+np.arange(len(connected['q']))*c.dt
    for j,color in enumerate(('#1769aa','#8d3c95')):
        axes[1].plot(tb,connected['v'][:,j]*m.MASS[j],color=color,label=f'{"AB"[j]} canonical response momentum')
        axes[1].plot(tb,disconnected['v'][:,j]*m.MASS[j],color=color,ls=':',alpha=.55)
    axes[1].set_ylabel('p = M dq/dt (model units)')
    axes[1].set_xlabel('Exterior time / antipodal transit time')
    axes[1].legend(loc='upper right',fontsize=8)
    for ax in axes:
        ax.axvspan(-.5,.5,color='0.8',alpha=.2)
        ax.axvline(-.5,color='0.4',lw=.8,ls='--')
        ax.grid(alpha=.2);ax.set_xlim(-1.5,4.)
    fig.text(.02,.025,'Shading: source support. D = 3, time advance = 1.5, source amplitude = 1. Dotted momenta: disconnected control.\nResponse coordinates are not gravitational mouth centre-of-mass recoil.',fontsize=8)
    path=p.ROOT/'docs/figures/mty_packet_response.png';path.parent.mkdir(exist_ok=True)
    fig.savefig(path,dpi=160,bbox_inches='tight');plt.close(fig)
    print(path)


if __name__=='__main__':plot()
