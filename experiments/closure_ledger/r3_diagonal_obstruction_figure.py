"""Plot only the archived diagonal obstruction and its numerical error estimates."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from experiments.closure_ledger import r3_diagonal_obstruction_probe as p


def main():
    rows=json.loads((p.RUN/'scan.json').read_text())['rows'];score=json.loads((p.RUN/'scan_result.json').read_text())
    import gmpy2 as gp
    with gp.context(gp.get_context(),precision=336):
        differences=[float(abs(gp.mpfr(r['primary']['lam'])-gp.mpfr(r['check']['lam']))) for r in rows]
    phi=np.array([r['phi'] for r in rows])/(2*np.pi)
    lam=np.array([float(r['check']['lam']) for r in rows]);scale=10.**np.floor(np.log10(max(abs(lam))))
    plt.rcParams['svg.hashsalt']='r3-diagonal-obstruction'
    fig,ax=plt.subplots(1,2,figsize=(11,4.6),layout='constrained')
    fig.suptitle('Registered result: numerically unresolved (seed-locality gate failed)',fontsize=11)
    ax[0].plot(phi,lam/scale,color='#205b82',lw=1.4,label='336-bit velocity formulation')
    ax[0].scatter(phi[::4],np.array([float(r['primary']['lam']) for r in rows])[::4]/scale,s=15,facecolors='none',edgecolors='#ba652e',label='224-bit momentum formulation')
    ax[0].axhline(0,color='#888',lw=.6);ax[0].set(xlabel='Historical slice phase / $2\\pi$',ylabel=f'Obstruction / {scale:.0e}',title='Measured diagonal-loop obstruction');ax[0].legend(fontsize=8)
    ax[1].semilogy(phi,differences,'.',ms=3,label='Formulation difference',color='#205b82')
    ax[1].semilogy(phi,[float(r['error_estimate']) for r in rows],color='#ba652e',label='Residual-propagated error estimate')
    ax[1].axhline(float(score['F']),color='#555',ls='--',label='Registered resolution floor')
    ax[1].set(xlabel='Historical slice phase / $2\\pi$',ylabel='Absolute coordinate obstruction',title='Error checks and fixed resolution');ax[1].legend(fontsize=8)
    for a in ax:a.grid(alpha=.2)
    for a in ax:a.legend(loc='upper center',bbox_to_anchor=(.5,-.18),fontsize=8)
    path=p.ROOT/'docs/figures/r3_diagonal_obstruction.svg';fig.savefig(path,metadata={'Date':None});plt.close(fig)
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')


if __name__=='__main__':main()
