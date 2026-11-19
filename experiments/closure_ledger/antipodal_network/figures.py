import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def main(results, output):
    rows=json.loads((results/'transport.json').read_text())
    roots=json.loads((results/'refinement.json').read_text())
    fig,ax=plt.subplots(1,2,figsize=(8.6,3.3),layout='constrained')
    colors=['#939598','#19868a','#223f74']
    for row,col in zip([rows[0],rows[1],rows[-1]],colors):
        t,a,b,et=np.array(row['history']).T
        ax[0].plot(t,b,color=col,label=f"{row['cover_tetrahedra']:,} tetrahedra",lw=2)
    ax[0].axhline(0,color='#414141',lw=.7)
    ax[0].set(title='Net work out through mouth B',xlabel='Model time',ylabel='Fraction of initial energy')
    ax[0].legend(frameon=False,fontsize=8)
    x=np.arange(3)
    ax[1].bar(x-.18,[r['coarse_zero_clock_cross_reaction'] for r in roots],width=.36,color='#939598',label='Zero clock offset')
    ax[1].bar(x+.18,[r['coarse_clock_cross_reaction'] for r in roots],width=.36,color='#19868a',label='Coarse dark setting')
    ax[1].set(title='Holding the calibration fixed',xticks=x,xticklabels=['3,072','24,576','196,608'],xlabel='Cover tetrahedra',ylabel='Cross-reaction magnitude')
    ax[1].legend(frameon=False,fontsize=8)
    for a in ax:
        a.spines[['right','top']].set_visible(False)
        a.grid(axis='y',alpha=.2)
    fig.savefig(output, dpi=160)
    plt.close(fig)
if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Plot the archived or replayed network results")
    parser.add_argument("results", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    main(args.results, args.output)
