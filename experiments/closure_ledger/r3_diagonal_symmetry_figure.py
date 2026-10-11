"""Render the archived action measurements; no integrations or fitted result changes."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from experiments.closure_ledger import r3_diagonal_symmetry_probe as p
from geometrodynamics.waves import r3_diagonal_symmetry as ds


def main():
    rows=json.loads((p.RUN/'actions.json').read_text())['records'];nodes=p.data()[:,:6]
    f,_=ds.fourier_fit(nodes,26);theta=np.linspace(0,2*np.pi,1000);curve=f(theta)
    plt.rcParams['svg.hashsalt']='r3-diagonal-symmetry'
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    ax=axes[0];ax.plot(curve[:,2],curve[:,4],color='#205b82',label='archived component')
    chain=np.array(rows[0]['chain'])
    ax.scatter(chain[:4,2],chain[:4,4],c=range(4),cmap='viridis',s=45,zorder=3)
    for j,z in enumerate(chain[:4]):
        ax.annotate(f'$H^{j}$',(z[2],z[4]),xytext=(5,5),textcoords='offset points')
    cz=chain[0]
    for j in (1,2):
        cz=ds.spatial(cz)
        ax.scatter(cz[2],cz[4],marker='s',facecolors='none',edgecolors='#a35b35',s=60,zorder=3,
                   label='cyclic axis images' if j==1 else None)
        ax.annotate(f'$C^{j}$',(cz[2],cz[4]),xytext=(5,5),textcoords='offset points',color='#a35b35')
    ax.set(xlabel='$x_1$',ylabel='$x_2$',title='Half-clock and cyclic-axis images',
           xlim=(-.105,.105),ylim=(-.105,.12),aspect='equal')
    ax.legend(fontsize=8,loc='upper center');ax.grid(alpha=.2)
    ax=axes[1];a=[];b=[]
    for row in rows:
        chain=np.array(row['chain']);a.append(ds.angle(chain[0])/(2*np.pi)%1)
        b.append((ds.angle(chain[1])-ds.angle(chain[0]))/(2*np.pi)%1)
    ax.scatter(a,b,color='#205b82',s=25,label='measured azimuth increment')
    ax.axhline(.25,color='#666',ls='--',label='rotation number 1/4')
    ax.set(xlabel='Initial geometric azimuth / $2\\pi$',ylabel='Half-step azimuth increment / $2\\pi$',title='Geometric azimuth is not a uniform angle')
    ax.grid(alpha=.2);ax.legend(fontsize=8)
    path=p.ROOT/'docs/figures/r3_diagonal_symmetry.svg';fig.savefig(path,metadata={'Date':None})
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__=='__main__':main()
