"""Render immutable first-arrival evidence without rerunning propagation."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from experiments.closure_ledger import aperture_transfer_probe as p


def plot():
    result=json.loads((p.RUN/'result.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(11,4.3))
    for b,color in [(.9,'#d8791d'),(1.,'#1769aa'),(1.1,'#8d3c95')]:
        name=f'b{b:g}_a0.4_w12_fine';r=p.read(p.RUN/(name+'.npz.b64'))
        axes[0].plot(r['time'],r['outgoing'][:,1],label=f'b = {b:g}',color=color)
    axes[0].set(xlim=(.65,1.35),xlabel='Exterior time / round antipodal transit',ylabel='Outgoing B characteristic derivative',title='Resolved first-arrival waveform: a = 0.4, w = 12')
    axes[0].legend()
    for angle,color in [(.4,'#1769aa'),(.6,'#d8791d')]:
        for w,style in [(8,'-'),(12,'--')]:
            bs=[.9,1,1.1]
            values=[100*result['diagnostics'][f'b{b:g}_a{angle:g}_w{w}_fine']['capture_fraction'] for b in bs]
            axes[1].plot(bs,values,style,marker='o',color=color,label=f'a = {angle:g}, w = {w}')
    axes[1].set(xlabel='Hopf fiber scale b',ylabel='Captured energy / incident energy (%)',title='Finite aperture and bandwidth change capture')
    axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    fig.text(.02,.015,'Prescribed Berger S3 with distributed scalar ports. Frozen verdict: unresolved; source-label correction disclosed in report.',fontsize=8)
    fig.tight_layout(rect=(0,.06,1,1));path=p.ROOT/'docs/figures/aperture_transfer.png';fig.savefig(path,dpi=160);plt.close(fig)
    print(path)


if __name__=='__main__':plot()
