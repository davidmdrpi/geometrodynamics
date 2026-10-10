"""Prospectively frozen finite-aperture transfer experiment and raw replay."""
import argparse
import base64
from dataclasses import asdict, replace
import hashlib
import io
import json
from pathlib import Path
import platform
import subprocess
import numpy as np
import scipy
from geometrodynamics.transaction import aperture_transfer as a

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'experiments/closure_ledger/runs/20261009_aperture_transfer'
SOURCES=('geometrodynamics/transaction/aperture_transfer.py',
         'experiments/closure_ledger/aperture_transfer_probe.py','docs/aperture_transfer_prereg.md')
ARRAYS=('time','incoming','outgoing','energy','final_q','final_v')


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sources():return {p:digest(ROOT/p) for p in SOURCES}


def schedule():
    cases=[]
    for squash in (.9,1.,1.1):
        for aperture in (.4,.6):
            for carrier in (8.,12.):
                for tag,L,dt in [('coarse',40,1/512),('fine',56,1/1024)]:
                    c=a.Config(squash=squash,aperture=aperture,carrier=carrier,lmax=L,dt=dt)
                    cases.append((f'b{squash:g}_a{aperture:g}_w{carrier:g}_{tag}',c))
    c=a.Config(squash=1.1,lmax=56,dt=1/1024)
    controls={'disconnected':replace(c,connected=False),'reverse_source':replace(c,source=1),
              'time_refine':replace(c,dt=1/2048),'mode_refine':replace(c,lmax=72),
              'extended':replace(c,stop=2.25),'minimal':replace(c,xi=0),
              'round_minimal':replace(c,squash=1,xi=0)}
    for axis in ('horizontal','fiber'):
        for d in (.25,.5):controls[f'{axis}_{d:g}']=replace(c,axis=axis,offset=d)
    cases.extend(controls.items())
    return cases


def save(path,r):
    if path.exists():raise FileExistsError(path)
    buf=io.BytesIO();meta={k:v for k,v in r.items() if k not in ARRAYS}
    np.savez_compressed(buf,metadata=json.dumps(meta,allow_nan=False),**{k:r[k] for k in ARRAYS})
    path.write_text(base64.b64encode(buf.getvalue()).decode()+'\n')


def read(path):
    with np.load(io.BytesIO(base64.b64decode(path.read_text().strip(),validate=True)),allow_pickle=False) as z:
        r=json.loads(str(z['metadata']));r.update({k:z[k].copy() for k in ARRAYS})
    return r


def difference(x,y,swap=False):
    tx=x['time'];ty=y['time'];vy=y['outgoing'][:,::-1] if swap else y['outgoing']
    yy=np.column_stack([np.interp(tx,ty,vy[:,j]) for j in range(2)])
    # Absolute L2 difference normalized by source L2, including weak B output.
    return float(np.linalg.norm(x['outgoing']-yy)/np.linalg.norm(x['incoming']))


def assess(records):
    cases=schedule()
    if list(records)!=[n for n,c in cases]:raise ValueError('wrong inventory/order')
    ds={}
    for name,c in cases:
        if records[name]['config']!=asdict(c):raise ValueError('configuration mismatch')
        ds[name]=a.diagnose(records[name])
    diffs={n:difference(records[n],records[n.replace('_coarse','_fine')]) for n,c in cases if n.endswith('_coarse')}
    base='b1.1_a0.4_w12_fine'
    for variant in ('time_refine','mode_refine','extended'):
        diffs[variant]=difference(records[base],records[variant])
    diffs['reverse_source']=difference(records[base],records['reverse_source'],True)
    # No outcome-dependent peak selection: all twelve primary fine cases count.
    fine=[n for n,c in cases if n.endswith('_fine')]
    ratios={n:ds[n]['capture_fraction']/ds[n.replace(n.split('_')[0],'b1',1)]['capture_fraction'] for n in fine}
    displacement={n:difference(records[base],records[n]) for n in records if n.startswith(('horizontal_','fiber_'))}
    spectrum={};unitarity=0.;reciprocity=0.
    omega=np.pi*(np.arange(1,19)+.37)
    for name,c in cases:
        if name not in fine:continue
        S=a.scattering(c,omega)
        spectrum[name]=dict(omega=omega.tolist(),real=S.real.tolist(),imag=S.imag.tolist())
        unitarity=max(unitarity,float(np.max(abs(S.conj().transpose(0,2,1)@S-np.eye(2)))))
        reciprocity=max(reciprocity,float(np.max(abs(S-S.transpose(0,2,1)))))
    numeric=all(d['valid'] for d in ds.values()) and all(v<(.0000001 if n=='reverse_source' else .03) for n,v in diffs.items()) and unitarity<1e-10 and reciprocity<1e-10
    checks=dict(numerical_validity=bool(numeric),
                finite_capture=all(ds[n]['capture_fraction']>1e-4 for n in fine),
                perturbed_retention=all(v>.1 for v in ratios.values()),
                shifted_return=all(ds[n]['advanced_return_fraction']>1e-5 for n in fine),
                disconnected_control=ds['disconnected']['capture_fraction']<1e-12,
                causal_control=all(ds[n]['causal_return_fraction']<1e-12 for n in fine),
                position_sensitivity=all(displacement[axis+'_0.5']>.01 for axis in ('horizontal','fiber')))
    label='FINITE_APERTURE_TRANSFER_SUPPORTED' if all(checks.values()) else 'FINITE_APERTURE_TRANSFER_FAILED' if numeric else 'NUMERICALLY_UNRESOLVED'
    return dict(label=label,checks=checks,diagnostics=ds,comparisons=diffs,retention=ratios,
                displacement=displacement,spectrum=spectrum,unitarity_error=unitarity,reciprocity_error=reciprocity,
                gr_support='NOT_ESTABLISHED',excised_mouth_matching='NOT_ESTABLISHED',
                gravitational_recoil='NOT_ESTABLISHED',closed_feedback_history='NOT_ESTABLISHED')


def produce(freeze):
    if RUN.exists():raise FileExistsError('append-only evidence exists')
    src=sources()
    for p,h in src.items():
        if hashlib.sha256(subprocess.check_output(['git','show',freeze+':'+p],cwd=ROOT)).hexdigest()!=h:raise ValueError('freeze mismatch')
    RUN.mkdir(parents=True);records={}
    provenance=dict(freeze=freeze,sources=src,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    (RUN/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    for name,c in schedule():
        r=a.simulate(c);save(RUN/(name+'.npz.b64'),r);records[name]=r
        print(name,a.diagnose(r),flush=True)
    if sources()!=src:raise ValueError('sources changed')
    result=assess(records);(RUN/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    names=[n+'.npz.b64' for n,c in schedule()]+['provenance.json','result.json']
    (RUN/'manifest.json').write_text(json.dumps({n:digest(RUN/n) for n in names},indent=2)+'\n')
    print(result['label'],flush=True)


def replay(directory,manifest_sha):
    if digest(directory/'manifest.json')!=manifest_sha:raise ValueError('manifest mismatch')
    manifest=json.loads((directory/'manifest.json').read_text())
    names=[n+'.npz.b64' for n,c in schedule()]+['provenance.json','result.json']
    if set(manifest)!=set(names):raise ValueError('inventory mismatch')
    for n in names:
        if digest(directory/n)!=manifest[n]:raise ValueError('archive mismatch')
    provenance=json.loads((directory/'provenance.json').read_text())
    if provenance['sources']!=sources():raise ValueError('source mismatch')
    fresh=assess({n:read(directory/(n+'.npz.b64')) for n,c in schedule()})
    saved=json.loads((directory/'result.json').read_text())
    if fresh['label']!=saved['label'] or fresh['checks']!=saved['checks']:raise ValueError('decision mismatch')
    return fresh


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--freeze',required=True)
    produce(parser.parse_args().freeze)
