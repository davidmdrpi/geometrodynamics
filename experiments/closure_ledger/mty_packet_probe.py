"""Registered reduced MTY packet experiment: immutable production and replay."""
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
from geometrodynamics.transaction import mty_packet as m

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20261008_mty_packet'
SOURCES = ('geometrodynamics/transaction/mty_packet.py',
           'geometrodynamics/transaction/network.py',
           'experiments/closure_ledger/mty_packet_probe.py',
           'docs/mty_packet_prereg.md')
DURATIONS = (0., 2.25, 2.75, 3., 3.25)


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def sources(): return {p:digest(ROOT/p) for p in SOURCES}


def schedule():
    out = []
    for D, connected in [(d,True) for d in DURATIONS]+[(3.,False)]:
        for amp in (1.,2.):
            for res, dt in [('coarse',1/128),('fine',1/256)]:
                name=f'D{D:g}_A{amp:g}_{"on" if connected else "off"}_{res}'
                out.append((name,m.Config(D,amp,dt=dt,connected=connected)))
    for amp in (1.,2.):
        c=m.Config(3.,amp,dt=1/256)
        for name, config in [('extent',replace(c,start=-24.,stop=96.)),
                             ('seed',replace(c,seed=1)),
                             ('midpoint',replace(c,scheme='midpoint')),
                             ('basis',replace(c,basis_B=-1))]:
            out.append((f'D3_A{amp:g}_{name}',config))
    return out


def save(path, record):
    if path.exists(): raise FileExistsError(path)
    arrays={key:np.asarray(record[key]) for key in ('iterate','q','v') if key in record}
    meta={key:value for key,value in record.items() if key not in ('iterate','q','v','incoming','outgoing')}
    buf=io.BytesIO();np.savez_compressed(buf,metadata=json.dumps(meta,allow_nan=False),**arrays)
    path.write_text(base64.b64encode(buf.getvalue()).decode()+'\n')


def read(path):
    with np.load(io.BytesIO(base64.b64decode(path.read_text())),allow_pickle=False) as z:
        r=json.loads(str(z['metadata']))
        for key in ('iterate','q','v'):
            if key in z: r[key]=z[key].copy()
    if 'iterate' in r:
        config=m.Config(**r['config']);t=m.grid(config)
        inc=np.zeros((len(t),2,3));inc[:,:,:2]=r['iterate'][:,:,:2]
        inc[:,0,2]=m.packet(t,config.amplitude)
        r['incoming']=inc;r['outgoing']=(r['v'][:-1]+r['v'][1:])[:,:,None]/2-inc
    return r


def state(record):
    p=np.asarray(record['v'])*m.MASS;q=np.asarray(record['q'])
    out=np.concatenate([q,p],axis=1).copy()
    out[:,[1,3]]*=record['config']['basis_B']
    return out


def compare(a,b):
    """Relative L2 canonical state difference on a's full boundary-time grid."""
    ca,cb=(m.Config(**r['config']) for r in (a,b))
    if cb.start>ca.start or cb.stop<ca.stop: raise ValueError('comparison lacks coverage')
    ta=ca.start+np.arange(len(a['q']))*ca.dt
    tb=cb.start+np.arange(len(b['q']))*cb.dt
    av=state(a);bv=state(b)
    projected=np.array([np.interp(ta,tb,bv[:,j]) for j in range(4)]).T
    denom=np.linalg.norm(av)
    if denom<=1e-14: raise ValueError('degenerate response')
    return float(np.linalg.norm(projected-av)/denom)


def assess(records):
    expected=schedule()
    if list(records)!=[name for name,_ in expected]: raise ValueError('wrong case inventory/order')
    diagnostics={}
    for name,c in expected:
        if records[name]['config']!=asdict(c): raise ValueError('case configuration changed')
        diagnostics[name]=m.diagnose(records[name])
    numeric=all(d['valid'] for d in diagnostics.values())
    comparisons={};mechanism={}
    if all('q' in r for r in records.values()):
        for name,c in expected:
            if name.endswith('_coarse'):
                comparisons[name]=compare(records[name],records[name.replace('_coarse','_fine')])
        for amp in (1.,2.):
            stem=f'D3_A{amp:g}';base=records[stem+'_on_fine']
            for variant in ('extent','seed','midpoint','basis'):
                comparisons[stem+'_'+variant]=compare(base,records[stem+'_'+variant])
            mechanism[stem+'_disconnection_effect']=compare(base,records[stem+'_off_fine'])
        for D in (2.75,3.,3.25):
            a=state(records[f'D{D:g}_A1_on_fine'])
            b=state(records[f'D{D:g}_A2_on_fine'])/2
            mechanism[f'D{D:g}_nonlinear_shape_change']=float(np.linalg.norm(b-a)/np.linalg.norm(a))
    numeric=numeric and bool(comparisons) and all(v <= (1e-5 if k.endswith(('_seed','_basis')) else .02) for k,v in comparisons.items())
    causal=[diagnostics[f'D0_A{amp:g}_on_fine'].get('before_source_energy_fraction',1.) for amp in (1.,2.)]
    disconnected=[diagnostics[f'D3_A{amp:g}_off_fine'].get('before_source_energy_fraction',1.) for amp in (1.,2.)]
    advanced=[diagnostics[f'D{D:g}_A{amp:g}_on_fine'].get('before_source_energy_fraction',0.) for D in (2.75,3.,3.25) for amp in (1.,2.)]
    checks=dict(numerical_validity=bool(numeric),causal_control=bool(max(causal)<1e-9),
                disconnected_control=bool(max(disconnected)<1e-9),
                advanced_response=bool(min(advanced)>1e-5),
                dynamical_coupling=bool(mechanism and all(v>.01 for k,v in mechanism.items() if k.endswith('disconnection_effect'))),
                nonlinear_response=bool(mechanism and all(v>1e-5 for k,v in mechanism.items() if k.endswith('shape_change'))))
    label=('REDUCED_MTY_PACKET_SCATTERING_SUPPORTED' if all(checks.values()) else
           'REDUCED_MTY_PACKET_SCATTERING_FAILED' if numeric else 'NUMERICALLY_UNRESOLVED')
    return dict(label=label,checks=checks,diagnostics=diagnostics,comparisons=comparisons,
                mechanism=mechanism,gr_traversable_support='NOT_ESTABLISHED',
                gravitating_mouth_recoil='NOT_ESTABLISHED',unique_history='NOT_ESTABLISHED',
                action_discreteness='NOT_ESTABLISHED')


def produce(freeze):
    if RUN.exists(): raise FileExistsError('append-only run exists')
    src=sources()
    for path,sha in src.items():
        prior=subprocess.check_output(['git','show',freeze+':'+path],cwd=ROOT)
        if hashlib.sha256(prior).hexdigest()!=sha: raise ValueError('not the frozen implementation')
    RUN.mkdir(parents=True);records={}
    meta=dict(freeze=freeze,sources=src,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    (RUN/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    for name,config in schedule():
        r=m.solve(config);save(RUN/(name+'.npz.b64'),r)
        records[name]=read(RUN/(name+'.npz.b64'))
        d=m.diagnose(records[name]);print(name,r['status'],d,flush=True)
    if sources()!=src: raise ValueError('sources changed during production')
    result=assess(records)
    (RUN/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    names=[name+'.npz.b64' for name,_ in schedule()]+['provenance.json','result.json']
    (RUN/'manifest.json').write_text(json.dumps({name:digest(RUN/name) for name in names},indent=2)+'\n')
    print(result['label'],flush=True)


def replay(directory, manifest_sha):
    if not manifest_sha or digest(directory/'manifest.json')!=manifest_sha: raise ValueError('manifest fingerprint mismatch')
    manifest=json.loads((directory/'manifest.json').read_text())
    names=[name+'.npz.b64' for name,_ in schedule()]+['provenance.json','result.json']
    if set(manifest)!=set(names): raise ValueError('archive inventory')
    for name in names:
        if digest(directory/name)!=manifest[name]: raise ValueError('evidence fingerprint mismatch')
    provenance=json.loads((directory/'provenance.json').read_text())
    if provenance['sources']!=sources(): raise ValueError('measured sources changed')
    fresh=assess({name:read(directory/(name+'.npz.b64')) for name,_ in schedule()})
    saved=json.loads((directory/'result.json').read_text())
    if fresh['label']!=saved['label'] or fresh['checks']!=saved['checks']: raise ValueError('decisions changed')
    return fresh


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=['run','replay']);ap.add_argument('--freeze');ap.add_argument('--manifest-sha')
    a=ap.parse_args()
    if a.mode=='run':
        if not a.freeze: ap.error('--freeze required')
        produce(a.freeze)
    else: print(json.dumps(replay(RUN,a.manifest_sha),indent=2))
