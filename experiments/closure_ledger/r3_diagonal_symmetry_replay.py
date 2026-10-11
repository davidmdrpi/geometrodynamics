"""Authenticate the diagonal-loop archive, rescore it, optionally replay all H edges."""
import argparse
import hashlib
import json
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from experiments.closure_ledger import r3_diagonal_symmetry_probe as probe
from geometrodynamics.waves import r3_diagonal_symmetry as ds

FREEZE = 'efbd82553b69c0ddcfc465406ff30c9f3bb8d92e'
MANIFEST_SHA256 = '60325afd032f044231e7bcca9ad3dc46e6aebf2303c7781b44f704e8b2019a94'


def authenticate(directory):
    directory=Path(directory)
    raw=(directory/'manifest.json').read_bytes()
    if hashlib.sha256(raw).hexdigest()!=MANIFEST_SHA256:
        raise ValueError('manifest fingerprint mismatch')
    manifest=json.loads(raw)
    if manifest['freeze_commit']!=FREEZE or set(manifest['sha256'])!={'actions.json','result.json'}:
        raise ValueError('manifest provenance mismatch')
    loaded={}
    for name,fingerprint in manifest['sha256'].items():
        raw=(directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=fingerprint:
            raise ValueError('archive fingerprint mismatch: '+name)
        loaded[name]=json.loads(raw)
    actions=loaded['actions.json']
    if actions['sources']!=probe.hashes():
        raise ValueError('source fingerprint mismatch')
    checkpoints={f"{r['sample']:02d}.json":hashlib.sha256((json.dumps(
        dict(sources=actions['sources'],record=r),indent=2,allow_nan=False)+'\n').encode()).hexdigest()
        for r in actions['records']}
    if checkpoints!=manifest['checkpoint_sha256']:
        raise ValueError('checkpoint reconstruction mismatch')
    return loaded


def validate(records):
    v=probe.data()
    indices=np.floor(np.arange(24)*len(v)/24).astype(int)
    if len(records)!=24:
        raise ValueError('incomplete sample grid')
    for j,rec in enumerate(records):
        if rec['sample']!=j or rec['index']!=indices[j]:
            raise ValueError('wrong sample ordering')
        for key,shape in [('node',(6,)),('paired_node',(6,)),('chain',(5,6)),
                          ('times',(4,)),('constraints',(4,)),('cycle',(6,)),('reflection',(6,))]:
            a=np.asarray(rec[key])
            if a.shape!=shape or not np.isfinite(a).all():
                raise ValueError('invalid '+key)
        if not np.array_equal(rec['node'],v[indices[j],:6]) or not np.array_equal(rec['paired_node'],v[indices[j],6:]):
            raise ValueError('input node mismatch')
        if not np.array_equal(rec['chain'][0],rec['node']):
            raise ValueError('chain start differs')
        if (rec['radau_half'] is not None)!=(j%4==0):
            raise ValueError('wrong Radau grid')
        for key in ('full_return','half_cycle','matrix_half','radau_half'):
            a=rec[key]
            if a is not None and (np.asarray(a['z']).shape!=(6,) or not np.isfinite(np.r_[a['z'],a['time'],a['constraint']]).all()):
                raise ValueError('invalid map output')
    return probe.score(records,v[:,:6])


def compare(a,b):
    if isinstance(a,dict):
        if not isinstance(b,dict) or a.keys()!=b.keys(): raise ValueError('score keys differ')
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b): raise ValueError('score length differs')
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,(float,np.floating)):
        # Different LAPACK versions can move interpolation residuals by roundoff.
        if not np.isfinite(b) or not np.isclose(a,b,rtol=1e-8,atol=5e-14): raise ValueError('score number differs')
    elif a!=b: raise ValueError('score label differs')


def edge(args):
    z,target=args
    result=ds.clock_map(z,matrix=True)
    return float(np.max(abs(result['z']-target)))


def replay(directory=None,full=False):
    loaded=authenticate(directory or probe.RUN)
    records=loaded['actions.json']['records']
    result=validate(records);compare(result,loaded['result.json'])
    output=dict(replay='VERIFIED',label=result['label'],samples=len(records))
    if full:
        jobs=[(r['chain'][i],r['chain'][i+1]) for r in records for i in range(4)]
        with Pool(4) as pool: errors=pool.map(edge,jobs)
        if not np.isfinite(errors).all() or max(errors)>1e-10:
            raise ValueError('full matrix edge replay failed')
        output.update(full_matrix_edges=len(errors),max_full_matrix_edge_error=max(errors))
    return output


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--full',action='store_true')
    args=parser.parse_args();print(json.dumps(replay(full=args.full),indent=2))
