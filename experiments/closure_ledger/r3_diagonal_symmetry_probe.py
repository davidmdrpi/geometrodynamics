"""Registered group-action test on the archived #319 loop; checkpoint each job."""
import argparse
import hashlib
import json
from multiprocessing import Pool
from pathlib import Path
import platform
import numpy as np
import scipy
from geometrodynamics.waves import r3_diagonal_symmetry as ds

ROOT = Path(__file__).resolve().parents[2]
INPUT = 'experiments/closure_ledger/runs/20261001_r3_family/stage_F.json'
RUN = ROOT/'experiments/closure_ledger/runs/20261011_r3_diagonal_symmetry'
SOURCES = [INPUT, 'geometrodynamics/waves/r3_diagonal_symmetry.py',
           'geometrodynamics/waves/r3_family.py','geometrodynamics/waves/r3_extension.py',
           'geometrodynamics/waves/nonlinear_supported_tt.py',
           'experiments/closure_ledger/r3_diagonal_symmetry_probe.py',
           'docs/r3_diagonal_symmetry_prereg.md', 'tests/test_r3_diagonal_symmetry.py']


def hashes():
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False,
                                   default=lambda a: a.tolist() if isinstance(a,np.ndarray) else a.item())+'\n')
    temporary.replace(path)


def data():
    record = json.loads((ROOT/INPUT).read_text())
    # All 110 continuation points; start is separate, not a duplicated landing.
    return np.array([p['v'] for p in record['points']])


def job(args):
    j, index, v = args
    z = v[:6]
    chain, times, constraints = [z], [], []
    for _ in range(4):
        h = ds.clock_map(chain[-1])
        chain.append(h['z']); times.append(h['time']); constraints.append(h['constraint'])
    p = ds.clock_map(z,half=False)
    cz = ds.spatial(z)
    hc = ds.clock_map(cz)
    matrix = ds.clock_map(z,matrix=True)
    alternate = ds.clock_map(z,matrix=True,method='Radau') if j%4==0 else None
    rec = dict(sample=j,index=index,node=z,paired_node=v[6:],chain=chain,times=times,
               constraints=constraints,full_return=p,cycle=cz,half_cycle=hc,
               matrix_half=matrix,radau_half=alternate,reflection=ds.spatial(z,ds.REFLECTION))
    return rec


def score(records, nodes):
    f20, condition = ds.fourier_fit(nodes,20)
    f26, condition26 = ds.fourier_fit(nodes,26)
    even,_ = ds.fourier_fit(nodes[::2],20)
    odd,_ = ds.fourier_fit(nodes[1::2],20)
    grid = 2*np.pi*np.arange(1024)/1024
    interpolation = max(np.max(abs(f20(ds.angle(nodes))-nodes)),
                        np.max(abs(even(ds.angle(nodes[1::2]))-nodes[1::2])),
                        np.max(abs(odd(ds.angle(nodes[::2]))-nodes[::2])),
                        np.max(abs(f20(grid)-f26(grid))))
    membership_gate = max(20*interpolation,1e-9)
    metrics = {k:0. for k in ('membership_H','membership_C','square','fourth','pair',
                             'commutator','matrix_difference','radau_difference','constraint')}
    turns, deviations, times = [], [], []
    separations = []
    for rec in records:
        chain = np.array(rec['chain']); z=chain[0]
        for name, values in [('membership_H',chain[1:]),('membership_C',np.array([rec['cycle']]))]:
            metrics[name]=max(metrics[name],float(np.max(abs(values-f20(ds.angle(values))))))
        measurements = dict(square=np.max(abs(chain[2]-rec['full_return']['z'])),
                            fourth=np.max(abs(chain[4]-z)),
                            pair=np.max(abs(np.array(rec['full_return']['z'])-rec['paired_node'])),
                            commutator=np.max(abs(np.array(rec['half_cycle']['z'])-ds.spatial(chain[1]))),
                            matrix_difference=np.max(abs(chain[1]-rec['matrix_half']['z'])),
                            constraint=max(rec['constraints']+[rec[k]['constraint'] for k in ('full_return','half_cycle','matrix_half')]))
        if rec['radau_half'] is not None:
            measurements['radau_difference']=np.max(abs(chain[1]-rec['radau_half']['z']))
            measurements['constraint']=max(measurements['constraint'],rec['radau_half']['constraint'])
        for k,v in measurements.items(): metrics[k]=max(metrics[k],float(v))
        increments=np.mod(np.diff(ds.angle(chain)),2*np.pi)/(2*np.pi)
        turn=float(increments.sum());turns.append(turn)
        branch=round(turn)/4
        deviations.append(float(max(abs(increments-branch))))
        separations.append(min(np.max(abs(chain[1]-z)),np.max(abs(chain[2]-z)),np.max(abs(np.array(rec['cycle'])-z))))
        times.append(sum(rec['times']))
    branches = sorted(set(round(x) for x in turns))
    jvalues=ds.chirality(nodes)
    reflection_j=ds.chirality(np.array([r['reflection'] for r in records]))
    gates = dict(interpolation_resolved=bool(interpolation<=5e-9 and membership_gate<=1e-7),
                 same_component=bool(max(metrics['membership_H'],metrics['membership_C'])<=membership_gate),
                 order_four=bool(metrics['fourth']<=1e-7 and min(separations)>1e-2),
                 rotation_branch=bool(len(branches)==1 and branches[0] in (1,3) and max(abs(np.array(turns)-branches[0]))<=1e-6 and max(deviations)<.10),
                 map_identities=bool(max(metrics[k] for k in ('square','pair','commutator'))<=1e-10),
                 independent_maps=bool(max(metrics['matrix_difference'],metrics['radau_difference'])<=1e-10),
                 constraints=bool(metrics['constraint']<=1e-10),
                 reflection_exchanges_chirality=bool((np.max(jvalues)<-1e-3 and np.min(reflection_j)>1e-3) or (np.min(jvalues)>1e-3 and np.max(reflection_j)<-1e-3)))
    if not gates['interpolation_resolved'] or not gates['independent_maps'] or not gates['constraints']:
        label='NUMERICALLY_UNRESOLVED'
    elif all(gates.values()): label='ORDER_TWELVE_SELECTION_SUPPORTED'
    else: label='PROPOSED_GROUP_ACTION_FAILED'
    return dict(label=label,gates=gates,metrics=metrics,interpolation_error=float(interpolation),
                membership_gate=float(membership_gate),fit_conditions=[condition,condition26],
                oriented_half_rotation=branches[0]/4 if len(branches)==1 else None,
                four_step_windings=turns,max_azimuth_increment_deviation=max(deviations),
                min_nontrivial_separation=min(separations),chirality_range=[min(jvalues),max(jvalues)],
                reflected_chirality_range=[min(reflection_j),max(reflection_j)],
                four_half_time_range=[min(times),max(times)],sample_count=len(records),
                radau_count=sum(r['radau_half'] is not None for r in records))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();v=data();indices=np.floor(np.arange(24)*len(v)/24).astype(int)
    fingerprint=hashes(); records=[];todo=[]
    for j,index in enumerate(indices):
        path=RUN/'checkpoints'/f'{j:02d}.json'
        if path.exists():
            saved=json.loads(path.read_text())
            if saved['sources']!=fingerprint or saved['record']['index']!=index:
                raise ValueError('checkpoint provenance mismatch')
            records.append(saved['record'])
        else:todo.append((j,int(index),v[index]))
    with Pool(args.workers) as pool:
        for rec in pool.imap_unordered(job,todo):
            write(RUN/'checkpoints'/f"{rec['sample']:02d}.json",dict(sources=fingerprint,record=rec))
            records.append(rec);print('checkpoint',rec['sample'],flush=True)
    records.sort(key=lambda r:r['sample'])
    archive=dict(sources=fingerprint,environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),records=records)
    write(RUN/'actions.json',archive)
    result=score(records,v[:,:6]);write(RUN/'result.json',result)
    print(json.dumps(result,indent=2,default=lambda a:a.item()))


if __name__=='__main__':main()
