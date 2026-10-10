"""Frozen phase family: append-only production, hash authentication and replay."""
import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import numpy as np
import scipy
from geometrodynamics.transaction import aperture_phase as a
from experiments.closure_ledger.aperture_transfer_probe import save, read, digest, difference

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'experiments/closure_ledger/runs/20261010_aperture_phase'
SOURCES=('geometrodynamics/transaction/aperture_phase.py',
         'experiments/closure_ledger/aperture_phase_probe.py',
         'docs/aperture_phase_prereg.md',
         'geometrodynamics/transaction/aperture_transfer.py',
         'experiments/closure_ledger/aperture_transfer_probe.py')


def utc(): return datetime.now(timezone.utc).isoformat()
def sources(): return {p:digest(ROOT/p) for p in SOURCES}
def name(b,w,f,tag): return f'b{b:g}_w{w}_f{f:g}_{tag}'


def schedule():
    cases=[]
    for footprint,bs in [(4.8,(.8,.9,.98,1.,1.02,1.1,1.2)),(7.2,(.8,1.))]:
        for b in bs:
            for w in (12,24,48):
                scale=w//12
                for tag,L,dt in [('coarse',40*scale,1/(512*scale**2)),
                                 ('fine',56*scale,1/(1024*scale**2))]:
                    cases.append((name(b,w,footprint,tag),a.Config(b,w,footprint,L,dt)))
    c=a.Config(.8,48,4.8,224,1/16384)
    cases.extend([('time_refine',replace(c,dt=1/32768)),
                  ('mode_refine',replace(c,lmax=288)),
                  ('extended',replace(c,stop=2.25))])
    return cases


def fit_phase(ds,footprint,tag):
    keys=[name(.8,w,footprint,tag) for w in (12,24,48)]
    phi=np.array([ds[n]['phase']['proxy'] for n in keys])
    retention=np.array([ds[n]['capture']/ds[name(1.,w,footprint,tag)]['capture']
                        for n,w in zip(keys,(12,24,48))])
    if np.any(retention<=0):return dict(beta=None,max_log_residual=None,pass_law=False)
    x=np.log(phi);y=np.log(retention)
    slope,intercept=np.polyfit(x,y,1)
    residual=float(np.max(abs(y-(intercept+slope*x))))
    return dict(beta=float(-slope),max_log_residual=residual,
                pass_law=bool(.5<=-slope<=1.5 and residual<=.25))


def score(ds,comparisons):
    expected=[n for n,c in schedule()]
    if list(ds)!=expected:raise ValueError('diagnostic inventory mismatch')
    fine=[n for n,c in schedule() if n.endswith('_fine')]
    configs=dict(schedule())
    retention={n:ds[n]['capture']/ds[name(1.,configs[n].carrier,configs[n].footprint,'fine')]['capture']
               for n in fine}
    fits={str(f):{tag:fit_phase(ds,f,tag) for tag in ('coarse','fine')} for f in (4.8,7.2)}
    beta_convergence=all(v['fine']['beta'] is not None and v['coarse']['beta'] is not None
                         and abs(v['fine']['beta']-v['coarse']['beta'])<.1 for v in fits.values())
    numeric=all(d['valid'] for d in ds.values()) and all(x<.03 for x in comparisons.values()) and beta_convergence
    # Geometry/source-weighted phase spread, independent of captured waveform.
    phase_coverage=all(ds[name(.8,48,f,'fine')]['phase']['weighted_std']>=5
                       and ds[name(.8,48,f,'fine')]['phase']['weighted_std'] /
                       ds[name(.8,12,f,'fine')]['phase']['weighted_std']>=3 for f in (4.8,7.2))
    candidates=[n for n in fine if configs[n].squash!=1 and ds[n]['phase']['proxy']>=10]
    robust=all(retention[n]>=.1 for n in candidates)
    inverse_law=all(v['fine']['pass_law'] for v in fits.values())
    def label(passed):
        return ('NUMERICALLY_UNRESOLVED' if not numeric else 'PHASE_RANGE_INSUFFICIENT'
                if not phase_coverage else 'SUPPORTED_IN_DECLARED_FAMILY' if passed
                else 'FAILED_IN_DECLARED_FAMILY')
    return dict(numerical_validity=bool(numeric),phase_coverage=bool(phase_coverage),
                retention_label=label(robust),inverse_phase_label=label(inverse_law),
                retention=retention,fits=fits,comparisons=comparisons,
                beta_convergence=bool(beta_convergence),diagnostics=ds,
                r3_propagation='NOT_TESTED',physical_scale_extrapolation='NOT_ESTABLISHED',
                closed_feedback='NOT_ESTABLISHED')


def assess(records):
    cases=schedule()
    if list(records)!=[n for n,c in cases]:raise ValueError('record inventory mismatch')
    ds={}
    for n,c in cases:
        if records[n]['config']!=asdict(c):raise ValueError('configuration mismatch')
        ds[n]=a.diagnose(records[n])
    comparisons={n:difference(records[n],records[n.replace('_coarse','_fine')])
                 for n,c in cases if n.endswith('_coarse')}
    base=name(.8,48,4.8,'fine')
    for variant in ('time_refine','mode_refine','extended'):
        comparisons[variant]=difference(records[base],records[variant])
    return score(ds,comparisons)


def produce(freeze):
    if RUN.exists():raise FileExistsError('append-only run already exists')
    src=sources()
    for p,h in src.items():
        actual=hashlib.sha256(subprocess.check_output(['git','show',freeze+':'+p],cwd=ROOT)).hexdigest()
        if actual!=h:raise ValueError('freeze mismatch')
    RUN.mkdir(parents=True)
    provenance=dict(freeze=freeze,sources=src,started_utc=utc(),python=platform.python_version(),
                    numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform())
    (RUN/'started.json').write_text(json.dumps(provenance,indent=2)+'\n')
    records={};ds={};timings={}
    for n,c in schedule():
        start=utc();r=a.simulate(c);finish=utc()
        save(RUN/(n+'.npz.b64'),r);records[n]=r;ds[n]=a.diagnose(r)
        timings[n]=dict(started_utc=start,finished_utc=finish)
        print(n,json.dumps(ds[n]),flush=True)
    comparisons={n:difference(records[n],records[n.replace('_coarse','_fine')])
                 for n,c in schedule() if n.endswith('_coarse')}
    for v in ('time_refine','mode_refine','extended'):
        comparisons[v]=difference(records[name(.8,48,4.8,'fine')],records[v])
    result=score(ds,comparisons)
    if sources()!=src:raise ValueError('sources changed')
    (RUN/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    provenance.update(finished_utc=utc(),case_times=timings)
    (RUN/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    files=[n+'.npz.b64' for n,c in schedule()]+['started.json','provenance.json','result.json']
    (RUN/'manifest.json').write_text(json.dumps({p:digest(RUN/p) for p in files},indent=2)+'\n')
    print(result['retention_label'],result['inverse_phase_label'],flush=True)


def replay(directory,manifest_sha):
    if digest(directory/'manifest.json')!=manifest_sha:raise ValueError('manifest mismatch')
    manifest=json.loads((directory/'manifest.json').read_text())
    files=[n+'.npz.b64' for n,c in schedule()]+['started.json','provenance.json','result.json']
    if set(manifest)!=set(files):raise ValueError('manifest inventory mismatch')
    for p in files:
        if digest(directory/p)!=manifest[p]:raise ValueError('archive mismatch')
    provenance=json.loads((directory/'provenance.json').read_text())
    if provenance['sources']!=sources():raise ValueError('source hash mismatch')
    result=assess({n:read(directory/(n+'.npz.b64')) for n,c in schedule()})
    saved=json.loads((directory/'result.json').read_text())
    for key in ('numerical_validity','phase_coverage','retention_label','inverse_phase_label','beta_convergence'):
        if result[key]!=saved[key]:raise ValueError('decision mismatch')
    # Do not accept agreement of labels alone as numeric reproducibility.
    for n in result['diagnostics']:
        for key in ('capture','remaining','reflected'):
            if abs(result['diagnostics'][n][key]-saved['diagnostics'][n][key])>1e-10:
                raise ValueError('measurement mismatch')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze');parser.add_argument('--manifest-sha')
    args=parser.parse_args()
    if bool(args.freeze)==bool(args.manifest_sha):parser.error('choose production or authenticated replay')
    if args.freeze:produce(args.freeze)
    else:
        r=replay(RUN,args.manifest_sha)
        print(json.dumps({k:r[k] for k in ('numerical_validity','phase_coverage','retention_label','inverse_phase_label','fits')},indent=2))
