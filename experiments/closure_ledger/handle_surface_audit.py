"""Post-hoc audit of #308's immutable initial/final data; no evolution replay."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import handle_evolution as e
from geometrodynamics.waves import handle_surfaces as surfaces
from . import handle_evolution_probe as prior
from .evidence_archive import read_bytes
from .restore_handle_evidence import restore, RAW_SHA256

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'experiments/closure_ledger/runs/20260926_handle_surface_audit'
SOURCES = ('geometrodynamics/waves/handle_surfaces.py',
           'experiments/closure_ledger/handle_surface_audit.py')


def read_inputs(directory=None):
    path = restore() if directory is None else restore(directory)
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != RAW_SHA256:
        raise ValueError('evolution evidence changed')
    old = gzip.decompress(read_bytes(prior.OLD/'refinement.json.gz'))
    if hashlib.sha256(old).hexdigest() != prior.parent.REFINEMENT_HASH:
        raise ValueError('initial evidence changed')
    initial = json.loads(old)
    records = [next(r for r in initial['solutions']
                    if (r['L'],r['eta'],r['n_initial']) == (5.5,eta,513))
               for eta in (0.,.3)]
    return json.loads(raw), records, initial['profiles']['5.5']


def measure(y, n):
    dx = 2*e.L/n
    dy, _ = e.derivatives(y, dx)
    direct = surfaces.expansions(y, dy)
    area = surfaces.area_expansions(y, dx)
    results = {}
    for frame, pair in direct.items():
        labels = surfaces.classify(pair)
        seam_margin = float(-np.max(pair[:,0]))
        radius = y[1] if frame == 'Einstein' else y[1]/np.sqrt(1-(y[4]**2+y[5]**2)/6)
        results[frame] = dict(seam_expansions=pair[:,0].tolist(),
            central_bulk_expansions=pair[:,n//2].tolist(),
            seam_radius=float(radius[0]), central_bulk_radius=float(radius[n//2]),
            seam_future_trapped=bool(labels['future_trapped'][0]),
            seam_past_trapped=bool(labels['past_trapped'][0]),
            seam_near_marginal=bool(labels['near_marginal'][0]),
            seam_negative_margin=seam_margin,
            grid_sphere_counts={k:int(v.sum()) for k,v in labels.items()},
            area_route_max_difference=float(np.max(abs(pair-area[frame]))))
    plus = surfaces.section_flux(y, 0, 1)
    minus = surfaces.section_flux(y, 0, -1)
    return dict(frames=results, cut_faces=dict(positive_normal=plus,negative_normal=minus,
        sum=plus+minus,independent_surfaces=False))


def audit(directory=None):
    data, records, profiles = read_inputs(directory)
    if [g['eta'] for g in data['groups']] != [0.,.3]:
        raise ValueError('unexpected amplitude schedule')
    rows = []
    for group, record in zip(data['groups'], records):
        if [r['N'] for r in group['runs']] != [512,1024,2048]:
            raise ValueError('unexpected grid schedule')
        for run in group['runs']:
            n = run['N']
            _, y0 = e.prepared(record, profiles, n)
            yf = np.asarray(run['final_fields'])
            if yf.shape != (8,n):
                raise ValueError('missing final fields')
            rows.append(dict(eta=group['eta'],N=n,initial=measure(y0,n),
                             final=measure(yf,n)))
    convergence = []
    for eta in (0.,.3):
        group = [r for r in rows if r['eta']==eta]
        for frame in ('Einstein','Jordan'):
            vals = [np.asarray(r['final']['frames'][frame]['seam_expansions']) for r in group]
            differences = [float(np.max(abs(a-b))) for a,b in zip(vals,vals[1:])]
            margin = group[-1]['final']['frames'][frame]['seam_negative_margin']
            convergence.append(dict(eta=eta,frame=frame,seam_differences=differences,
                fine_negative_margin=margin,
                margin_over_medium_fine=margin/max(differences[-1],np.finfo(float).eps)))
    return dict(source_commit='7a9aba6e8ac2c4520075c77f1a14662864f94da3',
        evolution_raw_sha256=RAW_SHA256, initial_raw_sha256=prior.parent.REFINEMENT_HASH,
        audit_source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
        method='post-hoc initial/final snapshot audit; no evolution rerun',
        times=[0.,e.T],sign_tolerance=1e-8,rows=rows,convergence=convergence,
        status=dict(FUTURE_TRAPPED_NECK_AT_FINAL_TIME=all(
            r['final']['frames'][f]['seam_future_trapped'] for r in rows for f in ('Einstein','Jordan')),
            TWO_INDEPENDENT_GRAVITATING_OBJECTS='NOT_TESTED',
            EXTERIOR_BULK_RECIPROCAL_MOMENTUM_TRANSFER='NOT_TESTED',
            DISCRETE_RECIPROCAL_MOMENTUM_EXCHANGE='NOT_TESTED'),
        limitations=['No full field snapshots at intermediate diagnostic times: onset and trapping at each crossing not measured here.',
                     'Counts are grid sections in spherical symmetry, not object or horizon counts.',
                     'A trapped sphere does not identify an event horizon or prove two separate black holes.',
                     'The two cut fluxes use the same section with opposite normals.'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=DEFAULT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    target = args.output_dir/'audit.json'
    target.write_text(json.dumps({'status':{},'error':'incomplete'})+'\n')
    try:
        result = audit()
        target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(json.dumps({'status':result['status'],'convergence':result['convergence']},indent=2))
    except Exception as error:
        target.write_text(json.dumps({'status':{},'error':str(error)})+'\n')
        raise


if __name__ == '__main__':
    main()
