"""Preregistered finite-packet first-transit measurements and replay.

Raw evidence consists of explicit modal solutions from two integrators,
coordinate spatial checks, and quadrature matrices. All claims are rebuilt.
Partial replay validates recorded evidence; --full independently remeasures it.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves import finite_packet as p
from . import finite_packet_spatial as spatial
from . import esu_floquet_probe as old

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20260927_finite_packet'
FREEZE = '2eb1fecb8d9884e768bdc1b4c3bc895f5aac09ce'
BASELINE = '99c94fee6dc6f5158c1ccd8e791ccd51f52e3b5a'
SOURCES = ('geometrodynamics/waves/finite_packet.py', 'geometrodynamics/waves/esu_floquet.py',
           'experiments/closure_ledger/finite_packet_probe.py',
           'experiments/closure_ledger/finite_packet_spatial.py',
           'docs/finite_packet_first_transit_prereg.md')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_modes(run, method, M):
    for j, ix in enumerate(np.array_split(np.arange(79), 8)):
        np.savez_compressed(run/f'modes_{method}_{j}.npz', degrees=p.DEGREES[ix],
                            times=p.TIMES, matrices=M[:, ix])


def load_modes(run, method):
    blocks = []
    for j, ix in enumerate(np.array_split(np.arange(79), 8)):
        with np.load(run/f'modes_{method}_{j}.npz', allow_pickle=False) as a:
            if not np.array_equal(a['degrees'], p.DEGREES[ix]) or not np.array_equal(a['times'], p.TIMES):
                raise ValueError('changed preparation grid')
            b = a['matrices']
            if b.shape != (len(p.TIMES), len(ix), 2, 2) or not np.isfinite(b).all():
                raise ValueError('invalid modal solution')
            blocks.append(b)
    return np.concatenate(blocks, axis=1)


def relative_difference(a, b):
    return float(np.max(np.linalg.norm(a-b, axis=(-2,-1))/np.maximum(1.,np.linalg.norm(b,axis=(-2,-1)))))


def inspect_packet(M, grams, c, phase, even, kind):
    state = p.packet_state(M, c, phase)
    fields, powers = p.observables(state, grams, kind)
    j = int(np.flatnonzero(p.TIMES == np.pi)[0])
    parity = (-1.)**p.DEGREES
    target = -parity[:, None]*state[0]
    error = np.linalg.norm(state[j]-target)/np.linalg.norm(target)
    overlap = np.sum(state[j]*target)/(np.linalg.norm(state[j])*np.linalg.norm(target))
    Etarget = -parity*fields['weyl'][0]
    e_error = np.linalg.norm(fields['weyl'][j]-Etarget)/np.linalg.norm(Etarget)
    def fraction(field, at, destination):
        q = powers[field]
        numerator = q[destination][at]+(q['north' if destination == 'south' else 'south'][at] if even else 0.)
        return float(numerator/q['full'][at])
    initial_metric = fraction('metric', 0, 'north')
    initial_weyl = fraction('weyl', 0, 'north')
    final_weyl = fraction('weyl', j, 'south')
    mean_ratio = (powers['weyl']['south'][j]/p.volume('south'))/(powers['weyl']['belt'][j]/p.volume('belt'))
    window = np.flatnonzero((p.TIMES >= np.pi-.15-1e-14)&(p.TIMES <= np.pi+.15+1e-14))
    peak = window[np.argmax(powers['weyl']['south'][window])]
    # Raw squared powers, not energy and not a claim about a sharp causal front.
    early = p.TIMES <= np.pi-2*p.CAP
    checks = dict(initial_localization=bool(initial_metric >= .5),
                  antipodal_state=bool(error <= .1),
                  weyl_retention=bool(final_weyl/initial_weyl >= .9),
                  target_over_belt=bool(mean_ratio >= 10),
                  arrival_time=bool(abs(p.TIMES[peak]-np.pi) <= .05))
    out = dict(state_error=float(error), state_overlap=float(overlap), weyl_error=float(e_error),
               initial_metric_cap_fraction=initial_metric, initial_weyl_cap_fraction=initial_weyl,
               final_weyl_cap_fraction=final_weyl, weyl_cap_retention=final_weyl/initial_weyl,
               target_over_belt_mean_weyl_power=float(mean_ratio),
               peak_time=float(p.TIMES[peak]), peak_offset=float(p.TIMES[peak]-np.pi),
               initial_south_metric_fraction=float(powers['metric']['south'][0]/powers['metric']['full'][0]),
               initial_south_weyl_fraction=float(powers['weyl']['south'][0]/powers['weyl']['full'][0]),
               early_south_weyl_over_initial_north=float(np.max(powers['weyl']['south'][early])/powers['weyl']['north'][0]),
               initial_weyl_norm=float(np.linalg.norm(fields['weyl'][0])),
               final_weyl_norm=float(np.linalg.norm(fields['weyl'][j])),
               final_stress_norm=float(np.linalg.norm(fields['stress'][j])),
               physical_criteria=checks)
    return out, powers


def score(run=RUN, include_curves=False):
    M, independent = (load_modes(run, method) for method in ('DOP853','Radau'))
    norms = p.normalization()
    grams = p.grams(norms)
    coarse = p.grams(norms,256)
    factors, coarse_factors = p.power_factors(norms), p.power_factors(norms,256)
    # Recompute quadratures even in partial replay; never trust stored matrices.
    with np.load(run/'quadrature.npz', allow_pickle=False) as a:
        quadrature_record_error = float(np.max(abs(a['norms']-norms)/norms))
        for name in p.REGIONS:
            quadrature_record_error = max(quadrature_record_error,
                float(np.max(abs(a['fine_'+name]-grams[name]))),
                float(np.max(abs(a['coarse_'+name]-coarse[name]))))
    spatial_record = json.loads((run/'spatial.json').read_text())
    if set(spatial_record) != {'2','3','5'} or not all(np.isfinite(v) and v>=0 for v in spatial_record.values()):
        raise ValueError('invalid spatial check')
    gram_error = float(np.max(abs(grams['full']-np.eye(79))))
    difference = relative_difference(independent, M)
    scalar = {}
    for method in ('DOP853','Radau'):
        with np.load(run/f'scalar_{method}.npz', allow_pickle=False) as a:
            scalar[method] = (a['matrices'],a['transfer'])
    SM, transfer = scalar['DOP853']
    count = np.count_nonzero(p.TIMES <= np.pi)
    for sm, tr in scalar.values():
        if sm.shape != (count,4,4) or tr.shape != (count,2,4) or not np.isfinite(sm).all() or not np.isfinite(tr).all():
            raise ValueError('invalid scalar evidence')
    scalar_error = max(relative_difference(scalar['Radau'][0],SM),
                       relative_difference(scalar['Radau'][1],transfer))
    s_gain = float(np.max(np.linalg.svd(SM,compute_uv=False)))
    metric_gain = float(np.max(np.linalg.norm(transfer,axis=2))*3/(np.sqrt(2)*np.pi))
    hom_gain = float(np.exp(np.sqrt(2)*np.pi))
    budget = dict(n2_spectral_radius=float(np.max(abs(np.linalg.eigvals(SM[-1])))),
                  n2_max_scaled_state_gain=s_gain, n2_max_metric_transfer_gain=metric_gain,
                  n2_seed_bound_for_point01_metric=.01/metric_gain,
                  homogeneous_gain=hom_gain, homogeneous_seed_bound_for_point01_drift=.01/hom_gain,
                  seed_scan=[dict(seed=s, homogeneous_drift_bound=s*hom_gain,
                                  n2_metric_bound=s*metric_gain) for s in (0.,1e-8,1e-6,1e-4)],
                  independent_difference=scalar_error)
    rows, controls, curves = [], [], {}
    quad_error, scaling_error, free_error = 0., 0., 0.
    for even in (False, True):
        for center,width in p.WINDOWS:
            c = p.coefficients(center,width,norms,even)
            for phase in p.PHASES:
                key = f'{"paired" if even else "cover"}_{center}_{phase:.8f}'
                record, powers = inspect_packet(M,factors,c,phase,even,'supported')
                fields, low = p.observables(p.packet_state(M,c,phase),coarse_factors)
                for field in powers:
                    denominator = max(float(np.max(powers[field]['full'])),1e-300)
                    for region in powers[field]:
                        quad_error = max(quad_error,float(np.max(abs(powers[field][region]-low[field][region])))/denominator)
                        if include_curves:
                            curves[key+'_'+field+'_'+region] = powers[field][region]
                # Explicit linear scaling and zero controls for every preparation.
                state = p.packet_state(M,c,phase)
                # Field amplitudes are linear at all times. Check squared
                # scaling with direct reconstruction at initial, middle, final
                # samples; this is a consistency control, not a nonlinear run.
                sample = np.array([0,len(p.TIMES)//2,len(p.TIMES)-1])
                for amplitude in (0.,)+p.AMPLITUDES:
                    for field in powers:
                        for region in powers[field]:
                            scaled = np.sum(((amplitude*fields[field][sample])@coarse_factors[region])**2,axis=1)
                            expected = amplitude**2*low[field][region][sample]
                            denom = max(amplitude**2*float(np.max(powers[field]['full'])),1e-300)
                            scaling_error = max(scaling_error,float(np.max(abs(scaled-expected)))/denom)
                record.update(preparation=key, even_only=even, center=center,width=width,phase=float(phase),
                              amplitude_scan=[dict(epsilon=e, final_weyl_norm=e*record['final_weyl_norm'],
                                                   final_stress_norm=e*record['final_stress_norm']) for e in p.AMPLITUDES])
                rows.append(record)
                for kind in ('free','bare'):
                    control, _ = inspect_packet(p.analytic_modes(kind),factors,c,phase,even,kind)
                    if kind == 'free':
                        free_error = max(free_error,control['state_error'])
                    controls.append(dict(preparation=key,kind=kind,state_error=control['state_error'],
                                         weyl_error=control['weyl_error'],peak_offset=control['peak_offset']))
    gates = dict(G1=bool(max(spatial_record.values())<1e-9 and gram_error<1e-8),
                 G2=bool(difference<1e-7 and scalar_error<1e-7),
                 G3=bool(quad_error<1e-7 and quadrature_record_error<1e-12),
                 G4=bool(free_error<1e-9 and scaling_error<1e-10))
    for row in rows:
        row['verdict'] = ('UNRESOLVED' if not all(gates.values()) else
                          ('PAIRED_RECURRENCE' if row['even_only'] else 'LOCALIZED_FIRST_TRANSIT')
                          if all(row['physical_criteria'].values()) else 'NOT_ESTABLISHED')
    result = dict(gates=gates, diagnostics=dict(coordinate_checks=spatial_record,gram_error=gram_error,
                    integrator_difference=difference,quadrature_difference=quad_error,
                    quadrature_record_error=quadrature_record_error,free_state_error=free_error,
                    amplitude_squared_error=scaling_error),packets=rows,controls=controls,instability_budget=budget)
    return (result,curves) if include_curves else result


def measure(run=RUN):
    run.mkdir(parents=True,exist_ok=True)
    for method in ('DOP853','Radau'):
        save_modes(run,method,p.modes(method))
        sm,tr=p.scalar_budget(method)
        np.savez_compressed(run/f'scalar_{method}.npz',matrices=sm,transfer=tr)
    (run/'spatial.json').write_text(json.dumps(spatial.run(),indent=2)+'\n')
    save_quadrature(run)


def save_quadrature(run):
    norms=p.normalization()
    data=dict(norms=norms)
    for prefix,nodes in (('fine',512),('coarse',256)):
        data.update({prefix+'_'+k:v for k,v in p.grams(norms,nodes).items()})
    np.savez_compressed(run/'quadrature.npz',**data)


def publish_record(run=RUN):
    result,curves=score(run,True)
    np.savez_compressed(run/'powers.npz',times=p.TIMES,**curves)
    raw = sorted(list(run.glob('*.npz'))+[run/'spatial.json'])
    record=dict(freeze=FREEZE,baseline=BASELINE,sources={s:sha(ROOT/s) for s in SOURCES},
                raw_files={f.name:sha(f) for f in raw},result=result)
    (run/'packet.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    return record


def replay(run=RUN,full=False):
    try:
        archive=json.loads((run/'packet.json').read_text())
        if archive['freeze']!=FREEZE or archive['baseline']!=BASELINE:
            return False
        if archive['sources']!={s:sha(ROOT/s) for s in SOURCES}:
            return False
        required={f'modes_{method}_{j}.npz' for method in ('DOP853','Radau') for j in range(8)}
        required.update({'scalar_DOP853.npz','scalar_Radau.npz','quadrature.npz','powers.npz','spatial.json'})
        if set(archive['raw_files'])!=required or any(sha(run/name)!=digest for name,digest in archive['raw_files'].items()):
            return False
        derived,curves=score(run,True)
        if not old.close(archive['result'],derived,1e-10):
            return False
        with np.load(run/'powers.npz',allow_pickle=False) as a:
            if set(a.files)!=set(curves)|{'times'} or not np.array_equal(a['times'],p.TIMES):
                return False
            if any(not np.allclose(a[k],v,rtol=1e-10,atol=1e-13) for k,v in curves.items()):
                return False
        if full:
            import tempfile
            with tempfile.TemporaryDirectory() as tmp:
                fresh=Path(tmp)
                measure(fresh)
                actual=score(fresh)
                if old.decisions(actual)!=old.decisions(derived) or not old.close(derived,actual,1e-7):
                    return False
                for method in ('DOP853','Radau'):
                    if relative_difference(load_modes(run,method),load_modes(fresh,method))>=1e-7:
                        return False
        return True
    except (ValueError,KeyError,TypeError,OSError,IndexError):
        return False


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--run-dir',type=Path,default=RUN)
    ap.add_argument('--replay',action='store_true')
    ap.add_argument('--full',action='store_true')
    ap.add_argument('--assemble',action='store_true',help='assemble already measured raw evidence')
    args=ap.parse_args()
    if args.replay:
        passed=replay(args.run_dir,args.full)
        print('packet replay:',passed)
        raise SystemExit(0 if passed else 1)
    args.run_dir.mkdir(parents=True,exist_ok=True)
    (args.run_dir/'packet.json').write_text(json.dumps(dict(result=dict(verdict='UNRESOLVED'),error='incomplete'))+'\n')
    if not args.assemble:
        measure(args.run_dir)
    result=publish_record(args.run_dir)['result']
    print(json.dumps(dict(gates=result['gates'],diagnostics=result['diagnostics'],
        packets=[{k:r[k] for k in ('preparation','verdict','physical_criteria','state_error','peak_offset')} for r in result['packets']],
        instability_budget=result['instability_budget']),indent=2))


if __name__=='__main__':
    main()
