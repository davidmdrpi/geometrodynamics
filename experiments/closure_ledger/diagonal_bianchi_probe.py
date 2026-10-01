"""Prospective diagonal continuation and closure; frozen at 0c51f9b."""
import argparse
import base64
import gzip
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from geometrodynamics.waves import diagonal_bianchi as d

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20260929_diagonal_bianchi'
SOURCES = ['docs/diagonal_bianchi_prereg.md', 'geometrodynamics/waves/diagonal_bianchi.py',
           'geometrodynamics/waves/nonlinear_supported_tt.py',
           'experiments/closure_ledger/diagonal_bianchi_probe.py']
FACTOR = 2**.25


def serial(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    raise TypeError(type(x).__name__)


def source_hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def encode(data):
    return base64.b64encode(gzip.compress(json.dumps(data, default=serial, allow_nan=False,
                                                   separators=(',', ':')).encode(), mtime=0))+b'\n'


def decode(path):
    return json.loads(gzip.decompress(base64.b64decode(Path(path).read_bytes())))


def write(path, data):
    with path.open('xb') as stream:
        stream.write(encode(data))


def array(value, shape):
    a = np.asarray(value)
    if a.shape != shape or not np.isfinite(a).all():
        raise ValueError('missing, nonfinite or wrong-shaped evidence')
    return a


def history_metrics(rec, initial):
    start = array(rec['initial'], (29,))
    end = array(rec['returned'], (29,))
    samples = array(rec['samples'], (257, 29))
    times = array(rec['times'], (257,))
    expected = np.asarray(initial).copy(); expected[2] = 0.
    if not np.allclose(start, expected, rtol=1e-11, atol=1e-12):
        raise ValueError('full-system preparation changed')
    if not np.allclose(samples[[0, -1]], [start, end], rtol=1e-11, atol=1e-12):
        raise ValueError('full history endpoints changed')
    if not 0 < rec['return_time'] < 8 or not np.allclose(times, np.linspace(0., rec['return_time'], 257), atol=1e-12, rtol=0):
        raise ValueError('invalid return time schedule')
    constraints = max(float(np.max(abs(d.full.constraints(y)['residual']))) for y in samples)
    chart = min(min(y[0], d.full.ingredients(y)[0], np.linalg.eigvalsh(d.full.unpack(y)[4]).min()) for y in samples)
    return dict(constraint=constraints, chart=chart, section=abs(end[3]), negative_clock=bool(end[7] < 0))


def circle_metrics(rec):
    if 'error' in rec:
        return dict(accepted=False, reason=rec['error'])
    n = len(rec['K'])
    if n not in (63, 95): raise ValueError('wrong circle grid')
    K = array(rec['K'], (n, 6)); image = array(rec['image'], (n, 6))
    omega = float(rec['omega'])
    if not np.isfinite([omega, rec['a']]).all(): raise ValueError('nonfinite circle parameters')
    array(rec['return_times'], (n,))
    mid = 2*np.pi*(np.arange(n)+.5)/n
    points = array(rec['offgrid_points'], (n, 6))
    if not np.allclose(points, d.spectral(K, mid), rtol=0, atol=1e-12): raise ValueError('changed off-grid points')
    array(rec['offgrid_times'], (n,))
    off = array(rec['offgrid_image'], (n, 6))
    T, _ = d.shift_matrix(n, omega)
    grid_error = float(np.max(abs(image-T@K)))
    off_error = float(np.max(abs(off-d.spectral(K, mid+omega))))
    c = np.fft.fft(K-K.mean(0), axis=0)/n
    tail = float(np.max(abs(c[abs(np.fft.fftfreq(n, 1/n)) >= n//2-2])))
    first = np.fft.fft(K[:, 2])[1]/n
    amplitude_error = float(max(abs(first.real-rec['a']/2), abs(first.imag)))
    points = array(rec['full_points'], (4, 6))
    if not np.allclose(points, d.spectral(K, np.arange(4)*np.pi/2), rtol=0, atol=1e-12): raise ValueError('changed validation points')
    phase = array(rec['full_phase_image'], (4, 6)); times = array(rec['full_phase_times'], (4,))
    if len(rec['full_histories']) != 4: raise ValueError('missing full validation')
    diagnostics = []
    for z, target, t, hist in zip(points, phase, times, rec['full_histories']):
        met = history_metrics(hist, d.full_initial(z))
        met['formulation_error'] = float(max(np.max(abs(d.canonical(hist['returned'])-target)), abs(hist['return_time']-t)))
        diagnostics.append(met)
    action = abs(d.action(K)); area = d.circular_area(K)
    accepted = (action > 0 and area < 0 and 0 < omega < 2*np.pi and grid_error <= 1e-9 and off_error <= 1e-9
                and tail <= 1e-9 and amplitude_error <= 1e-10
                and all(h['formulation_error'] <= 1e-8 and h['constraint'] <= 1e-9 and h['chart'] > 0
                        and h['section'] <= 1e-10 and h['negative_clock'] for h in diagnostics))
    return dict(accepted=bool(accepted), action=action, circular_area=area, omega=omega, rho=1+omega/(2*np.pi),
                grid_error=grid_error, offgrid_error=off_error, fourier_tail=tail, amplitude_error=amplitude_error,
                full_checks=diagnostics)


def closure_metrics(rec, minimum_distance):
    if 'error' in rec: return dict(verified=False, reason=rec['error'])
    nodes = array(rec['nodes'], (2, 6)); image = array(rec['image'], (2, 6))
    residual = float(np.max(abs(image-nodes[::-1])))
    separation = float(np.linalg.norm(image[0]-nodes[0]))
    distance = float(np.min(np.linalg.norm(nodes-d.STAR, axis=1)))
    candidate = residual <= 1e-9 and separation > 1e-5 and distance > minimum_distance/4
    checks = []
    if candidate and 'validation' in rec:
        for method in ('DOP853', 'Radau'):
            records = rec['validation'][method]
            if len(records) != 2: raise ValueError('missing closure segment')
            start = d.full_initial(nodes[0])
            for i, hist in enumerate(records):
                if hist['method'] != method: raise ValueError('wrong closure validation method')
                met = history_metrics(hist, start)
                met['node_error'] = float(np.max(abs(d.canonical(hist['returned'])-nodes[1-i])))
                checks.append(met)
                start = np.array(hist['returned'])
            initial, final = d.full_initial(nodes[0]), np.asarray(records[-1]['returned']).copy()
            final[2] = initial[2]
            checks[-1]['full_closure'] = float(np.max(abs(final-initial)))
    verified = candidate and len(checks) == 4 and all(
        x['node_error'] <= 1e-8 and x.get('full_closure', 0.) <= 1e-8 and x['constraint'] <= 1e-9
        and x['chart'] > 0 and x['section'] <= 1e-10 and x['negative_clock'] for x in checks)
    return dict(candidate=bool(candidate), verified=bool(verified), residual=residual,
                separation=separation, minimum_distance=distance, full_checks=checks)


def analyze(raw):
    if raw['freeze'] != d.FREEZE: raise ValueError('wrong freeze')
    metrics = []; good = []; expected = .004; previous = None; retry = False
    for i, rec in enumerate(raw['circles']):
        if abs(rec['a']-expected) > 1e-12: raise ValueError('changed continuation schedule')
        m = circle_metrics(rec); metrics.append(m)
        if m['accepted']:
            good.append(i); previous=rec['a']; expected=previous*FACTOR; retry=False
        elif previous is not None and not retry:
            expected=previous*2**.125; retry=True
        elif i != len(raw['circles'])-1:
            raise ValueError('continued after registered stopping condition')
    pairs = [(i, j) for i, j in zip(good, good[1:]) if
             (raw['circles'][i]['omega']-np.pi)*(raw['circles'][j]['omega']-np.pi) < 0]
    crossing = 'NO_CROSSING_ON_ACCEPTED_LADDER'; refinements=[]
    if pairs:
        crossing='CROSSING_UNRESOLVED'
        if len(raw['refinements']) == 2:
            i,j=pairs[0]
            for old, rec in zip((raw['circles'][i],raw['circles'][j]),raw['refinements']):
                m=circle_metrics(rec)
                delta=max(abs(rec['omega']-old['omega']),1e-9)
                m['uncertainty']=delta
                m['resolved_side']=bool(m['accepted'] and delta<=1e-7 and abs(rec['omega']-np.pi)>10*delta)
                if rec['a'] != old['a'] or len(rec['K']) != 95: raise ValueError('wrong refinement endpoint')
                refinements.append(m)
            if all(m['resolved_side'] for m in refinements) and np.prod([r['omega']-np.pi for r in raw['refinements']]) < 0:
                crossing='CROSSING_BRACKETED_NUMERICALLY'
    if not good: return dict(crossing=crossing, closure='NOT_ATTEMPTED', circles=metrics, accepted=0, termination=raw['termination'])
    nearest=min(good,key=lambda i:abs(raw['circles'][i]['omega']-np.pi))
    attempt=bool(pairs) or abs(raw['circles'][nearest]['omega']-np.pi)<.03
    K=np.asarray(raw['circles'][nearest]['K'])
    distance=float(np.min(np.linalg.norm(K-d.STAR,axis=1)))
    if len(raw['shooting']) != (8 if attempt else 0): raise ValueError('wrong shooting schedule')
    shooting=[]
    for i,rec in enumerate(raw['shooting']):
        theta=2*np.pi*i/8
        expected_nodes=d.spectral(K,np.array([theta,theta+np.pi]))
        if not np.allclose(rec['seed_nodes'],expected_nodes,rtol=0,atol=1e-12): raise ValueError('changed shooting seed')
        shooting.append(closure_metrics(rec,distance))
    closure=('NONTRIVIAL_TWO_RETURN_ORBIT_VERIFIED_NUMERICALLY' if any(m['verified'] for m in shooting)
             else 'CLOSURE_UNRESOLVED' if attempt else 'NOT_ATTEMPTED')
    return dict(crossing=crossing,closure=closure,termination=raw['termination'],accepted=len(good),
                circles=metrics,refinements=refinements,shooting=shooting,
                closest_index=nearest,closest_omega=raw['circles'][nearest]['omega'],
                action_selection='NOT_ESTABLISHED')


def measure(output):
    output.mkdir(parents=True,exist_ok=False)
    sources=source_hashes()
    (output/'provenance.json').write_text(json.dumps(dict(freeze=d.FREEZE,sources=sources),indent=2)+'\n')
    raw=dict(freeze=d.FREEZE,sources=sources,circles=[],refinements=[],shooting=[],termination='AMPLITUDE_CAP')
    K,omega=d.linear_seed(.004)
    a=.004; previous=None; retry=False; accepted=[]
    while a<=.2*(1+1e-12):
        print('Circle',a,flush=True)
        seed=K if previous is None else d.STAR+(K-d.STAR)*(a/previous)
        rec=d.circle(a,seed,omega)
        try: rec=d.validate_circle(rec)
        except (ValueError,ArithmeticError) as exc: rec['error']='validation: '+str(exc)
        metrics=circle_metrics(rec)
        write(output/f'circle_{len(raw["circles"]):02d}.json.gz.b64',rec)
        raw['circles'].append(rec)
        print(json.dumps(dict(a=a,**{k:metrics[k] for k in metrics if k not in ['full_checks']})),flush=True)
        if metrics['accepted']:
            accepted.append(rec); K=np.asarray(rec['K']);omega=rec['omega'];previous=a;retry=False
            if len(accepted)>1 and (accepted[-2]['omega']-np.pi)*(omega-np.pi)<0:
                raw['termination']='BRACKET_FOUND'
                for old in accepted[-2:]:
                    seed=d.spectral(np.asarray(old['K']),2*np.pi*np.arange(95)/95)
                    refined=d.circle(old['a'],seed,old['omega'])
                    try:refined=d.validate_circle(refined)
                    except (ValueError,ArithmeticError) as exc:refined['error']='validation: '+str(exc)
                    write(output/f'refine_{len(raw["refinements"])}.json.gz.b64',refined)
                    raw['refinements'].append(refined)
                break
            a*=FACTOR
        elif previous is not None and not retry:
            a=previous*2**.125;retry=True
        else:
            raw['termination']='CONTINUATION_FAILED';break
    if accepted:
        nearest=min(accepted,key=lambda r:abs(r['omega']-np.pi))
        if raw['termination']=='BRACKET_FOUND' or abs(nearest['omega']-np.pi)<.03:
            K=np.asarray(nearest['K']); distance=float(np.min(np.linalg.norm(K-d.STAR,axis=1)))
            for i in range(8):
                print('Shooting seed',i,flush=True)
                theta=2*np.pi*i/8
                rec=d.shoot(d.spectral(K,np.array([theta,theta+np.pi])))
                metrics=closure_metrics(rec,distance)
                if metrics.get('candidate'):
                    rec['validation']={}
                    try:
                        for method in ('DOP853','Radau'):
                            first=d.full_return(rec['nodes'][0],method)
                            second=d.full_return(rec['nodes'][1],method,first['returned'])
                            rec['validation'][method]=[first,second]
                    except (ValueError,ArithmeticError) as exc:rec['error']='closure validation: '+str(exc)
                write(output/f'shoot_{i}.json.gz.b64',rec);raw['shooting'].append(rec)
                print(json.dumps(closure_metrics(rec,distance)),flush=True)
    if sources!=source_hashes():raise RuntimeError('sources changed during measurement')
    report=analyze(raw)
    report.update(freeze=d.FREEZE,sources=sources,environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.glob('*.b64'))}
    manifest['provenance.json']=hashlib.sha256((output/'provenance.json').read_bytes()).hexdigest()
    (output/'report.json').write_text(json.dumps(report,default=serial,indent=2,allow_nan=False)+'\n')
    manifest['report.json']=hashlib.sha256((output/'report.json').read_bytes()).hexdigest()
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['crossing','closure','termination','accepted']},indent=2))
    return report


def load_run(directory):
    directory=Path(directory)
    report=json.loads((directory/'report.json').read_text())
    return dict(freeze=report['freeze'],sources=report['sources'],termination=report['termination'],
                circles=[decode(p) for p in sorted(directory.glob('circle_*.b64'))],
                refinements=[decode(p) for p in sorted(directory.glob('refine_*.b64'))],
                shooting=[decode(p) for p in sorted(directory.glob('shoot_*.b64'))])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=RUN)
    args=parser.parse_args()
    measure(args.output)
