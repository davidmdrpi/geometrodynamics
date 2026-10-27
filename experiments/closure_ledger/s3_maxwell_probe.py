"""Reproducible S3 Maxwell focusing and finite receiver worldtube controls."""
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.s3_maxwell_encounter import (
    DipolePulse, Worldtube, encounter, surface_rates, stored,
    ConductingReceiver,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'experiments/closure_ledger/runs/20261027_s3_maxwell'


def save(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/(name+'.json')).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    phase_review()
    rows = []
    pulse = DipolePulse(64, .2)
    for radius in [.15, .3, .45]:
        for offset in [0., .12, .36]:
            for deformation in [0., .2]:
                tube = Worldtube(radius, offset, deformation, 24)
                row = dict(radius=radius, offset=offset, deformation=deformation,
                           cutoff=64, width=.2, angular_order=24, time_order=128)
                row['whole'] = encounter(pulse, tube, np.pi-1, np.pi+1, 128, 24)
                row['partial'] = encounter(pulse, tube, np.pi-.8, np.pi-.12, 96, 24)
                row['focus_stored_energy'] = float(stored(pulse, tube, [np.pi], 24)[0, 0])
                rows.append(row)
                print(json.dumps({k: row[k] for k in ['radius', 'offset', 'deformation', 'focus_stored_energy']}), flush=True)
    save('encounters', rows)

    refinement = []
    for width in [.12, .2]:
        for cutoff in [24, 48, 72]:
            for order in [16, 24, 36]:
                p = DipolePulse(cutoff, width)
                tube = Worldtube(.3, .12, .2, order)
                result = encounter(p, tube, np.pi-.8, np.pi-.12, 4*order, order)
                e, b, _ = p.fields(np.pi, p.prepare(np.array([[-1., 0, 0, 0]])))
                refinement.append(dict(width=width, cutoff=cutoff, angular_order=order,
                    time_order=4*order, radial_order=order, result=result,
                    focal_density=float(.5*np.sum(e*e+b*b))))
    save('refinement', refinement)

    # Full field antipodal transport, not a scalar-template phase fit.
    rng = np.random.default_rng(349)
    points = rng.normal(size=(257, 4)); points /= np.linalg.norm(points, axis=1)[:, None]
    transport = []
    for phase in [0., np.pi/4, np.pi/2]:
        p = DipolePulse(64, .2, phase)
        for t in [.17, .74, 1.31]:
            e, b, _ = p.fields(t+np.pi, p.prepare(points))
            ei, bi, _ = p.fields(t, p.prepare(-points))
            q, momentum = p.state(np.pi)
            parity = (-1.)**p.n
            defects = {}
            for eta in [-1, 1]:
                defects[str(eta)] = float(.5*np.sum((p.n*(q-eta*parity*p.q0))**2
                                                    +(momentum+eta*parity*p.p0)**2)/p.energy)
            transport.append(dict(phase=phase, time=t, max_field_refocus_error=float(max(
                np.max(abs(e+ei)), np.max(abs(b+bi)))), quotient_raw_defect_energy_fraction=defects))
    save('transport', transport)

    # Preserve a nonzero transient field-momentum curve: omitting it creates false recoil.
    tube = Worldtube(.3, .12, .2, 24)
    times = np.linspace(np.pi-.8, np.pi+.8, 161)
    cache = pulse.prepare(tube.points)
    rates = [surface_rates(pulse, tube, t, cache) for t in times]
    charges = stored(pulse, tube, times, 24)
    trace = dict(time=times.tolist(), stored_energy=charges[:, 0].tolist(),
                 stored_killing_momentum=charges[:, 1].tolist(),
                 rates={key: [r[key] for r in rates] for key in rates[0]})
    save('trace', trace)

    # Independent reflection of the receiver about x1=0 changes signed momentum.
    symmetry = []
    for offset in [-.12, .12]:
        result = encounter(pulse, Worldtube(.3, offset, .2, 24), np.pi-.8, np.pi-.12, 96, 24)
        symmetry.append(dict(offset=offset, result=result))
    save('offset_symmetry', symmetry)
    summary = dict(
        geometry='unit round S3 x time, prescribed metric, source-free TE dipole Maxwell',
        prepared_energy=1., modal_frequency='n=k+1, n>=2',
        total_encounter_cases=len(rows), refinement_cases=len(refinement),
        max_refocus_field_error=max(r['max_field_refocus_error'] for r in transport),
        max_whole_energy_residual=max(abs(r['whole']['energy_balance_residual']) for r in rows),
        max_whole_transparent_impulse=max(abs(r['whole']['transparent_matter_impulse']) for r in rows),
        max_partial_transparent_impulse=max(abs(r['partial']['transparent_matter_impulse']) for r in rows),
        max_partial_incident_only_impulse=max(abs(r['partial']['incident_only_impulse']) for r in rows),
        explicit_nonclaims=['No sourced birth law', 'No throat scattering or finite-mouth absorption',
                           'No curvature backreaction', 'No attraction or collapse established'],
        strict_router_modified=False, s2_dispersion_obstruction_applies_to_s3=False)
    reflected = reflection_controls()
    summary['reflection_cases'] = len(reflected)
    summary['reflection_scope'] = 'Ideal fixed conductor, inward radius load; zero vector impulse by symmetry'
    save('summary', summary)
    figure(rows, refinement, trace)
    print(json.dumps(summary, indent=2), flush=True)


def phase_review():
    """Reproduce the latest review's uniform-phase correction on its S2 operator."""
    from geometrodynamics.waves.pair_birth import PairBirth, Parameters
    from geometrodynamics.waves.quotient_flash import project_history
    model = PairBirth(Parameters(cutoff=48, mouth_angle=.3))
    rows = []
    for k in range(4):
        time = (2*k+1)*np.pi
        q, p = model.flash(model.initial(v=(np.pi-.6)/time))
        q = q.sum(0)
        for name, omega in [('s2_maxwell', model.omega), ('s2_conformal', model.ell+.5)]:
            weights = abs(omega*q)**2
            d = np.sum(weights*(-1.)**(model.ell+1)*np.exp(1j*omega*time))/weights.sum()
            theta = -np.angle(d)/2
            result = project_history(omega, (-1.)**(model.ell+1), q*np.cos(theta),
                                     -omega*q*np.sin(theta), time)
            rows.append(dict(k=k, operator=name, theta_degrees=float(theta*180/np.pi),
                minimum_correction_fraction=float(result['correction_energy']/result['initial_energy']),
                analytic_minimum=float((1-abs(d))/2),
                scope='Optimized prepared phase; not a derived birth law; S2, not S3'))
    save('s2_phase_review', rows)


def upstream_endpoint_review():
    # Reuse the byte-pinned reproducer in a fresh process; do not alter its old artifacts.
    from experiments.closure_ledger import quotient_endpoint_review as audit
    audit.PIN = 'a464530e5fc9a1b257582ee45e0bea3d92bcee6c'
    audit.OUT = OUT
    audit.main()


def reflection_controls():
    rows = []
    for radius in [.15, .3, .45]:
        for cells in [256, 512, 1024, 2048]:
            receiver = ConductingReceiver(radius=radius, cells=cells)
            row = receiver.encounter(np.pi-1, np.pi+1, 256)
            if cells == 2048:
                refined = receiver.encounter(np.pi-1, np.pi+1, 384)
                row['quadrature_change'] = refined['radius_impulse']-row['radius_impulse']
            rows.append(row)
    save('reflection', rows)
    return rows


def figure(rows, refinement, trace):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    t = np.array(trace['time'])-np.pi
    axes[0, 0].plot(t, trace['stored_energy'])
    axes[0, 0].set(xlabel='Time relative to antipodal focus', ylabel='Stored field energy / initial energy',
                   title='Finite receiver: a=0.3, offset=0.12, deformation=0.2')
    axes[0, 1].plot(t, trace['stored_killing_momentum'], label='Stored Killing momentum')
    axes[0, 1].plot(t, trace['rates']['momentum_flux'], label='Outward momentum flux')
    axes[0, 1].set(xlabel='Time relative to antipodal focus', title='Transient field momentum is not mouth recoil')
    axes[0, 1].legend(fontsize=8)
    for a in [.15, .3, .45]:
        sub = [r for r in rows if r['radius'] == a and r['deformation'] == 0]
        axes[1, 0].plot([r['offset'] for r in sub], [r['focus_stored_energy'] for r in sub], 'o-', label=f'a={a}')
    axes[1, 0].set(xlabel='Receiver offset (radians)', ylabel='Energy fraction at focus',
                   title='Finite size and alignment matter')
    axes[1, 0].legend(fontsize=8)
    for width in [.12, .2]:
        sub = [r for r in refinement if r['width'] == width and r['angular_order'] == 36]
        axes[1, 1].plot([r['cutoff'] for r in sub], [r['focal_density'] for r in sub], 'o-', label=f'width={width}')
    axes[1, 1].set(xlabel='Highest Maxwell frequency n', ylabel='Energy density at focus',
                   title='Fixed smooth pulse: spectral convergence')
    axes[1, 1].legend(fontsize=8)
    fig.suptitle('S³ Maxwell encounter reference — transparent observation surfaces, no throat force', fontsize=13)
    target = ROOT/'docs/images/history_throat/s3_maxwell_encounter.png'
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=160)
    plt.close(fig)


def manifest():
    paths = sorted(p for p in OUT.glob('*.json') if p.name != 'manifest.json')
    paths += [Path(__file__), ROOT/'geometrodynamics/waves/s3_maxwell_encounter.py',
              ROOT/'tests/test_s3_maxwell_encounter.py', ROOT/'docs/s3_maxwell_encounter.md',
              ROOT/'docs/images/history_throat/s3_maxwell_encounter.png',
              ROOT/'experiments/closure_ledger/quotient_endpoint_review.py',
              ROOT/'geometrodynamics/waves/pair_birth.py', ROOT/'geometrodynamics/waves/quotient_flash.py']
    save('manifest', dict(files=[dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]))


if __name__ == '__main__':
    import sys
    if '--upstream-review' in sys.argv:
        upstream_endpoint_review()
    else:
        main()
