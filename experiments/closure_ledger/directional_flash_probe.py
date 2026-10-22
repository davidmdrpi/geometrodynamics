"""Actual vector birth-flash projection and strict event-availability controls.

Run from repository root with OPENBLAS_NUM_THREADS=1:
python -m experiments.closure_ledger.directional_flash_probe
No trajectories are integrated and no force is assigned at empty foci.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.directional_flash_response import (
    coasting_focus_audit, electric, electric_integral, power_audit,
    rim_geometry, rim_modes,
)
from geometrodynamics.waves.pair_birth import PairBirth, Parameters

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'experiments/closure_ledger/runs/20261022_directional_flash'


def geometry(theta):
    return np.array([np.sin(theta), 0., np.cos(theta)]), np.array([np.cos(theta), 0., -np.sin(theta)])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    files = []

    def save(name, data):
        path = OUT/name
        path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
        files.append(path)

    projections, windows, convergence, powers = [], [], [], []
    for angle in [.16, .3, .5]:
        model = PairBirth(Parameters(mouth_angle=angle, cutoff=48))
        state = model.initial()
        qc, pc = model.flash(state)
        q, p = qc.sum(axis=0), pc.sum(axis=0)
        for n in [1, 2, 3, 4]:
            t = n*np.pi
            _, pt = model.free_wave(q, p, t)
            powers.append(dict(angle=angle, return_number=n,
                               **power_audit(model, pt, -state.x, -model.velocity(state))))
            for offset_ratio in [-.5, 0., .5]:
                center, first = geometry(angle+n*np.pi+offset_ratio*angle)
                points, er, ep, phi, frame = rim_geometry(center, angle, first=first)
                for source, qq, pp in [('A', qc[0], pc[0]), ('B', qc[1], pc[1]), ('coherent_pair', q, p)]:
                    _, momentum = model.free_wave(qq, pp, t)
                    center_field = (electric(model, momentum, center)@frame.T)[0]
                    ring_field = electric(model, momentum, points)
                    projections.append(dict(angle=angle, return_number=n, offset_ratio=offset_ratio,
                                            source=source, center_electric=center_field.tolist(),
                                            **rim_modes(ring_field, er, ep, phi)))
                    if offset_ratio == 0.:
                        for half_width in [.06, .12, .24, .48]:
                            integral = (electric_integral(model, qq, pp, center, t-half_width, t+half_width)@frame.T)[0]
                            windows.append(dict(angle=angle, return_number=n, source=source,
                                                half_width=half_width, integral_electric=integral.tolist()))
    save('directional_projections.json', projections)
    save('fixed_probe_windows.json', windows)
    save('reciprocal_power.json', powers)

    for cutoff in [16, 24, 32, 48, 64]:
        model = PairBirth(Parameters(cutoff=cutoff))
        state = model.initial(); q, p = model.flash(state)
        q, p = q.sum(axis=0), p.sum(axis=0)
        center, first = geometry(.3+np.pi)
        points, er, ep, phi, frame = rim_geometry(center, .3, first=first)
        pt = model.free_wave(q, p, np.pi)[1]
        convergence.append(dict(cutoff=cutoff, center_electric=(electric(model, pt, center)@frame.T)[0].tolist(),
                                **rim_modes(electric(model, pt, points), er, ep, phi)))
    save('cutoff_convergence.json', convergence)
    # Independent angular quadrature check at fixed spectral cutoff.
    angular = []
    for samples in [32, 64, 128, 256]:
        points, er, ep, phi, frame = rim_geometry(center, .3, samples=samples, first=first)
        angular.append(dict(samples=samples, **rim_modes(electric(model, pt, points), er, ep, phi)))
    save('rim_convergence.json', angular)

    # Use the SAME aperture eligibility function as PR #347; no mid-flight hits.
    routes = [coasting_focus_audit(a, v) for a in [.05, .16, .3, .5]
              for v in [.01, .03, .05, .1, .2, 1/3, .5, .65, .8, .95]]
    save('first_flight_events.json', routes)

    model = PairBirth(Parameters(cutoff=48))
    state = model.initial(); qc, pc = model.flash(state)
    center, frame_first = geometry(.3+np.pi)
    times = np.linspace(0, 4*np.pi, 801)
    traces = []
    for t in times:
        values = [float(electric(model, model.free_wave(qq, pp, t)[1], center)[0]@frame_first)
                  for qq, pp in zip(qc, pc)]
        traces.append(dict(time=float(t), source_A=values[0], source_B=values[1], total=sum(values)))
    save('fixed_antipode_waveform.json', traces)

    # At geometrically eligible first captures, show only the incident reference
    # field of the captured source. Its value is not a BAM event impulse.
    candidates = []
    for route in routes:
        model = PairBirth(Parameters(cutoff=48, mouth_angle=route['angle']))
        state = model.initial(v=route['speed']); qc, pc = model.flash(state)
        for event in route['events']:
            if event['target'] is None:
                continue
            target, source, t = event['target'], event['source'], event['time']
            sign = 1 if target == 0 else -1
            point, along = geometry(sign*(route['angle']+route['speed']*t))
            outgoing = sign*along
            momentum = model.free_wave(qc[source], pc[source], t)[1]
            force = model.charges[target]*electric(model, momentum, point)[0]
            candidates.append(dict(angle=route['angle'], speed=route['speed'], **event,
                                   reference_qE_along_velocity=float(force@outgoing),
                                   interpretation='Unglued Maxwell reference, not a boundary-matched event force or impulse'))
    save('eligible_reference_fields.json', candidates)

    ref = next(r for r in projections if r['angle']==.3 and r['return_number']==1
               and r['offset_ratio']==0 and r['source']=='coherent_pair')
    baseline = next(r for r in routes if r['angle']==.3 and r['speed']==1/3)
    summary = dict(
        interpretation='Vector Maxwell reference and antipodal eligibility audited separately; no BAM recoil trajectory',
        angular_radius=.3, birth_speed=1/3, footprint=.12, charge=.75, mass=1.,
        frequency_law='sqrt(l(l+1)); no extra imposed caustic phase',
        directional_cases=len(projections), window_cases=len(windows), first_flight_cases=len(routes),
        eligible_source_events=len(candidates),
        standard_first_flight_captures=sum(e['target'] is not None for e in baseline['events']),
        first_antipode_reference=ref,
        maximum_power_residual=max(abs(r['residual']) for r in powers),
        cutoff_48_64_center_difference=float(np.linalg.norm(np.array(convergence[-1]['center_electric'])-convergence[-2]['center_electric'])),
        cutoff_48_64_rim_difference=float(np.linalg.norm(np.array(convergence[-1]['transported_mean'])-convergence[-2]['transported_mean'])),
        rim_128_256_difference=float(np.linalg.norm(np.array(angular[-1]['transported_mean'])-angular[-2]['transported_mean'])),
        preserved_negative_result='No demonstrated sustained approach, collapse or bound trajectory; scalar candidate remains negative.',
    )
    save('summary.json', summary)
    figure(projections, windows, routes, traces)
    files.append(ROOT/'docs/images/history_throat/directional_flash.png')
    files += [Path(__file__), ROOT/'geometrodynamics/waves/directional_flash_response.py',
              ROOT/'geometrodynamics/waves/pair_birth.py', ROOT/'geometrodynamics/waves/antipodal_throat_transport.py',
              ROOT/'tests/test_directional_flash_response.py', ROOT/'docs/directional_flash_response.md']
    save('manifest.json', dict(files=[dict(path=str(p.relative_to(ROOT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                                     for p in files]))
    print(json.dumps(summary, indent=2))


def figure(projections, windows, routes, traces):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax = axes[0, 0]
    for key, label in [('source_A', 'Source A'), ('source_B', 'Source B'), ('total', 'Coherent pair')]:
        ax.plot([r['time']/np.pi for r in traces], [r[key] for r in traces], label=label)
    ax.axhline(0, color='.7', lw=.6); ax.set(xlabel='Time / π', ylabel='Tangential E', title='(a) Fixed probe at A’s antipode; a = 0.3')
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    for n in [1, 2, 3, 4]:
        rows = [r for r in windows if r['angle']==.3 and r['return_number']==n and r['source']=='coherent_pair']
        ax.plot([r['half_width'] for r in rows], [r['integral_electric'][0] for r in rows], 'o-', label=f'n = {n}')
    ax.axhline(0, color='.7', lw=.6); ax.set(xlabel='Diagnostic window half-width', ylabel='∫ E dt at each alternating focus', title='(b) Window integrals are not event kicks')
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    for source in ['A', 'B', 'coherent_pair']:
        rows = [r for r in projections if r['angle']==.3 and r['return_number']==1 and r['source']==source]
        ax.plot([r['offset_ratio'] for r in rows], [r['transported_mean'][0] for r in rows], 'o-', label=source.replace('_', ' '))
    ax.axhline(0, color='.7', lw=.6); ax.set(xlabel='Probe offset / aperture radius', ylabel='Transported rim mean E', title='(c) Directional content, not neck-normal stress')
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    for capture, marker, color, label in [(False, 'x', '.6', 'No capture'), (True, 'o', '#168457', 'Capture')]:
        rows = [r for r in routes if any(e['target'] is not None for e in r['events']) == capture]
        ax.scatter([r['speed'] for r in rows], [r['angle'] for r in rows], marker=marker, color=color, label=label)
    ax.set(xlabel='Prescribed coasting speed', ylabel='Angular aperture radius a', yticks=[.05, .16, .3, .5], title='(d) Strict antipodal eligibility before first contact')
    ax.legend(fontsize=8, loc='center', bbox_to_anchor=(.77, .73))
    fig.suptitle('Directional birth-flash diagnostics — reference fields, no recoil evolution', fontsize=14)
    path = ROOT/'docs/images/history_throat/directional_flash.png'
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    main()
