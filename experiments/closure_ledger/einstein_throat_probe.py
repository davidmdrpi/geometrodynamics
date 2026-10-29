"""Reproduce conditional Einstein evolution and finite-port group-delay audits."""
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves.maxwell_throat_response import Bridge
from geometrodynamics.waves.s3_maxwell_encounter import DipolePulse
from geometrodynamics.waves.einstein_throat_evolution import (
    initial_acceleration, central_history, integrate_central, optical_length, group_delays,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'experiments/closure_ledger/runs/20261028_einstein_throat'


def save(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/(name+'.json')).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    histories, checks, supports = [], [], []
    for angle in [.2, .3, .45]:
        for fraction in [0., .5, .9, .99]:
            b = Bridge(angle=angle, charge_fraction=fraction)
            exact = central_history(b, np.linspace(0, .97*np.pi, 101))
            coarse = integrate_central(b, exact['tau'], rtol=2e-10)
            fine = integrate_central(b, exact['tau'], rtol=3e-12)
            r, u = exact['radius'], exact['velocity']
            histories.append(dict(angle=angle, charge_fraction=fraction, mass=b.mass,
                geometric_charge=b.geometric_charge, initial_neck=b.b,
                neck_initial_acceleration=float(initial_acceleration(b, 0)),
                rim_normal_geodesic_initial_radius_acceleration=float(initial_acceleration(b, b.limit)),
                inner_radius=b.geometric_charge**2/b.b, endpoint_proper_time=float(np.pi*b.mass),
                endpoint='curvature singularity' if fraction == 0 else 'ideal RN inner Cauchy horizon',
                history={key: value.tolist() for key, value in exact.items()}))
            checks.append(dict(angle=angle, charge_fraction=fraction,
                max_radius_error_over_initial_radius=float(np.max(abs(fine['radius']-r))/b.b),
                max_velocity_error=float(np.max(abs(fine['velocity']-u))),
                max_refinement_radius_change_over_initial_radius=float(np.max(abs(fine['radius']-coarse['radius']))/b.b),
                max_exact_energy_residual=float(np.max(abs(u*u+1-2*b.mass/r+b.geometric_charge**2/r**2)))))
            p_hold = (b.support(0)[1]+b.geometric_charge**2/b.b**4)/(8*np.pi)
            supports.append(dict(angle=angle, charge_fraction=fraction, support_density=0.,
                holding_radial_pressure=float(p_hold),
                pressure_multipliers=[0., .5, 1., 1.5, 2.],
                initial_neck_accelerations=[float(initial_acceleration(b, 0, f*p_hold)) for f in [0., .5, 1., 1.5, 2.]],
                scope='Initial acceleration sensitivity at identical Hamiltonian/momentum data; not five evolved matter models'))
    save('histories', histories)
    save('evolution_checks', checks)
    save('pressure_sensitivity', supports)
    delays, delay_checks = [], []
    pulse = DipolePulse(48, .2)
    energy = .5*((pulse.n*pulse.q0)**2+pulse.p0**2)
    for depth in [0., .4, .8, .9, .95, .99]:
        b = Bridge(lapse_depth=depth)
        values = []
        for n in pulse.n:
            value = group_delays(b, float(n), step=1e-4)
            value['smith_eigenvalues'] = np.linalg.eigvalsh(value.pop('smith')).tolist()
            values.append(value)
        t = sum(e*v['transmission'] for e, v in zip(energy, values))
        rd = sum(e*v['reflection']*v['reflection_delay'] for e, v in zip(energy, values))
        td = sum(e*v['transmission']*v['transmission_delay'] for e, v in zip(energy, values))
        delays.append(dict(lapse_depth=depth, neck_lapse=1-depth, optical_length=optical_length(b),
            transmitted_energy_weight=float(t), reflected_energy_weight=float(sum(e*v['reflection'] for e, v in zip(energy, values))),
            transmitted_weighted_group_delay=float(td/t), reflected_weighted_group_delay=float(rd/(1-t)),
            incident_weighted_smith_delay=float(rd+td), modes=values,
            spectrum_energy_weights=energy.tolist(), scope='Modal group-delay diagnostics, not a time-domain pulse arrival or router gate'))
        print('completed delay depth', depth, flush=True)
        for omega in [5., 10., 20.]:
            coarse = group_delays(b, omega, step=2e-4)
            fine = group_delays(b, omega, step=5e-5, rtol=3e-12)
            delay_checks.append(dict(lapse_depth=depth, omega=omega,
                mean_delay_change=abs(fine['mean_delay']-coarse['mean_delay']),
                transmission_delay_change=abs(fine['transmission_delay']-coarse['transmission_delay']),
                hermiticity_error=fine['hermiticity_error']))
    save('delays', delays)
    save('delay_checks', delay_checks)
    summary = dict(evolution_cases=len(histories), pressure_cases=len(supports),
        delay_modes=sum(len(x['modes']) for x in delays), delay_refinement_cases=len(delay_checks),
        max_radius_error_over_initial_radius=max(x['max_radius_error_over_initial_radius'] for x in checks),
        max_velocity_error=max(x['max_velocity_error'] for x in checks),
        max_delay_refinement_change=max(x['mean_delay_change'] for x in delay_checks),
        conditional_electrovacuum_radius_evolution_solved=True,
        bam_mouth_center_acceleration_solved=False, nonsingular_annihilation_solved=False,
        router_modified=False, evolved_support_constitutive_law=False,
        charged_history_continued_beyond_cauchy_horizon=False)
    save('summary', summary)
    figure(histories, delays)
    print(json.dumps(summary, indent=2), flush=True)


def figure(histories, delays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    for case in histories:
        if case['angle'] != .3:
            continue
        h = case['history']
        axes[0].plot(np.array(h['tau'])/case['endpoint_proper_time'],
                     np.array(h['radius'])/case['initial_neck'], label=f"q/qmax={case['charge_fraction']}")
        axes[0].plot(1, case['inner_radius']/case['initial_neck'], 'o', color=axes[0].lines[-1].get_color(), ms=4)
    axes[0].set(xlabel='Central proper time / (pi M)', ylabel='Central areal radius / initial radius',
                title='Unsupported electrovacuum control (a=0.3)')
    axes[0].legend(fontsize=8)
    axes[0].text(.03, .04, 'Neutral: singular endpoint\nCharged: stop at inner horizon; no stable bounce claimed',
                 transform=axes[0].transAxes, fontsize=8)
    for key, label in [('optical_length', 'Optical length'),
                       ('transmitted_weighted_group_delay', 'Transmitted group delay'),
                       ('reflected_weighted_group_delay', 'Reflected group delay')]:
        axes[1].plot([x['neck_lapse'] for x in delays], [x[key] for x in delays], 'o-', label=label)
    axes[1].set(xscale='log', xlabel='Supported static neck lapse', ylabel='Time (R=c=1)',
                title='Finite-port modal delay diagnostics')
    axes[1].invert_xaxis()
    axes[1].legend(fontsize=8)
    fig.suptitle('Different Lorentzian completions: collapse control and supported scattering')
    p = ROOT/'docs/images/history_throat/einstein_throat_evolution.png'
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=160)
    plt.close(fig)


def manifest():
    files = sorted(p for p in OUT.glob('*.json') if p.name != 'manifest.json')
    files += [Path(__file__), ROOT/'geometrodynamics/waves/einstein_throat_evolution.py',
        ROOT/'geometrodynamics/waves/maxwell_throat_response.py', ROOT/'geometrodynamics/waves/gr_handle_checks.py',
        ROOT/'geometrodynamics/waves/s3_maxwell_encounter.py', ROOT/'tests/test_einstein_throat_evolution.py',
        ROOT/'docs/einstein_throat_evolution.md', ROOT/'docs/images/history_throat/einstein_throat_evolution.png']
    save('manifest', dict(files=[dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]))


if __name__ == '__main__':
    main()
