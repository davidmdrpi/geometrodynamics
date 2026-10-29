"""Geometry-derived Maxwell two-ports, fixed-flux tension, and shape loads."""
import hashlib
import json
from pathlib import Path
import numpy as np

from geometrodynamics.waves.maxwell_throat_response import (
    Bridge, solve_port, force_matrices, scattering_force_derivative, static_flux,
)
from geometrodynamics.waves.s3_maxwell_encounter import DipolePulse

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'experiments/closure_ledger/runs/20261028_maxwell_throat'


def save(name, value):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/(name+'.json')).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def complex_matrix(matrix):
    return np.stack([matrix.real, matrix.imag], axis=-1).tolist()


def serialize(bridge, value):
    return dict(angle=bridge.angle, lapse_depth=bridge.lapse_depth,
        charge_fraction=bridge.charge_fraction, neck_radius=bridge.b,
        length=2*bridge.half_length, omega=value['omega'], ell=value['ell'],
        transmission=value['transmission'], reflection=value['reflection'],
        unitarity_error=value['unitarity_error'], scattering=complex_matrix(value['scattering']),
        waist_force=complex_matrix(value['waist_force']), stretch_force=complex_matrix(value['stretch_force']),
        stored_energy=complex_matrix(value['stored_energy']),
        waist_eigenvalues=value['waist_eigenvalues'].tolist(),
        stretch_eigenvalues=value['stretch_eigenvalues'].tolist())


def main():
    rows = []
    for angle in [.2, .3, .45]:
        for depth in [0., .4, .8]:
            bridge = Bridge(angle, depth)
            for ell in [1, 2]:
                for omega in [5., 10., 20., 40., 60., 100.]:
                    rows.append(serialize(bridge, force_matrices(bridge, omega, ell)))
            print('completed neutral', angle, depth, flush=True)
    for fraction in [.5, .9, .99]:
        for depth in [0., .8]:
            bridge = Bridge(lapse_depth=depth, charge_fraction=fraction)
            for omega in [5., 10., 20., 40., 60.]:
                rows.append(serialize(bridge, force_matrices(bridge, omega)))
    save('scattering', rows)

    weighted = []
    pulse = DipolePulse(48, .2)
    energies = .5*((pulse.n*pulse.q0)**2+pulse.p0**2)
    for fraction in [0., .5, .9]:
        for depth in [0., .4, .8]:
            bridge = Bridge(lapse_depth=depth, charge_fraction=fraction)
            values = [force_matrices(bridge, float(n), order=192) for n in pulse.n]
            transmitted = sum(e*v['transmission'] for e, v in zip(energies, values))
            waist = sum(e*v['waist_force'] for e, v in zip(energies, values))
            stretch = sum(e*v['stretch_force'] for e, v in zip(energies, values))
            static = static_flux(bridge, np.sqrt(4*np.pi)*bridge.geometric_charge)
            weighted.append(dict(charge_fraction=fraction, lapse_depth=depth,
                spectrum_mode_numbers=pulse.n.tolist(), spectrum_energy_weights=energies.tolist(),
                transmission=float(transmitted), waist_force_per_incident_power=complex_matrix(waist),
                stretch_force_per_incident_power=complex_matrix(stretch),
                static_matched_flux=static,
                zero_stretch_load_incident_power=float(-static['stretch_force']/stretch[0, 0].real),
                scope='Prescribed phase-averaged port drive using S3 modal energy weights; not a sourced closed-universe encounter'))
            print('completed weighted', fraction, depth, transmitted, flush=True)
    save('weighted_drive', weighted)

    verification = []
    for fraction, depth, omega in [(0., 0., 10.), (0., 0., 60.), (0., .8, 10.), (.9, .4, 20.)]:
        b = Bridge(lapse_depth=depth, charge_fraction=fraction)
        fine = force_matrices(b, omega, order=384, rtol=2e-12)
        coarse = force_matrices(b, omega, order=192, rtol=2e-10)
        for parameter in ['waist', 'stretch']:
            for step in [2e-5, 5e-6]:
                derivative = scattering_force_derivative(b, omega, parameter=parameter, step=step)
                key = parameter+'_force'
                verification.append(dict(charge_fraction=fraction, lapse_depth=depth, omega=omega,
                    parameter=parameter, step=step,
                    phase_derivative_error=float(np.max(abs(fine[key]-derivative))),
                    quadrature_ode_change=float(np.max(abs(fine[key]-coarse[key]))),
                    fine_unitarity_error=fine['unitarity_error']))
    save('verification', verification)

    supports = []
    for fraction in [0., .5, .9, .99]:
        for depth in [0., .4, .8]:
            b = Bridge(lapse_depth=depth, charge_fraction=fraction)
            x = np.linspace(-b.limit, b.limit, 65)
            p = b.profile(x)
            stress = b.support(x)
            em = b.geometric_charge**2/p['radius']**4
            support = stress-np.array([em, -em, em])
            supports.append(dict(charge_fraction=fraction, lapse_depth=depth,
                geometric_charge=b.geometric_charge, mass=b.mass, neck_radius=b.b,
                length=2*b.half_length, s=p['s'].tolist(), radius=p['radius'].tolist(),
                lapse=p['lapse'].tolist(), total_stress_8pi=stress.tolist(),
                support_stress_8pi=support.tolist(),
                neck_nec_8pi=float(stress[0, 32]+stress[1, 32]),
                max_support_density_8pi=float(np.max(abs(support[0])))))
    save('support', supports)
    summary = dict(scattering_cases=len(rows), weighted_drive_cases=len(weighted),
        verification_cases=len(verification), support_cases=len(supports),
        max_unitarity_error=max(r['unitarity_error'] for r in rows),
        minimum_waist_force_eigenvalue=min(min(r['waist_eigenvalues']) for r in rows),
        minimum_stretch_force_eigenvalue=min(min(r['stretch_eigenvalues']) for r in rows),
        max_phase_force_error=max(r['phase_derivative_error'] for r in verification),
        max_quadrature_ode_change=max(r['quadrature_ode_change'] for r in verification),
        actual_support_constitutive_law_derived=False,
        mouth_translation_or_collapse_derived=False,
        findings=['Maxwell geometry determines a finite-port scattering map for each stated positive lapse',
                  'Fixed-lapse nonnegative waist variations have nonnegative radiative and static-flux loads',
                  'Fixed topological flux supplies intrinsic-length tension',
                  'Charged Hamiltonian constraint is solved before including finite static flux',
                  'Positive-lapse static completions require additional null-energy-violating stress',
                  'Time-reversing seams require real-linear matching, not a complex gain'])
    save('summary', summary)
    figure(rows, weighted)
    print(json.dumps(summary, indent=2), flush=True)


def figure(rows, weighted):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    for depth in [0., .4, .8]:
        sub = [r for r in rows if r['angle'] == .3 and r['ell'] == 1
               and r['charge_fraction'] == 0 and r['lapse_depth'] == depth]
        axes[0].semilogy([r['omega'] for r in sub], [r['transmission'] for r in sub], 'o-', label=f'depth={depth}')
    axes[0].set(xlabel='Frequency (R=c=1)', ylabel='Finite-port transmitted power fraction',
                title='Same spatial bridge, different lapse')
    axes[0].legend(fontsize=8)
    sub = [r for r in rows if r['angle'] == .3 and r['ell'] == 1 and r['charge_fraction'] == 0 and r['lapse_depth'] == 0]
    w = [r['omega'] for r in sub]
    axes[1].fill_between(w, [r['waist_eigenvalues'][0] for r in sub],
                        [r['waist_eigenvalues'][1] for r in sub], alpha=.3, label='All coherent two-port inputs')
    axes[1].plot(w, [r['waist_force'][0][0][0] for r in sub], 'o-', label='Left incidence')
    axes[1].axhline(0, color='black', lw=.7)
    axes[1].set(xlabel='Frequency', ylabel='Waist load / incident power', title='Positive load opens the waist')
    axes[1].legend(fontsize=8)
    sub = [r for r in weighted if r['lapse_depth'] == 0 and r['charge_fraction'] > 0]
    q = [r['charge_fraction'] for r in sub]
    axes[2].plot(q, [r['static_matched_flux']['waist_force'] for r in sub], 'o-', label='Waist opening load')
    axes[2].plot(q, [r['static_matched_flux']['stretch_force'] for r in sub], 'o-', label='Intrinsic-length load')
    axes[2].axhline(0, color='black', lw=.7)
    axes[2].set(xlabel='Charge / matched extremal limit', ylabel='Generalized load (G=c=R=1)',
                title='Static flux: opens radius, shortens length')
    axes[2].legend(fontsize=8)
    fig.suptitle('Maxwell throat action — scattering and loads; support dynamics remains unspecified')
    p = ROOT/'docs/images/history_throat/maxwell_throat_response.png'
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=160)
    plt.close(fig)


def manifest():
    files = sorted(p for p in OUT.glob('*.json') if p.name != 'manifest.json')
    files += [Path(__file__), ROOT/'geometrodynamics/waves/maxwell_throat_response.py',
        ROOT/'geometrodynamics/waves/physical_throat.py', ROOT/'geometrodynamics/waves/gr_handle_checks.py',
        ROOT/'geometrodynamics/waves/s3_maxwell_encounter.py', ROOT/'tests/test_maxwell_throat_response.py',
        ROOT/'docs/maxwell_throat_response.md', ROOT/'docs/images/history_throat/maxwell_throat_response.png']
    save('manifest', dict(files=[dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]))


if __name__ == '__main__':
    main()
