#!/usr/bin/env python
"""Reproduce the frozen antipodal-kick comparison and convergence checks."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.waves.antipodal_kick import AntipodalKick, Parameters, harmonics


def save_trace(path, a):
    scalar_keys = ['t', 'q', 'pq', 'energy', 'shape_energy', 'ubar']
    vector_keys = ['x', 'angular', 'wave_angular', 'torque']
    columns = [a[k] for k in scalar_keys]
    names = list(scalar_keys)
    for k in vector_keys:
        columns.extend(a[k].T)
        names.extend(f'{k}_{i}' for i in 'xyz')
    np.savetxt(path, np.column_stack(columns), delimiter=',',
               header=','.join(names), comments='', fmt='%.12e')


def bare_cap_energy(model, state, time, cap=.6):
    """Quadrature over the target cap, independently of the attachment footprint.

    Gauss-Legendre in cos(theta_about_x), uniform azimuth. Rotation derivatives
    give the full intrinsic gradient norm, including both angular directions.
    """
    s = state.copy()
    model.drift(s, time)
    z, weights = np.polynomial.legendre.leggauss(48)
    z = (1+np.cos(cap))/2 + z*(1-np.cos(cap))/2
    weights *= (1-np.cos(cap))/2
    phi = np.arange(96)*2*np.pi/96
    z, phi = np.meshgrid(z, phi, indexing='ij')
    points = np.column_stack([z.ravel(), (np.sqrt(1-z*z)*np.cos(phi)).ravel(),
                              (np.sqrt(1-z*z)*np.sin(phi)).ravel()])
    coeff = np.column_stack([s.a, s.p]+[d@s.a for d in model.generators])
    values = []
    for chunk in np.array_split(points, 12):
        values.append((harmonics(model.ell, model.order, chunk)@coeff).real)
    values = np.concatenate(values)
    density = .5*(values[:, 1]**2+.25*values[:, 0]**2
                   +np.sum(values[:, 2:]**2, axis=1))
    return float(np.sum(density.reshape(48, 96)*weights[:, None])*2*np.pi/96)


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    par = Parameters()
    summary = {'parameters': asdict(par), 'dt': .002, 'horizon': np.pi+1.2,
               'focus_window': [np.pi-.6, np.pi+.6], 'cases': {}, 'convergence': {},
               'runtime': {'python': platform.python_version(), 'numpy': np.__version__,
                           'scipy': scipy.__version__}}
    protocol = Path(__file__).resolve().parents[2]/'docs/antipodal_kick_protocol.md'
    summary['protocol_sha256'] = hashlib.sha256(protocol.read_bytes()).hexdigest()
    render = {}
    for name, eta, dephased, coupling, energy in [
        ('symmetric', 0., False, 1., .1), ('directional', .7, False, 1., .1),
        ('mirror', -.7, False, 1., .1), ('dephased', .7, True, 1., .1),
        ('uncoupled', .7, False, 0., .1), ('half_amplitude', .7, False, 1., .025),
    ]:
        model = AntipodalKick(replace(par, coupling=coupling, energy=energy))
        s = model.initial(eta, dephased)
        arrays, stats = model.evolve(s)
        stats['bare_target_cap_energy_at_focus'] = bare_cap_energy(model, s, np.pi)
        stats['bare_target_cap_energy_at_midtransit'] = bare_cap_energy(model, s, np.pi/2)
        # Antipodal parity maps a source-cap integral to a target-cap integral.
        source = s.copy()
        source.a *= (-1.)**model.ell
        source.p *= (-1.)**model.ell
        stats['initial_source_cap_energy'] = bare_cap_energy(model, source, 0.)
        summary['cases'][name] = stats
        save_trace(output/f'{name}.csv', arrays)
        if name in ['symmetric', 'directional', 'dephased']:
            # Output coefficients only at a sparse movie cadence. Full-precision
            # traces and deterministic replay cover the finer diagnostic cadence.
            idx = np.unique(np.r_[np.arange(0, len(arrays['t']), 4),
                                   np.argmin(abs(arrays['t']-np.pi)), len(arrays['t'])-1])
            for key in ['t', 'a', 'x', 'q']:
                render[f'{name}_{key}'] = arrays[key][idx]
        print(name, json.dumps(stats), flush=True)
    np.savez_compressed(output/'render_states.npz', **render)
    for label, cutoff, dt in [('L24', 24, .002), ('L40', 40, .002), ('dt_half', 32, .001)]:
        model = AntipodalKick(replace(par, cutoff=cutoff))
        arrays, stats = model.evolve(dt=dt)
        summary['convergence'][label] = stats
        save_trace(output/f'convergence_{label}.csv', arrays)
        print(label, json.dumps(stats), flush=True)
    cases = summary['cases']
    baseline = np.array(cases['directional']['final_impulse'])
    changes = {key: float(np.linalg.norm(np.array(v['final_impulse'])-baseline)
                          /max(np.linalg.norm(baseline), 1e-7))
               for key, v in summary['convergence'].items()}
    all_stats = list(cases.values())+list(summary['convergence'].values())
    gates = {
        'energy': all(v['max_energy_relative'] < 1e-4 for v in all_stats),
        'angular_momentum': all(v['max_angular_relative'] < 1e-8 for v in all_stats),
        'space_convergence': np.linalg.norm(np.array(summary['convergence']['L40']['final_impulse'])-baseline)
                < max(.02*np.linalg.norm(baseline), 1e-7),
        'time_convergence': np.linalg.norm(np.array(summary['convergence']['dt_half']['final_impulse'])-baseline)
                < max(.02*np.linalg.norm(baseline), 1e-7),
        'symmetric_null': cases['symmetric']['max_attachment_angular'] < 1e-8,
        'mirror_reversal': np.linalg.norm(np.array(cases['mirror']['final_impulse'])+baseline) < 1e-7,
        'uncoupled_null': cases['uncoupled']['max_attachment_angular'] < 1e-14,
        'localized_activity': cases['directional']['focus_activity_fraction'] >= .75,
        'coherence_enhancement': np.linalg.norm(cases['directional']['focus_impulse'])
                >= 2*np.linalg.norm(cases['dephased']['focus_impulse']),
        'nonzero_focus_impulse': np.linalg.norm(cases['directional']['focus_impulse']) > 1e-5,
    }
    gates = {k: bool(v) for k, v in gates.items()}
    numerical_keys = ['energy', 'angular_momentum', 'space_convergence', 'time_convergence',
                      'symmetric_null', 'mirror_reversal', 'uncoupled_null']
    summary.update(relative_final_impulse_changes=changes, gates=gates,
                   verdict='SUPPORTED_FOR_THIS_SCALAR_FAMILY' if all(gates.values()) else
                   ('UNRESOLVED_NUMERICALLY' if not all(gates[k] for k in numerical_keys) else 'NOT_SUPPORTED'))
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    files = sorted(output.glob('*.csv'))+[output/'render_states.npz', output/'summary.json']
    (output/'sha256.json').write_text(json.dumps({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                for p in files}, indent=2)+'\n')
    print(summary['verdict'], gates, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('experiments/closure_ledger/runs/20261016_antipodal_kick'))
    run(parser.parse_args().output)
