"""Reproduce the coherent recoil mechanism audit; no old trajectories rerun."""
from dataclasses import asdict, replace
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from geometrodynamics.waves.coherent_mouth_recoil import (
    PathParameters, closed_paths, finite_mass_scatter, focus_overlap,
    junction, path_sum, scatter,
)

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO/'experiments/closure_ledger/runs/20261022_coherent_recoil'
IMAGE = REPO/'docs/images/history_throat/coherent_recoil.png'


def save(name, data):
    (ROOT/name).write_text(json.dumps(data, separators=(',', ':'), allow_nan=False)+'\n')


def sizing_controls():
    # Import the repository geometry rather than treating a chosen radius
    # ratio as the vacuum bridge prediction. This remains initial-data geometry.
    from geometrodynamics.waves.physical_throat import VacuumThroat
    rows=[]
    for a in [.025,.05,.1,.16,.3,.5]:
        geometry=VacuumThroat(a)
        ratio=geometry.neck_radius()/geometry.mouth_f()
        r,t=junction(ratio**2)[0]
        rows.append(dict(angular_mouth_radius=a,neck_over_mouth=ratio,
            geometric_length=geometry.length(),transmitted_energy_fraction=float(t*t),
            maximum_balanced_interference_impulse_per_energy=float(2*abs(r*t))))
    save('geometric_sizes.json',dict(rows=rows,
        interpretation='VacuumThroat scalar-flat initial-data geometry plus the assumed equal-speed area impedance. Not a traversable dynamic throat solution.'))
    return rows


def write_manifest():
    import scipy
    save('environment.json',dict(python=platform.python_version(),numpy=np.__version__,
                                 scipy=scipy.__version__,matplotlib=matplotlib.__version__))
    paths = list(ROOT.glob('*.json'))+[IMAGE,Path(__file__),ROOT/'README.md',
        REPO/'docs/coherent_mouth_recoil.md',
        REPO/'geometrodynamics/waves/coherent_mouth_recoil.py',
        REPO/'geometrodynamics/waves/physical_throat.py',
        REPO/'tests/test_coherent_mouth_recoil.py']
    hashes={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths if p.name!='sha256.json'}
    save('sha256.json',hashes)


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    sizing_controls()
    radii = np.array([.05, .1, .2, np.sqrt(2)-1, .6, .8, 1.])
    lengths = np.array([1/8, 1/4, 1/2, 1., 2.])
    frequencies = np.linspace(0, 8, 65)
    caps = ['before_caustic', 'after_caustic']; handles = ['torus', 'klein']
    shape = (len(caps), len(handles), len(radii), len(lengths), len(frequencies))
    pressure = np.zeros(shape); exposure = np.zeros(shape)
    max_unitarity = 0.; max_walk_error = 0.
    for ci, cap in enumerate(caps):
        for hi, handle in enumerate(handles):
            for ri, radius in enumerate(radii):
                for li, length in enumerate(lengths):
                    for wi, omega in enumerate(frequencies):
                        result = path_sum(omega, radius**2, np.pi*length, survival=.8,
                                          handle=handle, capture_phase=cap, terms=1)
                        pressure[ci, hi, ri, li, wi] = result['impulse_per_exposure']
                        exposure[ci, hi, ri, li, wi] = result['exposure']
                        max_unitarity = max(max_unitarity, result['unitary_error'])
                    check = path_sum(1.137, radius**2, np.pi*length, survival=.8,
                                     handle=handle, capture_phase=cap, terms=256)
                    max_walk_error = max(max_walk_error, check['walk_sum_relative_error'])
    axes = dict(capture_phase=caps, handle=handles, neck_over_mouth=radii.tolist(),
                length_over_pi=lengths.tolist(), omega=frequencies.tolist())
    save('spectrum.json', dict(axes=axes, axis_order=list(axes), survival=.8,
                              source=[1, 0, -1, 0], impulse_per_exposure=pressure.tolist()))
    save('exposure.json', dict(axes=axes, axis_order=list(axes), exposure=exposure.tolist(),
                              interpretation='Externally fed resolvent, not a finite-seed amplification'))
    print('frequency grid complete', flush=True)

    closed = []; worst_energy = 0.
    for cap in caps:
        for handle in handles:
            for radius in radii:
                for length in lengths:
                    for fraction, phase in [(0., 0.), (.1, np.pi), (.5, 0.), (.5, np.pi)]:
                        p = PathParameters(throat_bins=round(24*length), impedance_ratio=radius**2,
                            capture_phase=cap, handle=handle, throat_seed_fraction=fraction,
                            throat_seed_phase=phase, pulse_bins=2, periods=32)
                        coherent = closed_paths(p)
                        incoherent = closed_paths(p, coherent=False)
                        worst_energy = max(worst_energy, coherent['max_energy_error'],
                                           incoherent['max_energy_error'])
                        closed.append(dict(parameters=asdict(p), coherent={k:v for k,v in coherent.items() if k!='sample'},
                                           incoherent={k:v for k,v in incoherent.items() if k!='sample'}))
    save('closed_paths.json', closed)
    print('closed finite-seed grid complete', flush=True)

    # Check exact refinement of the same pulse/delays, not a changed pulse width.
    p = PathParameters(throat_bins=6, throat_seed_fraction=.5, pulse_bins=1, periods=32)
    coarse = closed_paths(p)
    fine = closed_paths(replace(p, surface_bins=48, throat_bins=12))
    refinement = float(np.max(np.abs(np.array(coarse['impulse'])-fine['impulse'])))
    # Broader finite pulse sensitivity, including non-simple delay ratios.
    width_controls = []
    for bins in [5, 7, 11, 24, 48]:
        for width in [1, 2, 4]:
            for handle in handles:
                p = PathParameters(throat_bins=bins, pulse_bins=width, handle=handle,
                                   throat_seed_fraction=.5, periods=32)
                result = closed_paths(p)
                worst_energy=max(worst_energy,result['max_energy_error'])
                width_controls.append(dict(parameters=asdict(p), result=result))
    save('pulse_controls.json', width_controls)

    duration_controls = []
    for periods in [16, 32, 64, 128, 256]:
        for label,p in [('surface_seed_exception', PathParameters(throat_bins=24,
                impedance_ratio=.64, capture_phase='before_caustic', pulse_bins=2)),
                ('balanced_opposite_phase', PathParameters(throat_seed_fraction=.5,pulse_bins=2))]:
            p=replace(p,periods=periods)
            result=closed_paths(p)
            worst_energy=max(worst_energy,result['max_energy_error'])
            duration_controls.append(dict(label=label,parameters=asdict(p),result=result))
    save('duration_controls.json',duration_controls)

    clock_controls=[]
    for offset in [0., .2825, .70625, 2., np.pi+.565]:
        for handle in handles:
            for omega in np.linspace(.013,8.013,65):
                result=path_sum(omega,(np.sqrt(2)-1)**2,.565,survival=.8,
                                clock_offset=offset,handle=handle,terms=256)
                max_walk_error=max(max_walk_error,result['walk_sum_relative_error'])
                max_unitarity=max(max_unitarity,result['unitary_error'])
                clock_controls.append(dict(offset=offset,handle=handle,omega=float(omega),
                    impulse_per_exposure=result['impulse_per_exposure'],exposure=result['exposure'],
                    throat_delays=result['throat_delays'],advanced_boundary=result['advanced_boundary']))
    save('clock_controls.json',clock_controls)

    mass_rows = []; energy_error = 0.; momentum_error = 0.
    for velocity in [-.3, 0., .3]:
        for energy in [.001, .0063, .01, .1, 1.]:
            for radius in radii:
                for phi in [0., np.pi/2, np.pi]:
                    a = np.sqrt(energy/2)*np.array([1, np.exp(1j*phi)])
                    result = finite_mass_scatter(a, radius**2, mass=1,
                                                momentum=velocity/np.sqrt(1-velocity**2))
                    energy_error = max(energy_error, abs(result['energy_residual']))
                    momentum_error = max(momentum_error, abs(result['momentum_residual']))
                    mass_rows.append(dict(velocity=velocity, energy_over_mass=energy,
                        neck_over_mouth=float(radius), relative_phase=phi,
                        impulse=result['impulse'], outgoing_energy=float(np.sum(np.abs(result['outgoing'])**2)),
                        rest_energy_scale=result['rest_energy_scale'],
                        energy_residual=result['energy_residual'], momentum_residual=result['momentum_residual']))
    save('finite_mass.json', mass_rows)

    summary = dict(spectral_cases=int(np.prod(shape)), closed_cases=len(closed),
        pulse_controls=len(width_controls), duration_controls=len(duration_controls),
        clock_boundary_cases=len(clock_controls), finite_mass_cases=len(mass_rows),
        max_closed_energy_error=worst_energy, max_scattering_unitarity_error=max_unitarity,
        max_walk_sum_relative_error=max_walk_error, delay_refinement_impulse_error=refinement,
        max_collision_energy_error=energy_error, max_collision_momentum_error=momentum_error,
        optimum_balanced_neck_over_mouth=float(np.sqrt(2)-1),
        spectrum_negative_fraction={cap:{h:float(np.mean(pressure[ci,hi]<-1e-10))
            for hi,h in enumerate(handles)} for ci,cap in enumerate(caps)},
        closed_negative_counts={str(f):sum(sum(x['coherent']['impulse']) < -1e-10
              for x in closed if x['parameters']['throat_seed_fraction']==f) for f in [0., .1, .5]},
        closed_sample_counts={str(f):sum(x['parameters']['throat_seed_fraction']==f for x in closed)
                              for f in [0., .1, .5]},
        interpretation='Negative PORT impulse, not demonstrated S3 separation attraction. No moving-mouth history or GR field closure.')
    save('summary.json', summary)

    # Standalone scientific figure, not an illustrative BAM geometry animation.
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axs = plt.subplots(2, 3, figsize=(14, 8.3), constrained_layout=True)
    fig.suptitle('Wave-exchange recoil: conservative mechanisms and their limits', fontsize=16)
    rs = np.linspace(.001, 1, 240); phases = np.linspace(0, 2*np.pi, 240)
    reflection = ((1-rs**2)/(1+rs**2))**2
    values = 2*np.sqrt(reflection*(1-reflection))[None,:]*np.cos(phases[:,None])
    im = axs[0,0].imshow(values, origin='lower', aspect='auto', extent=[0,1,0,2],
                         cmap='RdBu_r', vmin=-1, vmax=1)
    axs[0,0].axvline(np.sqrt(2)-1,color='k',lw=.8,ls='--')
    axs[0,0].set(xlabel='Neck / mouth radius', ylabel='Relative phase / π',
                 title='a  Balanced coherent input: Δp / E')
    fig.colorbar(im,ax=axs[0,0],label='Local port impulse')
    for ph,label in [(0,'In phase'),(np.pi,'Opposite phase')]:
        axs[0,1].plot(rs,2*np.sqrt(reflection*(1-reflection))*np.cos(ph),label=label)
    axs[0,1].axhline(0,color='k',lw=.8,ls='--',label='Phase average')
    axs[0,1].set(xlabel='Neck / mouth radius',ylabel='Δp / E',title='b  Partial mixing, not a sign-only force')
    axs[0,1].legend(fontsize=8)
    for hi,handle in enumerate(handles):
        axs[0,2].plot(frequencies,pressure[1,hi,3,1],label=handle.title())
    axs[0,2].axhline(0,color='k',lw=.6)
    axs[0,2].set(xlabel='Frequency × sphere radius / c',ylabel='Impulse / incident exposure',
                 title='c  All-path susceptibility; externally fed')
    axs[0,2].legend(fontsize=8)
    for frac,ph,label in [(0,0,'Surface seed only'),(.5,0,'Balanced, in phase'),(.5,np.pi,'Balanced, opposite phase')]:
        p = PathParameters(throat_seed_fraction=frac, throat_seed_phase=ph, periods=32)
        data=np.array(closed_paths(p,sample=True)['sample'])
        axs[1,0].plot(data[:,0]/np.pi,data[:,1]+data[:,2],label=label)
    axs[1,0].axhline(0,color='k',lw=.6)
    axs[1,0].set(xlabel='Elapsed time / π',ylabel='Cumulative port impulse / seed energy',
                 title='d  Finite seed, fixed supports (no mouth motion)')
    axs[1,0].legend(fontsize=8)
    for v in [0,.3]:
        rows=[x for x in mass_rows if x['velocity']==v and
              abs(x['neck_over_mouth']-(np.sqrt(2)-1))<1e-12 and x['relative_phase']==np.pi]
        axs[1,1].plot([x['energy_over_mass'] for x in rows],[-x['impulse']/x['energy_over_mass'] for x in rows],
                       'o-',label=f'Initial v={v}c')
    axs[1,1].set(xscale='log',xlabel='Incident wave energy / mouth rest energy',
                 ylabel='−Δp / incident energy',title='e  One finite-mass encounter, exact E/p balance')
    axs[1,1].legend(fontsize=8)
    x=np.linspace(0,4,150)
    axs[1,2].plot(x,focus_overlap(x),color='#346a8b')
    axs[1,2].set(xlabel='Mouth radius / focal width',ylabel='Captured energy fraction',ylim=(0,1.05),
                 title='f  Gaussian aperture overlap (assumed profile)')
    fig.supxlabel('Port stress ≠ global orbital recoil. Klein sign, impedance law and coherent preparation are explicit assumptions.',fontsize=10)
    IMAGE.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(IMAGE,dpi=145);plt.close(fig)
    write_manifest()
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
