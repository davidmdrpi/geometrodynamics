"""PR #348 review controls: source degeneracy, flash weighting and dimensions.

The flash supplies a modal ENERGY spectrum. Its projection into incoming neck
modes is not known: the three explicitly labelled source vectors are controls,
not a claim to reproduce the incident field of an actual birth trajectory.
"""
import numpy as np

from geometrodynamics.waves.coherent_mouth_recoil import (
    electric_stress_cross, junction, path_sum, scatter, size_impedance,
)
from geometrodynamics.waves.pair_birth import PairBirth, Parameters


SOURCES = {
    'one_sided': np.array([1.,0.,0.,0.],complex),
    'biased_pair': np.array([1.,0.,-.35*np.exp(.37j),0.])/np.sqrt(1+.35**2),
    'antisymmetric_control': np.array([1.,0.,-1.,0.])/np.sqrt(2),
}


def flash_spectrum(angle=.3, cutoff=32):
    model=PairBirth(Parameters(mouth_angle=angle,cutoff=cutoff))
    state=model.initial(v=1/3)
    q,p=model.flash(state)
    # Sum source copies BEFORE forming energy: retain their cross terms.
    q,p=q.sum(axis=0),p.sum(axis=0)
    energy=.5*(np.abs(p)**2+np.abs(model.omega*q)**2)
    shell=np.array([energy[model.ell==ell].sum() for ell in range(1,cutoff+1)])
    return dict(angular_mouth_radius=angle,cutoff=cutoff,footprint=model.par.footprint,
                birth_speed=1/3,charge=model.par.charge,
                omega=np.sqrt(np.arange(1,cutoff+1)*(np.arange(1,cutoff+1)+1)).tolist(),
                shell_energy=shell.tolist(),weights=(shell/shell.sum()).tolist(),
                total_energy=float(shell.sum()),neck_length=float(np.sin(2*angle)),
                interpretation='Actual PairBirth transverse Maxwell flash energy weights; assumed port sources and exact-focus caustic convention remain separate controls.')


def review_controls(save):
    source_controls=[]
    for name,source in SOURCES.items():
        for h in ['torus','klein']:
            for omega in [.7,1.137,2.43,5.16]:
                values={}
                for cap in ['before_caustic','after_caustic']:
                    r=path_sum(omega,.17,.565,survival=.8,source=source,handle=h,capture_phase=cap)
                    values[cap]={k:r[k] for k in ['exposure','impulse_per_exposure','impulse_per_source_energy']}
                source_controls.append(dict(source=name,handle=h,omega=omega,result=values))
    save('source_controls.json',source_controls)

    # Analytic reduced vector from the review, evaluated over a complete throat
    # phase period. Equal fractions are exact phase-shift equivalence, not a
    # topology preference inferred from a finite frequency grid.
    phases=2*np.pi*(np.arange(16384)+.5)/16384
    uniform=[]
    for ratio in [.05,.2,np.sqrt(2)-1]:
        r,t=junction(ratio**2)[0]
        for h,chi in [('torus',1),('klein',-1)]:
            dt=-.8*chi*np.exp(1j*phases)
            a=np.stack([1+r*dt,t*dt],axis=1)
            b=a@junction(ratio**2).T
            dp=(np.abs(a[:,0])**2-np.abs(a[:,1])**2
                +np.abs(b[:,0])**2-np.abs(b[:,1])**2)/np.sum(np.abs(a)**2,axis=1)
            uniform.append(dict(neck_over_mouth=float(ratio),spatial_dimension=3,
                                handle=h,negative_fraction=float(np.mean(dp < 0)),samples=len(phases)))
    save('uniform_throat_phase.json',uniform)

    spectra=[]; weighted=[]
    for angle in [.05,.16,.3,.5]:
        spec=flash_spectrum(angle);spectra.append(spec)
        weights=np.array(spec['weights'])
        for dimension in [2,3]:
            optimum=(3-2*np.sqrt(2))**(1/(dimension-1))
            for size_name,ratio in [('balanced',optimum),('vacuum_geometry',np.sin(angle)**2)]:
                z=size_impedance(ratio,dimension)
                for name,source in SOURCES.items():
                    for cap in ['before_caustic','after_caustic']:
                        for handle in ['torus','klein']:
                            results=[path_sum(omega,z,spec['neck_length'],survival=.8,source=source,
                                handle=handle,capture_phase=cap,terms=1) for omega in spec['omega']]
                            exp=np.array([r['exposure'] for r in results])
                            raw=np.array([r['impulse_per_source_energy'] for r in results])
                            normalized=np.array([r['impulse_per_exposure'] for r in results])
                            weighted.append(dict(angle=angle,dimension=dimension,size_model=size_name,
                                neck_over_mouth=float(ratio),length=spec['neck_length'],source=name,
                                capture_phase=cap,handle=handle,
                                weighted_normalized_ratio=float(weights@normalized),
                                weighted_exposure=float(weights@exp),
                                weighted_impulse_per_source_energy=float(weights@raw),
                                exposure_weighted_pressure=float(weights@raw/(weights@exp))))
    save('birth_flash_spectra.json',spectra)
    save('birth_weighted_pressure.json',weighted)

    # Geometry and charged cross term on a local 2D tangent patch. Axisymmetric
    # incoming radial traction sums to zero tangential force. A uniform external
    # electric field crossed with q/(2*pi*r) e_r integrates instead to q E.
    theta=2*np.pi*(np.arange(1024)+.5)/1024
    radial=np.stack([np.cos(theta),np.sin(theta)],axis=1)
    axisymmetric=radial.mean(axis=0)
    field=np.array([.4,-.2]);charge=.75;radius=.13
    force=sum((electric_stress_cross(charge/(2*np.pi*radius)*n,field)@n)
              for n in radial)*(2*np.pi*radius/len(theta))
    geometry=dict(axisymmetric_tangential_resultant=axisymmetric.tolist(),
        charged_cross_force=force.tolist(),expected_qE=(charge*field).tolist(),
        residual=float(np.linalg.norm(force-charge*field)),
        scope='Local Maxwell-stress identity; no global background field, aperture-mode scattering, or advanced recoil history derived.')
    save('stress_projection_controls.json',geometry)
    summary=dict(source_cases=len(source_controls),birth_weighted_cases=len(weighted),
        symmetric_max_caustic_difference=max(abs(x['result']['before_caustic']['impulse_per_exposure']
            -x['result']['after_caustic']['impulse_per_exposure']) for x in source_controls if x['source']=='antisymmetric_control'),
        one_sided_max_caustic_difference=max(abs(x['result']['before_caustic']['impulse_per_exposure']
            -x['result']['after_caustic']['impulse_per_exposure']) for x in source_controls if x['source']=='one_sided'),
        charged_cross_force_residual=geometry['residual'],
        centred_axisymmetric_tangential_force='zero by symmetry',
        balanced_optimum_radius_ratio={'surface_2d':float(3-2*np.sqrt(2)),'space_3d':float(np.sqrt(2)-1)})
    save('review_summary.json',summary)
    return summary
