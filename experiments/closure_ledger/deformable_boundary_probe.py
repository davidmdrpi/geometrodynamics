"""PR #349 continuation: deformable-wave benchmark and time-oriented gluing.

All models and controls are specified below. The benchmark's continuous local
action is NOT installed in the antipodal-only BAM router. Each completed case
is checkpointed; rerunning resumes cases with identical source/parameter hashes.
"""
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.alice_boundary import (
    METRIC,faraday,frequency_map,gluing_jacobian,solve_history,stress,transport_faraday,
    time_deck,time_reversing_initial,
)
from geometrodynamics.waves.deformable_mode_conversion import (
    DeformableDefect,Parameters,step_scattering,
)
from geometrodynamics.waves.pair_birth import PairBirth,Parameters as BirthParameters

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'experiments/closure_ledger/runs/20261023_deformable_boundary'
MODEL=ROOT/'geometrodynamics/waves/deformable_mode_conversion.py'


def save(name,value):
    path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def run_case(name,*,pair=False,dx=.1,dt=.02,end=80.,amplitude=.04,channel=1,
             coupling=1.,offset=0.,reverse_check=False):
    par=Parameters(length=400. if pair else 200.,dx=dx,twisted=not pair,coupling=coupling)
    spec=dict(name=name,pair=pair,parameters=asdict(par),dt=dt,end=end,amplitude=amplitude,
              channel=channel,offset=offset,reverse_check=reverse_check,
              source_sha256=hashlib.sha256(MODEL.read_bytes()).hexdigest())
    signature=hashlib.sha256(json.dumps(spec,sort_keys=True).encode()).hexdigest()
    record=OUT/(name+'.json')
    if record.exists():
        old=json.loads(record.read_text())
        if old['signature']==signature:
            state=OUT/(name+'.npz')
            if not state.exists() or hashlib.sha256(state.read_bytes()).hexdigest()!=old['state_sha256']:
                raise RuntimeError('Missing or corrupt state archive: '+name)
            print('Reusing hash-verified case',name,flush=True)
            return old
        raise RuntimeError('Existing case has different source or parameters: '+name)
    model=DeformableDefect(par);u,v=model.initial(pair=pair,separation=20.,speed=.02 if pair else 0.)
    background=model.energy(u,v)
    if pair:
        for center,direction in [(-3.5,-1),(3.5,1)]:
            model.add_packet(u,v,center=center+offset,width=2.,amplitude=amplitude,direction=direction,channel=channel)
    else:
        model.add_packet(u,v,center=25.+offset,width=5.,amplitude=amplitude,channel=channel)
    u0,v0=u.copy(),v.copy();e0=model.energy(u,v);p0=model.momenta(u,v)
    rows=[];snapshots=[];times=[]
    stride=max(1,round(.5/dt));snapshot_stride=max(1,round(5./dt));steps=round(end/dt)
    for step in range(steps+1):
        if step%stride==0 or step==steps:
            momenta=model.momenta(u,v);energy=model.energy(u,v)
            rows.append(dict(time=step*dt,centers=model.centers(u),energy=energy,
                             energy_error=energy-e0,momenta=momenta.tolist(),
                             momentum_error=float(np.sum(momenta-p0))))
        if step%snapshot_stride==0 or step==steps:
            snapshots.append(u[:,::max(1,round(.4/dx))].copy());times.append(step*dt)
        if step<steps:u,v=model.step(u,v,dt)
    reversal=None
    if reverse_check:
        ur,vr=u.copy(),-v.copy()
        for _ in range(steps):ur,vr=model.step(ur,vr,dt)
        reversal=dict(max_field_error=float(np.max(abs(ur-u0))),max_velocity_error=float(np.max(abs(vr+v0))))
    seed=e0-background
    result=dict(spec=spec,signature=signature,background_energy=background,seed_energy=seed,
                maximum_energy_error=max(abs(r['energy_error']) for r in rows),
                maximum_momentum_error=max(abs(r['momentum_error']) for r in rows),
                reversal=reversal,rows=rows)
    state_path=OUT/(name+'.npz')
    np.savez_compressed(state_path,x=model.x,u=u,v=v,u_initial=u0,v_initial=v0,
                        snapshot_x=model.x[::max(1,round(.4/dx))],snapshot_t=times,snapshots=snapshots)
    result['state_sha256']=hashlib.sha256(state_path.read_bytes()).hexdigest()
    save(name+'.json',result)
    print(json.dumps(dict(case=name,seed=seed,energy_error=result['maximum_energy_error'],
                          momentum_error=result['maximum_momentum_error'],final=rows[-1],reversal=reversal)),flush=True)
    return result


def boundary_controls():
    rows=[]
    e=np.array([.3,-.5,.2]);b=np.array([.1,.4,.7]);f=faraday(e,b)
    for angle in [0.,.3,np.pi/2,np.pi]:
        r=np.array([[np.cos(angle),-np.sin(angle),0.],[np.sin(angle),np.cos(angle),0.],[0.,0.,-1.]])
        for time_orientation in [-1,1]:
            for clock in [-2.,0.,.7]:
                j=gluing_jacobian(r,time_orientation,[.2,0.,0.],[-.15,.1,0.])
                tf=transport_faraday(f,j)
                k=np.linalg.solve(j,np.array([1.,1.,0.,0.]))
                u,v=frequency_map(r,1.7,clock,time_orientation)
                # A known consistent source, an unrelated source, and the
                # homogeneous problem distinguish closure from source selection.
                target=np.array([.2+.4j,-.3+.1j,.7-.2j])
                sources={'compatible':target-u@target-v@target.conj(),
                         'unselected':np.array([1.,.2j,.3]),'homogeneous':np.zeros(3)}
                outcomes={}
                for name,source in sources.items():
                    solved=solve_history(u,v,source)
                    outcomes[name]={key:solved[key] for key in ['status','rank','nullity','residual']}
                rows.append(dict(angle=angle,time_orientation=time_orientation,clock_offset=clock,
                                 metric_error=float(np.max(abs(j.T@METRIC@j-METRIC))),
                                 stress_error=float(np.max(abs(stress(tf)-j.T@stress(f)@j))),
                                 local_energy=float(stress(tf)[0,0]),oriented_frequency=float(k[0]),
                                 closure=outcomes))
    save('alice_boundary_controls.json',rows)
    return rows


def gap_screen():
    rows=[]
    for angle in [.16,.3,.5]:
        model=PairBirth(BirthParameters(cutoff=96,mouth_angle=angle))
        q,p=model.flash(model.initial());q,p=q.sum(axis=0),p.sum(axis=0)
        energies=.5*(abs(p)**2+abs(model.omega*q)**2)
        for label,ratio in [('vacuum_spatial_slice',np.sin(angle)**2),('vacuum_4D_spatial_slice',np.sin(angle)),
                            ('rim_balanced',3-2*np.sqrt(2)),
                            ('area_balanced',np.sqrt(2)-1),('wide_neck',1.)]:
            radius=np.sin(angle)*ratio
            for cross_section,eigenvalue in [('circle_m1',1.),('sphere_l1',2.),('three_sphere_l1',3.)]:
                cutoff=np.sqrt(eigenvalue)/radius
                for doppler in [1.,2.]:
                    fraction=float(np.sum(energies[model.omega*doppler>cutoff])/np.sum(energies))
                    rows.append(dict(angle=angle,size_model=label,neck_over_mouth=float(ratio),neck_radius=float(radius),
                                     cross_section=cross_section,gap=float(cutoff),doppler_factor=doppler,
                                     flash_energy_above_gap=fraction,
                                     interpretation='Optimistic product-waveguide propagation screen; not a transmission coefficient, tensor GR reduction or traversability result'))
    save('directional_gap_screen.json',rows)
    return rows


def nonlinear_time_quotient():
    cell=DeformableDefect(Parameters(length=200.,dx=.1,twisted=True))
    uc,vc=cell.initial();cell.add_packet(uc,vc,center=25.,width=5.,amplitude=.04)
    u0,v0=time_reversing_initial(uc,vc)
    cover=DeformableDefect(Parameters(length=400.,dx=.1))
    up,vp=u0.copy(),v0.copy();um,vm=u0.copy(),v0.copy()
    e0=cover.energy(u0,v0);rows=[];fields_saved=[]
    for step in range(4001):
        if step%25==0:
            rows.append(dict(time=step*.02,
                             field_closure_error=float(np.max(abs(um-time_deck(up)))),
                             velocity_closure_error=float(np.max(abs(vm+time_deck(vp)))),
                             cover_energy_error=cover.energy(up,vp)-e0,
                             future_cell_centers=cell.centers(up[:,:cell.n]),
                             past_cell_centers=cell.centers(um[:,:cell.n])))
        if step%250==0:fields_saved.append(up[:,::4].copy())
        if step<4000:
            up,vp=cover.step(up,vp,.02)
            um,vm=cover.step(um,vm,-.02)
    state=OUT/'nonlinear_time_quotient.npz'
    np.savez_compressed(state,x=cover.x,u_initial=u0,v_initial=v0,u_future=up,v_future=vp,
                        u_past=um,v_past=vm,snapshots=fields_saved,snapshot_x=cover.x[::4])
    result=dict(identification='(t,x,u) ~ (-t,x+200,D u), D(phi,a,b)=(-phi,b,a)',
                spacetime='Flat 1+1 Lorentzian quotient; matter backreaction on metric not solved',
                initial_cover_energy=e0,state_sha256=hashlib.sha256(state.read_bytes()).hexdigest(),
                max_field_closure_error=max(r['field_closure_error'] for r in rows),
                max_velocity_closure_error=max(r['velocity_closure_error'] for r in rows),rows=rows)
    save('nonlinear_time_quotient.json',result)
    print('Nonlinear time quotient closure',result['max_field_closure_error'],result['max_velocity_closure_error'],flush=True)
    return result


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    boundary=boundary_controls();gaps=gap_screen()
    scattering=[dict(omega=w,coupling=g,**step_scattering(w,coupling=g))
                for g in [-1.,0.,1.] for w in [2.3,2.5,2.8,3.5,5.]]
    save('step_scattering.json',scattering)
    cases=[]
    schedule=[('single_heavy',dict()),
              ('single_heavy_fine',dict(dx=.05,dt=.01)),
              ('single_heavy_finer',dict(dx=.025,dt=.005)),
              ('single_light',dict(channel=2)),
              ('single_disconnected',dict(coupling=0.)),
              ('single_offset',dict(offset=7.)),
              ('single_reversed',dict(end=60.,reverse_check=True)),
              ('single_closed_return',dict(end=400.)),
              ('pair_no_wave',dict(pair=True,amplitude=0.,end=240.)),
              ('pair_heavy',dict(pair=True,amplitude=.2,end=240.)),
              ('pair_heavy_fine',dict(pair=True,amplitude=.2,dx=.05,dt=.01,end=240.)),
              ('pair_heavy_time_fine',dict(pair=True,amplitude=.2,dt=.01,end=240.)),
              ('pair_light',dict(pair=True,amplitude=.2,channel=2,end=120.)),
              ('pair_offset',dict(pair=True,amplitude=.2,offset=1.,end=120.))]
    save('schedule.json',[dict(name=name,**kw) for name,kw in schedule])
    for name,kwargs in schedule:cases.append(run_case(name,**kwargs))
    quotient=nonlinear_time_quotient()
    summary=dict(cases=len(cases),scattering_cases=len(scattering),boundary_cases=len(boundary),gap_cases=len(gaps),
                 max_metric_error=max(r['metric_error'] for r in boundary),
                 max_stress_error=max(r['stress_error'] for r in boundary),
                 max_scattering_flux_error=max(abs(r['flux_residual']) for r in scattering),
                 max_nonlinear_time_closure_error=quotient['max_field_closure_error'],
                 benchmark='Reciprocal deformable-defect analogue; continuous coupling, not BAM',
                 bam_status='Boundary-matched antipodal moving-mouth GR solution remains unresolved')
    save('summary.json',summary)
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
