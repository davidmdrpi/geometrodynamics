"""Targeted replication of review 6402240493; no BAM force law is changed."""
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from geometrodynamics.waves.deformable_mode_conversion import DeformableDefect,Parameters
from geometrodynamics.waves.directional_flash_response import coasting_focus_audit

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'experiments/closure_ledger/runs/20261025_review'


def save(name,value):
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/(name+'.json')).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def run(name,length,end,speed=.02,amplitude=.2,seed_ratio=None,dx=.1):
    model=DeformableDefect(Parameters(length=length,dx=dx))
    u0,v0=model.initial(pair=True,separation=20.,speed=speed)
    background=model.energy(u0,v0)
    def prepare(amp):
        u,v=u0.copy(),v0.copy()
        for center,direction in [(-3.5,-1),(3.5,1)]:
            model.add_packet(u,v,center=center,width=2.,amplitude=amp,direction=direction)
        return u,v
    if seed_ratio is not None:
        amplitude=brentq(lambda amp:model.energy(*prepare(amp))-background-seed_ratio*8/3,0.,.3,xtol=1e-14)
    u,v=prepare(amplitude);energy=model.energy(u,v);rows=[]
    dt=.2*dx;steps=round(end/dt);stride=round(.5/dt)
    for step in range(steps+1):
        if step%stride==0:
            roots=model.centers(u)
            rows.append(dict(time=step*dt,centers=roots,
                energy_error=model.energy(u,v)-energy,
                psi_max=float(np.max(abs(u[1:])))))
        if step<steps:u,v=model.step(u,v,dt)
    first=next((r['time'] for r in rows if len(r['centers'])==2 and r['centers'][1]-r['centers'][0]<2),None)
    lost=next((r['time'] for r in rows if len(r['centers'])!=2),None)
    late=[r for r in rows if r['time']>=end-15 and len(r['centers'])==2]
    speed_final=float(np.polyfit([r['time'] for r in late],[r['centers'][1] for r in late],1)[0]) if len(late)>4 else None
    result=dict(name=name,length=length,end=end,dx=dx,dt=dt,initial_speed=speed,
        amplitude=amplitude,seed_energy=energy-background,seed_over_continuum_rest=(energy-background)/(8/3),
        final_fitted_right_velocity=speed_final,first_separation_below_2=first,first_lost_pair_zeros=lost,
        maximum_relative_recorded_energy_error=max(abs(r['energy_error']) for r in rows)/energy,rows=rows)
    save(name,result)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
    return result


def geometry():
    a=.3;rows=[]
    for k in range(4):
        n=2*k+1;resonance=(np.pi-2*a)/(n*np.pi)
        lower=(np.pi-3*a)/(n*np.pi)
        for factor in [.0,.25,.75]:
            speed=resonance-factor*(resonance-lower)
            end=(np.pi-2*a)/speed;t=np.linspace(0,end,101)
            qa=a+speed*t;qb=-a-speed*(end-t)
            mismatch=float(np.max(abs(np.exp(1j*(qa+np.pi))-np.exp(1j*qb))))
            audit=coasting_focus_audit(a,speed)
            rows.append(dict(k=k,speed=speed,center_resonance=resonance,
                capture_band_lower=lower,contact_time=end,worldline_quotient_mismatch=mismatch,
                events=audit['events']))
    save('quotient_geometry',rows)
    return rows


def main():
    geometry()
    carryovers()
    cases=[run('pair_return_80',80.,140.),run('pair_return_60',60.,120.),
           run('pair_no_wave_80',80.,140.,amplitude=0.),
           run('pair_bam_scale',400.,100.,speed=.3,seed_ratio=.006),
           run('pair_bam_no_wave',400.,100.,speed=.3,amplitude=0.)]
    seed,control=cases[-2:]
    dv=seed['final_fitted_right_velocity']-control['final_fitted_right_velocity']
    p=lambda v:(4/3)*v/np.sqrt(1-v*v)
    dp=p(seed['final_fitted_right_velocity'])-p(control['final_fitted_right_velocity'])
    save('summary',dict(cases=[{k:v for k,v in r.items() if k!='rows'} for r in cases],
        bam_scale_velocity_change=dv,bam_scale_inferred_momentum_change=dp,
        momentum_estimator='Continuum kink M gamma v, fitted late position; not direct field momentum',
        bam_solution_established=False))
    refinement()


def refinement():
    fine=run('pair_bam_scale_fine',400.,100.,speed=.3,seed_ratio=.006,dx=.05)
    control=run('pair_bam_no_wave_fine',400.,100.,speed=.3,amplitude=0.,dx=.05)
    p=lambda v:(4/3)*v/np.sqrt(1-v*v)
    save('refinement',dict(cases=[{k:v for k,v in r.items() if k!='rows'} for r in [fine,control]],
        velocity_change=fine['final_fitted_right_velocity']-control['final_fitted_right_velocity'],
        inferred_momentum_change=p(fine['final_fitted_right_velocity'])-p(control['final_fitted_right_velocity'])))


def carryovers():
    from geometrodynamics.waves.pair_birth import PairBirth,Parameters as BirthParameters
    model=PairBirth(BirthParameters(cutoff=48,mouth_angle=.3))
    def flash_energy(speed):
        q,p=model.flash(model.initial(v=speed));q,p=q.sum(0),p.sum(0)
        return float(.5*np.sum(abs(p)**2+abs(model.omega*q)**2))
    kick=1/3
    coast=brentq(lambda v:2/np.sqrt(1-v*v)+flash_energy(v)-2/np.sqrt(1-kick*kick),0.,kick)
    projections=json.loads((ROOT/'experiments/closure_ledger/runs/20261022_directional_flash/directional_projections.json').read_text())
    sources=[r for r in projections if r['angle']==.3 and r['return_number']==1 and r['offset_ratio']==0]
    save('carryovers',dict(kick_speed=kick,flash_paid_speed=coast,flash_energy=flash_energy(coast),
        paid_first_flight=coasting_focus_audit(.3,coast),source_projections=sources,
        note='Flash-paid scalar energy balance uses existing PairBirth coefficients; not a replay of the #345 multicyle run'))


def manifest():
    files=sorted(p for p in OUT.glob('*.json') if p.name!='manifest.json')
    files += [Path(__file__),ROOT/'docs/pr349_review_assessment.md',
              ROOT/'tests/test_pr349_review.py',ROOT/'geometrodynamics/waves/pair_birth.py',
              ROOT/'geometrodynamics/waves/deformable_mode_conversion.py',
              ROOT/'geometrodynamics/waves/directional_flash_response.py']
    save('manifest',dict(files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]))


if __name__=='__main__':main()
