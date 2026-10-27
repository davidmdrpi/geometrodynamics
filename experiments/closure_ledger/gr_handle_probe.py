"""Reproducible GR orientation, support and Maxwell-horizon controls for #349."""
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves.gr_handle_checks import (
    deck_matrix,radial_metric,killing_time,schwarzschild_profile,
    stress_components,maxwell_scattering,
)

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'experiments/closure_ledger/runs/20261024_gr_handle'


def save(name,data):
    (OUT/(name+'.json')).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[]
    for epsilon in [-1,1]:
        for offset in [-.7,0.,.7]:
            for point in [np.array([.3,1.2]),np.array([-.5,-1.4]),np.array([.6,.1])]:
                j=deck_matrix(epsilon,offset);mapped=j@point
                rows.append(dict(epsilon=epsilon,offset=offset,point=point.tolist(),mapped=mapped.tolist(),
                    metric_error=float(np.max(abs(j.T@radial_metric(*mapped)@j-radial_metric(*point)))),
                    square_error=float(np.max(abs(j@j-np.eye(2)))),
                    future_tangent_time_component=float((j@np.array([1.,.4]))[0]),
                    killing_time_source=float(killing_time(*point)) if abs(point[0])<abs(point[1]) else None,
                    killing_time_destination=float(killing_time(*mapped)) if abs(point[0])<abs(point[1]) else None))
    save('orientation',rows)
    supports=[]
    for angle in [.16,.3,.5]:
        b=np.sin(angle)**3
        for ratio in [1.0001,1.1,2.,10.]:
            r=b*ratio
            supports.append(dict(angle=angle,neck=b,radius=r,
                vacuum=schwarzschild_profile(r,b).tolist(),
                forced_unit_lapse=schwarzschild_profile(r,b,False).tolist()))
        supports.append(dict(angle=angle,neck=b,radius=b,
            throat_finite_lapse=stress_components(b,0.,1/(2*b)).tolist(),
            ellis_throat=stress_components(b,0.,1/b).tolist()))
    save('support',supports)
    scatter=[]
    for ell in [1,2]:
        for omega in [.2,.5,1.,2.]:
            for right in [120.,240.,480.]:
                case=maxwell_scattering(omega,ell,right=right)
                scatter.append(case)
                save('scattering',scatter)
                print(json.dumps({k:case[k] for k in ['ell','omega_b','right','reflection','horizon_flux','flux_error']}),flush=True)
    convergence=[]
    for ell in [1,2]:
        for omega in [.2,.5,1.,2.]:
            r=[x for x in scatter if x['ell']==ell and x['omega_b']==omega]
            convergence.append(dict(ell=ell,omega_b=omega,
                last_change=abs(r[2]['horizon_flux']-r[1]['horizon_flux']),
                previous_change=abs(r[1]['horizon_flux']-r[0]['horizon_flux'])))
    summary=dict(cases=len(scatter),max_flux_error=max(x['flux_error'] for x in scatter),
        max_metric_error=max(x['metric_error'] for x in rows),
        endpoint_convergence=convergence,bam_solution_established=False,
        interpretation='T enters the future horizon; no second-exterior transmission, pair recoil, or nonlinear Einstein evolution computed')
    save('summary',summary)
    print(json.dumps(summary,indent=2))


def manifest():
    files=sorted(p for p in OUT.glob('*.json') if p.name!='manifest.json')
    files += [Path(__file__),ROOT/'geometrodynamics/waves/gr_handle_checks.py',
              ROOT/'tests/test_gr_handle_checks.py',ROOT/'docs/gr_handle_orientation.md']
    save('manifest',dict(files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]))


if __name__=='__main__':main()
