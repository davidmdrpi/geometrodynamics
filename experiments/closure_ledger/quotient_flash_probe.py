"""Follow-up review: ray closure bands versus actual Maxwell phase data."""
import hashlib
import json
from pathlib import Path
import numpy as np
from geometrodynamics.waves.quotient_flash import closure_band,evolve,project_history,antipodal_overlap,driven_zero_initial,energy
from geometrodynamics.waves.pair_birth import PairBirth,Parameters
from geometrodynamics.waves.directional_flash_response import electric_integral
from geometrodynamics.waves.antipodal_throat_transport import AntipodalTransport,inertial_mouth

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'experiments/closure_ledger/runs/20261025_quotient_flash'


def save(name,data):
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/(name+'.json')).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def main():
    a=.3;bands=[closure_band(a,k) for k in range(4)];save('bands',bands)
    rows=[];source_checks=[];windows=[]
    for cutoff in [24,48,64]:
        model=PairBirth(Parameters(cutoff=cutoff,mouth_angle=a))
        for band in bands:
            for label,speed in [('endpoint',band['upper']),('interior',(band['lower']+band['upper'])/2)]:
                duration=(np.pi-2*a)/speed
                qc,pc=model.flash(model.initial(v=speed));q,p=qc.sum(0),pc.sum(0)
                parity=(-1.)**(model.ell+1)
                result=project_history(model.omega,parity,q,p,duration)
                e=result['initial_energy']
                rows.append(dict(cutoff=cutoff,k=band['k'],case=label,speed=speed,contact_time=duration,
                    source_A_antipodal_template=antipodal_overlap(model.omega,model.ell,qc[0],pc[0],duration),
                    conformal_template_control=antipodal_overlap(model.ell+.5,model.ell,qc[0],pc[0],(2*band['k']+1)*np.pi),
                    initial_flash_energy=e,projection_correction_fraction=result['correction_energy']/e,
                    raw_endpoint_defect_norm=np.sqrt(result['raw_endpoint_defect_energy']/e),
                    projected_endpoint_defect_norm=np.sqrt(result['projected_endpoint_defect_energy']/e),
                    projected_energy_fraction=result['projected_energy']/e,
                    interpretation='Prescribed free-field jump matching; projection changes the emitted Q/P data, not a physical absorption calculation'))
                def source(t):
                    theta=a+speed*t
                    x=np.array([[np.sin(theta),0,np.cos(theta)],[-np.sin(theta),0,np.cos(theta)]])
                    v=speed*np.array([[np.cos(theta),0,-np.sin(theta)],[-np.cos(theta),0,-np.sin(theta)]])
                    return model.source(x,v).sum(0)
                defects=[]
                scale=np.linalg.norm(source(0.))
                for t in np.linspace(0,duration,9):
                    defects.append(np.linalg.norm(source(duration-t)-parity*source(t))/scale)
                source_checks.append(dict(cutoff=cutoff,k=band['k'],case=label,max_current_parity_error=float(max(defects))))
        if cutoff==48:
            q,p=model.flash(model.initial())
            for n in [1,2,3,4]:
                theta=a+n*np.pi;point=np.array([np.sin(theta),0,np.cos(theta)])
                tangent=np.array([np.cos(theta),0,-np.sin(theta)])
                value=electric_integral(model,q[0],p[0],point,n*np.pi-np.pi/2,n*np.pi+np.pi/2)[0]@tangent
                windows.append(dict(n=n,start=n*np.pi-np.pi/2,stop=n*np.pi+np.pi/2,integral=float(value)))
    save('maxwell_closure',rows);save('current_parity',source_checks);save('whole_passage_windows',windows)
    mouths=[inertial_mouth(a),inertial_mouth(a+np.pi)]
    _,fa=mouths[0].frame(0);_,fb=mouths[1].frame(0);A=fb@(-np.eye(3))@fa.T
    rims=[]
    for handle in ['klein','torus']:
        H=AntipodalTransport(mouths,handle=handle).rim_map(0,1);loop=A@H
        rims.append(dict(handle=handle,antipodal=A.tolist(),handle_map=H.tolist(),loop=loop.tolist(),
            determinant=float(np.linalg.det(loop)),eigenvalues=np.linalg.eigvalsh(loop).tolist(),
            involution_residual=float(np.max(abs(loop@loop-np.eye(2))))))
    save('rim_loops',rims)
    summary=dict(bands=bands,max_current_parity_error=max(r['max_current_parity_error'] for r in source_checks),
        max_projected_endpoint_error=max(r['projected_endpoint_defect_norm'] for r in rows),
        cutoff_48=[r for r in rows if r['cutoff']==48],windows=windows,
        bam_solution_established=False,exact_original_maxwell_flash_closure=False)
    save('summary',summary);print(json.dumps(summary,indent=2))
    driven()


def driven():
    rows=[];a=.3
    for cutoff in [24,48]:
        m=PairBirth(Parameters(cutoff=cutoff,mouth_angle=a));parity=(-1.)**(m.ell+1)
        band=closure_band(a,1)
        for label,u in [('endpoint',band['upper']),('interior',(band['lower']+band['upper'])/2)]:
            T=(np.pi-2*a)/u
            def source(t):
                theta=a+u*t
                x=np.array([[np.sin(theta),0,np.cos(theta)],[-np.sin(theta),0,np.cos(theta)]])
                v=u*np.array([[np.cos(theta),0,-np.sin(theta)],[-np.cos(theta),0,-np.sin(theta)]])
                return m.source(x,v).sum(0)
            previous=None
            for order in [64,128,256,512]:
                dq,dp=driven_zero_initial(m.omega,source,T,order)
                q0,p0=.5*parity*dq,-.5*parity*dp
                qt,pt=evolve(m.omega,q0,p0,T);qt+=dq;pt+=dp
                scale=np.sqrt(max(energy(m.omega,q0,p0),1e-30))
                defect=np.sqrt(energy(m.omega,qt-parity*q0,pt+parity*p0))/scale
                change=None if previous is None else np.sqrt(energy(m.omega,q0-previous[0],p0-previous[1]))/scale
                interior=[]
                if order==512:
                    for t in [T*.2,T*.5,T*.8]:
                        def at(t):
                            q,p=evolve(m.omega,q0,p0,t);a,b=driven_zero_initial(m.omega,source,t,order)
                            return q+a,p+b
                        q,p=at(t);qr,pr=at(T-t)
                        interior.append(float(np.sqrt(energy(m.omega,q-parity*qr,p+parity*pr))/scale))
                rows.append(dict(cutoff=cutoff,case=label,speed=u,order=order,initial_energy=scale**2,
                    relative_endpoint_defect=float(defect),relative_initial_change=change,
                    relative_interior_defects=interior,
                    scope='Symmetrized forced Maxwell reference with prescribed continuous currents and selected homogeneous data; no strict-boundary matching or mouth backreaction'))
                previous=q0,p0
    save('driven_reference',rows)


def manifest():
    files=sorted(p for p in OUT.glob('*.json') if p.name!='manifest.json')
    files += [Path(__file__),ROOT/'geometrodynamics/waves/quotient_flash.py',
        ROOT/'experiments/closure_ledger/quotient_endpoint_review.py',
        ROOT/'geometrodynamics/waves/pair_birth.py',ROOT/'geometrodynamics/waves/directional_flash_response.py',
        ROOT/'geometrodynamics/waves/antipodal_throat_transport.py',ROOT/'tests/test_quotient_flash.py',
        ROOT/'docs/quotient_flash_review.md']
    save('manifest',dict(files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]))


if __name__=='__main__':main()
