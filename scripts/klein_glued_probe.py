"""Test #344's incident data on the actual glued metric, without forcing collapse.

Uses existing cotangent FEM with the selected constant +1/4 term. Contact uses
only #338's neck in rim units. Element energies partition exactly into sphere
and neck, including the appropriate shares at common rim vertices.
"""
import argparse,json
from pathlib import Path
import numpy as np
from geometrodynamics.viz.klein_handles import build_surface,HandleSpec
from klein_release_model import KleinRelease


def element_energy(surf,area,KT,u,v):
    tri=surf.tris
    return (.5*np.einsum('ti,tij,tj->t',u[tri],KT,u[tri])
            +area/6*np.sum(v[tri]**2+.25*u[tri]**2,axis=1))


def run(n_sphere=6000,glued=True,end=2.6,output_dt=.01):
    reference=KleinRelease(detach=False)
    handles=[HandleSpec((1,0,0),(0,-1,0),True)] if glued else []
    surf=build_surface(handles=handles,n_sphere=n_sphere,n_ring=32,ds_target=.01)
    K,M=surf.operators();area,KT=surf.local_stiffness()
    mu=surf.points[:,0]/np.linalg.norm(surf.points,axis=1)
    a,p=reference.incident(0)
    Y=reference.basis(mu);u=np.where(surf.on_sphere,Y@a,0.);v=np.where(surf.on_sphere,Y@p,0.)
    accel=lambda q:-(K@q)/M-.25*q
    vec=np.random.default_rng(0).standard_normal(len(u))
    for _ in range(80):
        nxt=-accel(vec);lam=np.linalg.norm(nxt)/np.linalg.norm(vec);vec=nxt/np.linalg.norm(nxt)
    sub=max(1,int(np.ceil(output_dt/min(.002,1.8/np.sqrt(1.05*lam)))));dt=output_dt/sub
    times=np.linspace(0,end,round(end/output_dt)+1)
    neck=surf.tri_tube>=0
    near_b=surf.on_sphere&(mu<-np.cos(1.))
    sphere_b=(~neck)&near_b[surf.tris].all(axis=1)
    pole=int(np.argmax(np.where(surf.on_sphere,mu,-2)))
    if glued:
        ids=surf.tube_ids[0];s=surf.s[ids[:,0]];rim=np.sin(surf.a)
        r0=np.sqrt(s*s+surf.neck.b**2)/rim
        taper=(1-(s/surf.neck.S)**2)**2
    records=[];ring_records=[];acc=accel(u)
    eps=.025;half=eps*np.sqrt(-np.log(1-np.sqrt(.5)))
    for index,t in enumerate(times):
        E=element_energy(surf,area,KT,u,v)
        row=[t,float(u[pole]),float(E.sum()),float(E[neck].sum()),float(E[~neck].sum()),float(E[sphere_b].sum())]
        if glued:
            drive=u[ids].mean(axis=1)/reference.peak
            r=r0-.18*drive*taper
            positive=drive*taper>1e-12
            critical=float(np.min((r0[positive]-half)/(drive[positive]*taper[positive]))) if positive.any() else np.inf
            row += [float(drive[len(s)//2]),float(r.min()),critical]
            ring_records.append(drive)
        else:row += [np.nan,np.nan,np.nan]
        records.append(row)
        if index==len(times)-1:break
        for _ in range(sub):
            v+=.5*dt*acc;u+=dt*v;acc=accel(u);v+=.5*dt*acc
    rec=np.array(records);E0=rec[0,2]
    report=dict(glued=glued,requested_sphere_vertices=n_sphere,vertices=surf.n_vertices,triangles=len(surf.tris),dt=dt,
                horizon=end,initial_energy=E0,max_relative_energy_drift=float(np.max(abs(rec[:,2]-E0))/E0),
                max_partition_residual=float(np.max(abs(rec[:,2]-rec[:,3]-rec[:,4]))),
                max_neck_energy_fraction=float(rec[:,3].max()/E0),neck_peak_time=float(times[np.argmax(rec[:,3])]),
                max_B_side_energy_fraction=float(rec[:,5].max()/E0),B_peak_time=float(times[np.argmax(rec[:,5])]),
                pole_peak_over_reference=float(rec[:,1].max()/reference.peak))
    if glued:
        report.update(waist_peak_over_reference=float(rec[:,6].max()),waist_peak_time=float(times[np.argmax(rec[:,6])]),
                      minimum_radius_in_rim_units=float(rec[:,7].min()),contour_radius_threshold=half,
                      contact_detected=bool(np.any(rec[:,7]<=half)),minimum_sampled_gain_for_contact=float(rec[:,8].min()),
                      gain=.18,neck=dict(rim_radius=rim,waist_sphere_units=surf.neck.b,waist_rim_units=surf.neck.b/rim,
                      arclength=2*surf.neck.S,axial_halfheight_rim_units=surf.neck.b*np.arcsinh(surf.neck.S/surf.neck.b)/rim),
                      topology_change_performed=False,release_energy_injected=0.)
    return rec,report


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path('docs/images/klein_release/glued_probe'));args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    results={};data={}
    for label,n,glued in [('glued_6000',6000,True),('glued_12000',12000,True),('bare_12000',12000,False)]:
        rec,r=run(n,glued);results[label]=r;data[label]=rec
        np.savetxt(args.output/(label+'.csv'),rec,delimiter=',',header='t,pole,total_energy,neck_energy,sphere_energy,B_energy,waist_over_reference,min_rim_radius,critical_gain',comments='',fmt='%.12g')
        print(label,json.dumps(r),flush=True)
    results['verdict']='NO_CONTACT_AT_GAIN_0_18_ON_GLUED_FIXED_GEOMETRY'
    results['scope']='Fixed intrinsic metric; no evolving neck geometry or topology change. Null does not rule out a dynamical thin-neck instability.'
    results['energy_partition']='Per-element kinetic, gradient and +1/4 potential energy; shared rim mass belongs fractionally to each side.'
    (args.output/'report.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
    render(data,results,args.output/'comparison.png')


def render(data,results,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('dark_background')
    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained');fig.patch.set_facecolor('#0c1622')
    m=KleinRelease(detach=False)
    for label,color in [('glued_6000','#efb96e'),('glued_12000','#55d7cc')]:
        v=data[label];axes[0,0].plot(v[:,0],v[:,6],color=color,label=label.replace('_',' '))
        axes[0,1].plot(v[:,0],v[:,7],color=color)
    t=data['glued_12000'][:,0]
    axes[0,0].plot(t,m.transmitted(0,t),':',color='#ec91b5',label='isolated storyboard readout')
    axes[0,0].set(title='The glued waist receives less drive',ylabel='waist field / reference pole peak');axes[0,0].legend(fontsize=9)
    axes[0,1].axhline(m.neck.level_radius,color='#ec91b5',linestyle='--',label='half-density contact threshold')
    axes[0,1].set(title='One neck: #338 radius and intrinsic metric',ylabel='minimum radius / rim radius');axes[0,1].legend(fontsize=8)
    v=data['glued_12000'];E=v[0,2]
    axes[1,0].plot(t,v[:,3]/E,color='#55d7cc',label='inside neck')
    axes[1,0].plot(t,v[:,4]/E,color='#efb96e',label='sphere side')
    axes[1,0].set(title='Energy enters the throat; it is not duplicated',ylabel='fraction of initial field energy');axes[1,0].legend(fontsize=9)
    axes[1,1].plot(t,v[:,5]/E,color='#55d7cc',label='glued surface')
    bare=data['bare_12000'];axes[1,1].plot(bare[:,0],bare[:,5]/bare[0,2],color='#ec91b5',label='bare sphere')
    axes[1,1].set(title='The throat provides a shortcut to B',ylabel='energy fraction near B');axes[1,1].legend(fontsize=9)
    for ax in axes.flat:ax.set_xlabel('t');ax.grid(alpha=.15);ax.set_facecolor('#101e2c')
    fig.suptitle('PR #344 check: the wave enters the throat, but default contact does not occur\nFixed glued geometry; gain 0.18 unchanged; no cut or release imposed',fontsize=14)
    fig.savefig(path,dpi=120);plt.close(fig)

if __name__=='__main__':main()
