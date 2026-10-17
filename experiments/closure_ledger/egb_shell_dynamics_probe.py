"""Execute the frozen EGB shell protocol; retain failed gates, never tune them."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from geometrodynamics.bulk.egb_shell_dynamics import RoundEGBShell, junction_tensor, casimir_scale

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'experiments/closure_ledger/runs/20261012_egb_shell_dynamics'


def derivative(fun, x, step):
    return (-fun(x+2*step)+8*fun(x+step)-8*fun(x-step)+fun(x-2*step))/(12*step)


def independent_checks():
    rng = np.random.default_rng(330)
    tensor_errors, flux_errors = [], []
    for _ in range(32):
        m = RoundEGBShell(alpha=rng.uniform(-.24,-.035))
        r,v,acc = rng.uniform(.15,.75), rng.uniform(-.4,.4), rng.uniform(-.8,.8)
        f,fp,delta,q = m.geometry(r,v)
        T = junction_tensor((acc-fp*v*v/(2*f))/delta,q,acc/r,(1+v*v)/r**2,m.alpha)
        reference = np.diag([m.density(r,v)]+[m.pressure(r,v,acc)]*3)
        tensor_errors.append(float(np.max(abs(T-reference))/(1+np.max(abs(reference)))))
        sr = derivative(lambda x:m.density(x,v),r,1e-5)
        sv = derivative(lambda x:m.density(r,x),v,1e-5)
        lhs = sr*v+sv*acc+3*v/r*(m.density(r,v)+m.pressure(r,v,acc))
        flux_errors.append(float(abs(lhs-m.flux(r,v))/(1+abs(lhs))))
    m = RoundEGBShell()
    return dict(max_tensor_relative_error=max(tensor_errors),
                max_flux_normalized_error=max(flux_errors),
                omitted_flux_residual_at_r031_v012=float(m.flux(.31,.12)))


def scan():
    bounds = RoundEGBShell().static_bounds()
    alphas = np.linspace(bounds['alpha_lower']+1e-4,bounds['alpha_WEC']-1e-4,101)
    etas = np.linspace(0,1,101)
    omega2, stable, dec, jacobian_error = [], [], [], []
    for alpha in alphas:
        row, signs, dom = [], [], []
        for eta in etas:
            m = RoundEGBShell(alpha=float(alpha),eta=float(eta))
            w2 = m.linear_data()['omega2']
            jac = derivative(lambda x:m.acceleration(x,0.),m.r0,1e-5*m.r0)
            jacobian_error.append(abs(w2+jac)/(1+abs(w2)))
            s,p = m.equilibrium
            row.append(w2); signs.append(w2>0 and s>=0 and s+p>=0)
            dom.append(s>=abs(p))
        omega2.append(row); stable.append(signs); dec.append(dom)
    stable, dec = np.array(stable), np.array(dec)
    return dict(alpha_over_R2=alphas,eta=etas,omega2_R2=np.array(omega2),
                stable_WEC=stable,DEC=dec), dict(
                    points=int(stable.size),stable_WEC_points=int(stable.sum()),
                    stable_DEC_points=int((stable&dec).sum()),
                    max_jacobian_normalized_error=float(max(jacobian_error)))


def trajectory_metrics(m,t,y,fraction):
    r,v,work,exchange = y
    s = np.array([m.density(x,z) for x,z in zip(r,v)])
    p = np.array([m.material_pressure(x,z) for x,z in zip(r,v)])
    D = np.array([m.acceleration_coefficient(x,z) for x,z in zip(r,v)])
    energy = s*r**3
    balance = energy-energy[0]+work-exchange
    w2 = m.linear_data()['omega2']
    out = dict(fraction=fraction,duration_over_R=float(t[-1]/m.radius),
               minimum_density_R=float(s.min()*m.radius),
               minimum_NEC_R=float((s+p).min()*m.radius),
               minimum_DEC_margin_R=float((s-abs(p)).min()*m.radius),
               minimum_abs_D=float(abs(D).min()),D_remains_negative=bool(np.all(D<0)),
               minimum_r_over_R=float(r.min()/m.radius),maximum_r_over_R=float(r.max()/m.radius),
               maximum_work_balance_error=float(np.max(abs(balance))/(abs(energy[0])+1e-30)),
               maximum_abs_work=float(np.max(abs(work))),maximum_abs_flux_integral=float(np.max(abs(exchange))))
    if w2>0:
        omega=np.sqrt(w2)
        norm = np.hypot((r-m.r0)/m.r0,v/(m.r0*omega))
        # Downward crossings of the equilibrium, linearly interpolated.
        k = np.flatnonzero((r[:-1]>=m.r0)&(r[1:]<m.r0))
        crossings = t[k]+(t[k+1]-t[k])*(r[k]-m.r0)/(r[k]-r[k+1])
        if len(crossings)<3:
            raise RuntimeError('insufficient crossings')
        period = float(np.mean(np.diff(crossings)))
        predicted = 2*np.pi/omega
        out.update(maximum_coordinate_norm_over_initial=float(norm.max()/abs(fraction)),
                   period_over_R=period/m.radius,linear_period_over_R=float(predicted/m.radius),
                   relative_period_error=abs(period/predicted-1),
                   relative_frequency_shift=float(predicted/period-1),
                   phase_drift_radians=float((2*np.pi/period-omega)*t[-1]),
                   crossing_count=len(crossings))
    else:
        predicted = fraction*m.r0*np.cosh(np.sqrt(-w2)*t)
        out.update(final_displacement_growth=float((r[-1]-m.r0)/(fraction*m.r0)),
                   maximum_relative_cosh_error=float(np.max(abs((r-m.r0)-predicted)/predicted)))
    return out


def measure():
    checks = independent_checks()
    samples, grid_report = scan()
    m = RoundEGBShell()
    duration = 20*2*np.pi/np.sqrt(m.linear_data()['omega2'])
    records=[]; trajectories=[]
    for name,fraction,tolerance in [('stable',1e-3,2e-10),('small',1e-4,2e-10),('replay',1e-3,2e-12)]:
        t,y=m.evolve(fraction,duration,rtol=tolerance)
        samples[name+'_t']=t; samples[name+'_y']=y
        trajectories.append(y)
        records.append(trajectory_metrics(m,t,y,fraction))
    y,z=trajectories[0],trajectories[2]
    replay_error=float(np.max(np.hypot((y[0]-z[0])/m.r0,(y[1]-z[1])/(m.r0*np.sqrt(m.linear_data()['omega2'])))))
    unstable=RoundEGBShell(eta=.3)
    t,y=unstable.evolve(1e-5,5/np.sqrt(-unstable.linear_data()['omega2']),samples=1001)
    samples['unstable_t']=t; samples['unstable_y']=y
    unstable_record=trajectory_metrics(unstable,t,y,1e-5)
    gates=dict(tensor_contraction=checks['max_tensor_relative_error']<1e-10,
               flux_identity=checks['max_flux_normalized_error']<1e-7,
               source_free_countercontrol_fails=abs(checks['omitted_flux_residual_at_r031_v012'])>1e-3,
               acceleration_jacobian=grid_report['max_jacobian_normalized_error']<1e-6,
               stable_norm=all(x['maximum_coordinate_norm_over_initial']<1.1 for x in records),
               period_agreement=all(x['relative_period_error']<1e-3 for x in records),
               amplitude_convergence=records[1]['relative_period_error']<records[0]['relative_period_error'],
               integrator_replay=replay_error<1e-6,
               unstable_control=(unstable_record['final_displacement_growth']>50 and unstable_record['maximum_relative_cosh_error']<.02),
               work_balance=all(x['maximum_work_balance_error']<1e-8 for x in records+[unstable_record]),
               stable_energy_conditions=all(x['minimum_density_R']>0 and x['minimum_NEC_R']>0 and x['minimum_DEC_margin_R']>0 for x in records))
    gates={key:bool(value) for key,value in gates.items()}
    if not all(gates.values()): label='NUMERICALLY_UNRESOLVED'
    elif grid_report['stable_WEC_points']==0: label='NO_RADIAL_WINDOW'
    else: label='OPEN_RADIAL_WINDOW_COUPLED_ADMISSION_OPEN'
    sources=['geometrodynamics/bulk/egb_shell_dynamics.py',
             'experiments/closure_ledger/egb_shell_dynamics_probe.py',
             'tests/test_egb_shell_dynamics.py','docs/egb_shell_dynamics_protocol.md']
    result=dict(label=label,numerical_gates=gates,independent_checks=checks,grid=grid_report,
                parameters=dict(R=1.,a=.3,alpha=-.06,eta=.1),linear=m.linear_data(),
                static_bounds=m.static_bounds(),equilibrium=dict(s=float(m.equilibrium[0]),p=float(m.equilibrium[1])),
                stable_trajectories=records,unstable_trajectory=unstable_record,replay_coordinate_error=replay_error,
                corrective_kick_budget=0.,
                admission=dict(radial_stability_window=grid_report['stable_WEC_points']>0,
                               exterior_TT_principal_part=m.bulk_C>0,
                               canonical_kinetic_norm='UNCOMPUTED: D<0 is not a canonical energy calculation',
                               bulk_matter_action='UNCOMPUTED',coupled_characteristics='UNCOMPUTED',
                               nonspherical_stability='UNCOMPUTED',smooth_finite_thickness='UNCOMPUTED',
                               focal_activation='UNCOMPUTED'),
                casimir=dict(D=5,assumption='one-scale Einstein curvature balance; no magnetic or near-extremal enhancement',
                             rows=[dict(N=N,L_over_b=L,b_over_lP=casimir_scale(5,N,L)) for N in [1,1000,1000000] for L in [1,10,100]]),
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    return result,samples


def write_plot(path, result, samples):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['svg.hashsalt']='egb-shell-dynamics'
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    ax=axes[0]
    alpha=samples['alpha_over_R2'];r=RoundEGBShell().r0
    threshold=-2/(3*(1-4*alpha+8*alpha/r**2))
    ax.fill_between(alpha,0,threshold,color='#c8dbe8',label='WEC + radial stability')
    cut=result['static_bounds']['alpha_DEC']
    ax.fill_between(alpha,0,threshold,where=alpha<=cut,color='#72b8a3',label='also DEC')
    ax.plot(alpha,threshold,color='#294e70',label='neutral boundary')
    ax.scatter([-.06,-.06],[.1,.3],c=['#17664e','#ac3d34'],zorder=5)
    ax.set(xlabel=r'$\alpha/R^2$',ylabel=r'$\eta=dp/d\sigma$',ylim=(0,.42),
           title='Restricted radial stability, a = 0.3')
    ax.legend(loc='upper left',fontsize=8)
    ax=axes[1]
    period=result['stable_trajectories'][0]['linear_period_over_R']
    t=samples['stable_t'];y=samples['stable_y']
    ax.plot(t/period,(y[0]-r)/(r*1e-3),color='#17664e',linewidth=1)
    ax.set(xlabel='Proper time / linear period',ylabel='Radial displacement / initial displacement',
           title='Nonlinear stable control, 20 periods',xlim=(0,20),ylim=(-1.15,1.15))
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Fixed supported exterior; canonical energy and coupled stability uncomputed',fontsize=10)
    fig.savefig(path,metadata={'Date':None})
    plt.close(fig)
    path=Path(path)
    if path.suffix=='.svg':
        path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=RUN)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    result,samples=measure()
    archive=args.output/'samples.npz'
    np.savez_compressed(archive,**samples)
    result['samples_archive']=dict(file=archive.name,sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                                  schema='grid arrays (alpha,eta); trajectory y rows: r, dr/dtau, pressure work, integrated exterior flux')
    write_plot(args.output/'radial_window.svg',result,samples)
    (args.output/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['label','numerical_gates','independent_checks','grid','linear','stable_trajectories','unstable_trajectory']},indent=2))
    if result['label']=='NUMERICALLY_UNRESOLVED': raise SystemExit(1)
