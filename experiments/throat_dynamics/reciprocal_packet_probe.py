"""Produce 2D packet, leading Einstein mass and null-junction compatibility data."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import simpson
from geometrodynamics.tangherlini.characteristic_packet import solve_pair


def summarize(run,N,ell):
    row=dict(N=N,ell=ell,
             incoming_symplectic=run['incoming_symplectic'],
             outgoing_symplectic=run['outgoing_symplectic'],
             observer_peak_abs_u=float(np.max(np.abs(run['observer_u']))))
    if ell==0:
        seam=run['seam']
        row['seam_end']={key:float(val[-1]) for key,val in seam.items()}
        for side in ('incoming_mass','outgoing_mass'):
            row[side]={key:val for key,val in run[side].items()
                       if key not in ('delta_m','other_path')}
        omega3=2*np.pi**2
        V,U=run['V'],run['out_U']
        a=run['field'][-1,:,0]
        ap=np.gradient(a,V,edge_order=2)
        bp=np.gradient(run['outgoing'][:,0,0],U,edge_order=2)
        Ein=omega3*simpson(V*ap**2,x=V)
        Eout=omega3*simpson(-U*bp**2,x=U)
        row['seam_killing_energy']=dict(incoming=float(Ein),outgoing=float(Eout),
                                       relative_difference=float(abs(Eout/Ein-1)))
        Eobs=float(simpson(run['observer_outward_killing_flux'],x=run['observer_t']))
        dM=float(run['observer_delta_M'][-1]-run['observer_delta_M'][0])
        row['observer_window']=dict(t_start=float(run['observer_t'][0]),
                                   t_end=float(run['observer_t'][-1]),
                                   outward_killing_energy=Eobs,delta_M=dM,
                                   relative_balance=abs(Eobs+dM)/abs(Eobs))
    return row


def run(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    rows=[];convergence=[];finest={}
    for ell in (0,1):
        previous=None
        for N in (64,128,256,512):
            current=solve_pair(N,ell)
            rows.append(summarize(current,N,ell))
            if previous is not None:
                convergence.append(dict(ell=ell,N=N,**{
                    key:float(np.linalg.norm(current[key][::2,::2]-previous[key])/
                              np.linalg.norm(current[key][::2,::2]))
                    for key in ('field','outgoing')}))
            if N==64:
                checkpoint={key:current[key] for key in ('U','V','field','out_U','out_V','outgoing')}
                if ell==0:
                    checkpoint.update(incoming_delta_m=current['incoming_mass']['delta_m'],
                                      outgoing_delta_m=current['outgoing_mass']['delta_m'])
                np.savez_compressed(output/f'fields_N64_l{ell}.npz',**checkpoint)
            previous=current
        finest[ell]=current
    fine=finest[0]
    base=solve_pair(128)
    boosts=[]
    for L in (.5,1.,2.,4.):
        r=solve_pair(128,boost=L)
        boosts.append(dict(boost=L,
            observer_relative_L2=float(np.linalg.norm(r['observer_u']-base['observer_u'])/np.linalg.norm(base['observer_u'])),
            clock_max_difference=float(np.max(np.abs(r['observer_t']-base['observer_t']))),
            mass_relative_L2=float(np.linalg.norm(r['outgoing_mass']['delta_m']-base['outgoing_mass']['delta_m'])/
                                   np.linalg.norm(base['outgoing_mass']['delta_m']))))
    scaling=[]
    for factor in (.5,2.):
        r=solve_pair(128,amplitude=.01*factor)
        scaling.append(dict(amplitude_factor=factor,
            field_relative_error=float(np.linalg.norm(r['outgoing']/factor-base['outgoing'])/np.linalg.norm(base['outgoing'])),
            mass_relative_error=float(np.linalg.norm(r['incoming_mass']['delta_m']/factor**2-base['incoming_mass']['delta_m'])/
                                      np.linalg.norm(base['incoming_mass']['delta_m']))))
    orders=[]
    for ell in (0,1):
        seq=[x for x in convergence if x['ell']==ell]
        for a,b in zip(seq[:-1],seq[1:]):
            orders.append(dict(ell=ell,N=b['N'],**{key:float(np.log2(a[key]/b[key])) for key in ('field','outgoing')}))
    finest_rows=[x for x in rows if x['N']==512]
    numerical=dict(
        field_second_order=all(1.7<r[key]<2.5 for r in orders for key in ('field','outgoing')),
        mass_path_conservation=all(finest_rows[0][side]['relative_path_discrepancy']<1e-3 for side in ('incoming_mass','outgoing_mass')),
        symplectic_conservation=all(r[side]['relative_residual']<1e-4 for r in finest_rows for side in ('incoming_symplectic','outgoing_symplectic')),
        passive_boost_covariance=max(x['observer_relative_L2'] for x in boosts)<1e-10 and max(x['mass_relative_L2'] for x in boosts)<1e-10,
        perturbative_amplitude_orders=max(x['mass_relative_error'] for x in scaling)<1e-10 and max(x['field_relative_error'] for x in scaling)<1e-10)
    seq=[r['seam_end']['induced_radius_mismatch'] for r in rows if r['ell']==0]
    mismatch_resolved=abs(seq[-1])>100*abs(seq[-1]-seq[-2]) and seq[-1]<0
    label=('RECIPROCAL_FIXED_NULL_MATCHING_FAILS_AT_ORDER_A2' if mismatch_resolved
           else 'JUNCTION_RESPONSE_UNRESOLVED') if all(numerical.values()) else 'NUMERICAL_VALIDATION_FAILED'
    columns=[fine['V'],fine['field'][-1,:,0],fine['observer_t'],fine['observer_u'],
             fine['observer_delta_M'],fine['observer_outward_killing_flux']]
    names=['incoming_V','horizon_u','observer_Killing_t','observer_u','observer_delta_M','observer_outward_Killing_flux']
    for key in ('incoming_delta_r','outgoing_delta_r','induced_radius_mismatch','incoming_delta_m','outgoing_delta_m'):
        names.append(key);columns.append(fine['seam'][key])
    columns.append(finest[1]['observer_u']);names.append('observer_u_l1_test_field')
    np.savetxt(output/'waveforms_N512.csv',np.column_stack(columns),delimiter=',',header=','.join(names),comments='')
    root=Path(__file__).resolve().parents[2]
    paths=['geometrodynamics/tangherlini/characteristic_packet.py',
           'experiments/throat_dynamics/reciprocal_packet_probe.py',
           'tests/test_reciprocal_characteristic_packet.py',
           'docs/reciprocal_packet_protocol.md','docs/reciprocal_packet_evolution.md']
    result=dict(label=label,numerical_gates=numerical,rows=rows,convergence_differences=convergence,
                measured_orders=orders,boost_controls=boosts,amplitude_controls=scaling,
                mismatch_resolved=mismatch_resolved,
                scope=dict(field='2D linear characteristic evolution, l=0 and l=1',
                           gravity='O(A^2) spherical Einstein mass and null-area response, l=0 only',
                           seam='prescribed reciprocal field transfer; perturbed induced metric fails matching',
                           checkpoint='full 2D fields at N=64; finest N=512 observer/seam arrays in CSV',
                           uncomputed=['nonlinear coupled metric evolution','surface material response',
                                       'finite-thickness focus resolution','global causal embedding','quantum event']),
                parameters=dict(b=1.,c=1.,kappa5=1.,A=.01,observer_r=2.,Omega3=2*np.pi**2),
                sources={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths},
                artifacts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir()) if p.suffix in ('.npz','.csv')},
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:result[key] for key in ('label','numerical_gates','measured_orders')},indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('experiments/throat_dynamics/runs/20261014_reciprocal_packet'))
    args=parser.parse_args();run(args.output)
