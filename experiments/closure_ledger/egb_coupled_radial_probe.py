"""Reproduce the constrained radial spectrum, norms, and independent controls."""
import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
import numpy as np
import scipy
from geometrodynamics.bulk.egb_coupled_radial import CoupledRadial
from experiments.closure_ledger.egb_coupled_action_check import check_action_hessian

EVEN_BRACKETS=[(-6.5,-6.2),(-2.6,-2.4),(1.8,2.1),(6.2,6.6),(30.,32.),(67.,69.),(116.,118.),(177.,179.)]


def run():
    m=CoupledRadial();even=[]
    for bracket in EVEN_BRACKETS:
        value=m.root(bracket);tight=m.root(bracket,tolerance=2e-13)
        norm=m.physical_norm(value);profile,z=m.mode(value)
        gauges=[m.symplectic_norm(value,p) for p in (1,3)]
        junctions=m.junction_residuals(value)
        even.append(dict(omega2_R2=value,shooting_tight_difference=abs(value-tight),
                         growth_rate_R=math.sqrt(-value) if value<0 else None,
                         frequency_R=math.sqrt(value) if value>0 else None,
                         physical_norm_no_S3_over_kappa=norm,
                         norm_spectral=m.spectral_norm(value),norm_gauges=gauges,
                         shell_fractional_displacement=z,bulk_psi_at_mouth=float(profile(m.a)[0]),
                         charge_residual=m.charge_residual(value),junction_residuals=junctions))
    # Odd physical branch has zero trace at both mouth and equator. The h=0
    # coordinate mode is analyzed separately and is not in this determinant.
    grid=np.linspace(-m.e+.001,180,361);values=[m.determinant(v,'odd') for v in grid]
    odd=[]
    for lo,hi,fl,fh in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
        if fl*fh<0:
            value=m.root((lo,hi),'odd')
            odd.append(dict(omega2_R2=value,physical_norm_no_S3_over_kappa=m.physical_norm(value,'odd'),
                            growth_rate_R=math.sqrt(-value) if value<0 else None))
    fem=[m.finite_element_spectrum(n) for n in (80,160,320)]
    errors=[[abs(row['low_eigenvalues'][i]-even[i]['omega2_R2']) for i in range(4)] for row in fem]
    action=[check_action_hessian(eps) for eps in (1e-3,5e-4)]
    gauges_ok=all(abs(row['physical_norm_no_S3_over_kappa']-v['total'])/(1+abs(row['physical_norm_no_S3_over_kappa']))<1e-7 for row in even for v in row['norm_gauges'])
    gates=dict(two_gauge_norm_agreement=gauges_ok,
               spectral_norm_agreement=all(abs(v['norm_spectral']-v['physical_norm_no_S3_over_kappa'])/(1+abs(v['norm_spectral']))<1e-7 for v in even),
               tighter_shooting=max(v['shooting_tight_difference'] for v in even)<1e-7,
               junctions=max(abs(v) for row in even for v in row['junction_residuals'].values())<1e-7,
               total_charge=max(abs(v['charge_residual']) for v in even)<1e-7,
               action_hessian=action[-1]['relative_error']<2e-7,
               fem_refinement=all(errors[n+1][i]<.3*errors[n][i] for n in (0,1) for i in range(4)) and max(errors[-1])<4e-4,
               physical_inertia=all(v['physical_negative_directions']==2 for v in fem),
               finite_pencil_real=all(v['maximum_imaginary_part']<1e-8 for v in fem))
    failures=[]
    if any(v['omega2_R2']>0 and v['physical_norm_no_S3_over_kappa']<0 for v in even):
        failures.append('oscillatory shell mode has negative physical kinetic norm')
    if any(v['omega2_R2']<0 for v in even+odd):
        failures.append('negative squared frequencies in the computed radial spectrum')
    root=Path(__file__).resolve().parents[2]
    sources=['geometrodynamics/bulk/egb_coupled_radial.py','experiments/closure_ledger/egb_coupled_radial_probe.py','experiments/closure_ledger/egb_coupled_action_check.py','docs/egb_coupled_radial_protocol.md']
    result=dict(label=('COUPLED_RADIAL_ADMISSION_FAILED' if failures else 'NO_FAILURE_FOUND_IN_COMPUTED_MODES') if all(gates.values()) else 'NUMERICALLY_UNRESOLVED',
                parameters=dict(R=m.R,a=m.a,alpha=m.shell.alpha,eta=m.shell.eta,bulk_lambda=0),
                scope='linear spherical perturbations of the #331 material action; even and odd sectors',
                normalization=dict(even='z=delta r_shell/r0=1',odd='half-bulk integral sin(u)^3 psi^2 du=1',norm='coefficient of Qdot^2 in S2=(Omega3/2kappa5) integral(K Qdot^2-K omega^2 Q^2)dt'),
                constants=dict(C=m.C,D=m.D,e=m.e,x=m.x,d=m.d,b=m.b,j=m.j),
                auxiliary_mass_is_physical_norm=False,
                even_modes=even,odd_modes_in_scanned_range=odd,odd_scan_range=[float(grid[0]),float(grid[-1])],
                finite_elements=fem,finite_element_errors_first_four=errors,full_action_hessian=action,
                numerical_gates=gates,old_prescribed_exterior_omega2_R2=m.shell.linear_data()['omega2']*m.R**2,
                physical_negative_kinetic_directions_even=2,
                failure_reasons=failures,
                exceptional_modes=dict(constant_scalar_shift='quotiented at fixed total charge',odd_h_zero='cos(u), z=0, pure coordinate mode with zero symplectic norm; removed'),
                limitations=['no nonlinear evolution','no nonspherical stability claim','no finite-thickness regularity','no no-go for different material actions'],
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                sources={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources})
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path('experiments/closure_ledger/runs/20261013_egb_coupled_radial/result.json'))
    args=parser.parse_args();result=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ['label','numerical_gates','finite_element_errors_first_four']},indent=2))
    for row in result['even_modes'][:4]:print(row['omega2_R2'],row['physical_norm_no_S3_over_kappa'])
    print('odd low:',result['odd_modes_in_scanned_range'][:2])
