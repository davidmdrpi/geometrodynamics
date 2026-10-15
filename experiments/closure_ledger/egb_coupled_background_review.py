"""Post-review controls: shell-free support and the odd small-mouth limit.

Preserves the original #333 producer/protocol/archive. Numbers attributed to
reviewer exploration are independently recomputed here under stated norms.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import eval_gegenbauer
from geometrodynamics.bulk.egb_coupled_radial import CoupledRadial
from geometrodynamics.bulk.egb_shell_action import ScalarShellAction
from geometrodynamics.bulk.egb_shell_dynamics import RoundEGBShell


def no_shell_mode(ell, C=.76, R=1.):
    if not isinstance(ell,int) or ell<0 or C<=0 or R<=0:
        raise ValueError('require integer ell>=0, C>0, R>0')
    spatial=ell*(ell+3);h=spatial-4
    # Full-sphere weight: Gegenbauer orthogonality, alpha=3/2.
    squared_norm=(ell+1)*(ell+2)/(ell+1.5)
    f=lambda u:eval_gegenbauer(ell,1.5,math.cos(u))/math.sqrt(squared_norm)
    fp=lambda u:0. if ell==0 else -3*math.sin(u)*eval_gegenbauer(ell-1,2.5,math.cos(u))/math.sqrt(squared_norm)
    I=quad(lambda u:math.sin(u)**3*f(u)**2,0.,math.pi,epsabs=1e-12)[0]
    direct=3*C*R**4*quad(lambda u:math.sin(u)**3*(fp(u)**2-4*f(u)**2),0.,math.pi,epsabs=1e-12)[0]
    lam=spatial-4-2/C
    return dict(ell=ell,omega2_R2=lam,h=h,gauge=(ell==1),growing_physical=(ell!=1 and lam<0),
                full_sphere_L2=I,kinetic_norm=3*C*R**4*h,kinetic_norm_quadrature=direct,
                growth_rate_R=math.sqrt(-lam) if ell!=1 and lam<0 else None)


def odd_mouth_control(angle):
    m=CoupledRadial(ScalarShellAction(RoundEGBShell(mouth_angle=angle)))
    lam=m.root((-m.e+1e-8,0.),'odd',tolerance=2e-12)
    f,_=m.mode(lam,'odd',tolerance=2e-12);h=lam+m.e
    # rho0=3(C+1)/R^2; background rho is constant, so delta rho is gauge invariant.
    ratio=m.C/(m.C+1)
    measured=math.sqrt(quad(lambda u:math.sin(u)**3*(ratio*h*f(u)[0])**2,angle,m.mid,epsabs=1e-13)[0])
    return dict(a=angle,omega2_R2=lam,h=h,growth_rate_R=math.sqrt(-lam),
                even_lowest_omega2_R2=m.root((-6.8,-6.),tolerance=2e-12),
                norm_half_bulk_L2_unit=6*m.C*m.R**4*h,
                density_contrast_L2_half_bulk=ratio*h,density_contrast_quadrature=measured,
                norm_gradient_quadrature=m.physical_norm(lam,'odd'))


def run():
    root=Path(__file__).resolve().parents[2]
    old_path='experiments/closure_ledger/runs/20261013_egb_coupled_radial/result.json'
    old=json.loads((root/old_path).read_text())
    baseline=[no_shell_mode(l) for l in range(7)]
    mouths=[odd_mouth_control(a) for a in [.3,.2,.1,.05,.02]]
    old_modes=old['even_modes']+old['odd_modes_in_scanned_range']
    growing=sum(v['omega2_R2']<0 for v in old_modes)
    background_growing=sum(v['growing_physical'] for v in baseline)
    negative_oscillations=sum(v['omega2_R2']>0 and v['physical_norm_no_S3_over_kappa']<0 for v in old_modes)
    sources=['experiments/closure_ledger/egb_coupled_background_review.py','geometrodynamics/bulk/egb_coupled_radial.py']
    gates=dict(no_shell_norm=max(abs(v['kinetic_norm']-v['kinetic_norm_quadrature']) for v in baseline)<1e-9,
               no_shell_L2=max(abs(v['full_sphere_L2']-1) for v in baseline)<1e-12,
               odd_density=max(abs(v['density_contrast_L2_half_bulk']-v['density_contrast_quadrature']) for v in mouths)<1e-10,
               odd_norm=max(abs(v['norm_half_bulk_L2_unit']-v['norm_gradient_quadrature']) for v in mouths)<2e-8,
               positive_physical_odd_h=all(v['h']>0 for v in mouths),
               odd_density_decreases=all(mouths[i+1]['density_contrast_L2_half_bulk']<mouths[i]['density_contrast_L2_half_bulk'] for i in range(4)))
    shell_failures=[]
    if old['physical_negative_kinetic_directions_even']>1:shell_failures.append('one additional negative even kinetic direction')
    if negative_oscillations:shell_failures.append('real-frequency negative-norm shell oscillation')
    if growing>background_growing:shell_failures.append('two additional growing modes in the computed benchmark spectrum')
    return dict(label='BACKGROUND_AND_SHELL_FAILURES_SEPARATED' if all(gates.values()) else 'NUMERICALLY_UNRESOLVED',
                original_label=old['label'],original_archive_sha256=hashlib.sha256((root/old_path).read_bytes()).hexdigest(),
                original_archive=old_path,baseline=baseline,odd_mouth_controls=mouths,numerical_gates=gates,
                normalization=dict(no_shell='unit full-sphere weighted L2; G=3 C R^4 integral_0^pi sin^3(u)(fprime^2-4f^2)',
                                   odd='unit half-bulk weighted L2, as in original #333',density='L2 of delta rho/rho0 per unit mode amplitude'),
                comparison=dict(no_shell_negative_even_directions=1,with_shell_negative_even_directions=old['physical_negative_kinetic_directions_even'],
                                no_shell_growing_physical=background_growing,with_shell_growing_found=growing,
                                no_shell_negative_norm_oscillations=0,with_shell_negative_norm_oscillations_found=negative_oscillations),
                failure_reasons=dict(background=['homogeneous Eddington growth for every C>0','one closed-universe conformal negative kinetic direction; baseline signature, not a new focus defect'],shell=shell_failures),
                qualification='growing-mode counts compare the benchmark C=.76; higher no-shell harmonics can grow for smaller C; finite with-shell count is not an infinite-spectrum certificate',
                small_mouth_scope='fixed R=1, alpha=-.06, eta=.1; each a has its own equilibrium-anchored surface action, not evolution of one fixed surface material',
                sources={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources},
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path('experiments/closure_ledger/runs/20261014_egb_coupled_review/result.json'))
    args=parser.parse_args();result=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ['label','comparison','numerical_gates','odd_mouth_controls']},indent=2))
