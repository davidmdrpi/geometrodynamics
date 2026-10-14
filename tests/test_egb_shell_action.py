"""Variational material checks and physical-tangent rejection, not a ghost claim."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.bulk.egb_shell_action import ScalarShellAction,second_derivative
from geometrodynamics.bulk.egb_shell_dynamics import RoundEGBShell
from geometrodynamics.bulk import negative_egb as ne
from experiments.closure_ledger.egb_shell_dynamics_probe import derivative


def test_static_action_matching_and_exact_affine_eos():
    for eta in [.05,.1,.3,1.]:
        a=ScalarShellAction(shell=RoundEGBShell(eta=eta));m=a.shell
        assert a.bulk_stress(a.X0)==pytest.approx((ne.exterior_density(m.alpha,m.radius),ne.exterior_pressure(m.alpha,m.radius)))
        assert a.surface_stress(a.X0)==pytest.approx(m.equilibrium)
        for ratio in [.8,1.,1.2]:
            s,p=a.surface_stress(a.X0*ratio);s0,p0=m.equilibrium
            assert p-p0==pytest.approx(eta*(s-s0),abs=1e-12)
        assert a.frozen_current_residual(m.r0,0.,0.)==0


def test_material_hessians_from_lagrangian_differencing():
    for lam in [-1.,0.,1.]:
        a=ScalarShellAction(bulk_lambda=lam);q=math.sqrt(a.q2);c=a.bare_coefficients()
        for kind in ['bulk','surface']:
            ht=second_derivative(lambda z:a.material_lagrangian(kind,z),q,5e-4*q)
            hs=-second_derivative(lambda z:a.material_lagrangian(kind,q,z),0.,5e-4*q)
            assert ht==pytest.approx(c[kind+'_time'],rel=1e-5,abs=1e-6)
            assert hs==pytest.approx(c[kind+'_space'],rel=1e-5,abs=1e-6)
        assert c['surface_sound2']==pytest.approx(a.shell.eta)


def test_identical_static_data_do_not_fix_bulk_kinetic_sign():
    good=ScalarShellAction();bad=ScalarShellAction(bulk_lambda=-1.)
    assert good.bulk_stress(good.X0)==bad.bulk_stress(bad.X0)
    assert good.bulk_PX(good.X0)==bad.bulk_PX(bad.X0)==1
    assert good.surface_stress(good.X0)==bad.surface_stress(bad.X0)
    assert good.bare_coefficients()['bulk_time']==1
    assert bad.bare_coefficients()['bulk_time']==-1
    # The elliptic/ghost control is rejected, not supplied as another healthy model.
    assert bad.bare_coefficients()['bulk_sound2'] is None


def test_shift_stationary_local_energy_has_the_same_bare_kinetic_hessian():
    for lam in [-1.,0.,1.]:
        a=ScalarShellAction(bulk_lambda=lam);q=math.sqrt(a.q2)
        for kind in ['bulk','surface']:
            energy=lambda z:a.local_shift_stationary_energy(kind,z)
            assert derivative(energy,q,1e-4*q)==pytest.approx(0.,abs=1e-8)
            hessian=second_derivative(energy,q,5e-4*q)
            assert hessian==pytest.approx(a.bare_coefficients()[kind+'_time'],rel=1e-5,abs=1e-6)


def test_linear_density_constraint_excludes_frozen_shell_motion():
    for alpha,eta in [(-.06,.1),(-.1,.05),(-.2,.01)]:
        a=ScalarShellAction(shell=RoundEGBShell(alpha=alpha,eta=eta));m=a.shell;r=m.r0
        slope=a.linear_obstruction()['geometric_density_slope']
        material_slope=derivative(lambda x:a.surface_stress(a.induced_Y(x,0.))[0],r,1e-4*r)
        geometric_slope=derivative(lambda x:m.density(x,0.),r,1e-4*r)
        assert material_slope==pytest.approx(0,abs=1e-7)
        assert geometric_slope==pytest.approx(slope,rel=1e-8)
        assert slope<0
        eps=1e-6
        assert a.frozen_density_residual(r*(1+eps),0)/(slope*r*eps)==pytest.approx(1.,rel=1e-5)


def test_global_shift_charge_and_local_scalar_boundary_agree():
    a=ScalarShellAction();r=.31;v=.12;acc=-.2
    dQdr=derivative(lambda x:a.frozen_total_charge(x,v),r,1e-5)
    dQdv=derivative(lambda z:a.frozen_total_charge(r,z),v,1e-5)
    assert dQdr*v+dQdv*acc==pytest.approx(r**3*a.frozen_current_residual(r,v,acc),rel=1e-8)
    m=a.shell;r=m.r0;linear=a.linear_obstruction()
    dQ=derivative(a.frozen_total_charge,r,1e-5)
    assert dQ==pytest.approx(linear['total_charge_slope_no_S3_area'],rel=1e-8)
    assert a.frozen_current_residual(r,1e-7,0.)/1e-7==pytest.approx(linear['current_residual_per_velocity'],rel=1e-9)


def test_scalar_ward_identity_off_shell_and_required_trace_response():
    a=ScalarShellAction();m=a.shell;r=.31;v=.12;acc=-.2
    Y=a.induced_Y(r,v);s,p=a.surface_stress(Y)
    dsdr=derivative(lambda x:a.surface_stress(a.induced_Y(x,v))[0],r,1e-5)
    dsdv=derivative(lambda z:a.surface_stress(a.induced_Y(r,z))[0],v,1e-5)
    lhs=dsdr*v+dsdv*acc+3*v/r*(s+p)-m.flux(r,v)
    chi_dot=math.sqrt(2*Y)
    assert lhs==pytest.approx(chi_dot*a.frozen_current_residual(r,v,acc),rel=1e-8)
    linear=a.linear_obstruction();s,p=m.equilibrium
    rate=linear['required_scalar_trace_rate_per_q_per_displacement']
    assert (s+p)/m.eta*rate==pytest.approx(linear['geometric_density_slope'])


def test_scaling_and_declared_domain():
    a=ScalarShellAction();k=2.3
    b=ScalarShellAction(shell=RoundEGBShell(radius=k,alpha=-.06*k*k))
    assert b.q2*k*k==pytest.approx(a.q2)
    assert b.potential*k*k==pytest.approx(a.potential)
    assert b.bare_coefficients()['surface_time']/k==pytest.approx(a.bare_coefficients()['surface_time'])
    assert b.linear_obstruction()['geometric_density_slope']*k*k==pytest.approx(a.linear_obstruction()['geometric_density_slope'])
    with pytest.raises(ValueError):ScalarShellAction(shell=RoundEGBShell(eta=0.))
    with pytest.raises(ValueError):ScalarShellAction(shell=RoundEGBShell(alpha=-.25))
    with pytest.raises(ValueError):a.surface_P(0.)


def test_archive_retains_failed_admission_and_provenance():
    root=Path(__file__).resolve().parents[1]
    report=json.loads((root/'experiments/closure_ledger/runs/20261013_egb_shell_kinetic/result.json').read_text())
    assert report['label']=='STATIC_ACTION_MATCH_FROZEN_MODE_REJECTED'
    assert all(report['numerical_gates'].values())
    assert report['physical_radial_kinetic_norm'] is None
    assert not report['admission']['original_frozen_radial_mode_satisfies_constraints']
    for name,digest in report['sources'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
    original=root/'experiments/closure_ledger/runs/20261012_egb_shell_dynamics/samples.npz'
    assert hashlib.sha256(original.read_bytes()).hexdigest()==report['original_trajectory_residual']['archive_sha256']
