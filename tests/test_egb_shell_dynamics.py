"""Tensor, geometric, conservation and stability checks for the restricted shell."""
import math
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.bulk.egb_shell_dynamics import RoundEGBShell, junction_tensor, casimir_scale
from experiments.closure_ledger.egb_shell_dynamics_probe import independent_checks, derivative, trajectory_metrics


def test_full_junction_and_bulk_flux():
    data=independent_checks()
    assert data['max_tensor_relative_error']<1e-10
    assert data['max_flux_normalized_error']<1e-7
    assert data['omitted_flux_residual_at_r031_v012']>.5


def test_ultrastatic_embedding_not_schwarzschild_lapse():
    m=RoundEGBShell();r=.31;v=.12;acc=.23
    f,fp,delta,q=m.geometry(r,v)
    # In the flat (t,l) plane, u=(sqrt(1+w^2),w), n=(w,sqrt(1+w^2)).
    w=v/math.sqrt(f);wdot=acc/math.sqrt(f)-fp*v*v/(2*f**1.5)
    gamma=math.sqrt(1+w*w)
    normal_acceleration=-w*(w*wdot/gamma)+gamma*wdot
    k0=(acc-fp*v*v/(2*f))/delta
    assert normal_acceleration==pytest.approx(k0,rel=1e-14)
    assert abs(k0-(acc+fp/2)/delta)>.1
    assert q==pytest.approx(gamma*math.sqrt(f)/r)


def test_einstein_limit_and_static_shell():
    m=RoundEGBShell(alpha=0);r=.41;v=.21;acc=-.17
    f,fp,d,q=m.geometry(r,v);k0=(acc-fp*v*v/(2*f))/d
    assert m.density(r,v)==pytest.approx(-6*q)
    assert m.pressure(r,v,acc)==pytest.approx(2*k0+4*q)
    T=junction_tensor(k0,q,acc/r,(1+v*v)/r**2,0.)
    assert T[0,0]==pytest.approx(-6*q)
    assert T[1,1]==pytest.approx(2*k0+4*q)
    assert m.equilibrium[0]<0


def test_exact_energy_boundaries_and_nonempty_stability_window():
    m=RoundEGBShell();b=m.static_bounds()
    assert RoundEGBShell(alpha=b['alpha_WEC']).equilibrium[0]==pytest.approx(0,abs=1e-13)
    d=RoundEGBShell(alpha=b['alpha_DEC']);s,p=d.equilibrium
    assert s==pytest.approx(p)
    assert m.equilibrium[0]>abs(m.equilibrium[1])
    lin=m.linear_data()
    assert lin['omega2']>0 and lin['D']<0 and lin['bulk_C']>0
    assert RoundEGBShell(eta=.3).linear_data()['omega2']<0
    assert RoundEGBShell(eta=lin['eta_critical']).linear_data()['omega2']==pytest.approx(0,abs=1e-13)
    for alpha in np.linspace(-.249,b['alpha_WEC'],10):
        model=RoundEGBShell(alpha=alpha,eta=0)
        assert model.linear_data()['D']<0
        assert model.linear_data()['omega2']>0


def test_linearization_and_static_sourced_conservation():
    for alpha in [-.24,-.06,-.032]:
        for eta in [0.,.1,1.]:
            m=RoundEGBShell(alpha=alpha,eta=eta);r=m.r0
            jac=derivative(lambda x:m.acceleration(x,0),r,1e-5*r)
            assert jac==pytest.approx(-m.linear_data()['omega2'],rel=1e-8,abs=1e-8)
            s,p=m.equilibrium
            sr=derivative(lambda x:m.density(x,0),r,1e-5*r)
            expected=-3*(s+p)/r+6*m.bulk_C/(m.radius**2*math.cos(m.mouth_angle))
            assert sr==pytest.approx(expected,rel=1e-9)


def test_radius_scaling():
    m=RoundEGBShell();scale=3.7
    n=RoundEGBShell(radius=scale,alpha=m.alpha*scale**2)
    r,v,acc=.31,.12,.23
    assert n.density(scale*r,v)*scale==pytest.approx(m.density(r,v))
    assert n.pressure(scale*r,v,acc/scale)*scale==pytest.approx(m.pressure(r,v,acc))
    assert n.flux(scale*r,v)*scale**2==pytest.approx(m.flux(r,v))
    assert n.linear_data()['omega2']*scale**2==pytest.approx(m.linear_data()['omega2'])
    assert n.acceleration(scale*r,v)*scale==pytest.approx(m.acceleration(r,v))


def test_nonlinear_stable_work_ledger_and_unstable_control():
    m=RoundEGBShell();period=2*np.pi/np.sqrt(m.linear_data()['omega2'])
    t,y=m.evolve(1e-3,4*period,samples=1001)
    record=trajectory_metrics(m,t,y,1e-3)
    assert record['relative_period_error']<1e-3
    assert record['maximum_coordinate_norm_over_initial']<1.1
    assert record['maximum_work_balance_error']<1e-8
    assert record['maximum_abs_flux_integral']>1e-6
    assert record['minimum_DEC_margin_R']>0
    u=RoundEGBShell(eta=.3)
    t,y=u.evolve(1e-5,5/np.sqrt(-u.linear_data()['omega2']),samples=501)
    record=trajectory_metrics(u,t,y,1e-5)
    assert record['final_displacement_growth']>50
    assert record['maximum_relative_cosh_error']<.02


def test_degeneracy_and_domain_are_not_regularized():
    m=RoundEGBShell();r=m.r0
    critical=-1/(8/r**2-4/m.radius**2)
    with pytest.raises(ValueError,match='degenerate'):
        RoundEGBShell(alpha=critical).acceleration(r,0)
    for radius in [0.,1.,-1.,np.nan]:
        with pytest.raises(ValueError):m.density(radius,0)
    with pytest.raises(ValueError):RoundEGBShell(radius=0)


def test_casimir_scale_is_conditional_and_dimension_dependent():
    assert casimir_scale(5,1000)==pytest.approx(10.)
    assert casimir_scale(4,1000)==pytest.approx(math.sqrt(1000))
    assert casimir_scale(5,1,10)<casimir_scale(5,1,1)
    b=casimir_scale(5,1e6,10)
    assert b*(10*b)**2==pytest.approx(1e6)
    with pytest.raises(ValueError):casimir_scale(2,1)


def test_archive_provenance_and_recomputed_trajectory():
    root=Path(__file__).resolve().parents[1]
    run=root/'experiments/closure_ledger/runs/20261012_egb_shell_dynamics'
    report=json.loads((run/'result.json').read_text())
    assert report['label']=='OPEN_RADIAL_WINDOW_COUPLED_ADMISSION_OPEN'
    assert all(report['numerical_gates'].values())
    for path,digest in report['sources'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest
    archive=run/report['samples_archive']['file']
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==report['samples_archive']['sha256']
    with np.load(archive,allow_pickle=False) as data:
        assert int(data['stable_WEC'].sum())==report['grid']['stable_WEC_points']
        assert int((data['stable_WEC']&data['DEC']).sum())==report['grid']['stable_DEC_points']
        measured=trajectory_metrics(RoundEGBShell(),data['stable_t'],data['stable_y'],1e-3)
        for key in ['maximum_work_balance_error','relative_period_error','minimum_DEC_margin_R']:
            assert measured[key]==pytest.approx(report['stable_trajectories'][0][key],abs=1e-12)
