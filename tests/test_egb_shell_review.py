"""Checks added after independent review; the original producer/archive stays frozen."""
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
import pytest

from geometrodynamics.bulk.egb_shell_dynamics import RoundEGBShell
from geometrodynamics.bulk import negative_egb as ne


def test_exact_bracket_identity_and_ordering():
    # Exact rational checks of the general algebra proved in the report.
    for x in [F(1,100), F(1,10), F(1,2), F(9,10), F(99,100)]:
        alpha_wec=-3*x/(4*(2+x))
        alpha_D=-x/(8-4*x)
        assert F(-1,4)<alpha_wec<alpha_D<0
        assert alpha_D-alpha_wec==x*(1-x)/(4-x*x)
        assert alpha_wec+F(1,4)==(1-x)/(2*(2+x))
        for alpha in [F(-3,10), F(-1,4), alpha_wec, alpha_D, F(0)]:
            P=1+4*alpha/3+8*alpha/(3*x)
            C=1+4*alpha
            E=1-4*alpha+8*alpha/x
            assert E==3*P-2*C
        E_at_D=1-4*alpha_D+8*alpha_D/x
        assert E_at_D==0


def test_negative_D_for_varied_angles_and_scales():
    for R in [1.,2.7]:
        for a in [.03,.3,.8,1.2,1.5]:
            original=RoundEGBShell(radius=R,mouth_angle=a)
            low,high=-R*R/4,original.static_bounds()['alpha_WEC']
            for fraction in [.01,.25,.5,.99]:
                alpha=low+(high-low)*fraction
                m=RoundEGBShell(radius=R,mouth_angle=a,alpha=alpha,eta=0.)
                r=m.r0;P=1+4*alpha/(3*R*R)+8*alpha/(3*r*r)
                linear=m.linear_data()
                assert P<0 and m.bulk_C>0
                assert linear['E']==pytest.approx(3*P-2*m.bulk_C,rel=1e-11)
                assert linear['D']<0 and linear['omega2']>0
            critical=-original.r0**2*R**2/(8*R**2-4*original.r0**2)
            at_D=RoundEGBShell(radius=R,mouth_angle=a,alpha=critical)
            assert at_D.acceleration_coefficient(at_D.r0,0.)==pytest.approx(0.,abs=1e-12)
            assert at_D.equilibrium[0]<0  # Degeneracy is outside shell WEC.
    reference=RoundEGBShell()
    critical=-reference.r0**2/(8-4*reference.r0**2)
    assert critical==pytest.approx(-.011414971304056791)


def test_bulk_stress_projection_and_bidirectional_flux():
    m=RoundEGBShell();r=.31
    rho=ne.exterior_density(m.alpha,m.radius)
    pressure=ne.exterior_pressure(m.alpha,m.radius)
    assert pressure/rho==pytest.approx(-.5681818181818182)
    for v in [-.12,.12]:
        f,_,delta,_=m.geometry(r,v)
        u=np.array([delta/math.sqrt(f),v/math.sqrt(f)])
        normal=np.array([u[1],u[0]])
        assert -u[0]**2+u[1]**2==pytest.approx(-1.)
        assert -normal[0]**2+normal[1]**2==pytest.approx(1.)
        assert -u[0]*normal[0]+u[1]*normal[1]==pytest.approx(0.)
        flux=2*u@np.diag([rho,pressure])@normal
        assert m.flux(r,v)==pytest.approx(flux,rel=1e-14)
        assert flux*v>0
        assert m.flux(r,-v)==pytest.approx(-flux,rel=1e-14)


def test_critical_exterior_does_not_remove_all_shell_inertia():
    r,v=.31,.12
    critical=RoundEGBShell(alpha=-.25)
    _,_,delta,_=critical.geometry(r,v)
    assert critical.bulk_C==0
    assert critical.flux(r,v)==0
    assert critical.acceleration_coefficient(r,v)==pytest.approx(16*critical.alpha*delta/r**2)
    assert critical.acceleration_coefficient(r,v)<0
    ratios=[]
    for C in [.76,.1,.01]:
        m=RoundEGBShell(alpha=(C-1)/4)
        ratios.append(m.flux(r,v)/m.bulk_C)
    assert max(ratios)-min(ratios)<1e-12


def test_archived_oscillation_requires_both_exchange_signs():
    root=Path(__file__).resolve().parents[1]
    path=root/'experiments/closure_ledger/runs/20261012_egb_shell_dynamics/samples.npz'
    m=RoundEGBShell()
    with np.load(path,allow_pickle=False) as archive:
        r,v=archive['stable_y'][:2]
        flux=np.array([m.flux(x,z) for x,z in zip(r,v)])
        assert flux.min()<0<flux.max()
        assert np.all(flux*v>=0)
        # This is the original scan margin, not a new retrospective admission cut.
        assert (1+4*archive['alpha_over_R2']).min()==pytest.approx(.0004)
