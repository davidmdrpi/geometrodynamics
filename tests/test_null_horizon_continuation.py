"""Independent equation and normalization checks, not a packet-evolution test."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
import sympy as s

from experiments.throat_dynamics.null_horizon_continuation_probe import (
    junction, kruskal_product, metric_factor, tortoise, trace_control, wave_coefficient,
)
from geometrodynamics.tangherlini.dynamics import master_potential


def test_metric_pullback_from_static_coordinates():
    t, r, b = s.symbols('t r b', positive=True)
    star = r+b/2*s.log((r-b)/(r+b))
    U, V = -b*s.exp(-(t-star)/b), b*s.exp((t+star)/b)
    F = (r+b)**2/r**2*s.exp(-2*r/b)
    f = 1-b*b/r**2
    assert s.simplify(-F*s.diff(U,t)*s.diff(V,t)+f)==0
    assert s.simplify(-F*s.diff(U,r)*s.diff(V,r)-1/f)==0
    assert s.simplify(s.diff(U,t)*s.diff(V,r)+s.diff(U,r)*s.diff(V,t))==0
    for radius in (.5,1.,2.):
        assert metric_factor(radius,radius)==pytest.approx(4/math.e**2)
        assert kruskal_product(radius,radius)==0


def test_scalar_operator_from_radial_divergence():
    r, b, L = s.symbols('r b L', positive=True)
    w = s.Function('w')(r)
    f = 1-b*b/r**2
    # Multiply Box(r^-3/2 w Y)=0 by f*r^3/2; time term is -w_tt.
    radial = s.expand(f*r**s.Rational(3,2)*(s.diff(r**3*f*s.diff(r**s.Rational(-3,2)*w,r),r)/r**3
                      -L/r**2*r**s.Rational(-3,2)*w))
    V = f*(L/r**2+3*f/(4*r**2)+3*s.diff(f,r)/(2*r))
    assert s.simplify(radial-(f*s.diff(f*s.diff(w,r),r)-V*w))==0
    for ell in (0,1,3):
        for radius in (.5,1.,2.):
            r0=1.8*radius
            assert wave_coefficient(r0,ell,radius)==pytest.approx(
                metric_factor(r0,radius)*master_potential(r0,ell,radius)/(4*(1-radius**2/r0**2)))
            assert wave_coefficient(radius,ell,radius)==pytest.approx(
                math.exp(-2)*(ell*(ell+2)+3)/radius**2)


def test_common_frame_jump_and_retained_side_condition():
    # Implicitly differentiate UV=H(r), independently of junction().
    r,b=s.symbols('r b',positive=True)
    H=-b*b*(r-b)/(r+b)*s.exp(2*r/b)
    assert s.simplify(1/s.diff(H,r).subs(r,b)+2/(b*s.E**2))==0
    for radius in (.5,1.,2.):
        for V in (radius,2*radius,3*radius):
            for alpha in (.4,1.,3.):
                linear=junction(V,-alpha*V,-alpha,b=radius)
                assert linear['jump_nr']==pytest.approx(0.,abs=1e-13)
                assert not linear['physical_stress_interpretation']
                # g' h=1 fixes h<0: output V has wrong sign for retained exterior.
                assert linear['continuous_metric_transverse_slope']<0
            c=2*radius**2
            control=junction(V,-c/V,c/V**2,-2*c/V**3,b=radius)
            assert control['physical_stress_interpretation']
            assert control['formal_mu_8piG']==pytest.approx(-6*V/radius**2)
            assert control['formal_p_8piG']==pytest.approx(2/V)


def test_linear_and_reciprocal_maps_transform_as_same_law():
    alpha,c=1.3,2.
    for L in (.5,1.,2.,4.):
        for V in (1.,2.,3.):
            assert -(alpha/L**2)*(L*V)==pytest.approx(-alpha*V/L)
            assert -c/(L*V)==pytest.approx((-c/V)/L)
    # Normal vector to UV=constant is (U,V), not (-1,1) at generic t.
    U,V=-.1,2.
    assert abs(U+V)>1
    # Fixed alpha at nonunit boost is a changed physical map.
    assert -alpha*(2*V)!=pytest.approx((-alpha*V)/2)


@pytest.mark.parametrize('kind,orientation',[('linear',-1),('reciprocal',1)])
def test_finite_energy_trace_covariance_and_symplectic_orientation(kind,orientation):
    base=trace_control(kind)
    for ell in (0,1):
        for L in (.5,1.,2.,4.):
            row=trace_control(kind,ell,L)
            assert row['killing_out']==pytest.approx(base['killing_in'],rel=1e-10)
            assert row['symplectic_out']==pytest.approx(orientation*base['symplectic_in'],rel=1e-10)
            assert row['raw_in']==pytest.approx(base['raw_in']/L,rel=1e-10)
            assert row['raw_out']==pytest.approx(base['raw_out']*L,rel=1e-10)
    assert base['symplectic_in']>0


def test_cutoff_delay_matches_integrated_radial_null_speed():
    from scipy.integrate import quad
    for eps in (.1,.01,.001,.0001):
        delay=2*(tortoise(2.)-tortoise(1.+eps))
        direct=2*quad(lambda r:1/(1-1/r**2),1.+eps,2.,epsabs=1e-10)[0]
        assert delay==pytest.approx(direct,rel=1e-10)
    d1=2*(tortoise(2.)-tortoise(1.00001))
    d2=2*(tortoise(2.)-tortoise(1.000001))
    assert d2-d1==pytest.approx(math.log(10),abs=2e-5)


def test_archive_preserves_protocol_and_does_not_claim_evolution():
    root=Path(__file__).resolve().parents[1]
    path=root/'experiments/throat_dynamics/runs/20261014_null_horizon_preflight/result.json'
    data=json.loads(path.read_text())
    assert data['label']=='SUBMITTED_PROTOCOL_REQUIRES_CORRECTION'
    assert all(data['checks'].values())
    assert 'finite-packet PDE evolution' in data['uncomputed']
    assert all(row['G3'].startswith('NOT_RUN') for row in data['laws'].values())
    for path,digest in data['sources'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest
