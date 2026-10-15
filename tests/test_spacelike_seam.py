"""Leading-order spacelike embedding and junction validation."""
import numpy as np
import pytest
import sympy as s
from geometrodynamics.tangherlini.characteristic_packet import solve_pair
from geometrodynamics.tangherlini.spacelike_seam import deflected_seam,field_derivative_source_control


def test_first_form_and_leading_conservation_symbolically():
    v,b,kappa=s.symbols('v b kappa',positive=True)
    I,B,C,D,source=s.symbols('I B C D S')
    L=s.log(v/b)*B-D;K=2/(b*s.E**2)
    def dv(e):
        return s.diff(e,v)+s.diff(e,I)*source+s.diff(e,B)*v*source+s.diff(e,C)*v**2*source+s.diff(e,D)*v*s.log(v/b)*source
    x=(-v*I+B+L)/(K*v);difference=-v*I+2*B-C/v
    y=v/b**2*(difference/K-v*x)
    assert s.simplify(dv(x)-b*b/v**2*dv(y))==0
    assert s.simplify(dv(x)+L/(K*v*v))==0
    H=2*b*L/v**2
    Sl=-6*v*s.sqrt(H)/(kappa*b*b)
    Sa=2/(kappa*v*s.sqrt(H))
    assert s.simplify(dv(Sl)-3*(-B/v)/b*Sa)==0
    # Leading normal-normal (Gauss) jump; O(A^2) terms are not included.
    mean_l=-B/(2*v*L*s.sqrt(H));mean_a=-B/(b*v*s.sqrt(H))
    jump_l=-2/(v*s.sqrt(H));jump_a=2*v*s.sqrt(H)/b**2
    assert s.simplify(mean_l*jump_a+mean_a*jump_l)==0


def test_unit_normal_and_acceleration_curvatures():
    p=solve_pair(128);d=deflected_seam(p['V'],p['field'][-1,:,0])
    mask=p['V']>=2.5;F=4/np.e**2
    H=d['H'][mask];root=np.sqrt(H)
    xp=d['xprime'][mask];yp=d['yprime'][mask];gp=d['gprime'][mask]
    for tU,tV,nU,nV in ((xp,np.ones_like(xp),-xp/root,1/root),
                       (gp,yp,gp/root,-yp/root)):
        assert np.allclose(-F*nU*nV,-1)
        assert np.allclose(-F/2*(tU*nV+tV*nU),0,atol=1e-14)
        assert np.all(nU>0) and np.all(nV>0)
    # K_ss=-n.acceleration/H; horizon Christoffel terms contribute only
    # subleading O(A), excluded from this independent leading calculation.
    kin=F*d['xsecond'][mask]/(2*H**1.5)
    kout=F*(gp*d['ysecond'][mask]-yp*d['gsecond'][mask])/(2*H**1.5)
    assert np.allclose(kin,d['K_length_minus'][mask],rtol=1e-12)
    assert np.allclose(kout,d['K_length_plus'][mask],rtol=1e-12)
    # Independently use n(r)/r, retaining the computed O(A^2) bulk radius response.
    V=d['V'];S=d['S']
    from scipy.integrate import cumulative_trapezoid as ct
    rinprime=-ct(S,V,initial=0);routprime=-ct(V*V*S,V,initial=0)/V**2
    K=2/np.e**2
    amin=(K*V*d['xprime']+rinprime-K*d['x'])/d['rootH']
    aplus=(routprime-K*d['gprime']*d['y']+K*d['g']*d['yprime'])/d['rootH']
    assert np.allclose(amin[mask],d['K_angle_minus'][mask],rtol=1e-11)
    assert np.allclose(aplus[mask],d['K_angle_plus'][mask],rtol=1e-11)


def test_matching_scaling_and_null_prefix_are_not_regularized():
    p=solve_pair(128);V=p['V'];a=p['field'][-1,:,0]
    base=deflected_seam(V,a);mask=V>=2.5
    assert max(abs(base['radius_match_residual']))<1e-12
    assert max(abs(base['temporal_match_residual']))<1e-12
    assert np.all(base['H'][mask]>0)
    assert np.all(base['S_length_leading'][mask]<0)
    for factor in (.5,2.):
        d=deflected_seam(V,factor*a)
        for key,power in [('H',2),('x',2),('y',2),('S_length_leading',1),
                          ('S_angle_leading',-1),('scalar_normal_jump_leading',2)]:
            assert np.allclose(d[key][mask],factor**power*base[key][mask],rtol=1e-12,atol=1e-14)
        assert np.allclose(d['rootH'][mask]*d['S_angle_leading'][mask],2/V[mask])
        # Both components recover #334's independently computed null source,
        # accounting for conversion from unit tangent to the common d/dV.
        from experiments.throat_dynamics.null_horizon_continuation_probe import junction
        for i in np.flatnonzero(mask)[::32]:
            null=junction(V[i],-1/V[i],1/V[i]**2,-2/V[i]**3)
            assert d['S_length_leading'][i]/d['rootH'][i]==pytest.approx(null['formal_mu_8piG'])
            assert d['rootH'][i]*d['S_angle_leading'][i]==pytest.approx(null['formal_p_8piG'])
    zero=deflected_seam(V,np.zeros_like(a))
    assert not np.any(zero['active'])
    assert np.all(np.isnan(zero['S_angle_leading']))


def test_scalar_jump_has_independent_2d_confirmation():
    errors=[]
    for N in (128,256,512):
        p=solve_pair(N);d=deflected_seam(p['V'],p['field'][-1,:,0])
        q=field_derivative_source_control(p)
        expected=d['scalar_normal_jump_leading'][-1]
        assert expected<-3e-5
        errors.append(abs(q['scalar_normal_jump'][-1]/expected-1))
    assert errors[-1]<1e-2
    assert all(a>3*b for a,b in zip(errors[:-1],errors[1:]))


def test_archive_hashes_and_limits():
    import hashlib,json
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    folder=root/'experiments/throat_dynamics/runs/20261014_spacelike_seam'
    data=json.loads((folder/'result.json').read_text())
    assert data['label']=='SPACELIKE_MATCH_REQUIRES_SURFACE_STRESS_AND_SCALAR_SOURCE'
    assert all(data['numerical_gates'].values())
    assert 'full subleading conservation' in data['uncomputed']
    for path,digest in data['sources'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest
    for path,digest in data['artifacts'].items():
        assert hashlib.sha256((folder/path).read_bytes()).hexdigest()==digest
