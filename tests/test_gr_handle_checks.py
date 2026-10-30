import numpy as np
import pytest
from geometrodynamics.waves.gr_handle_checks import (
    deck_matrix,radial_metric,kruskal_radius,killing_time,
    stress_components,schwarzschild_profile,maxwell_potential,maxwell_scattering,
)


@pytest.mark.parametrize('epsilon',[-1,1])
def test_isometry_is_not_the_same_as_time_orientation(epsilon):
    point=np.array([.3,1.2]);future=np.array([1.,.4])
    for offset in [-.7,0.,.7]:
        j=deck_matrix(epsilon,offset);mapped=j@point
        np.testing.assert_allclose(j.T@radial_metric(*mapped)@j,radial_metric(*point),atol=3e-15)
        np.testing.assert_allclose(kruskal_radius(*point),kruskal_radius(*mapped),atol=1e-15)
        assert np.sign((j@future)[0])==epsilon
        expected=-epsilon*killing_time(*point)+2*offset
        assert abs(killing_time(*mapped)-expected)<1e-14
        if epsilon==1 or offset==0:
            np.testing.assert_allclose(j@j,np.eye(2),atol=5e-16)
        else:assert np.max(abs(j@j-np.eye(2)))>1


def test_bridge_is_vacuum_only_with_its_vacuum_lapse():
    r=np.geomspace(1.0001,100,60)
    np.testing.assert_allclose(schwarzschild_profile(r),0,atol=4e-16)
    forced=schwarzschild_profile(r,vacuum_lapse=False)
    np.testing.assert_allclose(forced[0],0,atol=3e-16)
    np.testing.assert_allclose(forced[1],-1/r**3,atol=3e-16)
    np.testing.assert_allclose(forced[2],.5/r**3,atol=3e-16)


def test_supported_ellis_profile_and_flare_out_nec():
    x=np.linspace(-10,10,101);b=.7;r=np.sqrt(x*x+b*b)
    t=stress_components(r,x/r,b*b/r**3)
    np.testing.assert_allclose(t,np.array([-b*b/r**4,-b*b/r**4,b*b/r**4]),atol=2e-15)
    # At a smooth throat with finite lapse, changing lapse gradients cannot
    # remove rho+p_r = -2 r''/r. This is a local geometric requirement.
    for lapse_first in [-2.,0.,2.]:
        t=stress_components(b,0.,1/b,2.,lapse_first,3.)
        assert abs(t[0]+t[1]+2/b**2)<2e-15


def test_stress_against_independent_christoffel_ricci_construction():
    import sympy as sp
    t,s,theta,phi=sp.symbols('t s theta phi',real=True)
    coords=[t,s,theta,phi]
    # Nontrivial positive lapse and asymmetric radius exercise all derivatives.
    n=2+s+s*s;r=3+s+s*s/2
    g=sp.diag(-n*n,1,r*r,r*r*sp.sin(theta)**2);inv=g.inv()
    gamma=[[[sp.simplify(sum(inv[a,d]*(sp.diff(g[d,c],coords[b])+
             sp.diff(g[d,b],coords[c])-sp.diff(g[b,c],coords[d]))/2 for d in range(4)))
             for c in range(4)] for b in range(4)] for a in range(4)]
    ricci=sp.zeros(4)
    for a in range(4):
        for b in range(4):
            ricci[a,b]=sp.simplify(sum(sp.diff(gamma[c][a][b],coords[c])-
                sp.diff(gamma[c][a][c],coords[b])+sum(
                gamma[c][a][b]*gamma[d][c][d]-gamma[d][a][c]*gamma[c][b][d]
                for d in range(4)) for c in range(4)))
    scalar=sp.simplify(sum(inv[a,b]*ricci[a,b] for a in range(4) for b in range(4)))
    einstein=ricci-g*scalar/2
    at={s:sp.Rational(2,5),theta:sp.pi/3}
    independent=np.array([float((einstein[0,0]/n**2).subs(at)),
        float(einstein[1,1].subs(at)),float((einstein[2,2]/r**2).subs(at))])
    evaluate=lambda f:float(f.subs(at))
    reduced=stress_components(evaluate(r),evaluate(sp.diff(r,s)),evaluate(sp.diff(r,s,2)),
        evaluate(n),evaluate(sp.diff(n,s)),evaluate(sp.diff(n,s,2)))
    np.testing.assert_allclose(reduced,independent,atol=3e-16)


def test_radiative_maxwell_potential_and_zero_mode_rejection():
    r=np.geomspace(1.00001,1000,100);x=r+np.log(r-1)
    np.testing.assert_allclose(maxwell_potential(x,2),(1-1/r)*6/r**2,atol=4e-16)
    with pytest.raises(ValueError):maxwell_potential(0,0)


def test_horizon_flux_conservation_and_endpoint_refinement():
    a=maxwell_scattering(.5,1,right=120.)
    b=maxwell_scattering(.5,1,right=240.)
    assert b['flux_error']<2e-8
    assert abs(a['horizon_flux']-b['horizon_flux'])<2e-4
    assert .4<b['horizon_flux']<.6


def test_ordinary_exteriors_are_causally_disconnected_on_cover():
    rng=np.random.default_rng(149)
    for _ in range(100):
        xa=rng.uniform(.1,4);xb=-rng.uniform(.1,4)
        ta=rng.uniform(-xa,xa);tb=rng.uniform(xb,-xb)
        # A necessary condition for any causal path, even allowing angular
        # motion, is |Delta T|>=|Delta X|. Every exterior pair fails it.
        assert abs(tb-ta)<abs(xb-xa)
