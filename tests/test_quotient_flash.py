import numpy as np
import pytest
from geometrodynamics.waves.quotient_flash import closure_band,evolve,energy,endpoint_image,project_history,antipodal_overlap,driven_zero_initial
from geometrodynamics.waves.pair_birth import PairBirth,Parameters
from geometrodynamics.waves.antipodal_throat_transport import AntipodalTransport,inertial_mouth


def test_ray_band_edges_follow_independent_transit_times():
    a=.3;length=np.sin(2*a)
    for k in range(3):
        b=closure_band(a,k);edge=b['lower'];n=2*k+1
        if b['lower_closed']:assert edge==b['capture_lower']
        else:
            assert abs((np.pi-2*a)/edge-n*np.pi-length/np.sqrt(1-edge**2))<1e-12
            for delta in [-1e-6,1e-6]:
                u=edge+delta
                assert (((np.pi-2*a)/u < n*np.pi+length/np.sqrt(1-u*u))==(delta>0))


@pytest.mark.parametrize('duration',[np.pi,3*np.pi,9.8])
def test_projected_history_at_interior_times_and_orthogonal_energy(duration):
    rng=np.random.default_rng(512);omega=np.sqrt(np.arange(1,9)*np.arange(2,10));s=(-1.)**np.arange(1,9)
    q=rng.normal(size=8)+1j*rng.normal(size=8);p=rng.normal(size=8)+1j*rng.normal(size=8)
    out=project_history(omega,s,q,p,duration);qp,pp=out['q'],out['p']
    for t in np.linspace(0,duration,13):
        a,b=evolve(omega,qp,pp,t);c,d=evolve(omega,qp,pp,duration-t)
        np.testing.assert_allclose(a,s*c,atol=4e-14)
        np.testing.assert_allclose(b,-s*d,atol=4e-13)
    assert abs(out['initial_energy']-out['projected_energy']-out['correction_energy'])<1e-12
    again=project_history(omega,s,qp,pp,duration)
    assert again['correction_energy']<1e-25


def test_source_parity_is_measured_from_actual_maxwell_current():
    m=PairBirth(Parameters(cutoff=12));a=.3;u=.26;tc=(np.pi-2*a)/u
    def source(t):
        theta=a+u*t;x=np.array([[np.sin(theta),0,np.cos(theta)],[-np.sin(theta),0,np.cos(theta)]])
        v=u*np.array([[np.cos(theta),0,-np.sin(theta)],[-np.cos(theta),0,-np.sin(theta)]])
        return m.source(x,v).sum(0)
    for t in [0.,.7,tc/2]:
        np.testing.assert_allclose(source(tc-t),(-1.)**(m.ell+1)*source(t),atol=3e-16)


def test_conformal_phase_control_is_perfect_template_but_not_real_jump_closure():
    ell=np.arange(1,12);omega=ell+.5;q=np.exp(-ell*.2);p=np.zeros_like(q)
    assert antipodal_overlap(omega,ell,q,p,3*np.pi)['orthogonal_energy_fraction']<1e-14
    out=project_history(omega,(-1.)**(ell+1),q,p,3*np.pi)
    assert abs(out['correction_energy']/out['initial_energy']-.5)<1e-14


def test_rim_loop_parities_do_not_exclude_either_handle():
    mouths=[inertial_mouth(.3),inertial_mouth(.3+np.pi)]
    _,fa=mouths[0].frame(0);_,fb=mouths[1].frame(0)
    antipodal=fb@(-np.eye(3))@fa.T
    np.testing.assert_allclose(antipodal,np.diag([1.,-1.]),atol=2e-16)
    for handle,det in [('klein',-1.),('torus',1.)]:
        loop=antipodal@AntipodalTransport(mouths,handle=handle).rim_map(0,1)
        np.testing.assert_allclose(loop@loop,np.eye(2),atol=5e-16)
        assert abs(np.linalg.det(loop)-det)<5e-16


def test_driven_symmetric_history_against_independent_ode():
    from scipy.integrate import solve_ivp
    omega=np.array([1.3,2.7]);s=np.array([1.,-1.]);duration=4.3
    source=lambda t:np.array([.4+.2*(t-duration/2)**2,.3*(t-duration/2)])
    dq,dp=driven_zero_initial(omega,source,duration,64)
    q0,p0=.5*s*dq,-.5*s*dp
    initial=np.r_[q0,p0]
    sol=solve_ivp(lambda t,y:np.r_[y[2:],source(t)-omega**2*y[:2]],
                  [0,duration],initial,rtol=2e-11,atol=1e-12,dense_output=True)
    for t in np.linspace(0,duration,11):
        x=sol.sol(t);y=sol.sol(duration-t)
        np.testing.assert_allclose(x[:2],s*y[:2],atol=5e-11)
        np.testing.assert_allclose(x[2:],-s*y[2:],atol=5e-11)
