from dataclasses import replace
import numpy as np
import pytest
from geometrodynamics.transaction import aperture_transfer as a


def test_round_spectrum_and_analytic_antipodal_propagator():
    c=a.Config(lmax=12);k,W,_=a.operators(c)
    assert np.allclose(np.sqrt(k)/np.pi,np.round(np.sqrt(k)/np.pi))
    # At t=1 the conformal wave velocity changes by (-1)^(l+1).
    # The antipodal port differs by (-1)^l, hence exact minus-refocusing.
    assert np.linalg.norm(np.cos(np.sqrt(k))*W[:,0]+W[:,1])<1e-12


def test_cap_parseval_tail_and_quadrature():
    for radius in (.4,.6):
        x=a.cap_coefficients(radius,56)
        assert 0<1-x@x<1e-4
        assert np.max(abs(x-a.cap_coefficients(radius,56,384)))<1e-12


def test_rotation_character_and_port_norms():
    for l in (1,2,5,12):
        for axis in ('horizontal','fiber'):
            theta=.31
            assert np.sum(a.rotation_diagonal(l,theta,axis))==pytest.approx(np.sin((l+1)*theta)/np.sin(theta),abs=1e-12)
    c=a.Config(offset=.31,lmax=12);k,W,tail=a.operators(c)
    assert np.allclose(np.diag(W.T@W),c.gamma*(1-tail))


def test_low_rank_step_against_dense_equations():
    c=a.Config(lmax=3,dt=1/64,stop=.5)
    r=a.simulate(c);k,W,_=a.operators(c);q=np.zeros(len(k));v=q.copy()
    A=np.eye(len(k))+c.dt*c.dt*np.diag(k)/4+c.dt/2*W@W.T
    for i,inc in enumerate(r['incoming']):
        vm=np.linalg.solve(A,v-c.dt*k*q/2+c.dt*W@inc)
        q+=c.dt*vm;v=2*vm-v
        assert np.allclose(W.T@vm-inc,r['outgoing'][i],atol=1e-12)
    assert np.allclose(q,r['final_q'],atol=1e-12)
    E=r['energy'];flux=np.sum(r['incoming']**2-r['outgoing']**2,axis=1)
    assert np.max(abs(np.diff(E)-c.dt*flux))<1e-13


def test_two_port_scattering_is_reciprocal_and_unitary():
    c=a.Config(squash=1.1,offset=.25,lmax=12)
    S=a.scattering(c,np.pi*np.array([.37,2.37,5.37]))
    for s in S:
        assert np.max(abs(s-s.T))<1e-12
        assert np.max(abs(s.conj().T@s-np.eye(2)))<1e-12
    # Zero frequency cannot distinguish a pure time delay.
    assert np.allclose(a.scattering(c,[0])[0],-np.eye(2))


def test_shifted_return_retains_all_flux_and_never_wraps():
    r={'time':np.arange(5.),'outgoing':np.column_stack([np.zeros(5),np.arange(5.)])}
    t,b=a.shifted_return(r)
    assert np.array_equal(t,r['time']-1.375)
    assert np.array_equal(b,-r['outgoing'][:,1])
    assert b@b==r['outgoing'][:,1]@r['outgoing'][:,1]


@pytest.mark.parametrize('change',[{'squash':0},{'aperture':float('nan')},{'dt':0},{'lmax':2.5},{'offset':-1},{'axis':'invalid'}])
def test_invalid_config_rejected(change):
    with pytest.raises(ValueError):a.operators(replace(a.Config(),**change))
