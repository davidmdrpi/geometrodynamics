import numpy as np
import pytest

from geometrodynamics.waves.alice_boundary import (
    METRIC,compose_maps,faraday,fields,frequency_map,gluing_jacobian,
    solve_history,stress,transport_faraday,
    time_deck,time_reversing_initial,
)


@pytest.mark.parametrize('time_orientation',[-1,1])
@pytest.mark.parametrize('spatial',[np.eye(3),np.diag([-1.,1.,1.]),-np.eye(3)])
def test_tensor_gluing_preserves_metric_and_stress(time_orientation,spatial):
    j=gluing_jacobian(spatial,time_orientation,[.2,-.1,.05],[-.3,.04,.1])
    np.testing.assert_allclose(j.T@METRIC@j,METRIC,atol=5e-16)
    f=faraday([.2,.3,-.4],[.1,.7,.2])
    for charge in [-1,1]:
        new=transport_faraday(f,j,charge)
        np.testing.assert_allclose(stress(new),j.T@stress(f)@j,atol=8e-16)
        assert stress(new)[0,0]>0
    k=np.array([1.,.6,.8,0.]);mapped=np.linalg.solve(j,k)
    assert abs(mapped@METRIC@mapped)<16*np.finfo(float).eps*np.dot(mapped,mapped)
    assert np.sign(mapped[0])==time_orientation


def test_static_signs_follow_faraday_pullback():
    e=np.array([.2,.5,-.3]);b=np.array([.7,.1,.8]);r=np.diag([-1.,1.,1.])
    for time in [-1,1]:
        f=transport_faraday(faraday(e,b),gluing_jacobian(r,time))
        en,bn=fields(f)
        np.testing.assert_allclose(en,time*r.T@e)
        np.testing.assert_allclose(bn,np.linalg.det(r)*r.T@b)


@pytest.mark.parametrize('orientation',[-1,1])
def test_frequency_map_matches_real_time_wave(orientation):
    r=np.diag([-1.,1.,-1.]);a=np.array([.7+.2j,-.3+.6j,.1-.4j])
    omega,offset=1.7,.38
    u,v=frequency_map(r,omega,offset,orientation)
    mapped=u@a+v@a.conj()
    for t in np.linspace(-3,3,21):
        expected=orientation*r.T@np.real(a*np.exp(-1j*omega*(orientation*t+offset)))
        np.testing.assert_allclose(np.real(mapped*np.exp(-1j*omega*t)),expected,atol=1e-15)


def test_two_time_reversals_are_linear_with_correct_clock_composition():
    first=frequency_map(np.eye(3),2.,.4,-1)
    second=frequency_map(np.eye(3),2.,.9,-1)
    u,v=compose_maps(first,second)
    np.testing.assert_allclose(v,0,atol=1e-16)
    np.testing.assert_allclose(u,np.eye(3)*np.exp(2j*(.4-.9)),atol=2e-16)


def test_antilinear_closed_loop_distinguishes_compatibility_from_residual():
    phase=.71;v=np.array([[np.exp(1j*phase)]]);u=np.zeros((1,1),complex)
    allowed=np.array([1j*np.exp(.5j*phase)])
    good=solve_history(u,v,allowed)
    assert good['status']=='nonunique' and good['nullity']==1
    np.testing.assert_allclose(good['solution']-v@good['solution'].conj(),allowed,atol=5e-16)
    bad=solve_history(u,v,np.array([np.exp(.5j*phase)]))
    assert bad['status']=='incompatible' and bad['residual']>.99


def test_open_mixed_linear_antilinear_loop_matches_direct_iteration():
    u=np.array([[.1,.07j],[.03j,.2]])
    v=np.array([[.15j,.03],[.08,.1j]])
    b=np.array([.2+.5j,-.6+.1j]);a=np.zeros(2,complex)
    for _ in range(200):a=u@a+v@a.conj()+b
    result=solve_history(u,v,b)
    assert result['status']=='unique'
    np.testing.assert_allclose(result['solution'],a,atol=1e-15)


def test_nonlinear_time_reversing_quotient_by_independent_backward_evolution():
    from geometrodynamics.waves.deformable_mode_conversion import DeformableDefect,Parameters
    cell=DeformableDefect(Parameters(length=80.,dx=.2,twisted=True))
    u,v=cell.initial();cell.add_packet(u,v,center=8.,width=2.,amplitude=.08)
    u0,v0=time_reversing_initial(u,v)
    cover=DeformableDefect(Parameters(length=160.,dx=.2))
    np.testing.assert_array_equal(time_deck(u0),u0)
    np.testing.assert_array_equal(time_deck(v0),-v0)
    up,vp=u0.copy(),v0.copy();um,vm=u0.copy(),v0.copy()
    for _ in range(500):
        up,vp=cover.step(up,vp,.04)
        um,vm=cover.step(um,vm,-.04)
    np.testing.assert_allclose(um,time_deck(up),atol=3e-14)
    np.testing.assert_allclose(vm,-time_deck(vp),atol=3e-13)
    # A velocity with the ordinary spatial twist would be the wrong Cauchy data.
    wrong=np.concatenate([v,cell.deck(v)],axis=1)
    assert np.max(abs(wrong+time_deck(wrong)))>.1
