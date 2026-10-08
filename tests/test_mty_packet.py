import numpy as np
import pytest
from geometrodynamics.transaction import mty_packet as m


def test_clock_history_and_locally_forward_handle():
    c=m.Config(3.,1.)
    h=m.clock_history(3.)
    assert h['proper_duration']==1.5
    assert h['positive_work_per_rest_mass']==2
    ds=m.delays(c)
    assert ds[0,1]==-1.375
    assert ds[1,1]==1.625
    assert ds[0,1]+h['offset']==m.HANDLE_TIME
    assert ds[1,1]-h['offset']==m.HANDLE_TIME


def test_shift_does_not_wrap_and_advance_is_explicit():
    a=np.arange(5.)
    assert np.array_equal(m.shift(a,2.,1.),[0,0,0,1,2])
    assert np.array_equal(m.shift(a,-2.,1.),[2,3,4,0,0])
    with pytest.raises(ValueError): m.shift(a,.3,1.)


def test_filter_against_direct_two_state_trapezoidal_solve():
    dt=.01; f=np.column_stack([np.sin(np.arange(50)*.2),np.cos(np.arange(50)*.1)])
    q,v=m.linear_response(f,dt)
    for j in range(2):
        z=np.zeros(2)
        L=np.array([[0,1],[-m.STIFFNESS[j]/m.MASS[j],-3/m.MASS[j]]])
        for n in range(len(f)):
            z=np.linalg.solve(np.eye(2)-dt*L/2,(np.eye(2)+dt*L/2)@z+dt*np.array([0,f[n,j]/m.MASS[j]]))
            assert np.allclose(z,[q[n+1,j],v[n+1,j]],atol=1e-12,rtol=1e-10)


def test_discrete_quartic_force_accounts_for_potential_work():
    q=np.array([[.2,-.3],[.4,.1],[-.1,.25]])
    force=m.nonlinear_force(q,m.Config(0,1))
    work=force*np.diff(q,axis=0)
    assert np.allclose(work,np.diff(m.QUARTIC*q**4/4,axis=0),atol=1e-15)


def test_local_scattering_energy_and_momentum_from_action():
    c=m.Config(0,1,dt=1/128,start=-1,stop=1,nonlinear=False)
    t=m.grid(c); x=np.zeros((len(t),2,3)); x[:,0,0]=np.sin(3*t)
    _,r=m.response(x,c);q,v=r['q'],r['v']
    vm=(v[:-1]+v[1:])/2;qm=(q[:-1]+q[1:])/2
    power=np.sum(r['incoming']**2-r['outgoing']**2,axis=2)
    assert np.max(abs(np.diff(m.energy(q,v,False),axis=0)-c.dt*power))<1e-12
    force=np.sum(r['incoming']-r['outgoing'],axis=2)-m.STIFFNESS*qm
    assert np.max(abs(m.MASS*np.diff(v,axis=0)-c.dt*force))<1e-12
    assert np.max(abs(np.diff(q,axis=0)-c.dt*vm))<1e-12


def test_inventory_and_admissible_grid():
    from experiments.closure_ledger import mty_packet_probe as p
    cases=p.schedule()
    assert len(cases)==32 and len({n for n,_ in cases})==32
    for _,c in cases:
        ds=m.delays(c)/c.dt
        assert np.max(abs(ds-np.round(ds)))<1e-10
    with pytest.raises(ValueError): p.assess({})


@pytest.mark.parametrize('amplitude',[0.,-1.,float('nan'),float('inf')])
def test_invalid_amplitude_is_rejected(amplitude):
    with pytest.raises(ValueError):m.grid(m.Config(0,amplitude))
