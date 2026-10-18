"""Check the visual model's actual guarantees, not a claimed GR derivation."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from klein_release_model import KleinRelease,energy

@pytest.fixture(scope='module')
def model():return KleinRelease()

@pytest.mark.parametrize('kwargs',[dict(drive=0),dict(gain=.08)])
def test_wave_and_weak_coupling_nulls(kwargs):
    m=KleinRelease(**kwargs)
    assert m.tc is None and m.scale==0
    assert m.progress(3)==0 and abs(m.ledger(4)['work'])<1e-14


def test_causal_readout_and_contact(model):
    m=model;assert m.tc>m.focus
    x=np.linspace(-.6,.6,61);delays=(x+.6)/1.2*m.length
    np.testing.assert_allclose(m.transmitted(x,1.5),m.transmitted(-.6,1.5-delays),atol=1e-14)
    assert np.all(m.transmitted(x,-.1)==0)
    assert abs(m.minimum(m.tc)-m.neck.level_radius)<1e-11
    assert m.minimum(m.tc-.01)>m.neck.level_radius


def test_withdrawal_is_optional_and_continuous(model):
    m=model;n=KleinRelease(detach=False)
    assert n.tc==pytest.approx(m.tc) and n.scale==0
    x=np.linspace(-.6,.6,101)[None,:];r=np.linspace(-.1,.1,41)[:,None]
    assert np.max(m.density(x,r,m.tc+.1))>0
    assert np.max(m.density(x,r,m.tc+m.withdrawal+.01))==0
    assert np.max(n.density(x,r,m.tc+m.withdrawal+.01))>0
    changes=[]
    for dt in [1e-3,1e-4,1e-5]:changes.append(np.max(abs(m.density(x,r,m.tc+dt)-m.density(x,r,m.tc-dt))))
    assert changes[2]<.11*changes[1]<.012*changes[0]


def test_forced_wave_equation_and_continuous_state(model):
    m=model;dt=1e-6;t=m.tc+.17
    a,p=m.state(t);am,pm=m.state(t-dt);ap,pp=m.state(t+dt)
    u=(t-m.tc)/m.pulse_duration
    force=m.scale*m.h*30*u*u*(1-u)**2/m.pulse_duration
    np.testing.assert_allclose((ap-am)/(2*dt),p,atol=1e-8)
    np.testing.assert_allclose((pp-pm)/(2*dt)+m.omega**2*a,force,atol=1e-8)
    for t in [m.tc,m.tc+m.pulse_duration]:
        before=m.state(t-1e-7);after=m.state(t+1e-7)
        assert max(np.max(abs(x-y)) for x,y in zip(before,after))<1e-5


def test_budget_includes_interference_and_remains_funded(model):
    m=model;r=m.report()
    assert r['release_self_energy']!=pytest.approx(m.budget)
    assert r['release_self_energy']+r['release_interference']==pytest.approx(m.budget,abs=1e-13)
    assert r['minimum_reservoir']>-1e-12
    assert r['final_added_field_energy']==pytest.approx(m.budget,abs=1e-13)
    # Forced work may initially be negative: the signed ledger must preserve it.
    assert r['minimum_added_energy']<0 and r['maximum_reservoir']>m.budget


def test_return_is_a_standing_mode_superposition(model):
    m=model;t=m.tc+m.pulse_duration+.5
    a,p=m.state(t);b,q=m.state(t+2*np.pi)
    np.testing.assert_allclose(a,-b,atol=1e-13);np.testing.assert_allclose(p,-q,atol=1e-13)
    assert energy(a,p,m.omega)==pytest.approx(energy(b,q,m.omega),abs=1e-13)
    assert m.report()['cap_refocus_ratio']>10
