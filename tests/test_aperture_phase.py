from dataclasses import replace, asdict
import copy
import numpy as np
import pytest
from geometrodynamics.transaction import aperture_phase as a
from geometrodynamics.transaction import aperture_transfer as old
from experiments.closure_ledger import aperture_phase_probe as p


def test_degeneracy_compression_matches_full_gram_transfer():
    # Outside the production box: small cutoff, non-scheduled squash/footprint.
    c=a.Config(squash=1.07,footprint=5.64,lmax=8,dt=1/128,stop=.5)
    reference=old.simulate(old.Config(squash=c.squash,aperture=c.aperture,
                          carrier=12,lmax=8,dt=c.dt,stop=c.stop))
    result=a.simulate(c)
    assert np.max(abs(result['outgoing']-reference['outgoing']))<1e-12
    assert np.max(abs(result['energy']-reference['energy']))<1e-12


def test_independent_dense_midpoint_and_energy_ledger():
    c=a.Config(carrier=6,footprint=3.6,lmax=3,dt=1/64,stop=.75)
    r=a.simulate(c);k,W,_,_,_=a.operators(c);q=np.zeros(len(k));v=q.copy()
    M=np.eye(len(k))+c.dt**2*np.diag(k)/4+c.dt/2*W@W.T
    for i,inc in enumerate(r['incoming']):
        vm=np.linalg.solve(M,v-c.dt*k*q/2+c.dt*W@inc)
        q+=c.dt*vm;v=2*vm-v
        assert np.max(abs(W.T@vm-inc-r['outgoing'][i]))<1e-12
    assert np.linalg.norm(q-r['final_q'])<1e-12
    assert np.linalg.norm(v-r['final_v'])<1e-12
    flux=c.dt*np.sum(r['incoming']**2-r['outgoing']**2,axis=1)
    assert np.max(abs(np.diff(r['energy'])-flux))<1e-13


def test_analytic_source_transform_matches_quadrature():
    from scipy.integrate import quad
    c=a.Config(carrier=6,footprint=3.6)
    for omega in (0.,2.3,6*np.pi,35.):
        exact=quad(lambda t: float(a.packet(np.array([t]),c)[0])*np.cos(omega*t),
                   -c.halfwidth,c.halfwidth,epsabs=1e-13)[0]
        assert a.source_fourier(omega,c)==pytest.approx(exact,abs=1e-13)


def test_schedule_fixes_dimensionless_source_port_parameters_and_tail():
    cases=p.schedule();assert len(cases)==57
    assert len(set(n for n,c in cases))==57
    for n,c in cases:
        assert c.aperture*c.carrier==pytest.approx(c.footprint)
        assert c.gamma/c.carrier==pytest.approx(2/3)
        assert c.halfwidth*c.carrier==3
        coeff=old.cap_coefficients(c.aperture,c.lmax)
        assert abs(1-coeff@coeff)<1e-4
        assert np.max(abs(coeff-old.cap_coefficients(c.aperture,c.lmax,384)))<1e-12


def test_exact_round_phase_and_antipodal_propagator():
    c=a.Config(lmax=8);k,W,_,_,_=a.operators(c)
    assert np.linalg.norm(np.cos(np.sqrt(k))*W[:,0]+W[:,1])<1e-12
    phase=a.phase_summary(c)
    assert phase['weighted_std']<1e-13
    assert phase['weighted_free_coherence']==pytest.approx(1,abs=1e-14)


def synthetic_diagnostics(beta):
    ds={}
    for n,c in p.schedule():
        phi=np.pi*abs(c.squash**-2-1)*c.carrier/2
        ds[n]=dict(valid=True,capture=.5 if c.squash==1 else .5*max(phi,1)**(-beta),
                   phase=dict(proxy=phi,weighted_std=c.carrier/4))
    return ds


def test_preregistered_labels_can_fail_independently():
    inverse=p.score(synthetic_diagnostics(1),{'resolution':0.})
    assert inverse['inverse_phase_label']=='SUPPORTED_IN_DECLARED_FAMILY'
    assert inverse['retention_label']=='FAILED_IN_DECLARED_FAMILY'
    plateau=p.score(synthetic_diagnostics(0),{'resolution':0.})
    assert plateau['inverse_phase_label']=='FAILED_IN_DECLARED_FAMILY'
    assert plateau['retention_label']=='SUPPORTED_IN_DECLARED_FAMILY'
    unresolved=p.score(synthetic_diagnostics(1),{'resolution':.1})
    assert unresolved['retention_label']=='NUMERICALLY_UNRESOLVED'
    assert unresolved['inverse_phase_label']=='NUMERICALLY_UNRESOLVED'


def test_source_guards_and_reconstruction_reject_corruption():
    c=a.Config(carrier=6,footprint=3.6,lmax=3,dt=1/64,stop=.75)
    r=a.simulate(c);d=a.diagnose(r)
    assert d['energy_error']<1e-12 and d['port_error']<1e-12
    altered=copy.deepcopy(r);altered['outgoing'][20,1]+=.1
    assert not a.diagnose(altered)['valid']
    altered=copy.deepcopy(r);altered['incoming'][20,1]=1e-16
    with pytest.raises(ValueError,match='source'):a.diagnose(altered)
    altered=copy.deepcopy(r);altered['time'][0]+=.01
    with pytest.raises(ValueError,match='time'):a.diagnose(altered)
    altered=copy.deepcopy(r);altered['incoming'][20,0]+=2e-14
    with pytest.raises(ValueError,match='source'):a.diagnose(altered)
