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


from experiments.closure_ledger import aperture_transfer_probe as probe
from experiments.closure_ledger import aperture_transfer_replay as replay
import copy
import shutil


@pytest.fixture(scope='module')
def audited():
    return replay.replay()


def test_pinned_replay_keeps_frozen_failure_and_discloses_correction(audited):
    assert audited['original_label']=='NUMERICALLY_UNRESOLVED'
    assert not audited['original_checks']['numerical_validity']
    assert all(audited['audited']['checks'].values())
    assert audited['audited']['label']=='FINITE_APERTURE_TRANSFER_SUPPORTED_AFTER_DIAGNOSTIC_CORRECTION'
    assert len(audited['audited']['diagnostics'])==35
    assert audited['audited']['closed_feedback_history']=='NOT_ESTABLISHED'


def test_reverse_source_checks_receiver_not_prompt_reflection(audited):
    r=probe.read(probe.RUN/'reverse_source.npz.b64')
    old=replay.portable_diagnose(r);new=audited['audited']['diagnostics']['reverse_source']
    assert old['early_leak_fraction']>.3 and not old['valid']
    assert new['physical_receiver']=='A' and new['early_leak_fraction']<1e-15 and new['valid']
    assert new['capture_fraction']==pytest.approx(audited['audited']['diagnostics']['b1.1_a0.4_w12_fine']['capture_fraction'],abs=1e-12)


@pytest.mark.parametrize('field',['outgoing','energy','final_q','final_v'])
def test_recorded_dynamics_tampering_is_rejected(field):
    r=probe.read(probe.RUN/'b1.1_a0.4_w12_fine.npz.b64')
    r[field].flat[100]+=.01
    assert not replay.portable_diagnose(r)['valid']


def test_source_and_grid_tampering_are_rejected():
    r=probe.read(probe.RUN/'b1.1_a0.4_w12_fine.npz.b64')
    r['incoming'][0,0]+=.1
    with pytest.raises(ValueError,match='source'):replay.portable_diagnose(r)
    r=probe.read(probe.RUN/'b1.1_a0.4_w12_fine.npz.b64');r['time'][0]+=.01
    with pytest.raises(ValueError,match='time'):replay.portable_diagnose(r)


def test_missing_or_corrupt_evidence_rejected(tmp_path):
    shutil.copyfile(probe.RUN/'manifest.json',tmp_path/'manifest.json')
    with pytest.raises(FileNotFoundError):replay.replay(tmp_path)
    n=probe.schedule()[0][0]+'.npz.b64';(tmp_path/n).write_text('bad data')
    with pytest.raises(ValueError,match='archive'):replay.replay(tmp_path)


def test_wrong_schedule_rejected():
    with pytest.raises(ValueError,match='inventory'):probe.assess({})


def test_portable_source_accepts_scalar_math_without_mutating_archive():
    import math
    r=probe.read(probe.RUN/'b1.1_a0.4_w12_fine.npz.b64')
    r['incoming'][:,0]=[math.cos(2*math.pi*t)**4*math.cos(12*math.pi*t)
                        if abs(t)<.25 else 0. for t in r['time']]
    before=r['incoming'].copy()
    local,error=replay.portable_record(r)
    assert error<1e-14
    assert np.array_equal(r['incoming'],before)
    assert local['incoming'] is not r['incoming']
    assert replay.portable_diagnose(r)['valid']


def test_portable_source_accepts_ulp_change_but_rejects_material_change():
    r=probe.read(probe.RUN/'b1.1_a0.4_w12_fine.npz.b64')
    idx=np.argmax(abs(r['incoming'][:,0]))
    r['incoming'][idx,0]=np.nextafter(r['incoming'][idx,0],np.inf)
    assert replay.portable_diagnose(r)['valid']
    r['incoming'][idx,0]+=2e-14
    with pytest.raises(ValueError,match='source'):replay.portable_diagnose(r)


@pytest.mark.parametrize('kind',['inactive','support','nan','shape'])
def test_portable_source_keeps_exact_zero_and_shape_guards(kind):
    r=probe.read(probe.RUN/'b1.1_a0.4_w12_fine.npz.b64')
    if kind=='inactive':r['incoming'][100,1]=1e-16
    elif kind=='support':r['incoming'][-1,0]=1e-16
    elif kind=='nan':r['incoming'][100,0]=np.nan
    else:r['incoming']=r['incoming'][:-1]
    with pytest.raises(ValueError,match='source'):replay.portable_diagnose(r)
