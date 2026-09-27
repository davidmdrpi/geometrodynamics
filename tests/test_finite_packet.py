"""Physical identities, independent controls, and evidence failure paths."""
import copy
import json
import shutil
import numpy as np
import pytest
from geometrodynamics.waves import finite_packet as packet
from experiments.closure_ledger import finite_packet_probe as probe
from experiments.closure_ledger import finite_packet_spatial as spatial


def test_explicit_tensor_is_tt_eigenfunction_with_correct_weyl():
    assert max(spatial.check_degree(n) for n in (2,3,5)) < 1e-9


def test_harmonics_are_orthonormal_and_antipodal_pullback_has_tensor_parity():
    norms = packet.normalization()
    assert np.max(abs(packet.grams(norms)['full']-np.eye(79))) < 1e-8
    chi = np.array([.31,.88,1.22])
    left = packet.radial(chi)
    right = packet.radial(np.pi-chi)
    # The chi-angular component receives a minus from the pullback dchi.
    right[:,1,:] *= -1
    assert np.allclose(right,(-1.)**packet.DEGREES[:,None,None]*left,atol=1e-7,rtol=1e-10)


def test_even_projection_is_two_equal_lobes_not_one_cover_source():
    norms = packet.normalization()
    gram = packet.grams(norms)
    cover = packet.coefficients(24,6,norms)
    paired = packet.coefficients(24,6,norms,True)
    assert np.all(paired[packet.DEGREES%2==1]==0.)
    assert np.isclose(paired@gram['north']@paired,paired@gram['south']@paired,rtol=1e-10)
    assert cover@gram['north']@cover > 100*(cover@gram['south']@cover)


def test_free_dispersion_is_exact_for_arbitrary_phase_and_packet():
    norms=packet.normalization()
    c=packet.coefficients(12,3,norms)
    M=packet.analytic_modes('free')
    j=np.flatnonzero(packet.TIMES==np.pi)[0]
    for phase in (0.,np.pi/4,1.234):
        state=packet.packet_state(M,c,phase)
        target=-(-1.)**packet.DEGREES[:,None]*state[0]
        assert np.linalg.norm(state[j]-target)<1e-12


@pytest.fixture(scope='module')
def derived():
    return probe.score()


def test_archived_decisions_are_rebuilt_and_failure_retained(derived):
    archive=json.loads((probe.RUN/'packet.json').read_text())
    assert probe.old.close(archive['result'],derived,1e-10)
    assert all(derived['gates'].values())
    failures=[r for r in derived['packets'] if r['verdict']=='NOT_ESTABLISHED']
    assert len(failures)==2
    assert all(r['physical_criteria']['arrival_time'] is False for r in failures)
    assert all(r['verdict']!='LOCALIZED_FIRST_TRANSIT' for r in derived['packets'] if r['even_only'])


def test_replay_and_changed_verdict(tmp_path,derived):
    # Copy small manifests; link large immutable evidence to avoid duplication.
    for src in probe.RUN.iterdir():
        if src.name=='packet.json':
            shutil.copyfile(src,tmp_path/src.name)
        elif src.is_file():
            (tmp_path/src.name).symlink_to(src)
    assert probe.replay(tmp_path)
    archive=json.loads((tmp_path/'packet.json').read_text())
    archive['result']['packets'][0]['verdict']='NOT_ESTABLISHED'
    (tmp_path/'packet.json').write_text(json.dumps(archive))
    assert not probe.replay(tmp_path)


def test_changed_raw_file_is_rejected_before_claims(tmp_path):
    for src in probe.RUN.iterdir():
        if src.is_file() and src.name!='packet.json':
            (tmp_path/src.name).symlink_to(src)
    archive=json.loads((probe.RUN/'packet.json').read_text())
    archive['raw_files']['modes_DOP853_0.npz']='0'*64
    (tmp_path/'packet.json').write_text(json.dumps(archive))
    assert not probe.replay(tmp_path)


def test_coordinate_record_missing_a_degree_is_rejected(tmp_path):
    for src in probe.RUN.iterdir():
        if src.is_file() and src.name!='spatial.json':
            (tmp_path/src.name).symlink_to(src)
    (tmp_path/'spatial.json').write_text(json.dumps({'2':0.,'3':0.}))
    with pytest.raises(ValueError,match='spatial'):
        probe.score(tmp_path)


def test_scalar_budget_uses_transient_gain_not_just_floquet_radius(derived):
    budget=derived['instability_budget']
    assert budget['n2_max_scaled_state_gain']>budget['n2_spectral_radius']>1.24
    assert np.isclose(budget['homogeneous_gain'],np.exp(np.sqrt(2)*np.pi))
    assert budget['seed_scan'][0]['n2_metric_bound']==0.


def test_direct_reconstruction_preserves_positive_tiny_tails():
    with np.load(probe.RUN/'powers.npz') as a:
        for name in a.files:
            if name!='times':
                assert np.isfinite(a[name]).all() and np.min(a[name])>=0., name
        ratio=a['cover_24_0.00000000_weyl_south'][0]/a['cover_24_0.00000000_weyl_full'][0]
        assert 0 < ratio < 1e-16


def test_full_replay_rejects_joint_raw_and_decision_change(tmp_path,monkeypatch):
    for src in probe.RUN.iterdir():
        if src.is_file() and src.name not in ('packet.json','spatial.json'):
            (tmp_path/src.name).symlink_to(src)
    spatial_record=json.loads((probe.RUN/'spatial.json').read_text())
    spatial_record['2']=1.01e-9  # change smaller than full replay's numerical tol
    (tmp_path/'spatial.json').write_text(json.dumps(spatial_record))
    archive=json.loads((probe.RUN/'packet.json').read_text())
    archive['raw_files']['spatial.json']=probe.sha(tmp_path/'spatial.json')
    archive['result']=probe.score(tmp_path)
    (tmp_path/'packet.json').write_text(json.dumps(archive))
    assert not archive['result']['gates']['G1']
    assert probe.replay(tmp_path)  # internally consistent, explicitly partial
    def independently_measured_fixture(destination):
        for src in probe.RUN.iterdir():
            if src.is_file():
                (destination/src.name).symlink_to(src)
    monkeypatch.setattr(probe,'measure',independently_measured_fixture)
    assert not probe.replay(tmp_path,full=True)
