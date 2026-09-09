"""Review controls for the explicitly post-freeze uniqueness theorem."""

import copy
import json

import numpy as np
import pytest

from geometrodynamics.waves import scalar_esu_uniqueness as su
from experiments.closure_ledger import scalar_esu_uniqueness_probe as probe


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


def test_exact_expansion_rank_factorization_and_independent_routes(report):
    assert report["exact"]["all_zero"]
    assert set(report["exact"]["residuals"].values())=={"0"}


def test_generic_rank_three_does_not_imply_a_nonzero_isotropic_matrix(report):
    r=report["rank_controls"]
    assert r["min_generic_rank"]==r["max_rank"]==3
    assert r["trace_matched_identity_residual"]>1


def test_zero_B_is_a_real_exception_to_the_pointwise_alignment_claim():
    z=np.zeros(4)
    W=np.array([1.,2.,3.,4.])
    M=su.coefficient_conditions(1.,0.,0.,z,z,W)[2]
    assert np.array_equal(M,np.zeros((4,4)))
    assert np.linalg.norm(W)>0  # W cannot be a scalar multiple of B=0.
    assert np.linalg.norm(su.coefficient_conditions(1.,0.,0.,z,W,z)[2])>1


@pytest.mark.parametrize("lam",[-1.2,0.,.7])
def test_alignment_includes_zero_velocity_and_fixes_acceleration(lam):
    B=np.array([.2,.3,-.5,.7])
    V=lam*B
    W=2*lam*lam*B
    assert np.linalg.norm(su.coefficient_conditions(1.,0.,0.,B,V,W)[2])<1e-14
    direction=np.array([B[1],-B[0],0.,0.])
    assert np.linalg.norm(su.coefficient_conditions(1.,0.,0.,B,V,W+direction)[2])>.01


def test_three_conditions_keep_the_constant_quadratic_sphere_term():
    # A constant-on-sphere quadratic polynomial need not have C=0 before
    # imposing the rank theorem. This guards the order of the proof.
    x=np.array([.5,.5,.5,.5])
    mu=.7
    C=-mu
    assert C+x@(mu*np.eye(4))@x==pytest.approx(0.)
    assert np.linalg.matrix_rank(mu*np.eye(4))==4


@pytest.mark.parametrize("index",[0,1,2])
def test_offshell_expansion_and_full_stress_momentum_agree(report,index):
    r=report["pointwise"][index]
    assert r["wave_expansion_absolute"]/r["wave_scale"]<1e-10
    assert r["momentum_absolute"]<1e-10
    assert r["non_solution_wave_size"]>.1


def test_constant_B_candidate_has_a_pole_when_its_wave_equation_holds():
    b=.7
    for phase in (0.,.31,np.pi/2,np.pi):
        A=b*np.cos(phase)
        x0=-A/b
        assert abs(x0)<=1+1e-15
        assert abs(A+b*x0)<1e-15
    # This is only a control. The all-time theorem uses the exact identity
    # and the global regularity proof, not this list of phases.


@pytest.mark.parametrize("key",su.REQUIRED_CHECKS)
@pytest.mark.parametrize("mode",["missing","failed"])
def test_supplementary_gate_failure_overwrites_stale_success(report,tmp_path,monkeypatch,key,mode):
    broken=copy.deepcopy(report)
    if mode=="missing":del broken["checks"][key]
    else:broken["checks"][key]=False
    assert broken["checks_passed"]
    monkeypatch.setattr(probe,"run_probe",lambda:broken)
    (tmp_path/"uniqueness.json").write_text('{"stale":true}')
    (tmp_path/"uniqueness.md").write_text("STALE_SUCCESS")
    assert probe.main(["--output-dir",str(tmp_path)])==1
    result=json.loads((tmp_path/"uniqueness.json").read_text())
    assert result["verdict"]["uniqueness"]=="UNRESOLVED"
    assert key in result["verdict"]["failed_checks"]
    assert "STALE_SUCCESS" not in (tmp_path/"uniqueness.md").read_text()


def test_unknown_true_key_does_not_pass():
    assert su.verdict({"unrelated":True})["uniqueness"]=="UNRESOLVED"


def test_supplement_does_not_relabel_original_archive_or_trunk_gates(report,tmp_path,monkeypatch):
    original=tmp_path/"probe.json"
    original.write_text("ORIGINAL_ARCHIVE")
    monkeypatch.setattr(probe,"run_probe",lambda:copy.deepcopy(report))
    assert probe.main(["--output-dir",str(tmp_path)])==0
    assert original.read_text()=="ORIGINAL_ARCHIVE"
    v=report["verdict"]
    assert v["post_freeze"]
    assert v["BAM_support_selection"]==v["Phi_selection"]=="NOT_DERIVED"
    assert v["causality_gate"]=="OPEN"
