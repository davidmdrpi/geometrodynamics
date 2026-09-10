"""Global proof controls, direct geometric variation and unused constraints."""

import copy
import json

import numpy as np
import pytest

from geometrodynamics.waves import scalar_esu_support as ss
from experiments.closure_ledger import scalar_esu_support_probe as probe


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


def test_exact_local_and_background_certificates():
    certificate=ss.exact_certificate()
    assert certificate["all_zero"]
    assert set(certificate["residuals"].values())=={"0"}
    assert certificate["acceleration_determinant"]=="-2"


def test_isotropy_identity_against_the_inherited_full_stress(report):
    rows=report["spatial"]["inherited"]
    assert {r["degree"] for r in rows}=={1,3,5}
    assert len({r["points"] for r in rows})==2
    assert max(r["relative"] for r in rows)<1e-10
    assert min(r["max_anisotropic_stress"] for r in rows)>1e-3


def test_a_higher_order_odd_zero_is_not_mistaken_for_an_isotropic_field():
    # phi=x0^3 has zero gradient at every zero. A simple-node argument alone
    # would miss it, whereas the reciprocal/global proof still applies.
    rows=[]
    for x0 in (0.,.2,.4):
        dh=np.array([np.sqrt(1-x0*x0),0.,0.])
        hess_x=-x0*np.eye(3)
        phi=x0**3
        grad=3*x0*x0*dh
        hess=6*x0*np.outer(dh,dh)+3*x0*x0*hess_x
        rows.append(np.linalg.norm(ss.stf((2*np.outer(grad,grad)-phi*hess)/3)))
    assert rows[0]==0
    assert min(rows[1:])>1e-3


def test_nonconstant_reciprocal_controls_are_local_isotropy_not_history_witnesses(report):
    r=report["spatial"]["reciprocal"]
    assert r["max_absolute"]<1e-10
    assert r["nonconstant_range"]>.1
    assert r["h_lower_bound"]>0 and r["phi_lower_bound"]>0
    assert r["full_Einstein_solution"]=="NOT_ASSERTED"
    assert not r["pole_control"]["smooth"]


def test_even_support_satisfies_all_components_and_rejects_changed_parameters(report):
    assert len(report["background"])==6
    for r in report["background"]:
        assert r["Einstein_relative"]<1e-10
        assert r["wave_absolute"]<1e-10
        assert r["reversed_Lambda_residual"]>.1
        assert r["amplitude_1p1_residual"]>.1
        assert r["zero_history_enthalpy_residual"]==pytest.approx(2.)
        assert r["F_min"]==pytest.approx(1/(2*r["kappa"]))
        assert r["K_min"]==pytest.approx(1.)
        assert r["K_max"]==pytest.approx(4.)


def test_direct_geometry_recovers_the_round_curvature_and_offshell_trace():
    m=ss.HomogeneousSupport(radius=.7,kappa=.4)
    zero=(0.,0.,0.)
    d=ss.direct_geometry(m,.31,2,zero,zero,zero)
    assert d["R"]==pytest.approx(6/m.radius**2)
    assert np.max(np.abs(d["einstein"]-np.diag([3,-1,-1,-1])/m.radius**2))<1e-12
    d=ss.direct_geometry(m,.31,2,(.2,-.1,.3),(.11,.23,-.07),(-.3,.2,.1),epsilon=.08)
    assert abs(d["wave"])>1e-3  # the off-shell input must really be off shell
    assert abs(d["trace"]-d["phi"]*d["wave"])<1e-12


def test_stress_and_einstein_response_agree_with_direct_curvature_variation(report):
    assert len(report["variation"])==16
    for r in report["variation"]:
        assert max(r["stress_relative"])<1e-6
        assert max(r["Einstein_absolute"])/r["Einstein_scale"]<1e-6
        assert max(r["wave_absolute"])/r["wave_scale"]<1e-6
        assert r["off_shell_trace_absolute"]<1e-10
        assert r["response_trace_absolute"]<1e-10
    # Reported order is a real measured reduction, separate from the freeze's
    # residual-only acceptance gate. Steps are large enough to resolve it.
    ratios=[r["stress_error_ratio"] for r in report["variation"]]
    assert min(ratios)>3.5 and max(ratios)<4.5


@pytest.mark.parametrize("phase",[0.,np.pi/2,np.pi,3*np.pi/2])
def test_constraint_solve_and_accelerations_cross_both_turning_points(phase):
    m=ss.HomogeneousSupport(phase=phase)
    y=ss.initial_response(m,3)
    dy=ss.response_rhs(m,3,0.,y)
    r=ss.response_residuals(m,0.,3,y,(dy[1],dy[3]))
    assert all(np.isfinite(dy))
    for key in ("hamiltonian","momentum","spatial","anisotropic","KG","trace"):
        assert abs(r[key])<1e-10


def test_independent_evolution_propagates_both_unused_constraints(report):
    assert {(r["radius"],r["kappa"],r["degree"]) for r in report["evolution"]}=={
        (1.,1.,2),(1.,1.,3),(.7,.4,2)}
    for r in report["evolution"]:
        assert r["times"][-1]==pytest.approx(2*np.pi*r["radius"])
        for run in r["runs"]:
            assert max(run["relative"].values())<1e-8
        assert r["refinement_absolute"]/r["refinement_scale"]<1e-7


def test_wrong_initial_constraint_is_not_projected_away():
    m=ss.HomogeneousSupport()
    y=ss.initial_response(m,2)
    y[2]+=.01
    ts=np.array([0.,.01,.02])
    values=ss.integrate_response(m,2,ts,y)
    assert np.array_equal(values[0],y)
    r=ss.response_residuals(m,ts[-1],2,values[-1])
    assert abs(r["hamiltonian"])+abs(r["momentum"])>1e-3


def test_the_even_control_has_a_nonfluid_response_and_linear_cross_stress(report):
    for r in report["evolution"]:
        assert r["max_abs_Pi"]>1e-3
        assert r["omitted_Pi_Einstein_potential_residual"]>1e-3
        assert r["max_density_cross_term"]>1e-3
    assert report["verdict"]["homogeneous_even_control"]=="EXACT_BUT_OUTSIDE_ODD_SECTOR"
    assert report["verdict"]["BAM_support_selection"]=="NOT_DERIVED"
    assert report["verdict"]["Phi_selection"]=="NOT_DERIVED"
    assert report["verdict"]["causality_gate"]=="OPEN"


@pytest.mark.parametrize("name",ss.REQUIRED_CHECKS)
@pytest.mark.parametrize("mode",["missing","failed"])
def test_every_failed_gate_overwrites_stale_success_and_exits_nonzero(report,tmp_path,monkeypatch,name,mode):
    broken=copy.deepcopy(report)
    if mode=="missing":del broken["checks"][name]
    else:broken["checks"][name]=False
    # Stale top-level success is intentional, and must be overwritten.
    assert broken["checks_passed"]
    monkeypatch.setattr(probe,"run_probe",lambda **kwargs:broken)
    (tmp_path/"probe.json").write_text('{"verdict":"STALE_SUCCESS"}')
    (tmp_path/"probe.md").write_text("STALE_SUCCESS")
    assert probe.main(["--output-dir",str(tmp_path)])==1
    archived=json.loads((tmp_path/"probe.json").read_text())
    assert not archived["checks_passed"]
    assert name in archived["verdict"]["failed_checks"]
    assert all(archived["verdict"][key]=="UNRESOLVED" for key in ss.VERDICT_FIELDS)
    assert "STALE_SUCCESS" not in (tmp_path/"probe.md").read_text()


def test_unrelated_true_key_and_empty_checks_cannot_pass():
    for checks in ({},{"unrelated":True}):
        assert len(ss.failed_checks(checks))==len(ss.REQUIRED_CHECKS)
        assert ss.verdict(checks)["odd_sector_exact_ESU_support"]=="UNRESOLVED"


def test_passing_cli_writes_fresh_reports(report,tmp_path,monkeypatch):
    monkeypatch.setattr(probe,"run_probe",lambda **kwargs:copy.deepcopy(report))
    assert probe.main(["--output-dir",str(tmp_path)])==0
    assert json.loads((tmp_path/"probe.json").read_text())["checks_passed"]


@pytest.mark.parametrize("kwargs",[{"radius":0},{"kappa":-1},{"phase":np.nan}])
def test_invalid_background_parameters(kwargs):
    with pytest.raises(ValueError):ss.HomogeneousSupport(**kwargs)


def test_uncovered_degrees_do_not_silently_enter_the_response_solver():
    with pytest.raises(ValueError):ss.initial_response(ss.HomogeneousSupport(),1)
    with pytest.raises(ValueError):ss.integrate_response(ss.HomogeneousSupport(),0,[0.,1.])
