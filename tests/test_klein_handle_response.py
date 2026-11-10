"""Localized responding Klein handle. Freeze ``5d15963``."""

import numpy as np
import pytest

from geometrodynamics.bulk import klein_handle_response as kh
from experiments.closure_ledger import klein_handle_response_probe as probe

SMALL = dict(n_ext=7, port_a=0, port_b=3, n_tube=3)


@pytest.fixture(scope="module")
def net():
    return kh.Network(**SMALL)


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


def rng(offset=0):
    return np.random.default_rng(kh.SEED + offset)


# ---- construction ---------------------------------------------------------

def test_cover_is_a_fixed_point_free_double_cover(net):
    p = net.swap
    assert np.allclose(p @ p, np.eye(net.size))
    assert np.trace(p) == 0
    assert np.abs(p @ net.laplacian - net.laplacian @ p).max() < 1e-12
    assert np.allclose(net.laplacian.sum(axis=1), 0)


def test_tubes_cross_the_sheets(net):
    """Each tube starts on its own sheet and lands on the other."""
    last = net.n_tube - 1
    for s in (0, 1):
        assert (net.tube(s, last), net.ext(1 - s, net.port_b)) in net.edges
        assert (net.ext(s, net.port_a), net.tube(s, 0)) in net.edges
    assert np.array_equal(net.tube_sites[1], net.tube_sites[0] + net.n)


@pytest.mark.parametrize("components", [1, 2])
def test_jacobian_matches_finite_differences(net, components):
    pot = 0.3 * net.potential_profile(1) if components == 2 else None
    model = kh.Model(net, components=components, potential=pot)
    z = 0.3 * rng(1).standard_normal(model.dim)
    jac = model.jacobian(z)
    eps = 1e-6
    eye = np.eye(model.dim)
    fd = (model.vector_field(z + eps * eye)
          - model.vector_field(z - eps * eye)).T / (2 * eps)
    assert np.abs(jac - fd).max() < 1e-7


@pytest.mark.parametrize("components", [1, 2])
def test_vector_field_is_hamiltonian(net, components):
    """``W f = -grad H``, with the gradient by finite differences.

    With ``Omega(u, v) = u^T W v`` and ``W[x, p] = +1`` this is Hamilton's
    equations ``x' = dH/dp``, ``p' = -dH/dx``.
    """
    pot = 0.3 * net.potential_profile(1) if components == 2 else None
    model = kh.Model(net, components=components, potential=pot)
    z = 0.3 * rng(2).standard_normal(model.dim)
    eps = 1e-6
    eye = np.eye(model.dim)
    grad = (model.energy(z + eps * eye) - model.energy(z - eps * eye)) / (
        2 * eps)
    flow_field = model.vector_field(z)[0]
    assert np.abs(model.symplectic @ flow_field + grad).max() < 1e-7


def test_mouth_load_opens_the_waist(net):
    model = kh.Model(net, components=1)
    z = np.zeros(model.dim)
    z[net.tube_sites[0]] = 0.4
    force = model.vector_field(z)[0, model.field_dim + 2:]
    assert force[0] > 0 and force[1] == 0


def test_midpoint_is_symplectic_and_symmetric(net):
    model = kh.Model(net, components=1)
    z = 0.2 * rng(3).standard_normal(model.dim)
    there = kh.midpoint_step(model, z, kh.STEP)
    back = kh.midpoint_step(model, there, -kh.STEP)
    assert np.abs(back - z).max() < 1e-12
    frames = kh.tangent_flow(model, z, np.eye(model.dim), 20, sample=20)
    w = model.symplectic
    assert np.abs(frames[20].T @ w @ frames[20] - w).max() < 1e-10


# ---- H1 / H2 --------------------------------------------------------------

def test_localized_handle_keeps_a_lagrangian_half(net):
    r = kh.localized_report(net, rng(4))
    for t in r["twists"]:
        assert t["fixed_dim"] * 2 == t["field_dim"]
        assert t["omega_on_fixed"] < 1e-9
        assert t["frozen_seam_defect"] < 1e-9
    assert r["control"]["restricted_rank"] == r["control"]["fixed_dim"]


def test_sectors_are_the_ordinary_and_twisted_handle(net):
    r = kh.localized_report(net, rng(4))
    assert r["spectrum_error"]["even"] < 1e-9
    assert r["spectrum_error"]["odd"] < 1e-9
    assert r["even_odd_spectra_differ"] > 1e-3


def test_generic_data_violate_the_seam_relation(net):
    """Negative control: data that are not deck-fixed fail."""
    model = kh.Model(net, components=1, respond=False)
    z = np.zeros(model.dim)
    z[:model.field_dim] = rng(5).standard_normal(model.field_dim)
    assert kh.seam_defect(model, z, kh.deck(model, -1), steps=40) > 0.1


def test_sheet_wronskian_is_symmetric_and_not_conserved(net):
    r = kh.local_bracket_report(net, rng(6))
    assert r["sheet_data_rank"] == r["sheet_data_dim"]
    assert r["mirror_symmetry_error"] < 1e-9
    assert r["cover_wronskian"] < 1e-9
    assert r["relative_variation"] > 1e-6


# ---- H3 / H4 --------------------------------------------------------------

def test_response_preserves_the_structure(net):
    r = kh.response_structure_report(net, rng(7))
    assert r["reversibility_error"] < 1e-9
    assert all(a["seam_defect"] < 1e-9 for a in r["amplitudes"])
    assert r["tangent_rank"] * 2 == r["phase_space_dim"]
    assert r["tangent_lagrangian_error"] < 1e-9
    assert r["tangent_fd_error"] < 1e-6


def test_order_counting(net):
    r = kh.order_counting_report(net, rng(8))
    assert all(abs(s - 2) < 0.15 for s in r["q_slopes"])
    assert all(abs(s - 3) < 0.3 for s in r["deviation_slopes"])
    assert all(row["mouth_asymmetry"] > 0 for row in r["rows"])


def test_frozen_model_does_not_move_the_mouths(net):
    model = kh.Model(net, components=1, respond=False)
    z0 = kh.admissible_state(model, -1, 1, 0.5, rng(9))
    hist = kh.flow(model, z0, 30)[:, 0]
    assert np.abs(hist[:, model.field_dim:]).max() == 0.0


# ---- H5 - H8 --------------------------------------------------------------

@pytest.mark.parametrize("s_t,s_f", kh.CELLS)
def test_seam_cell(net, s_t, s_f):
    row = {(r["s_t"], r["s_f"]): r for r in kh.seam_cell_report(net)}[
        (s_t, s_f)]
    assert row["symplectic_sign_error"] < 1e-9
    assert row["charge_parity_error"] < 1e-9
    assert 2 * row["fixed_dim"] == row["phase_space_dim"]
    if s_t == -1:
        assert row["omega_on_fixed"] < 1e-9
    else:
        assert row["omega_rank_on_fixed"] == row["fixed_dim"]
    if s_t * s_f == 1:
        assert row["charge_positive"] == row["charge_negative"] == net.size
    else:
        assert row["charge_on_fixed"] < 1e-12
    assert row["fix_meet_i_fix"] == (0 if s_f == -1 else 2 * net.size)


def test_potential_parity(net):
    r = kh.potential_parity_report(net, rng(10))
    assert r["matching_seam_defect"] < 1e-9
    assert r["mismatched_seam_defect"] > 1e-4


def test_pairing(net):
    r = kh.pairing_report(net, rng(11))
    assert r["charge_minus_wronskian"] < 1e-12
    assert r["charge_drift"] < 1e-9
    assert r["sheet_relation_error"] < 1e-9
    assert max(r["mechanical_leak"]) < 1e-12
    assert min(r["mechanical_mouth_motion"]) > 1e-7
    assert all(abs(s - 1) < 0.1 for s in r["potential_slopes"])


def test_charge_is_not_conserved_sheet_by_sheet(net):
    """The seam carries charge between the sheets."""
    assert kh.pairing_report(net, rng(11))["sheet_charge_variation"] > 1e-6


def test_neutral_model_rejects_a_potential(net):
    with pytest.raises(ValueError):
        kh.Model(net, components=1, potential=net.potential_profile(1))


# ---- verdict and probe ----------------------------------------------------

def test_verdict_statuses():
    groups = {"localized_structure": ["a"], "local_bracket": ["b"],
              "neutral_sector_pairing": ["c"], "charged_sector_pairing": ["d"]}
    v = kh.verdict({"a": True, "b": False, "d": True}, groups)
    assert v["localized_structure"] == kh.CONFIRMED["localized_structure"]
    assert v["local_bracket"] == kh.UNRESOLVED_NUMERICALLY
    assert v["neutral_sector_pairing"] == kh.INCOMPLETE
    assert v["positive_response"] == kh.INCOMPLETE      # group absent
    assert v["numerical_validity"] == kh.INCOMPLETE
    assert v["bam_cell"] == "NOT_ESTABLISHED"
    full = {**groups, "positive_response": ["e"]}
    ok = kh.verdict(dict.fromkeys("abcde", True), full)
    assert ok["numerical_validity"] == "PASS"
    bad = kh.verdict({**dict.fromkeys("abcd", True), "e": False}, full)
    assert bad["numerical_validity"] == kh.UNRESOLVED_NUMERICALLY


def test_probe_passes_at_both_sizes(report):
    assert report["checks_passed"]
    assert len(report["checks"]) == 32
    for size in ("primary", "replicate"):
        assert all(report["checks_by_size"][size].values())
    v = report["verdict"]
    assert v["numerical_validity"] == "PASS"
    for name, value in {**kh.CONFIRMED, **kh.FIXED}.items():
        assert v[name] == value
    assert set(kh.VERDICT_FIELDS) <= set(v)


def test_probe_never_claims_the_bam_cell_or_a_norm(report):
    v = report["verdict"]
    assert v["bam_cell"] == "NOT_ESTABLISHED"
    assert v["positive_response"] == "NOT_DERIVED"
    assert v["spinor_matching"] == "NOT_DERIVED"
    assert v["event_law"] == "NOT_DERIVED"
