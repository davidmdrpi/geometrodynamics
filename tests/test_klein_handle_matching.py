"""The matching map. Freeze ``3d3ee38``."""

import numpy as np
import pytest

from geometrodynamics.bulk import klein_handle_matching as km
from geometrodynamics.bulk import klein_handle_response as kh
from experiments.closure_ledger import klein_handle_matching_probe as probe

SMALL = dict(n_ext=7, port_a=0, port_b=3, n_tube=3)


@pytest.fixture(scope="module")
def net():
    return kh.Network(**SMALL)


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


def rng(offset=0):
    return np.random.default_rng(km.SEED + offset)


# ---- Part I: the recorded maps -------------------------------------------

def test_antipode_is_a_fibre_translation():
    a = km.antipode_report(rng())
    assert a["is_fibre_translation_by_pi"] < 1e-12
    assert a["hopf_base_moved"] < 1e-12
    assert a["scalar_error"] < 1e-10
    assert {row["charge"] for row in a["rows"]} == set(range(-4, 5))


def test_charge_is_read_off_the_fibre_rotation():
    """A bidegree ``(p, r)`` monomial has charge ``p - r``."""
    z = km.random_sphere(50, rng(1))
    for p, r in ((2, 0), (1, 2), (0, 3), (2, 2)):
        rotated = km.evaluate(p, r, km.fibre_rotate(z, 0.4))
        assert np.allclose(rotated, np.exp(1j * (p - r) * 0.4)
                           * km.evaluate(p, r, z))


def test_the_two_spheres_are_different():
    t = km.two_spheres_report(rng(2))
    assert t["antipode_leaves_equator"] < 1e-12
    assert t["hopf_two_to_one_error"] < 1e-12
    assert t["hopf_formula_error"] < 1e-12
    assert t["sigma_leaves_equator"] > 0.5
    assert t["sigma_covers_base_antipode"] < 1e-12
    assert t["sigma_reverses_fibre"] < 1e-12


def test_antipode_moves_the_brane_mouth_but_not_the_hopf_base():
    """Same map, opposite behaviour on the two 2-spheres."""
    e = km.random_equator(60, rng(3))
    moved = np.abs(km.as_quaternion(km.antipode(e))
                   - km.as_quaternion(e)).max()
    assert moved > 1.0
    assert np.abs(km.hopf_base(km.antipode(e)) - km.hopf_base(e)).max() < 1e-12


def test_transport_is_the_repository_matrix():
    t = km.two_spheres_report(rng(4))
    assert t["sigma_matches_recorded_matrix"] < 1e-12
    assert t["sigma_squared_is_antipode"] < 1e-12
    assert abs(t["sigma_determinant"] - 1) < 1e-12


def test_transport_on_fibre_charged_functions():
    x = km.transport_report(rng(5))
    assert x["square_error"] < 1e-10
    assert x["fit_residual"] < 1e-10
    assert x["smallest_same_charge_residual"] > 0.1
    for row in x["rows"]:
        assert row["lands_in_opposite_charge"] < 1e-10


def test_doublet_matrix_is_j_up_to_sign():
    """``A(a z1 + b z2) = -conj(b) z1 + conj(a) z2``."""
    x = km.transport_report(rng(6))
    m = np.array(x["doublet_matrix"])
    assert np.allclose(m, [[0, -1], [1, 0]], atol=1e-10)
    a, b = 0.3 + 0.7j, -1.1 + 0.2j
    assert np.allclose(m @ np.conj([a, b]), [-np.conj(b), np.conj(a)])
    assert np.allclose(m @ np.conj(m), -np.eye(2))


def test_pr129_cell_has_no_charge(net):
    for twist in km.pr129_cell_report(net)["twists"]:
        assert twist["charge_on_fixed"] < 1e-12
        assert twist["omega_on_fixed"] < 1e-9


# ---- Part II: the multiplet ----------------------------------------------

@pytest.mark.parametrize("d,w", [(1, 1), (1, -1), (2, 1), (2, -1)])
def test_multiplet_vector_field_is_hamiltonian(net, d, w):
    model = km.Multiplet(net, components=d, loop_sign=w)
    z = 0.3 * rng(7).standard_normal(model.dim)
    eps, eye = 1e-6, np.eye(model.dim)
    grad = (model.energy(z + eps * eye) - model.energy(z - eps * eye)) / (
        2 * eps)
    assert np.abs(model.symplectic @ model.vector_field(z)[0] + grad
                  ).max() < 1e-7


def test_twisted_loop_has_holonomy_minus_one(net):
    """Exactly one edge carries the sign, and it lies on the handle loop."""
    plain = km.Multiplet(net, loop_sign=1).operator
    twisted = km.Multiplet(net, loop_sign=-1).operator
    changed = np.argwhere(np.abs(plain - twisted) > 1e-12)
    assert len(changed) == 2
    i, j = net.edges[-1]
    assert {tuple(c) for c in changed} == {(i, j), (j, i)}


def test_spatial_deck_squares_to_the_loop_sign(net):
    for w in (1, -1):
        model = km.Multiplet(net, components=1, loop_sign=w)
        v = rng(8).standard_normal((1, 1, model.n_sites))
        assert np.allclose(model.spatial_deck(model.spatial_deck(v)), w * v)
        k = model.operator
        deck = np.diag(model.sheet_sign) @ net.swap.T
        assert np.abs(deck.T @ k @ deck - k).max() < 1e-12


def test_consistency_rule(net):
    rows = {(r["components"], r["loop_sign"]): r
            for r in km.consistency_report(net, rng(9))}
    for r in rows.values():
        assert r["square_error"] < 1e-12
        assert r["energy_invariance_error"] < 1e-12
    assert rows[(1, -1)]["fixed_dim"] == 0
    assert rows[(2, 1)]["fixed_dim"] == 0
    assert 2 * rows[(1, 1)]["fixed_dim"] == rows[(1, 1)]["field_dim"]
    assert 2 * rows[(2, -1)]["fixed_dim"] == rows[(2, -1)]["field_dim"]


def test_j_commutes_with_spin_and_anticommutes_with_i(net):
    model = km.Multiplet(net, components=2, loop_sign=-1)
    deck = model.deck("J", -1)
    i_gen = model.generator("i")
    assert np.abs(deck @ i_gen + i_gen @ deck).max() < 1e-12
    for name in ("sigma_x", "sigma_y", "sigma_z"):
        g = model.generator(name)
        assert np.abs(deck @ g - g @ deck).max() < 1e-12


def test_which_charges_survive(net):
    rows = {(r["internal"], r["s_t"]): r for r in km.seam_charge_report(net)}
    full = 2 * net.size

    def alive(row, name):
        return row[name]["positive"] == row[name]["negative"] == full

    def dead(row, name):
        return row[name]["max"] < 1e-12

    assert alive(rows[("J", -1)], "i") and dead(rows[("J", -1)], "sigma_z")
    assert dead(rows[("J", 1)], "i") and alive(rows[("J", 1)], "sigma_z")
    assert dead(rows[("identity", -1)], "i")
    assert dead(rows[("identity", -1)], "sigma_x")
    assert alive(rows[("identity", 1)], "i")
    assert alive(rows[("identity", 1)], "sigma_y")


def test_transport_seam_dynamics(net):
    d = km.transport_seam_dynamics(net, rng(10))
    assert 2 * d["fixed_dim"] == d["phase_space_dim"]
    assert d["minus_j_is_quarter_rotated_j"] < 1e-12
    for row in d["rows"]:
        assert row["seam_defect"] < 1e-9
        assert row["total_spin"] < 1e-9
        assert all(v < 1e-9 for v in row["drift"].values())
        assert row["sheet_relation_error"] < 1e-9


def test_generic_doublet_data_violate_the_seam_relation(net):
    """Negative control."""
    model = km.Multiplet(net, components=2, loop_sign=-1)
    z = 0.2 * rng(11).standard_normal(model.dim)
    deck = model.deck("J", -1)
    fwd = km.flow(model, z, 40)[:, 0]
    bwd = km.flow(model, z, 40, -km.STEP)[:, 0]
    assert np.abs(fwd - bwd @ deck.T).max() > 0.05


# ---- Part III -------------------------------------------------------------

def test_pin_sign_and_flux_sign_are_independent(net):
    lab = km.label_independence_report(net, rng(12))
    assert lab["quarter_flips_pin_sign"] < 1e-12
    assert lab["quarter_preserves_flux"] < 1e-12
    assert lab["conjugation_flips_flux"] < 1e-12
    assert lab["conjugation_commutes_with_deck"] < 1e-12
    assert lab["rotated_seam_defect_opposite_pin"] < 1e-9
    assert lab["conjugated_seam_defect"] < 1e-9
    assert lab["rotated_seam_defect_same_pin"] > 1e-3


# ---- verdict and probe ----------------------------------------------------

def test_verdict_statuses():
    groups = {name: [name] for name in km.CONFIRMED}
    ok = km.verdict({name: True for name in groups}, groups)
    assert ok["numerical_validity"] == "PASS"
    assert ok["bam_cell"] == "DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION"
    bad = km.verdict({**{name: True for name in groups},
                      "pr129_cell": False}, groups)
    assert bad["pr129_cell"] == km.UNRESOLVED_NUMERICALLY
    partial = dict(groups)
    del partial["two_spheres"]
    missing = km.verdict({name: True for name in partial}, partial)
    assert missing["two_spheres"] == km.INCOMPLETE
    assert missing["numerical_validity"] == km.INCOMPLETE


def test_probe_passes(report):
    assert report["checks_passed"]
    assert len(report["checks"]) == 30
    v = report["verdict"]
    assert v["numerical_validity"] == "PASS"
    for name, value in {**km.CONFIRMED, **km.FIXED}.items():
        assert v[name] == value
    assert set(km.VERDICT_FIELDS) <= set(v)


def test_probe_does_not_decide_the_bam_cell(report):
    v = report["verdict"]
    assert v["bam_cell"] == "DEPENDS_ON_FIELD_TYPE_AND_TIME_EXTENSION"
    assert v["positive_response"] == "NOT_DERIVED"
    assert v["event_law"] == "NOT_DERIVED"
