"""Responding electric field on the Klein handle. Freeze ``4ab8bcd``."""

import numpy as np
import pytest

from geometrodynamics.bulk import klein_handle_electric as ke
from geometrodynamics.bulk import klein_handle_response as kh
from experiments.closure_ledger import klein_handle_electric_probe as probe

SMALL = dict(n_ext=7, port_a=0, port_b=3, n_tube=3)


@pytest.fixture(scope="module")
def model():
    return ke.Electric(kh.Network(**SMALL))


@pytest.fixture(scope="module")
def report():
    return probe.run_probe()


def rng(offset=0):
    return np.random.default_rng(ke.SEED + offset)


# ---- the Hamiltonian ------------------------------------------------------

def test_vector_field_is_hamiltonian(model):
    """``W f = -grad H`` with the gradient by finite differences."""
    z = 0.3 * rng(1).standard_normal(model.dim)
    eps, eye = 1e-6, np.eye(model.dim)
    grad = (model.energy(z + eps * eye) - model.energy(z - eps * eye)) / (
        2 * eps)
    assert np.abs(model.symplectic @ model.vector_field(z)[0] + grad
                  ).max() < 1e-7


def test_gauge_invariance(model):
    """Site phases with compensating link shifts change nothing physical."""
    z = 0.3 * rng(2).standard_normal(model.dim)
    chi = rng(3).standard_normal(model.n_sites)
    moved = model.gauge_transform(z, chi)
    assert abs(model.energy(moved)[0] - model.energy(z)[0]) < 1e-12
    assert np.abs(model.gauss(moved) - model.gauss(z)).max() < 1e-12
    there = ke.flow(model, z, 20)[-1]
    assert np.abs(ke.flow(model, moved, 20)[-1]
                  - model.gauge_transform(there, chi)).max() < 1e-10


def test_gauss_function_is_conserved_off_the_constraint_surface(model):
    """Conservation is a property of the flow, not of admissible data."""
    z = 0.2 * rng(4).standard_normal(model.dim)
    hist = ke.flow(model, z, 60)[:, 0]
    gauss = model.gauss(hist)
    assert np.abs(gauss[0]).max() > 1e-3
    assert np.abs(gauss - gauss[0]).max() < 1e-11


def test_uncoupled_field_reduces_to_round_thirteen(model):
    """With ``g = 0`` the charged field is the round-13 model."""
    free = ke.Electric(model.net, coupling=0.0)
    old = kh.Model(model.net, components=2)
    z_old = 0.2 * rng(5).standard_normal(old.dim)
    n = model.n_sites
    z = np.zeros(free.dim)
    z[:4 * n] = z_old[:4 * n]
    z[free.slices["q"]] = z_old[4 * n:4 * n + 2]
    z[free.slices["pq"]] = z_old[4 * n + 2:]
    a, b = ke.flow(free, z, 30)[-1, 0], kh.flow(old, z_old, 30)[-1, 0]
    assert np.abs(a[:4 * n] - b[:4 * n]).max() < 1e-11
    assert np.abs(a[free.slices["q"]] - b[4 * n:4 * n + 2]).max() < 1e-11


def test_handle_cycle_is_a_closed_loop(model):
    c = model.handle_cycle
    assert np.abs(model.divergence @ c).max() < 1e-12
    assert set(np.unique(c)) == {0.0, 1.0}
    assert np.all(c[model.tube_edges] == 1.0)
    assert np.abs(model.edge_swap @ c - c).max() < 1e-12


def test_seam_edges_bound_sheet_one(model):
    """Only the two seam edges join sheet 1 to sheet 2."""
    n = model.net.n
    crossing = [k for k, (i, j) in enumerate(model.net.edges)
                if (i < n) != (j < n)]
    assert sorted(crossing) == sorted(model.seam_edges)


# ---- E1 - E2 --------------------------------------------------------------

@pytest.mark.parametrize("s_t,s_f", ke.CELLS)
def test_deck_map_is_a_symmetry_of_the_one_hamiltonian(model, s_t, s_f):
    row = {(r["s_t"], r["s_f"]): r for r in ke.deck_report(model, rng(6))}[
        (s_t, s_f)]
    assert row["symplectic_sign_error"] < 1e-9
    assert row["energy_invariance_error"] < 1e-12
    assert row["gauss_covariance_error"] < 1e-12
    assert 2 * row["fixed_dim"] == row["phase_space_dim"]


def test_structure_survives_electric_response(model):
    st = ke.structure_report(model, rng(7))
    assert st["reversibility_error"] < 1e-9
    assert st["tangent_lagrangian_error"] < 1e-7
    assert 2 * st["tangent_rank"] == model.dim
    for cell in st["cells"]:
        for row in cell["amplitudes"]:
            assert row["symmetry_defect"] < 1e-9
            assert row["gauss_max"] < 1e-9
            assert row["total_charge"] < 1e-12


def test_admissible_state_is_on_the_gauss_surface(model):
    for s_t, s_f in ke.CELLS:
        z = ke.admissible_state(model, s_t, s_f, 0.3, rng(8), link_scale=0.1,
                                mouth=(0.02, 0.01))
        assert np.abs(model.deck(s_t, s_f) @ z - z).max() < 1e-12
        assert np.abs(model.gauss(z)).max() < 1e-12


def test_generic_data_violate_the_seam_relation(model):
    """Negative control."""
    z = 0.2 * rng(9).standard_normal(model.dim)
    d = model.deck(-1, -1)
    fwd, bwd = ke.two_sided(model, z, 40)
    assert np.abs(fwd[:, 0] - bwd[:, 0] @ d.T).max() > 0.05


def test_trapped_flux_admissibility(model):
    adm = ke.flux_admissibility_report(model)
    assert adm["cycle_dim"] == 3
    for cell in adm["cells"]:
        matched = cell["s_t"] * cell["s_f"] == 1
        assert cell["admissible_flux_dim"] == (2 if matched else 1)
        assert (cell["max_on_tube_edges"] > 1e-3) == matched


def test_handle_flux_is_not_deck_fixed_in_a_mixed_cell(model):
    z = np.zeros(model.dim)
    z[model.slices["E"]] = model.handle_cycle
    for cell in ((1, -1), (-1, 1)):
        assert np.abs(model.deck(*cell) @ z - z).max() > 1.0
    for cell in ((1, 1), (-1, -1)):
        assert np.abs(model.deck(*cell) @ z - z).max() < 1e-12


# ---- E3 - E7 --------------------------------------------------------------

@pytest.mark.parametrize("cell", [(-1, -1), (1, 1)])
def test_response_alone_generates_no_imaginary_part(model, cell):
    r = ke.spontaneous_report(model, rng(10), *cell)
    assert all(row["imaginary"] < 1e-12 and row["gauge"] < 1e-12
               for row in r["rows"])
    assert all(row["max_mouth"] > 1e-7 for row in r["rows"])


def test_conjugation_is_a_symmetry_of_the_flow(model):
    z = 0.2 * rng(11).standard_normal(model.dim)
    c = model.conjugation
    assert abs(model.energy(c @ z)[0] - model.energy(z)[0]) < 1e-12
    a = ke.flow(model, c @ z, 30)[-1, 0]
    b = c @ ke.flow(model, z, 30)[-1, 0]
    assert np.abs(a - b).max() < 1e-11


@pytest.mark.parametrize("cell", [(-1, -1), (1, 1)])
def test_trapped_flux_couples_the_sectors(model, cell):
    f = ke.flux_report(model, rng(12), *cell)
    assert all(abs(s - 1) < 0.1 for s in f["imaginary_slopes"])
    assert f["odd_error"] < 1e-9 and f["even_error"] < 1e-9
    assert f["imaginary_scale"] > 1e-6
    assert all(abs(s - 2) < 0.15 for s in f["flux_slopes_in_amplitude"])
    assert all(abs(s - 1) < 0.1 for s in f["flux_slopes_in_flux"])
    assert f["gauss_identity_error"] < 1e-9
    assert f["total_charge"] < 1e-12


def test_charge_separates_only_across_the_time_reversing_seam(model):
    """Observed after the freeze; see the write-up."""
    klein = ke.flux_report(model, rng(13), -1, -1)
    plain = ke.flux_report(model, rng(13), 1, 1)
    assert klein["sheet_charge_max"] > 1e-9
    assert klein["sheet_charge_antisymmetry_error"] < 1e-9
    assert plain["sheet_charge_max"] < 1e-12


# ---- verdict and probe ----------------------------------------------------

def test_verdict_statuses():
    groups = {name: [name] for name in ke.CONFIRMED}
    ok = ke.verdict({name: True for name in groups}, groups)
    assert ok["numerical_validity"] == "PASS"
    assert ok["flux_origin"] == "INITIAL_DATA_NOT_DERIVED"
    bad = ke.verdict({**{name: True for name in groups},
                      "flux_pairing": False}, groups)
    assert bad["flux_pairing"] == ke.UNRESOLVED_NUMERICALLY
    assert bad["numerical_validity"] == ke.UNRESOLVED_NUMERICALLY
    partial = dict(groups)
    del partial["pairing_carrier"]
    missing = ke.verdict({name: True for name in partial}, partial)
    assert missing["pairing_carrier"] == ke.INCOMPLETE
    assert missing["numerical_validity"] == ke.INCOMPLETE


def test_probe_passes_at_both_sizes(report):
    assert report["checks_passed"]
    assert len(report["checks"]) == 27
    for size in ("primary", "replicate"):
        assert all(report["checks_by_size"][size].values())
    v = report["verdict"]
    assert v["numerical_validity"] == "PASS"
    for name, value in {**ke.CONFIRMED, **ke.FIXED}.items():
        assert v[name] == value
    assert set(ke.VERDICT_FIELDS) <= set(v)


def test_probe_claims_no_origin_for_the_flux(report):
    v = report["verdict"]
    assert v["flux_origin"] == "INITIAL_DATA_NOT_DERIVED"
    assert v["bam_cell"] == "NOT_ESTABLISHED"
    assert v["positive_response"] == "NOT_DERIVED"
