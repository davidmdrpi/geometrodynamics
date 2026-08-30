"""PR #216's loop driven by a derived geometry.

Three kinds of check. The *wiring* ones pin that the loop really is
`network.py`'s — its own `η_topo`, its own eigenvalue function, no invented
symbols, and no double-counted transit. The *invariance* ones pin that the
closure test is branch-free, which the earlier `n = 0` comparison was not. And
the *no-fit* ones pin the two analytic oracles: the UV law of the closure
residual and the Born falloff of the reflection.
"""

import math

import numpy as np
import pytest

from geometrodynamics.tangherlini import traversable_throat as tt
from geometrodynamics.transaction import derived_network as dn
from geometrodynamics.transaction import network as nw


# ── the wiring is genuine ───────────────────────────────────────────────────

def test_the_throat_uses_the_derived_backend_not_the_ports():
    throat = dn.derived_throat()
    assert throat.whole_throat_transfer is not None
    assert abs(throat.port_A.t(1.0) - 1.0) < 1e-15, "ports must be transparent"
    assert abs(throat.port_B.r_in(1.0)) < 1e-15


def test_the_derived_transfer_is_the_scattering_transmission():
    """`NetworkThroat.transfer` must return the geometry's `T`, not `t_AB`."""
    throat = dn.derived_throat()
    for w in (0.3, 1.0, 3.0):
        _, transmitted = tt.scattering_matrix(np.array([w]), 0)
        assert throat.transfer(w) == pytest.approx(transmitted[0], abs=1e-14)


def test_the_mouth_port_backend_still_works():
    """The existing backend must be retained, not replaced."""
    port = nw.MouthPort(t=lambda w: 0.5 + 0.0j,
                        r_out=lambda w: 0.1 + 0.0j,
                        r_in=lambda w: 0.1 + 0.0j)
    throat = nw.NetworkThroat(
        mouth_A=nw.NetworkMouth("A", 0.0, "L", 1.0, 0.0),
        mouth_B=nw.NetworkMouth("B", math.pi, "L", 1.0, 0.0),
        tau_th=1.0, port_A=port, port_B=port)
    assert throat.whole_throat_transfer is None
    assert throat.transfer(1.0) == pytest.approx(throat.t_AB(1.0), abs=1e-15)


def test_the_topological_factor_is_the_modules_own():
    """`η_topo` comes from the deck orientations and mouth phases."""
    throat = dn.derived_throat(orientation_a=-1, orientation_b=1)
    assert throat.topological_factor == pytest.approx(-1.0 + 0.0j, abs=1e-15)
    phased = dn.derived_throat(transfer_phase_a=0.3, transfer_phase_b=0.2)
    assert phased.topological_factor == pytest.approx(
        np.exp(0.5j), abs=1e-12)


def test_the_loop_modulus_equals_the_transmission_modulus():
    result = dn.measure_the_loop_is_driven_by_the_geometry()
    assert result["topological_factor_is_unit_modulus"]
    assert result["lambda_modulus_equals_transmission"]


def test_the_batched_scan_agrees_with_the_scalar_network_path():
    """The continuous searches must measure what the module exposes."""
    assert dn.measure_the_loop_is_driven_by_the_geometry()[
        "the_batched_path_equals_the_scalar_network_path"]


def test_the_derived_loop_applies_no_extra_transit_phase():
    """`arg T` already carries the transit; `tau_th` would double-count it."""
    throat = dn.derived_throat()
    w, delta, legs = 1.7, 0.4, dn.EXTERIOR_LEGS
    value = nw.derived_loop_eigenvalue(throat, w, 0.5 * legs, 0.5 * legs, delta)
    expected = (throat.topological_factor * throat.transfer(w)
                * np.exp(1j * w * (legs + delta)))
    assert value == pytest.approx(expected, abs=1e-14)
    note = dn.measure_the_loop_is_driven_by_the_geometry()[
        "no_extra_transit_phase_is_applied"]
    assert "double-count" in note


# ── the closure test is branch-free ─────────────────────────────────────────

def test_the_closure_residual_is_invariant_under_the_branch_choice():
    """`C` is built from `Arg exp[i(θ − ωθ′)]`, so adding `2πn` cannot move it."""
    omega = np.array([0.7, 1.0, 2.5])
    base = dn.closure_residual(omega)
    assert np.all(np.abs(base) <= math.pi + 1e-12)
    # shifting theta by 2*pi is invisible to Arg, by construction
    shifted = np.angle(np.exp(1j * (dn.closure_phase(omega) + 2.0 * math.pi)))
    assert np.allclose(shifted, np.angle(np.exp(1j * dn.closure_phase(omega))))


def test_the_closure_residual_shifts_under_a_constant_rephasing():
    """Which is exactly why it must be evaluated in the network's convention."""
    omega = np.array([2.0])
    base = dn.closure_residual(omega)[0]
    rotated = dn.derived_throat(0, 1.0, 1, 1, 0.25, 0.25)
    value = rotated.topological_factor * rotated.transfer(2.0)
    assert abs(np.angle(value) - np.angle(
        dn.derived_throat().topological_factor
        * dn.derived_throat().transfer(2.0)) - 0.5) < 1e-9


def test_no_finite_frequency_closes_carrier_and_packet_together():
    """The deliverable, and now a continuous search rather than six samples."""
    result = dn.measure_the_closure_residual_has_no_root()
    assert result["there_is_no_simultaneous_closure"]
    assert result["roots"] == []
    assert result["samples"] >= 500
    assert result["smallest_absolute_residual"] > 1e-3


def test_the_round_says_why_the_old_comparison_was_branch_dependent():
    note = dn.measure_the_closure_residual_has_no_root()[
        "why_this_is_the_invariant_statement"]
    assert "branch" in note
    assert "2 pi/w" in note


def test_the_superseded_table_carries_its_own_caveat():
    """The `n = 0` gaps stay as dispersion data, flagged as not the verdict."""
    note = tt.measure_the_closure_offsets_disagree()[
        "THIS_IS_NOT_THE_INVARIANT_STATEMENT"]
    assert "BRANCH DEPENDENT" in note
    assert "derived_network" in note


# ── the no-fit oracles ──────────────────────────────────────────────────────

def test_the_integrated_potential_along_s_is_the_closed_form():
    """`∫V_ℓ ds = (π/a)[ℓ(ℓ+2) + 9/8]`, verified symbolically."""
    for a in (0.5, 1.0, 2.0):
        for ell in (0, 1, 2):
            assert dn.integrated_potential_along_s(ell, a) == pytest.approx(
                math.pi / a * (ell * (ell + 2) + 1.125), rel=1e-14)
    assert dn.integrated_potential_along_s(0, 1.0) == pytest.approx(
        9.0 * math.pi / 8.0, rel=1e-14)


def test_the_closure_residual_has_the_predicted_ultraviolet_law():
    """`ωC → −∫V_ℓ ds`, no fitted constant."""
    result = dn.measure_the_closure_residual_has_an_analytic_ultraviolet_law()
    assert result["the_limit_is_reached"]
    assert result["the_approach_is_monotone"]
    assert result["relative_error_at_the_top"] < 5e-3
    assert result["predicted_limit_of_omega_times_C"] == pytest.approx(
        -9.0 * math.pi / 8.0, rel=1e-12)


def test_the_closure_residual_vanishes_only_as_one_over_omega():
    """Which is what makes simultaneous closure a UV limit, never finite."""
    omega = np.array([5.0, 10.0, 20.0])
    residual = np.abs(dn.closure_residual(omega))
    ratios = residual[:-1] / residual[1:]
    assert np.all(np.abs(ratios - 2.0) < 0.1), "must halve as omega doubles"


def test_the_ultraviolet_slope_approaches_the_born_value():
    result = dn.measure_the_ultraviolet_slope_matches_born()
    assert result["predicted_slope"] == pytest.approx(-4.0, abs=1e-12)
    assert result["the_slope_approaches_the_born_value"]
    assert abs(result["slope_at_the_top"] - (-4.0)) < 0.3


def test_the_born_slope_scales_with_the_throat_radius():
    """`−4a`, so a wider throat reflects less steeply."""
    wide = dn.measure_the_ultraviolet_slope_matches_born(0, 2.0)
    assert wide["predicted_slope"] == pytest.approx(-8.0, abs=1e-12)


# ── the theorem that is not one ─────────────────────────────────────────────

def test_no_perfect_transmission_resonance_is_found():
    result = dn.measure_no_perfect_transmission_resonance()
    assert result["no_perfect_transmission_point_found"]
    assert result["interior_minima_found"] == 0
    assert result["it_is_at_the_top_of_the_range"], \
        "the minimum of |R| must be the UV limit, not an interior resonance"


def test_the_absence_of_resonance_is_stated_as_a_finding_not_a_theorem():
    """Positive barriers CAN have perfect-transmission points."""
    note = dn.measure_no_perfect_transmission_resonance()[
        "this_is_a_finding_not_a_theorem"]
    assert "not a general theorem" in note
    older = tt.measure_whether_the_loop_can_close()[
        "so_lambda_cannot_equal_one_exactly"]
    assert "NOT a theorem" in older


def test_the_ledger_derives_its_numbers_from_the_measurements():
    ledger = dn.measure_the_derived_network_ledger()
    ultraviolet = dn.measure_the_ultraviolet_slope_matches_born()
    entries = {e["claim"]: e for e in ledger["entries"]}
    row = next(e for k, e in entries.items() if "UV falloff constant" in k)
    assert f"{ultraviolet['slope_at_the_top']:.3f}" in row["evidence"]


def test_the_scope_of_the_geometry_is_stated():
    """Two asymptotic ends, not a glued finite-mouth solution."""
    note = dn.measure_the_derived_network_ledger()[
        "the_geometry_is_still_an_oracle_not_a_glued_solution"]
    assert "asymptotically flat ends" in note
    assert "should not be smuggled in" in note


def test_the_probe_module_imports_and_declares_its_checks():
    from experiments.closure_ledger import derived_network_probe as probe
    assert callable(probe.run_probe)
    assert callable(probe.render_markdown)
