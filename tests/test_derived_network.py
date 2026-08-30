"""PR #216's loop driven by a derived geometry.

Four kinds of check. The *wiring* ones pin that the loop really is
`network.py`'s — reached through the APIs that already existed, not a parallel
function beside them, with the double-counted transit made unconstructible.
The *topology* ones pin that `η_topo` is read off `ConjugatePair` rather than
chosen, which is what reversed this round's verdict. The *gauge* ones pin what
a rephasing can and cannot move: it can move the closure verdict, and it
cannot move `dΨ/dω = −ωθ″` or the total variation. And the *no-fit* ones pin
the two analytic oracles — the UV law of `Ψ` and the Born falloff of `|R|` —
plus the convergence study the phase needs and unitarity never supplies.
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


def test_every_existing_network_api_sees_the_derived_transfer():
    """The wiring point: dispatch lives in `t_AB`, which all of them call.

    An earlier draft dispatched only inside a new parallel function, so these
    five APIs silently read the *transparent* ports — `effective_green` in
    particular was not the `G₀/(1−Λ)` the round was discussing.
    """
    throat = dn.derived_throat()
    legs = 0.5 * dn.EXTERIOR_LEGS
    for w in (0.5, 1.0, 2.0):
        _, transmitted = tt.scattering_matrix(np.array([w]), 0)
        eta = throat.topological_factor
        assert throat.t_AB(w) == pytest.approx(transmitted[0], abs=1e-14)
        assert nw.traverse_throat(throat, w, 0.0).factor == pytest.approx(
            eta * transmitted[0], abs=1e-14)
        lam = nw.loop_eigenvalue(throat, w, legs, legs)
        assert lam == pytest.approx(
            nw.derived_loop_eigenvalue(throat, w, legs, legs, throat.delta_BA),
            abs=1e-14), "the two spellings must be one object"
        assert nw.effective_green(throat, w, legs, legs) == pytest.approx(
            1.0 / (1.0 - lam), abs=1e-14)
        assert abs(nw.projected_kernel(throat, w, legs)["greybody_magnitude"]
                   - abs(transmitted[0])) < 1e-14


def test_a_derived_throat_cannot_be_built_with_a_nonzero_transit_time():
    """The double-count is made unconstructible rather than merely avoided."""
    with pytest.raises(ValueError, match="double-count"):
        nw.NetworkThroat(
            mouth_A=nw.NetworkMouth("A", 0.0, "L", 1.0, 0.0),
            mouth_B=nw.NetworkMouth("B", math.pi, "L", 1.0, 0.0),
            tau_th=1.0, port_A=dn._transparent_port(),
            port_B=dn._transparent_port(),
            whole_throat_transfer=lambda w: 1.0 + 0.0j)


def test_the_ports_only_decompositions_refuse_a_derived_throat():
    """A smooth barrier has no echo train; answering from a transparent port
    would return a meaningless number silently."""
    throat = dn.derived_throat()
    with pytest.raises(NotImplementedError):
        throat.loop_expansion(1.0, 3)
    with pytest.raises(NotImplementedError):
        throat.r_AA(1.0)


def test_the_topological_factor_is_derived_from_the_declared_topology():
    """`η_topo` is read off `ConjugatePair`, not chosen.

    This is the test that would have caught the earlier draft: it asserts the
    *repository's* orientations reach the loop, so `η_topo = −1` for a scalar.
    """
    from geometrodynamics.embedding.topology import make_singlet_pair
    pair = make_singlet_pair()
    assert pair.mouth_a.orientation == -pair.mouth_b.orientation

    throat = dn.derived_throat()
    assert throat.mouth_A.orientation == pair.mouth_a.orientation
    assert throat.mouth_B.orientation == pair.mouth_b.orientation
    assert throat.topological_factor == pytest.approx(-1.0 + 0.0j, abs=1e-12)


def test_the_spinor_channel_carries_the_wrap_sign_and_the_scalar_does_not():
    """Four distinct operations, not one sign: the wrap parity is a spinor
    holonomy, so folding it into the orientation would be wrong for a scalar
    — and it flips the answer, so it cannot be done silently."""
    scalar = dn.derived_throat(channel="scalar")
    spinor = dn.derived_throat(channel="spinor")
    assert scalar.topological_factor.real == pytest.approx(-1.0, abs=1e-12)
    assert spinor.topological_factor.real == pytest.approx(+1.0, abs=1e-12)
    with pytest.raises(ValueError):
        dn.derived_throat(channel="tensor")


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


def test_the_loop_closes_at_a_finite_frequency_on_the_declared_topology():
    """The corrected deliverable. With `η_topo = −1` — which the repository
    prescribes rather than permits — one clock offset *does* serve both."""
    result = dn.measure_where_the_loop_closes()
    assert result["scalar_channel_closes"]
    assert len(result["scalar_roots"]) >= 1
    assert result["scalar_roots"][0] == pytest.approx(1.4617, abs=1e-3)


def test_the_two_topology_channels_give_different_answers():
    """Which is why the sign had to be derived instead of chosen."""
    result = dn.measure_where_the_loop_closes()
    assert result["the_two_channels_disagree"]
    rows = {r["channel"]: r for r in result["rows"]}
    assert rows["scalar"]["closes"] and not rows["spinor"]["closes"]
    assert rows["scalar"]["topological_factor"] == pytest.approx(-1.0)


def test_the_round_records_that_it_reversed_its_own_earlier_verdict():
    note = dn.measure_where_the_loop_closes()["what_this_reverses"]
    assert "eta_topo = +1" in note
    assert "chosen, not" in note


def test_the_search_looks_for_tangential_zeros_not_only_sign_changes():
    note = dn.measure_where_the_loop_closes()["how_it_is_searched"]
    assert "tangential" in note
    assert "minimiser" in note


def test_the_closure_verdict_is_gauge_dependent_and_says_so():
    """`Ψ` sweeps less than `2π`, so a constant rephasing can create or remove
    the root: neither answer is a property of the geometry alone."""
    result = dn.measure_what_survives_a_rephasing()
    assert result["the_swing_is_less_than_one_branch"]
    assert result["total_variation"] < 2.0 * math.pi
    assert result["so_the_verdict_is_gauge_dependent"]
    assert "not_invariant" in "".join(result.keys())


def test_the_gauge_invariant_content_is_the_derivative_and_the_variation():
    """`dΨ/dω = −ωθ″` has no constant in it, and converges second order."""
    result = dn.measure_what_survives_a_rephasing()
    assert result["the_derivative_identity_holds"]
    rows = result["derivative_convergence"]
    assert rows[-1]["relative_error"] < 1e-4
    assert rows[1]["ratio_to_previous"] > 3.0, "second order in the step"


def test_a_linear_reference_phase_cancels_but_a_constant_does_not():
    """`θ → θ + bω` leaves `θ − ωθ′` unchanged, while `θ → θ + c` shifts it by
    `c`. So the residual gauge freedom is exactly one constant, not a
    function — computed by rebuilding `Ψ` from a rephased transfer, not by
    rearranging the identity on paper."""
    omega = np.array([0.7, 1.4617, 3.0, 6.0])
    step = 1e-4

    def psi_from(rephase) -> np.ndarray:
        grid = np.concatenate([omega - step, omega, omega + step])
        phase = np.unwrap(np.angle(dn.loop_response(grid))) + rephase(grid)
        n = omega.size
        centre = phase[n:2 * n]
        derivative = (phase[2 * n:] - phase[:n]) / (2.0 * step)
        return centre - omega * derivative

    reference = psi_from(lambda w: np.zeros_like(w))
    assert np.allclose(reference, dn.closure_function(omega), atol=1e-9)

    for b in (0.5, -2.0):                      # linear: cancels identically
        assert np.allclose(psi_from(lambda w, b=b: b * w), reference, atol=1e-8)

    for c in (0.25, math.pi):                  # constant: shifts by exactly c
        assert np.allclose(psi_from(lambda w, c=c: np.full_like(w, c)),
                           reference + c, atol=1e-8)


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


def test_the_closure_function_has_the_predicted_ultraviolet_law():
    """`ωΨ → −∫V_ℓ ds`, no fitted constant."""
    result = dn.measure_the_closure_residual_has_an_analytic_ultraviolet_law()
    assert result["the_limit_is_reached"]
    assert result["the_approach_is_monotone"]
    assert result["relative_error_at_the_top"] < 5e-3
    assert result["predicted_limit_of_omega_times_C"] == pytest.approx(
        -9.0 * math.pi / 8.0, rel=1e-12)


def test_the_ultraviolet_law_constrains_the_tail_and_not_the_interior():
    """The `1/ω` decay was misread once as proving no finite root exists. It
    constrains the tail only, and the module has to say so — because on the
    declared topology there *is* a finite root."""
    result = dn.measure_the_closure_residual_has_an_analytic_ultraviolet_law()
    note = result["what_this_does_NOT_settle"]
    assert "asymptotic law constrains" in note
    assert "measure_where_the_loop_closes" in note
    assert dn.measure_where_the_loop_closes()["scalar_channel_closes"]


def test_the_closure_function_decays_to_the_topological_constant():
    """`Ψ → arg η_topo`, not to zero — and the offset halves as `ω` doubles.

    With `η_topo = −1` the ultraviolet limit is `π`, which is the *furthest*
    a phase can sit from a branch `2πn`. Subtracting the constant is what
    makes the `1/ω` law visible at all.
    """
    omega = np.array([5.0, 10.0, 20.0])
    constant = np.angle(dn.derived_throat().topological_factor)
    offset = np.abs(np.angle(np.exp(
        1j * (dn.closure_function(omega) - constant))))
    ratios = offset[:-1] / offset[1:]
    assert np.all(np.abs(ratios - 2.0) < 0.1), "must halve as omega doubles"
    assert abs(constant - math.pi) < 1e-12


def test_the_ultraviolet_law_removes_the_topological_constant():
    """Forgetting to would make `ωΨ` diverge, which is how the omission hid
    at `η_topo = +1`, where the constant is zero."""
    result = dn.measure_the_closure_residual_has_an_analytic_ultraviolet_law()
    assert result["topological_constant_removed"] == pytest.approx(math.pi)
    assert "diverges linearly" in result["why_the_constant_must_be_subtracted"]
    raw = np.array(result["psi"]) * np.array(result["omega"])
    assert abs(raw[-1]) > 10.0 * abs(
        result["predicted_limit_of_omega_times_C"]), (
        "the unsubtracted product must visibly diverge")


def test_the_ultraviolet_slope_approaches_the_born_value():
    result = dn.measure_the_ultraviolet_slope_matches_born()
    assert result["predicted_slope"] == pytest.approx(-4.0, abs=1e-12)
    assert result["the_slope_approaches_the_born_value"]
    assert abs(result["slope_at_the_top"] - (-4.0)) < 0.3


def test_the_born_slope_scales_with_the_throat_radius():
    """`−4a`: at fixed physical `ω` a **wider** throat falls off *more*
    steeply, since the exponent is `−4aω`. In the scale-free variable `aω`
    the curves coincide — `a` sets the scale, not the shape."""
    wide = dn.measure_the_ultraviolet_slope_matches_born(0, 2.0)
    narrow = dn.measure_the_ultraviolet_slope_matches_born(0, 1.0)
    assert wide["predicted_slope"] == pytest.approx(-8.0, abs=1e-12)
    assert wide["predicted_slope"] < narrow["predicted_slope"], (
        "larger a must give a steeper (more negative) slope in physical omega")


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


def test_the_absence_of_resonance_is_scoped_to_the_tested_band():
    """It says nothing below the low end or above the high end, and the
    module has to say so rather than implying all finite frequencies."""
    result = dn.measure_no_perfect_transmission_resonance()
    note = result["this_is_a_finding_not_a_theorem"]
    assert "ON THE TESTED BAND" in note
    assert "Nothing here rules one out" in note
    assert "tested band" in result["consequence_for_the_loop"]


# ── the phase-sensitive verdict has its own convergence study ───────────────

def test_the_closure_root_is_stable_under_refinement():
    """Unitarity at `1e-13` constrains moduli, not `arg T`. The root has to be
    re-measured against the matching edge, the spatial step, and the
    finite-difference step together."""
    result = dn.measure_the_closure_root_is_numerically_converged()
    assert result["the_root_is_stable"]
    assert result["worst_shift"] < 1e-3
    variants = {r["variant"] for r in result["rows"]}
    assert {"edge 150", "edge 300"} <= variants, "matching radius varied"
    assert {"steps 30000", "steps 120000"} <= variants, "spatial step varied"
    assert {"fd 1e-3", "fd 1e-5"} <= variants, "difference step varied"


def test_the_root_is_quoted_to_the_precision_the_spread_supports():
    result = dn.measure_the_closure_root_is_numerically_converged()
    assert result["quoted_root"] == pytest.approx(1.4617, abs=1e-3)
    assert "MODULI only" in result["why_unitarity_is_not_enough"]


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
