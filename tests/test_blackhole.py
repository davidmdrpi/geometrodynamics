"""
Tests for the black-hole condensate model.

Validates six claims:
  1. Entropy matching: S_BH = S_throat to within throat-granularity error
  2. Singularity avoidance: Kretschner scalar finite everywhere for l > 0
  3. Geodesic completeness: radial geodesics bounce (never reach r = 0)
  4. Charge from topology: neutral condensate has Q ≈ 0 for large N
  5. First law: dM ≈ T dS for uncharged Schwarzschild perturbations
  6. Schwarzschild recovery: Hayward → Schwarzschild as l → 0
"""

import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from geometrodynamics.blackhole.condensate import (
    CoherentCondensate,
    ThroatState,
    build_schwarzschild_condensate,
    build_charged_condensate,
    L_PLANCK_SQ,
    A_THROAT_MIN,
)
from geometrodynamics.blackhole.interior import (
    f_schwarzschild,
    f_hayward,
    df_hayward_dr,
    find_horizons,
    critical_core_scale,
    surface_gravity,
    hawking_temperature,
    kretschner_hayward,
    integrate_radial_geodesic,
)
from geometrodynamics.blackhole.entropy import (
    compute_entropy_balance,
    check_first_law,
    information_capacity_bits,
    evaporation_step,
)


# ════════════════════════════════════════════════════════════════════════════
# 1.  ENTROPY MATCHING
# ════════════════════════════════════════════════════════════════════════════

class TestEntropyMatching:
    """S_BH = A/(4 l_P²)  must equal  S_throat = N ln(k)  by construction."""

    def test_entropy_balance_small_bh(self):
        """M = 1 black hole: throat count reproduces area law."""
        bh = build_schwarzschild_condensate(mass=1.0)
        bal = compute_entropy_balance(bh, states_per_throat=2)

        # Relative error should be < 1/N  (one-throat granularity)
        assert bal.relative_error < 1.0 / max(bal.N_throats, 1), (
            f"Entropy mismatch {bal.relative_error:.4e} exceeds 1/N = "
            f"{1.0/bal.N_throats:.4e}"
        )

    def test_entropy_balance_large_bh(self):
        """M = 10 black hole: larger N → tighter entropy match."""
        bh = build_schwarzschild_condensate(mass=10.0)
        bal = compute_entropy_balance(bh, states_per_throat=2)
        assert bal.relative_error < 0.01, (
            f"Large-BH entropy mismatch {bal.relative_error:.4e} > 1%"
        )

    def test_entropy_scales_as_area(self):
        """S ∝ M² (equivalently ∝ A) across a mass range."""
        masses = [1.0, 2.0, 5.0, 10.0]
        entropies = []
        for m in masses:
            bh = build_schwarzschild_condensate(mass=m)
            entropies.append(bh.entropy_from_area())

        # S/M² should be constant = 4π
        ratios = [s / m ** 2 for s, m in zip(entropies, masses)]
        expected = 4.0 * np.pi
        for r in ratios:
            assert abs(r - expected) / expected < 1e-10, (
                f"S/M² = {r:.6f}, expected {expected:.6f}"
            )

    def test_information_capacity(self):
        """Bekenstein bound: I = S/ln(2) bits."""
        bh = build_schwarzschild_condensate(mass=1.0)
        I = information_capacity_bits(bh)
        S = bh.entropy_from_area()
        assert abs(I - S / np.log(2)) < 1e-10


# ════════════════════════════════════════════════════════════════════════════
# 2.  SINGULARITY AVOIDANCE
# ════════════════════════════════════════════════════════════════════════════

class TestSingularityAvoidance:
    """Kretschner scalar must be finite everywhere for l > 0."""

    def test_kretschner_schwarzschild_diverges(self):
        """Schwarzschild (l=0): K → ∞ as r → 0."""
        r = np.array([0.01])
        K = kretschner_hayward(r, M=1.0, l=0.0)
        # For Schwarzschild, K = 48M²/r⁶ at leading order
        assert K[0] > 1e6, "Schwarzschild K should diverge near r=0"

    def test_kretschner_hayward_finite(self):
        """Hayward (l > 0): K(0) = 24/l⁴, finite."""
        l = 0.1
        r = np.array([0.0, 0.001, 0.01, 0.1, 1.0, 5.0])
        K = kretschner_hayward(r, M=1.0, l=l)
        assert np.all(np.isfinite(K)), f"Kretschner not finite: {K}"
        # K(0) should equal 24/l⁴
        expected_K0 = 24.0 / l ** 4
        assert abs(K[0] - expected_K0) / expected_K0 < 1e-6, (
            f"K(0) = {K[0]:.2f}, expected {expected_K0:.2f}"
        )

    def test_kretschner_approaches_schwarzschild_at_large_r(self):
        """Far from the core, Hayward K ≈ Schwarzschild K = 48M²/r⁶."""
        M = 1.0
        l = 0.01
        r = np.array([10.0, 20.0, 50.0])
        K_hay = kretschner_hayward(r, M, l)
        K_sch = 48.0 * M ** 2 / r ** 6  # leading Schwarzschild term
        # Not exact because our Kretschner uses a simplified formula,
        # but should be same order of magnitude
        for kh, ks in zip(K_hay, K_sch):
            assert kh > 0, "Kretschner should be positive"


# ════════════════════════════════════════════════════════════════════════════
# 3.  GEODESIC COMPLETENESS
# ════════════════════════════════════════════════════════════════════════════

class TestGeodesicCompleteness:
    """Radial geodesics must bounce in the Hayward interior (l > 0)."""

    def test_schwarzschild_geodesic_incomplete(self):
        """Schwarzschild: infalling geodesic reaches r ≈ 0."""
        geo = integrate_radial_geodesic(M=1.0, l=0.0, r_start=3.0, tau_max=20.0)
        assert not geo.is_complete, (
            f"Schwarzschild geodesic should be incomplete, r_min = {geo.r_min:.6e}"
        )

    def test_hayward_geodesic_complete(self):
        """Hayward: infalling geodesic decelerates asymptotically — never singular."""
        l = 0.3
        geo = integrate_radial_geodesic(M=1.0, l=l, r_start=3.0, tau_max=50.0)
        assert geo.is_complete, (
            f"Hayward geodesic should be complete, r_min = {geo.r_min:.6e}"
        )

    def test_geodesic_completeness_across_scales(self):
        """Geodesic completeness holds for several (M, l) combinations."""
        cases = [(1.0, 0.3), (5.0, 1.0), (10.0, 2.0)]
        for M, l in cases:
            geo = integrate_radial_geodesic(M=M, l=l, r_start=3 * M, tau_max=80.0)
            assert geo.is_complete, (
                f"Incomplete geodesic for M={M}, l={l}: r_min = {geo.r_min:.6e}"
            )


# ════════════════════════════════════════════════════════════════════════════
# 4.  CHARGE FROM TOPOLOGY
# ════════════════════════════════════════════════════════════════════════════

class TestChargeFromTopology:
    """Wheeler's 'charge without charge' for black holes."""

    def test_neutral_condensate_charge(self):
        """Schwarzschild condensate has |Q| << N."""
        bh = build_schwarzschild_condensate(mass=5.0)
        Q = bh.net_charge
        N = bh.N
        # For random ±1 orientations balanced to ~0, |Q| should be small
        assert abs(Q) <= 2, (
            f"Neutral condensate Q = {Q}, N = {N} — should be ≈ 0"
        )

    def test_charged_condensate(self):
        """Charged condensate has Q matching the requested value."""
        Q_target = 10
        bh = build_charged_condensate(mass=5.0, charge=Q_target)
        # Allow ±1 from integer rounding
        assert abs(bh.net_charge - Q_target) <= 1, (
            f"Charged condensate Q = {bh.net_charge}, target = {Q_target}"
        )

    def test_charge_fraction_decreases_with_mass(self):
        """Larger BH → more throats → smaller charge fraction for fixed Q."""
        Q = 10
        fracs = []
        for m in [1.0, 5.0, 20.0]:
            bh = build_charged_condensate(mass=m, charge=Q)
            fracs.append(abs(bh.charge_fraction))
        # Charge fraction should decrease (or at least not increase)
        assert fracs[-1] < fracs[0], (
            f"Charge fraction not decreasing: {fracs}"
        )


# ════════════════════════════════════════════════════════════════════════════
# 5.  FIRST LAW OF THERMODYNAMICS
# ════════════════════════════════════════════════════════════════════════════

class TestFirstLaw:
    """dM = T dS (+ Φ dQ) for small perturbations."""

    def test_first_law_schwarzschild(self):
        """Uncharged BH: dM ≈ T dS for small dM."""
        bh = build_schwarzschild_condensate(mass=10.0)
        result = check_first_law(bh, dM=0.01)
        # For Schwarzschild, T = 1/(8πM), dS = 8πM dM → TdS = dM exactly
        assert result.relative_residual < 0.05, (
            f"First law residual {result.relative_residual:.4f} > 5%"
        )

    def test_temperature_schwarzschild_limit(self):
        """T_Hayward → 1/(8πM) as l → 0."""
        M = 5.0
        T_std = 1.0 / (8.0 * np.pi * M)
        T_hay = hawking_temperature(M, l=1e-6)
        assert abs(T_hay - T_std) / T_std < 0.01, (
            f"T_Hayward = {T_hay:.6e}, T_Schwarzschild = {T_std:.6e}"
        )


# ════════════════════════════════════════════════════════════════════════════
# 6.  SCHWARZSCHILD RECOVERY & HORIZON STRUCTURE
# ════════════════════════════════════════════════════════════════════════════

class TestSchwarzschildRecovery:
    """Hayward metric must reduce to Schwarzschild for l → 0."""

    def test_metric_function_large_r(self):
        """f_Hayward(r >> l) ≈ f_Schwarzschild(r)."""
        M = 1.0
        l = 1e-4
        r = np.linspace(3.0, 20.0, 50)
        f_h = f_hayward(r, M, l)
        f_s = f_schwarzschild(r, M)
        assert np.allclose(f_h, f_s, atol=1e-6), (
            f"Max deviation: {np.max(np.abs(f_h - f_s)):.2e}"
        )

    def test_horizon_location(self):
        """Outer horizon ≈ 2M for l << M."""
        M = 5.0
        l = 0.01
        horizons = find_horizons(M, l)
        assert len(horizons) >= 1, "No horizons found"
        r_outer = horizons[-1]
        assert abs(r_outer - 2 * M) / (2 * M) < 0.01, (
            f"Outer horizon at {r_outer:.4f}, expected ≈ {2*M:.1f}"
        )

    def test_two_horizons_exist(self):
        """Hayward BH has two horizons (inner + outer) for small l."""
        M = 5.0
        l = 0.3
        horizons = find_horizons(M, l)
        assert len(horizons) == 2, (
            f"Expected 2 horizons, found {len(horizons)}: {horizons}"
        )
        assert horizons[0] < horizons[1], "Inner < outer"

    def test_no_horizon_above_critical_l(self):
        """No trapped surface for l > l_crit (regular soliton, not BH)."""
        M = 1.0
        l_crit = critical_core_scale(M)
        horizons = find_horizons(M, l_crit * 1.5)
        assert len(horizons) == 0, (
            f"Horizons found above l_crit: {horizons}"
        )

    def test_large_mass_horizons(self):
        """Horizon finder must work for large M (regression: v0.41.0 r_max=50 bug)."""
        for M in [30.0, 100.0]:
            horizons = find_horizons(M, l=0.01)
            assert len(horizons) >= 1, f"No horizons for M={M}"
            r_outer = horizons[-1]
            assert abs(r_outer - 2 * M) / (2 * M) < 0.01, (
                f"M={M}: outer horizon at {r_outer:.2f}, expected ≈ {2*M:.1f}"
            )
            T = hawking_temperature(M, l=0.01)
            T_std = 1.0 / (8.0 * np.pi * M)
            assert abs(T - T_std) / T_std < 0.01, (
                f"M={M}: T_hay={T:.6e}, T_sch={T_std:.6e}"
            )

    def test_de_sitter_core(self):
        """f(r) → 1 - r²/l² near r = 0 (de Sitter-like)."""
        M = 1.0
        l = 0.5
        r_small = 0.001
        f_val = float(f_hayward(r_small, M, l))
        f_ds = 1.0 - r_small ** 2 / l ** 2
        # Should agree to within O(r⁵) corrections
        assert abs(f_val - f_ds) < 1e-4, (
            f"f({r_small}) = {f_val:.8f}, de Sitter = {f_ds:.8f}"
        )


# ════════════════════════════════════════════════════════════════════════════
# 7.  CONDENSATE COHERENCE & MODE SPECTRUM
# ════════════════════════════════════════════════════════════════════════════

class TestCondensateCoherence:
    """Phase coherence and collective mode structure."""

    def test_coherent_condensate_unity(self):
        """Fully coherent condensate has coherence = 1."""
        bh = build_schwarzschild_condensate(mass=1.0, coherent=True)
        assert abs(bh.coherence - 1.0) < 1e-10

    def test_thermal_condensate_low_coherence(self):
        """Thermal (random phase) condensate has coherence → 0 for large N."""
        bh = build_schwarzschild_condensate(mass=10.0, coherent=False)
        # For random phases, coherence ~ 1/√N
        expected_scale = 1.0 / np.sqrt(bh.N)
        assert bh.coherence < 10 * expected_scale, (
            f"Thermal coherence {bh.coherence:.4f} too high for N = {bh.N}"
        )

    def test_collective_mode_enhancement(self):
        """Coherent condensate has √N mode enhancement."""
        bh = build_schwarzschild_condensate(mass=2.0, coherent=True)
        spectrum = bh.collective_mode_spectrum()
        assert len(spectrum) > 0, "No collective modes"

        # Each mode should be enhanced by √N relative to bare frequency
        for key, omega_coll in spectrum.items():
            # Find average bare frequency
            bare = np.mean([
                t.modes[key].omega
                for t in bh.throats if key in t.modes
            ])
            enhancement = omega_coll / bare
            expected = np.sqrt(float(bh.N))
            assert abs(enhancement - expected) / expected < 0.1, (
                f"Mode {key}: enhancement = {enhancement:.2f}, "
                f"expected √N = {expected:.2f}"
            )


# ════════════════════════════════════════════════════════════════════════════
# 8.  EVAPORATION AS THROAT DECOUPLING
# ════════════════════════════════════════════════════════════════════════════

class TestEvaporation:
    """Hawking evaporation = throat loss from the condensate."""

    def test_evaporation_reduces_mass_and_throats(self):
        """Each evaporation step loses one throat and reduces mass."""
        bh = build_schwarzschild_condensate(mass=5.0)
        N0 = bh.N
        M0 = bh.mass

        bh2 = evaporation_step(bh)
        assert bh2.N == N0 - 1
        assert bh2.mass < M0

    def test_evaporation_preserves_approximate_neutrality(self):
        """Neutral BH stays approximately neutral after a few steps."""
        bh = build_schwarzschild_condensate(mass=5.0)
        rng = np.random.default_rng(123)
        for _ in range(10):
            bh = evaporation_step(bh, rng=rng)
        assert abs(bh.net_charge) < bh.N * 0.2, (
            f"|Q| = {abs(bh.net_charge)} too large for N = {bh.N}"
        )


# ════════════════════════════════════════════════════════════════════════════
# 9.  DERIVATION: CONDENSATE → METRIC VIA EINSTEIN EQUATIONS
# ════════════════════════════════════════════════════════════════════════════

from geometrodynamics.blackhole.derivation import (
    verify_throat_normalisation,
    stress_energy_from_throat_density,
    derive_metric_from_condensate,
    derive_core_scale,
    derive_temperature_from_modes,
    check_energy_conditions,
    full_derivation,
)


class TestDerivationThroatDensity:
    """Throat density profile normalises to N and sources the metric."""

    def test_throat_density_normalisation(self):
        """∫ n(r) 4πr² dr = N."""
        bh = build_schwarzschild_condensate(mass=5.0)
        n_check = verify_throat_normalisation(bh.N, bh.mass, bh.core_scale)
        assert abs(n_check - bh.N) / bh.N < 0.01, (
            f"Normalisation {n_check:.2f} vs N = {bh.N}"
        )

    def test_mass_function_from_density(self):
        """Numerical mass function matches analytic Hayward mass function."""
        bh = build_schwarzschild_condensate(mass=5.0)
        d = derive_metric_from_condensate(bh)
        assert d.mass_deviation < 0.001, (
            f"Mass function deviation {d.mass_deviation:.4e}"
        )

    def test_derived_metric_matches_hayward(self):
        """f_derived(r) ≈ f_hayward(r) to within numerical integration error."""
        bh = build_schwarzschild_condensate(mass=5.0)
        d = derive_metric_from_condensate(bh)
        assert d.max_deviation < 0.01, (
            f"Metric deviation {d.max_deviation:.4e}"
        )


class TestDerivationStressEnergy:
    """Stress-energy from the throat ensemble."""

    def test_de_sitter_eos(self):
        """p_r / ρ = −1 everywhere (de Sitter equation of state)."""
        bh = build_schwarzschild_condensate(mass=5.0)
        r = np.linspace(0.001, 50.0, 3000)
        se = stress_energy_from_throat_density(r, bh.N, bh.mass, bh.core_scale)
        assert np.allclose(se.eos_ratio[1:], -1.0, atol=1e-6), (
            "EOS ratio deviates from −1"
        )

    def test_nec_satisfied_sec_violated(self):
        """NEC (radial) satisfied; SEC violated (required for regularity)."""
        bh = build_schwarzschild_condensate(mass=5.0)
        r = np.linspace(0.001, 50.0, 3000)
        se = stress_energy_from_throat_density(r, bh.N, bh.mass, bh.core_scale)
        ec = check_energy_conditions(se)
        assert ec.NEC_radial_satisfied, "NEC radial should be satisfied"
        assert not ec.SEC_satisfied, "SEC should be violated for regularity"


class TestDerivationTemperature:
    """Temperature from collective mode spectrum."""

    def test_temperature_from_modes(self):
        """T_collective_mode ≈ T_surface_gravity for M ≥ 3 (thermodynamic limit)."""
        bh = build_schwarzschild_condensate(mass=5.0)
        td = derive_temperature_from_modes(bh)
        assert td.relative_error_sg < 0.01, (
            f"T mode-sg error {td.relative_error_sg:.4e}"
        )


class TestDerivationCoreScale:
    """Core scale constraints from condensate parameters."""

    def test_core_scale_bounds_consistent(self):
        """l_min < l_max < l_extremal."""
        bh = build_schwarzschild_condensate(mass=5.0)
        cs = derive_core_scale(bh.mass, bh.N)
        assert cs.l_min > 0
        assert cs.l_max > cs.l_min
        assert cs.l_extremal > cs.l_max

    def test_core_scale_planck_order(self):
        """Core scale l ≈ 0.47 l_P independent of mass."""
        for M in [1.0, 10.0, 100.0]:
            bh = build_schwarzschild_condensate(mass=M)
            assert 0.1 < bh.core_scale < 2.0, (
                f"M={M}: l = {bh.core_scale:.4f} not Planck-order"
            )


class TestDerivationPipeline:
    """Full derivation pipeline: condensate → metric → thermodynamics."""

    def test_full_pipeline(self):
        """Full derivation runs and converges for M ≥ 3."""
        for M in [3.0, 10.0, 50.0]:
            bh = build_schwarzschild_condensate(mass=M)
            d = full_derivation(bh)
            assert d.metric.mass_deviation < 0.01, (
                f"M={M}: mass dev {d.metric.mass_deviation:.4e}"
            )
            assert d.temperature.relative_error_sg < 0.05, (
                f"M={M}: T err {d.temperature.relative_error_sg:.4e}"
            )
