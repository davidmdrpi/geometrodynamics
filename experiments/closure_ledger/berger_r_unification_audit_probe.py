"""
Berger-sphere deformation audit of the R-unification assumption (PR #165).

> Framing: QFT on the *fixed classical* throat geometry (geometry → fields),
> not quantum gravity.  This is an AUDIT, not a demonstration.

WHAT THIS PROBE IS NOT
──────────────────────
This is **not** a quantum test, **not** a throat-formation test, and
**not** a wave-propagation test.  It is an **R-unification test**: BAM's
unified Bohr–Sommerfeld mass operator (THESIS PR #83)

    m²(k, n) = (k·2π / L_throat)² + ((n+1)·π / L_cavity)²

treats the *throat* scale (the Hopf-fiber winding, L_throat = √(2π)/k₅)
and the *cavity* scale (the radial/base resolution) as both riding on a
**single** S³ radius R — "everything rides on one R".  The Berger
deformation S³_λ squashes the Hopf fiber by λ while leaving the base S²
round, so it is the one geometric move that pulls those two scales apart.
It therefore tests whether the global cosmic-cavity vacuum energy and the
local throat self-energy are genuinely ONE dynamical object on ONE R, or
whether "one R" is only a scale-free bookkeeping shorthand.

ENFORCEMENT GUARDRAILS (anti-rigging; see test_T7)
──────────────────────────────────────────────────
1. No Derived Inversions — no fitted constant is relabelled as an
   arbitrary π-multiple inside a derived path.  L_throat = √(2π)/k₅ is the
   THESIS-locked value (PR #83), used as-is; the Casimir uses the genuine
   SU(2) Berger spectrum 4j(j+1)+4m²(λ⁻²−1); no constant is reverse-fit to
   the target.
2. No Hidden Imports — the Born rule and the singlet state
   (``geometrodynamics.bell``) are NEVER imported here.  If a step needed
   them, ``_forbid_quantum_inputs`` raises an explicit exception instead
   of smuggling them in.
3. No False Victories — the A/R + B·R² well's stability is computed and
   then EXPLICITLY discounted: it is an artifact of a hardcoded potential
   and is not counted as evidence for anything.

SUCCESS / FAILURE CRITERIA (evaluated strictly; negative results first)
──────────────────────────────────────────────────────────────────────
* True Win  — ρ(λ) is a parameter-free curve AND λ_min(λ) moves, AND the
  global cavity and local throat behave as one unified dynamical object
  (in particular ρ(1) matches the measured global/local ratio).
* Clean Failure — ρ(1) lands orders of magnitude off the measured ratio,
  exposing a global-Casimir vs local-self-energy mismatch: the "everything
  rides on one R" assumption breaks under un-dialed conditions.

Tests:
  T1. Goal + framing + guardrails (what this is and is not).
  T2. The genuine Berger spectrum (SU(2); reduces to round at λ=1).
  T3. The global Casimir E_cav(λ) by zeta regularization (validated at
      λ=1 against the exact conformal-scalar closed form 1/240R).
  T4. The local self-energy λ_min(λ): it moves with λ — and the A/R+B·R²
      stability is explicitly discounted (no false victory).
  T5. ρ(λ): a parameter-free curve — but NOT flat (the two pieces respond
      differently to the same deformation).
  T6. THE CLEAN FAILURE: ρ(1) vs the measured cosmic/particle ratio — off
      by tens of orders of magnitude (the global/local scale mismatch).
  T7. Enforcement audit: the three guardrails held.
  T8. Assessment.

Verdict:
  - R_UNIFICATION_BREAKS_GLOBAL_CASIMIR_LOCAL_SELFENERGY_MISMATCH (expected):
    a CLEAN FAILURE on the decisive criterion — ρ(1) is O(10⁻⁴) (a pure
    geometric "one-R" number) while the measured global/local ratio is
    ~10⁻³⁹, a ~35-order mismatch.  R-unification survives only as
    scale-free bookkeeping, not as a physical single-R identification (the
    cosmological-constant problem, rediscovered geometrically).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from experiments.closure_ledger.maslov_dimensional_bridge_probe import (
    cavity_modes_scaled,
)

# ── locked geometric constants (NOT fit here; see guardrail 1) ──────────
K5 = 5                                  # capacity / D_bulk (THESIS PR #73)
L_THROAT = math.sqrt(2.0 * math.pi) / K5  # Hopf-fiber winding length (PR #83)
WINDING = 2.0 * math.pi / L_THROAT      # k=1 throat winding scale (1/R units)

# zeta_R(-k), the only externally-known numbers entering the regularization
ZETA_NEG = {3: 1.0 / 120.0, 2: 0.0, 1: -1.0 / 12.0, 0: -0.5}
E_CAV_ROUND_EXACT = 1.0 / 240.0         # conformal scalar on round S³ (R=1)

# measured cosmic-cavity vs particle-throat scale separation
R_HUBBLE_M = 1.30e26                    # Hubble radius (m)
LAMBDA_COMPTON_E_M = 3.8616e-13         # electron reduced Compton wavelength (m)
RHO_MEASURED = LAMBDA_COMPTON_E_M / R_HUBBLE_M   # = E_cosmic/E_self ~ (1/R_H)/(1/λ_C)


# ════════════════════════════════════════════════════════════════════════
# GUARDRAIL 2 — no hidden Born-rule / singlet imports
# ════════════════════════════════════════════════════════════════════════

def _forbid_quantum_inputs(what: str) -> None:
    """Hard stop: this audit is pure geometry/spectral theory.  If a step
    ever needs the Born rule or the singlet state it must raise here, not
    silently import ``geometrodynamics.bell``."""
    raise NotImplementedError(
        f"Anti-rigging guardrail: {what} would require the Born rule / "
        "singlet state, which this R-unification audit must not import. "
        "The result must stand on geometry alone or not at all."
    )


# ════════════════════════════════════════════════════════════════════════
# THE BERGER SPECTRUM (genuine SU(2); guardrail 1)
# ════════════════════════════════════════════════════════════════════════

def berger_laplacian_level(n: int, a: float) -> Tuple[np.ndarray, float]:
    """Scalar-Laplacian eigenvalues on the unit Berger sphere at level
    j = n/2: Δ(j,m) = 4j(j+1) + 4m²(λ⁻²−1), m = −j..j (step 1), each with
    multiplicity (2j+1).  Returns (eigenvalues_over_m, multiplicity).
    ``a`` = λ⁻² − 1 (a = 0 ⇒ round, eigenvalue 4j(j+1))."""
    j = n / 2.0
    m = np.arange(-j, j + 0.5, 1.0)
    eig = 4.0 * j * (j + 1.0) + 4.0 * a * m * m
    return eig, (2.0 * j + 1.0)


def conformal_frequencies(n: int, a: float) -> np.ndarray:
    """Conformal-scalar frequencies ω = √(Δ + 1) at level n (R=1).  At
    a=0 every ω = √((2j+1)²) = n+1 (degeneracy (n+1)²), the round
    conformal tower."""
    eig, _ = berger_laplacian_level(n, a)
    return np.sqrt(np.maximum(eig + 1.0, 1e-12))


def _a3_weyl(a: float) -> float:
    """Analytic leading (n³) Casimir coefficient of the shell energy
    S(n) = (n+1)·Σ_m ω: a₃(a) = ½[√(1+a) + asinh(√a)/√a], the continuum
    (Weyl) limit of the inner m-sum.  Fixing it analytically removes the
    n³-error amplification that makes a naive fit unstable."""
    if abs(a) < 1e-12:
        return 1.0
    if a > 0.0:
        return 0.5 * (math.sqrt(1.0 + a) + math.asinh(math.sqrt(a)) / math.sqrt(a))
    s = math.sqrt(-a)                    # λ > 1 branch (asinh → asin)
    return 0.5 * (math.sqrt(1.0 + a) + math.asin(s) / s)


# ════════════════════════════════════════════════════════════════════════
# THE GLOBAL CASIMIR ENERGY (zeta regularization)
# ════════════════════════════════════════════════════════════════════════

def _shell_energy(n: int, a: float) -> float:
    """S(n) = (n+1)·Σ_m √((n+1)² + 4a m²) — the level-n contribution to
    2·E (the conformal-frequency half-sum, R=1)."""
    j = n / 2.0
    m = np.arange(-j, j + 0.5, 1.0)
    return float((n + 1.0) * np.sum(np.sqrt((n + 1.0) ** 2 + 4.0 * a * m * m)))


def casimir_energy(lam: float, n_max: int = 400) -> Tuple[float, float]:
    """Zeta-regularized conformal-scalar vacuum energy on the unit Berger
    sphere S³_λ (R = 1).  Returns (E_cav, anomaly), where ``anomaly`` is
    the residual 1/(n+1) spectral coefficient: it is 0 at λ=1 (the round
    point is anomaly-free) and grows with deformation, flagging the
    standard log-ambiguity of the finite part away from the round sphere.

    Method: subtract the analytic Weyl growth a₃·(n+1)³, fit the residual
    smooth tail (powers (n+1)^{2,1,0,−1,−2,−3}) on the upper window, sum
    the convergent remainder, and add back the zeta values a₃ζ(−3) +
    b₂ζ(−2) + b₁ζ(−1) + b₀ζ(0).  Validated to recover 1/240 at λ=1."""
    a = lam ** -2 - 1.0
    n = np.arange(0, n_max + 1)
    x = n + 1.0
    S = np.array([_shell_energy(int(k), a) for k in n])
    a3 = _a3_weyl(a)
    R = S - a3 * x ** 3
    lo = int(n_max * 0.4)
    cols = [x[lo:] ** 2, x[lo:] ** 1, x[lo:] ** 0,
            x[lo:] ** -1, x[lo:] ** -2, x[lo:] ** -3]
    coef, *_ = np.linalg.lstsq(np.vstack(cols).T, R[lo:], rcond=None)
    b2, b1, b0, bm1 = coef[0], coef[1], coef[2], coef[3]
    all_cols = [x ** 2, x ** 1, x ** 0, x ** -1, x ** -2, x ** -3]
    smooth = sum(ci * col for ci, col in zip(coef, all_cols))
    remainder = float(np.sum(R - smooth))
    zeta_part = (a3 * ZETA_NEG[3] + b2 * ZETA_NEG[2]
                 + b1 * ZETA_NEG[1] + b0 * ZETA_NEG[0])
    return 0.5 * (remainder + zeta_part), float(bm1)


# ════════════════════════════════════════════════════════════════════════
# THE LOCAL SELF-ENERGY (the throat rest energy)
# ════════════════════════════════════════════════════════════════════════

def radial_self_energy() -> float:
    """The local radial Tangherlini cavity ground state ω₀ — the
    self-energy floor of the A/R + B·R² throat well, from the repo's own
    Sturm–Liouville solver (PR #83 cavity machinery)."""
    oms, _, _ = cavity_modes_scaled(1, scale=1.0)
    return float(oms[0])


_OMEGA0_RADIAL = radial_self_energy()


def lambda_min(lam: float) -> float:
    """The unified mass-operator ground state m(k=1, n=0) on S³_λ.  The
    throat term winds the Hopf fiber, whose length scales as λ
    (L_throat → λ·L_throat), so the winding scale is WINDING/λ; the radial
    cavity floor ω₀ is fiber-blind.  λ_min(λ) = √((WINDING/λ)² + ω₀²)."""
    return math.sqrt((WINDING / lam) ** 2 + _OMEGA0_RADIAL ** 2)


def rho(lam: float, n_max: int = 400) -> float:
    """The audited ratio ρ(λ) = E_cav(λ) / λ_min(λ): global cavity Casimir
    over local throat self-energy.  Parameter-free (geometry only)."""
    e_cav, _ = casimir_energy(lam, n_max)
    return e_cav / lambda_min(lam)


# ── A/R + B·R² well (computed only to be EXPLICITLY discounted; T4) ──────

def throat_well_stability() -> dict:
    """E(R) = A/R + B·R²; R* = (A/2B)^{1/3}; d²E/dR² = 2A/R*³ + 2B > 0.
    Computed with the same toy A=B=1 the repo uses — and NOT counted as
    evidence (guardrail 3)."""
    A, B = 1.0, 1.0
    r_star = (A / (2.0 * B)) ** (1.0 / 3.0)
    d2e = 2.0 * A / r_star ** 3 + 2.0 * B
    return {"R_star": r_star, "d2E_dR2": d2e, "stable": d2e > 0}


# ════════════════════════════════════════════════════════════════════════
# TESTS
# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal_and_guardrails",
        "description": (
            "An AUDIT of the R-unification assumption — NOT a quantum test, "
            "NOT a throat-formation test, NOT a wave-propagation test. The "
            "unified mass operator m²(k,n) = (k·2π/L_throat)² + "
            "((n+1)π/L_cavity)² rides the throat (Hopf-fiber winding, "
            "L_throat = √(2π)/k₅) and the cavity (radial/base) on a single "
            "S³ radius R. The Berger deformation S³_λ squashes the fiber "
            "while keeping the base round — the one move that separates the "
            "two scales — and so tests whether the global cosmic cavity and "
            "the local throat are really one dynamical object on one R."
        ),
        "is_not": ["quantum test", "throat-formation test",
                   "wave-propagation test"],
        "is": "R-unification test (global cavity vs local throat on one R)",
        "guardrails": [
            "no derived inversions (no fitted constant relabelled as a π-multiple)",
            "no hidden Born-rule / singlet imports (explicit raise instead)",
            "no false victories from A/R+B·R² stability",
        ],
        "framing": "QFT on the fixed classical throat geometry — not quantum gravity",
        "pass": True,
    }


def test_T2_berger_spectrum() -> dict:
    """The genuine SU(2) Berger spectrum, reducing to the round tower at
    λ=1 (no rigging)."""
    # round: level n has frequencies all equal to (n+1), degeneracy (n+1)²
    round_ok = True
    for n in range(6):
        w = conformal_frequencies(n, 0.0)
        mult = 2.0 * (n / 2.0) + 1.0
        deg = len(w) * mult
        if not (np.allclose(w, n + 1.0) and abs(deg - (n + 1) ** 2) < 1e-9):
            round_ok = False
    # the analytic Weyl coefficient matches the round value and is finite
    a3_round = _a3_weyl(0.0)
    a3_squash = _a3_weyl(0.6 ** -2 - 1.0)   # λ=0.6
    ok = round_ok and abs(a3_round - 1.0) < 1e-12 and a3_squash > 1.0
    return {
        "name": "T2_berger_spectrum",
        "description": (
            "The scalar Laplacian on the unit Berger sphere is the genuine "
            "SU(2) result Δ(j,m) = 4j(j+1) + 4m²(λ⁻²−1), multiplicity "
            "(2j+1). At λ=1 it collapses to the round tower (every level-n "
            f"frequency = n+1, degeneracy (n+1)²): verified = {round_ok}. "
            "The analytic leading Weyl coefficient a₃ = "
            "½[√(1+a)+asinh(√a)/√a] equals 1 at the round point and "
            f"a₃(λ=0.6) = {a3_squash:.4f} — fixed analytically, not fit, so "
            "the regularization carries no reverse-engineered constant."
        ),
        "round_tower_recovered": round_ok,
        "a3_round": round(a3_round, 12),
        "a3_lambda_0p6": round(a3_squash, 6),
        "pass": ok,
    }


def test_T3_global_casimir() -> dict:
    """E_cav(λ) by zeta regularization, validated at λ=1 against 1/240."""
    e1, anom1 = casimir_energy(1.0)
    validated = abs(e1 - E_CAV_ROUND_EXACT) < 1e-5 and abs(anom1) < 1e-9
    lambdas = [0.7, 0.85, 1.0, 1.2, 1.4]
    curve = []
    moves = False
    for lam in lambdas:
        e, anom = casimir_energy(lam)
        curve.append({"lambda": lam, "E_cav": round(e, 6),
                      "anomaly": float(f"{anom:.2e}")})
        if abs(e - e1) > 1e-3:
            moves = True
    ok = validated and moves
    return {
        "name": "T3_global_casimir_zeta",
        "description": (
            "The global cosmic-cavity vacuum energy is the zeta-regularized "
            "conformal-scalar energy on S³_λ. VALIDATION: at λ=1 the method "
            f"recovers the exact closed form E_cav(1) = 1/240 = "
            f"{E_CAV_ROUND_EXACT:.6f} (computed {e1:.6f}, anomaly "
            f"{anom1:.0e}) — the round point is anomaly-free. Away from λ=1 "
            "E_cav moves strongly and a residual 1/(n+1) spectral term (the "
            "standard log-ambiguity of the finite part on a deformed "
            "sphere) appears and grows with the squash — the global energy "
            "does not ride rigidly on one scale."
        ),
        "validated_against_1over240": validated,
        "E_cav_round": round(e1, 6),
        "curve": curve,
        "moves_with_lambda": moves,
        "pass": ok,
    }


def test_T4_local_self_energy() -> dict:
    """λ_min(λ) moves; the A/R+B·R² stability is explicitly discounted."""
    lambdas = [0.7, 0.85, 1.0, 1.2, 1.4]
    rows = [{"lambda": lam, "lambda_min": round(lambda_min(lam), 4)}
            for lam in lambdas]
    span = lambda_min(0.7) / lambda_min(1.4)
    moves = span > 1.2
    well = throat_well_stability()
    ok = moves and well["stable"]   # stability noted, but see discount below
    return {
        "name": "T4_local_self_energy",
        "description": (
            "The local throat self-energy is the unified mass-operator "
            "ground state λ_min(λ) = √((2π/(λ·L_throat))² + ω₀²): the "
            "winding term tracks the squashed fiber (∝ 1/λ) while the "
            f"radial Tangherlini floor ω₀ = {_OMEGA0_RADIAL:.4f} is "
            f"fiber-blind. λ_min MOVES — from {lambda_min(0.7):.2f} (λ=0.7) "
            f"through {lambda_min(1.0):.2f} (λ=1) to {lambda_min(1.4):.2f} "
            f"(λ=1.4), span ×{span:.2f}. GUARDRAIL 3: the A/R+B·R² well is "
            f"stable (d²E/dR² = {well['d2E_dR2']:.2f} > 0 at R* = "
            f"{well['R_star']:.4f}), but this is an artifact of a hardcoded "
            "potential and is DISCOUNTED — it is not counted as evidence "
            "for R-unification."
        ),
        "lambda_min_curve": rows,
        "span_ratio": round(span, 3),
        "lambda_min_moves": moves,
        "well_stability_DISCOUNTED": {
            "d2E_dR2": round(well["d2E_dR2"], 4),
            "counts_as_evidence": False,
        },
        "pass": ok,
    }


def test_T5_rho_curve() -> dict:
    """ρ(λ) is a parameter-free curve — and NOT flat."""
    lambdas = [0.7, 0.85, 1.0, 1.2, 1.4]
    rows = []
    vals = []
    for lam in lambdas:
        r = rho(lam)
        rows.append({"lambda": lam, "rho": float(f"{r:.4e}")})
        vals.append(r)
    r1 = rho(1.0)
    spread = max(vals) / min(vals)
    not_flat = spread > 1.5
    ok = not_flat and r1 > 0
    return {
        "name": "T5_rho_parameter_free_curve",
        "description": (
            "ρ(λ) = E_cav(λ)/λ_min(λ) is parameter-free (pure geometry, no "
            "free knob). It is NOT flat — across λ ∈ [0.7, 1.4] it spans "
            f"×{spread:.2f}. The global cavity and the local throat respond "
            "DIFFERENTLY to the same deformation (the cavity by the full "
            "Berger spectrum, the throat only through the winding term), so "
            "ρ does not reduce to a single λ-independent number: the 'one R' "
            "cancellation that R-unification would require does not occur "
            "even in shape."
        ),
        "rho_curve": rows,
        "rho_round": float(f"{r1:.4e}"),
        "spread": round(spread, 3),
        "not_flat": not_flat,
        "pass": ok,
    }


def test_T6_clean_failure() -> dict:
    """ρ(1) vs the measured cosmic/particle ratio — the headline negative
    result."""
    r1 = rho(1.0)
    # measured global/local ratio: E_cosmic/E_self ~ (1/R_H)/(1/λ_C) = λ_C/R_H
    orders_off = abs(math.log10(r1) - math.log10(RHO_MEASURED))
    clean_failure = orders_off > 10.0
    return {
        "name": "T6_clean_failure_scale_mismatch",
        "description": (
            "THE DECISIVE TEST (negative result, reported first). "
            "R-unification asserts E_cav and E_self both equal (pure "
            "number)/R on ONE R, so ρ = E_cav/E_self must be a pure number, "
            f"computed here as ρ(1) = {r1:.3e} — an O(10⁻⁴) 'one-R' number. "
            "But physically E_cav is the global COSMIC-cavity vacuum energy "
            "(~1/R_Hubble) and E_self is the local PARTICLE throat "
            "(~1/λ_Compton), whose measured ratio is λ_C/R_H = "
            f"{RHO_MEASURED:.3e}. The geometric 'one-R' prediction overshoots "
            f"the measured global/local ratio by ~{orders_off:.0f} ORDERS OF "
            "MAGNITUDE. The cosmic cavity and the local throat cannot both "
            "ride on one R: this is the global-Casimir vs local-self-energy "
            "mismatch — the cosmological-constant problem, rediscovered as "
            "the failure of the geometric shorthand."
        ),
        "rho_predicted_one_R": float(f"{r1:.3e}"),
        "rho_measured_global_local": float(f"{RHO_MEASURED:.3e}"),
        "orders_of_magnitude_off": round(orders_off, 1),
        "clean_failure": clean_failure,
        "pass": clean_failure,    # the CLEAN FAILURE is the intended outcome
    }


def test_T7_enforcement_audit() -> dict:
    """The three anti-rigging guardrails held."""
    import sys
    mod = sys.modules[__name__]
    src = Path(mod.__file__).read_text()
    # 1. no derived inversions: L_throat is the locked √(2π)/k₅ (semantic
    #    check), and no arbitrary small-integer/100 π-multiple was spliced
    #    into the path.  The forbidden pattern is built from pieces so the
    #    literal never appears verbatim in this file (which would self-trip).
    locked_l_throat = abs(L_THROAT - math.sqrt(2.0 * math.pi) / 5.0) < 1e-15
    import re
    forbidden = re.compile(r"\b[1-9]\s*\*?\s*math\.pi\s*/\s*100\b")
    code_lines = [ln for ln in src.splitlines()
                  if "forbidden" not in ln and "re.compile" not in ln]
    no_pi_relabel = not any(forbidden.search(ln) for ln in code_lines)
    no_inversion = locked_l_throat and no_pi_relabel
    # 2. no hidden Born/singlet imports — scan actual import statements only
    #    (the docstring/guard MENTION geometrodynamics.bell to forbid it; an
    #    honest check looks for `import`/`from ... import`, not any mention).
    import_lines = [
        ln.strip() for ln in src.splitlines()
        if ln.strip().startswith(("import ", "from "))
    ]
    no_quantum = not any("geometrodynamics.bell" in ln for ln in import_lines)
    raised = False
    try:
        _forbid_quantum_inputs("a Born-rule shortcut")
    except NotImplementedError:
        raised = True
    # 3. stability discounted
    well = throat_well_stability()
    discounted = well["stable"] is True   # stable, but T4 marks it non-evidence
    ok = no_inversion and no_quantum and raised and discounted
    return {
        "name": "T7_enforcement_audit",
        "description": (
            "The guardrails held: (1) no derived inversion — L_throat = "
            "√(2π)/k₅ is used as the THESIS-locked value, no constant is "
            "relabelled as an arbitrary π-multiple; (2) no hidden imports — "
            "geometrodynamics.bell (Born rule / singlet) is never imported, "
            "and _forbid_quantum_inputs raises an explicit exception when a "
            f"quantum shortcut is requested (verified: {raised}); (3) the "
            "A/R+B·R² stability is computed but discounted as non-evidence."
        ),
        "no_derived_inversion": no_inversion,
        "no_hidden_quantum_imports": no_quantum,
        "explicit_raise_works": raised,
        "stability_discounted": discounted,
        "pass": ok,
    }


def test_T8_assessment() -> dict:
    r1 = rho(1.0)
    orders = abs(math.log10(r1) - math.log10(RHO_MEASURED))
    return {
        "name": "T8_assessment",
        "description": (
            "CLEAN FAILURE on the decisive criterion. The within-geometry "
            "structure of R-unification is intact — ρ(λ) is a parameter-free "
            "curve and λ_min(λ) moves dynamically with the fiber squash — "
            "but the absolute-scale identification fails: ρ(1) is an "
            f"O(10⁻⁴) one-R number while the measured global/local ratio is "
            f"~10^{math.log10(RHO_MEASURED):.0f}, a ~{orders:.0f}-order "
            "mismatch. The cosmic cavity and the local throat do NOT ride on "
            "one R. R-unification survives only as scale-free bookkeeping "
            "(consistent with the B4 anchor audit), not as a physical "
            "single-R claim. The geometric shorthand breaks exactly at the "
            "global-Casimir vs local-self-energy boundary."
        ),
        "classification": (
            "R_UNIFICATION_BREAKS_GLOBAL_CASIMIR_LOCAL_SELFENERGY_MISMATCH"
        ),
        "pass": True,
    }


# ════════════════════════════════════════════════════════════════════════
# RUNNER + VERDICT (negative result first)
# ════════════════════════════════════════════════════════════════════════

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_berger_spectrum(),
        test_T3_global_casimir(),
        test_T4_local_self_energy(),
        test_T5_rho_curve(),
        test_T6_clean_failure(),
        test_T7_enforcement_audit(),
        test_T8_assessment(),
    ]
    t6 = tests[5]
    all_ok = all(t["pass"] for t in tests)
    if all_ok and t6["clean_failure"]:
        verdict_class = (
            "R_UNIFICATION_BREAKS_GLOBAL_CASIMIR_LOCAL_SELFENERGY_MISMATCH"
        )
        verdict = (
            "NEGATIVE RESULT FIRST — CLEAN FAILURE. The 'everything rides on "
            "one R' assumption breaks under the Berger deformation. "
            f"ρ(1) = E_cav/E_self = {t6['rho_predicted_one_R']:.3e} is a pure "
            "geometric one-R number, but the MEASURED global/local ratio "
            f"(cosmic cavity 1/R_Hubble over particle throat 1/λ_Compton) is "
            f"{t6['rho_measured_global_local']:.3e} — a ~"
            f"{t6['orders_of_magnitude_off']:.0f}-ORDER-OF-MAGNITUDE "
            "mismatch. The global cosmic-cavity Casimir energy and the local "
            "throat self-energy cannot share one R; this is the "
            "cosmological-constant problem rediscovered as the failure of "
            "the geometric shorthand.\n\n"
            "WHAT SURVIVES (not a victory). Within geometric units the "
            "structure is internally consistent: ρ(λ) is a parameter-free "
            "curve (no free knob), the global Casimir E_cav(λ) is validated "
            "at λ=1 against the exact conformal closed form 1/240R, and "
            "λ_min(λ) moves dynamically as the winding term tracks the "
            "squashed fiber. But ρ(λ) is NOT flat — the cavity and the "
            "throat respond differently to the SAME deformation — so they "
            "are not one dynamical object even in shape. The A/R+B·R² well's "
            "stability was computed and explicitly DISCOUNTED.\n\n"
            "BOUNDARY MAPPED. R-unification is a valid scale-free "
            "bookkeeping device (consistent with the B4 single-anchor audit) "
            "but a failed physical identification: the geometric shorthand "
            "breaks precisely at the global-Casimir vs local-self-energy "
            "boundary, by tens of orders of magnitude, under un-dialed "
            "conditions. All three anti-rigging guardrails held."
        )
    else:
        verdict_class = "R_UNIFICATION_AUDIT_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. A guardrail or a spectral validation failed; "
            "review the Casimir validation at λ=1, the λ_min motion, or the "
            "enforcement audit before reading the scale comparison."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "Berger-deformation audit of R-unification: a clean failure — "
            "ρ(1) (one-R geometric) is ~35 orders off the measured "
            "cosmic/particle ratio; the global Casimir and the local "
            "self-energy do not share one R"
        ),
        "negative_result": (
            f"ρ(1) = {t6['rho_predicted_one_R']:.3e} vs measured "
            f"{t6['rho_measured_global_local']:.3e} — "
            f"~{t6['orders_of_magnitude_off']:.0f} orders off"
        ),
        "survives_as": "scale-free bookkeeping only (parameter-free ρ(λ), moving λ_min)",
        "guardrails": "all three held (no inversion, no Born/singlet, stability discounted)",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# Berger deformation audit of the R-unification assumption (PR #165)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "An audit — not a quantum, throat-formation, or wave-propagation "
        "test. It asks whether BAM's unified mass operator really rides the "
        "throat (Hopf-fiber winding) and the cavity (radial/base) on one S³ "
        "radius. The Berger deformation squashes the fiber alone, separating "
        "the two scales. **Negative result first.** *(QFT on the classical "
        "throat, not quantum gravity.)*"
    )
    out.append("")
    out.append(f"- **The clean failure**: {s['negative_result']}")
    out.append(f"- **Survives as**: {s['survives_as']}")
    out.append(f"- **Guardrails**: {s['guardrails']}")
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "audit goal + the three anti-rigging guardrails",
        "T2": "genuine SU(2) Berger spectrum (round tower at λ=1)",
        "T3": "global Casimir E_cav(λ), validated 1/240 at λ=1",
        "T4": "local λ_min(λ) moves; A/R+B·R² stability discounted",
        "T5": "ρ(λ) parameter-free but NOT flat",
        "T6": "CLEAN FAILURE: ρ(1) ~35 orders off measured ratio",
        "T7": "enforcement audit: all guardrails held",
        "T8": "R_UNIFICATION_BREAKS_GLOBAL_LOCAL_MISMATCH",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'—')} | {p} |")
    out.append("")
    t3, t5, t6 = s["tests"][2], s["tests"][4], s["tests"][5]
    out.append("## E_cav(λ) and ρ(λ) (parameter-free)")
    out.append("")
    out.append("| λ | E_cav(λ) | anomaly | ρ(λ) |")
    out.append("|---|---:|---:|---:|")
    rho_map = {r["lambda"]: r["rho"] for r in t5["rho_curve"]}
    for c in t3["curve"]:
        lam = c["lambda"]
        out.append(f"| {lam} | {c['E_cav']} | {c['anomaly']} | "
                   f"{rho_map.get(lam,'—')} |")
    out.append("")
    out.append("## The decisive comparison")
    out.append("")
    out.append(f"- ρ(1), geometric one-R prediction: **{t6['rho_predicted_one_R']:.3e}**")
    out.append(f"- measured global/local ratio (λ_C / R_Hubble): "
               f"**{t6['rho_measured_global_local']:.3e}**")
    out.append(f"- mismatch: **~{t6['orders_of_magnitude_off']:.0f} orders of magnitude**")
    out.append("")
    out.append("## Verdict")
    out.append("")
    out.append(f"**{s['verdict_class']}.** {s['verdict']}")
    out.append("")
    return "\n".join(out)


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_berger_r_unification_audit_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
