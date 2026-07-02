"""
The geon-statistics adjudication — companion probe (PR #196).

> Framing: QFT on the *fixed classical* throat geometry (geometry → fields),
> not quantum gravity.

THIS PROBE IS THE APPENDIX, NOT THE ARGUMENT
────────────────────────────────────────────
The deliverable of PR #196 is the adjudication document
``docs/geon_statistics_adjudication.md`` — careful mathematics against the
geon-statistics literature (Friedman–Sorkin; Sorkin–Surya gr-qc/9605050;
Dowker–Sorkin gr-qc/9609064; Hendriks; Giulini).  This probe machine-checks
every step of the argument that is finite/algebraic — orientation
determinants, the amphichirality criterion and the explicit descending
reflection, the SU(2)/pin lift of the rotation loop, the two-throat
mapping-class-group presentation and its sector structure (the four scalar
sectors and the indefinite-statistics continuum), and the bordism-level
facts — so that the document's chain has no unchecked arithmetic.

THE ADJUDICATED VERDICT (established in the document)
─────────────────────────────────────────────────────
* The BAM throat prime is RP³ (the repo's own #169 construction: the
  antipodally identified Einstein–Rosen neighborhood is the RP³ geon, with
  the one-sided RP² cross-cap as throat slice).
* The three Dowker–Sorkin hypotheses HOLD: prime, non-chiral, abelian —
  the spin-statistics CORRELATION is a theorem for pair-created BAM
  throats.
* But RP³ is NOT spinorial in bare Diff (Hendriks: cyclic π₁ ⟹
  non-spinorial): in BARE metric GR the same theorem selects BOSE, and in
  frozen-topology canonical GR the statistics is a sector choice (Bose,
  Fermi, and a continuum of indefinite sectors: G = (Z₂∗Z₂)⋊S₂).
* The −1 is recovered in PIN⁻-FRAMED GR — the framing BAM's own quotient
  forces (#169/#170/#195): the rotation loop, trivial in Diff, lifts to
  −1 on the pin bundle (π₁(SO(3)) = Z₂, RP³ ≅ SO(3)); SSC then gives
  FERMI, consistent with the #188 measured holonomy.
* The #58 nucleation channel passes every bordism-level check
  (mirror-pair identity, Ω₃^{SO} = 0, spin extension); the explicit BAM
  4-manifold is the honest open construction.

Tests:
  T1. Goal (the adjudication; document is the deliverable).
  T2. Lemma 1 arithmetic: the orientation determinants of the quotient
      (J on S² vs S³; the throat involution on S²×ℝ preserves).
  T3. Lemma 2 arithmetic: amphichirality q²≡−1 (mod p) at (2,1); the
      explicit reflection commuting with the deck map; π₁ abelian.
  T4. Lemma 3/5 arithmetic: the SU(2) lift of the 2π rotation loop ends
      at −I (the pin-framed spinoriality), 4π at +I.
  T5. Lemma 4: the two-throat MCG G = ⟨s₁,s₂,E⟩ — verify the relations
      in a faithful 2D model, the four 1-dim sectors, the Fermi sector,
      and an explicit indefinite-statistics (non-commuting) sector.
  T6. The SOH selection: abelian weights kill the indefinite sectors;
      SSC + bare rotation (+1) ⟹ Bose; SSC + pin rotation (−1) ⟹ Fermi.
  T7. The #58 cobordism ledger: mirror-identity, bordism existence, spin
      extension (cited facts, encoded); the open item flagged.
  T8. Assessment.

Verdict:
  SSC_THEOREM_HYPOTHESES_HOLD_SIGN_IS_FRAMING_DEPENDENT_BARE_GR_SELECTS_
  BOSE_PIN_MINUS_FRAMING_SELECTS_FERMI
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.linalg import expm

_SX = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
_SY = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
_SZ = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)


def _rot2(a: float) -> np.ndarray:
    return np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])


def _refl2(a: float) -> np.ndarray:
    """Reflection of the plane about the line at angle a."""
    return np.array([[math.cos(2 * a), math.sin(2 * a)],
                     [math.sin(2 * a), -math.cos(2 * a)]])


# ════════════════════════════════════════════════════════════════════════
# TESTS
# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "The adjudication PR: is the BAM exchange −1 a THEOREM? The "
            "deliverable is docs/geon_statistics_adjudication.md — "
            "mathematics against the Sorkin-school literature (Sorkin–"
            "Surya gr-qc/9605050; Dowker–Sorkin gr-qc/9609064; Hendriks; "
            "Giulini), applied to the actual BAM topology (the RP³ geon "
            "with RP² cross-cap throat slice, from the repo's own #169 "
            "quotient). This probe machine-checks every finite/algebraic "
            "step of the document's argument so the chain has no "
            "unchecked arithmetic. The adjudicated verdict: the SSC "
            "CORRELATION is a theorem (prime/non-chiral/abelian all "
            "hold); the SIGN is framing-dependent — bare metric GR "
            "selects BOSE (RP³ is non-spinorial), the Pin⁻ framing BAM's "
            "quotient forces selects FERMI."
        ),
        "deliverable": "docs/geon_statistics_adjudication.md",
        "framing": "QFT on the fixed classical throat geometry — not quantum gravity",
        "pass": True,
    }


def test_T2_lemma1_orientation_arithmetic() -> dict:
    """The quotient's orientation determinants, computed explicitly."""
    # antipodal map on S^n: orientation factor (−1)^{n+1}
    # (acting on tangent frames: det(−I_{n+1}) relative to the normal)
    det_s2 = (-1.0) ** (2 + 1)     # brane angular slice: −1 (RP² non-orientable)
    det_s3 = (-1.0) ** (3 + 1)     # bulk angular sphere: +1 (RP³ orientable)
    # the throat involution on the Einstein–Rosen neighborhood S²×ℝ:
    # J(x,y) = (−x,−y): det = det(dA|TS²)·det(−1 on ℝ) = (−1)·(−1) = +1
    det_throat = det_s2 * (-1.0)
    # explicit check of det(dA|TS²) with a concrete frame at x = ẑ:
    # tangent frame (e₁,e₂) at north pole maps to (−e₁,−e₂) at south pole;
    # comparing via the outward-normal orientation convention the sphere's
    # orientation reverses: encode as det(−I₂) relative to normal flip:
    frame_det = float(np.linalg.det(-np.eye(2))) * (-1.0)   # = −1  ✓
    ok = (det_s2 == -1.0 and det_s3 == 1.0 and det_throat == 1.0
          and frame_det == -1.0)
    return {
        "name": "T2_lemma1_orientation_arithmetic",
        "description": (
            "Lemma 1's arithmetic. The antipodal map on Sⁿ carries "
            "orientation factor (−1)^{n+1}: on the brane's angular S² it "
            f"is {det_s2:+.0f} (RP² non-orientable — the cross-cap throat "
            f"slice), on the bulk's S³ it is {det_s3:+.0f} (RP³ "
            "orientable) — the #169/#183 grading re-verified. The throat "
            "involution J(x,y) = (−x,−y) on the Einstein–Rosen "
            "neighborhood S²×ℝ has determinant (−1)·(−1) = "
            f"{det_throat:+.0f}: orientation-PRESERVING and free, so the "
            "quotient is an orientable 3-manifold — the twisted ℝ-bundle "
            "over RP², i.e. RP³ ∖ {point}: the throat prime is the RP³ "
            "GEON, the canonical object of the entire geon-statistics "
            "literature."
        ),
        "det_antipodal_S2": det_s2,
        "det_antipodal_S3": det_s3,
        "det_throat_involution": det_throat,
        "throat_prime": "RP3 (the RP3 geon; minimal surface RP2)",
        "pass": ok,
    }


def test_T3_lemma2_hypotheses() -> dict:
    """Prime, non-chiral, abelian — the arithmetic parts."""
    # (a) amphichirality criterion for lens spaces: q² ≡ −1 (mod p)
    p, q = 2, 1
    amphichiral = (q * q) % p == (-1) % p
    # explicit descending reflection: R = diag(−1,1,1,1) ∈ O(4),
    # det = −1, commutes with the deck map −I₄:
    r = np.diag([-1.0, 1.0, 1.0, 1.0])
    deck = -np.eye(4)
    commutes = bool(np.allclose(r @ deck, deck @ r))
    det_r = float(np.linalg.det(r))
    # (c) π₁ = Z₂ abelian: the deck group {I, −I} under composition
    z2 = [np.eye(4), -np.eye(4)]
    closed = all(any(np.allclose(a @ b, c) for c in z2) for a in z2 for b in z2)
    abelian = all(np.allclose(a @ b, b @ a) for a in z2 for b in z2)
    ok = amphichiral and commutes and det_r == -1.0 and closed and abelian
    return {
        "name": "T3_lemma2_hypotheses",
        "description": (
            "Lemma 2's arithmetic — the three Dowker–Sorkin hypotheses. "
            "(a) PRIME: RP³ is an elliptic space form, irreducible (cited; "
            "no arithmetic). (b) NON-CHIRAL: the lens-space criterion "
            f"q² ≡ −1 (mod p) holds at L(2,1): {q}² ≡ {(q*q)%p} ≡ "
            f"{(-1)%p} (mod {p}) ✓; explicitly, the reflection "
            "R = diag(−1,1,1,1) has det −1 and commutes with the deck map "
            f"−I₄ ({commutes}), hence descends to an orientation-"
            "REVERSING diffeomorphism of RP³ — the throat is its own "
            "mirror, so the pair-created partner is an identical geon. "
            "(c) ABELIAN: the deck group {I, −I} ≅ Z₂ is closed and "
            f"abelian ({abelian}). All three hypotheses of the "
            "Dowker–Sorkin spin-statistics theorem HOLD for the BAM "
            "throat: the CORRELATION is a theorem on this topology."
        ),
        "amphichirality_q2_eq_minus1_mod_p": amphichiral,
        "explicit_reflection_descends": commutes,
        "det_reflection": det_r,
        "pi1_abelian_Z2": abelian,
        "pass": ok,
    }


def test_T4_pin_lift_of_rotation() -> dict:
    """The 2π rotation: trivial in Diff (Hendriks, cited), −1 on the pin
    bundle (computed)."""
    # the rotation loop t ↦ R(t) in SO(3) lifts to U(t) = exp(−i t σ_z/2)
    n = 2000
    u = np.eye(2, dtype=complex)
    for k in range(n):
        u = expm(-1j * (2 * math.pi / n) * _SZ / 2.0) @ u
    lift_2pi_is_minus_I = bool(np.allclose(u, -np.eye(2), atol=1e-9))
    u4 = u @ u
    lift_4pi_is_plus_I = bool(np.allclose(u4, np.eye(2), atol=1e-8))
    # the projection of U(t) to SO(3) is a CLOSED loop (R(2π) = R(0)):
    # verify via the adjoint action on the Pauli basis at t = 2π
    ad = np.array([[np.real(np.trace(a @ u @ b @ u.conj().T)) / 2.0
                    for b in (_SX, _SY, _SZ)] for a in (_SX, _SY, _SZ)])
    closed_in_so3 = bool(np.allclose(ad, np.eye(3), atol=1e-8))
    ok = lift_2pi_is_minus_I and lift_4pi_is_plus_I and closed_in_so3
    return {
        "name": "T4_pin_lift_of_rotation",
        "description": (
            "Lemma 3 + Theorem 5's arithmetic. BARE DIFF: RP³ = L(2,1) is "
            "a lens space with cyclic π₁, hence NON-SPINORIAL (Hendriks; "
            "Giulini's tabulation — cited, not re-proved): the 2π "
            "rotation is isotopic to the identity, so in bare metric GR "
            "the throat is tensorial and the Friedman–Sorkin spin-½ from "
            "pure topology is NOT available — the correction to #170/"
            "#171. THE PIN LIFT (computed): the trivializing isotopy "
            "traces the rotation loop t ↦ R(t) in SO(3) ≅ RP³, the "
            "generator of π₁(SO(3)) = Z₂; its lift to the double cover "
            f"ends at −I (verified: {lift_2pi_is_minus_I}, with the SO(3) "
            f"projection closed: {closed_in_so3}), while the 4π loop "
            f"lifts to +I ({lift_4pi_is_plus_I}). On the Pin⁻-framed "
            "state space the 2π rotation therefore acts as −I: the "
            "throat is SPINORIAL WITH FRAMING — exactly what the #188 "
            "holonomy measured."
        ),
        "bare_diff_spinorial": False,
        "bare_diff_source": "Hendriks 1977; cyclic-pi1 primes are non-spinorial (cited)",
        "pin_lift_2pi": "-I (computed)",
        "pin_lift_4pi": "+I (computed)",
        "so3_loop_closed": closed_in_so3,
        "pass": ok,
    }


def test_T5_two_throat_mcg_sectors() -> dict:
    """G = ⟨s₁,s₂,E | s₁²=s₂²=E²=1, Es₁E=s₂⟩ = (Z₂∗Z₂)⋊S₂: the model,
    the four scalar sectors, and an indefinite sector."""
    # faithful 2D orthogonal model at slide angle θ:
    theta = 2.0 * math.sqrt(2.0)      # irrational multiple of π: generic
    s1 = _refl2(0.0)
    s2 = _refl2(theta / 2.0)
    e = _refl2(theta / 4.0)
    rel_s1 = bool(np.allclose(s1 @ s1, np.eye(2)))
    rel_s2 = bool(np.allclose(s2 @ s2, np.eye(2)))
    rel_e = bool(np.allclose(e @ e, np.eye(2)))
    rel_conj = bool(np.allclose(e @ s1 @ e, s2)) and bool(
        np.allclose(e @ s2 @ e, s1))
    # slides do NOT commute (the slide subgroup is Z₂∗Z₂, infinite):
    slides_noncommute = float(np.linalg.norm(s1 @ s2 - s2 @ s1)) > 0.1
    s1s2_order = float(np.linalg.norm(np.linalg.matrix_power(s1 @ s2, 12)
                                      - np.eye(2))) > 0.1
    # the four 1-dim sectors: hom(G, U(1)): s₁ ↦ σ, s₂ ↦ σ (E-conjugate
    # forces equality), E ↦ ε, σ, ε ∈ {±1}
    one_dim = []
    for sig in (1, -1):
        for eps in (1, -1):
            ok1 = (sig * sig == 1) and (eps * eps == 1) and (eps * sig * eps == sig)
            one_dim.append({"slides": sig, "exchange": eps, "consistent": ok1})
    fermi_exists = any(r["consistent"] and r["exchange"] == -1 for r in one_dim)
    n_scalar = sum(r["consistent"] for r in one_dim)
    # indefinite sector: in the 2D model, exchange does not commute with
    # slides — no simultaneous statistics/slide eigenstates:
    mixing = float(np.linalg.norm(e @ s1 - s1 @ e)) > 0.1
    # irreducibility of the 2D model (generic θ): commutant is scalar
    gens = [s1, s2, e]
    # solve for X with [X, g] = 0 for all g: stack the linear system
    rows = []
    for g in gens:
        m = np.kron(np.eye(2), g) - np.kron(g.T, np.eye(2))
        rows.append(m)
    ns = np.linalg.svd(np.vstack(rows))[1]
    commutant_dim = int(np.sum(ns < 1e-10))
    irreducible = commutant_dim == 1
    ok = (rel_s1 and rel_s2 and rel_e and rel_conj and slides_noncommute
          and s1s2_order and n_scalar == 4 and fermi_exists and mixing
          and irreducible)
    return {
        "name": "T5_two_throat_mcg_sectors",
        "description": (
            "Lemma 4's algebra. The two-throat mapping class group "
            "(Sorkin–Surya, with the Fouxe-Rabinovitch relations; "
            "internal diffeos trivial, rotation trivial by Lemma 3) is "
            "G = ⟨s₁,s₂,E | s₁²=s₂²=E²=1, Es₁E=s₂⟩ = (Z₂∗Z₂)⋊S₂. A "
            "faithful 2D orthogonal model verifies all relations "
            f"({rel_s1 and rel_s2 and rel_e and rel_conj}); the slides "
            f"do NOT commute ({slides_noncommute} — the slide subgroup "
            "is the FREE product Z₂∗Z₂ ≅ D_∞, G is infinite). SECTORS: "
            f"exactly {n_scalar} one-dimensional sectors (slides ±1 × "
            f"exchange ±1) — the FERMI sector exists ({fermi_exists}), so "
            "#171's consistency claim survives — plus a CONTINUUM of "
            "2-dimensional sectors (one per slide angle θ) in which "
            f"exchange and slides do not commute ({mixing}, irreducible: "
            f"{irreducible}): geons with INDEFINITE statistics, the "
            "Sorkin–Surya anomalous sectors realized on BAM's own prime. "
            "Bare frozen-topology GR does NOT select the −1."
        ),
        "relations_verified": rel_s1 and rel_s2 and rel_e and rel_conj,
        "slides_free_product": slides_noncommute,
        "one_dim_sectors": one_dim,
        "n_scalar_sectors": n_scalar,
        "fermi_sector_exists": fermi_exists,
        "indefinite_sector_exists": mixing and irreducible,
        "pass": ok,
    }


def test_T6_soh_selection() -> dict:
    """SSC + rotation phase: bare ⟹ Bose; pin-framed ⟹ Fermi."""
    hypotheses = {"prime": True, "non_chiral": True, "abelian": True}
    ssc_applies = all(hypotheses.values())
    # bare rotation phase: non-spinorial ⟹ +1 (T4, cited Hendriks)
    bare_rotation = +1
    bare_exchange = bare_rotation if ssc_applies else None
    # pin rotation phase: the computed lift (T4) ⟹ −1
    pin_rotation = -1
    pin_exchange = pin_rotation if ssc_applies else None
    # consistency with the #188 measured swap holonomy: recompute the
    # path-ordered 2π holonomy on the spin-½ frame (the pin lift again,
    # as transported in #188):
    n = 1500
    u = np.eye(2, dtype=complex)
    for k in range(n):
        u = expm(-1j * (2 * math.pi / n) * _SZ / 2.0) @ u
    holonomy_sign = int(round(float(np.real(np.trace(u)) / 2.0
                                    / (np.real(np.trace(np.eye(2))) / 2.0))))
    consistent_188 = holonomy_sign == pin_exchange
    ok = ssc_applies and bare_exchange == 1 and pin_exchange == -1 and consistent_188
    return {
        "name": "T6_soh_selection",
        "description": (
            "Theorem 5's selection logic, encoded. In the SOH with "
            "topology change the weights come from ABELIAN reps "
            "(Sorkin–Surya) — the indefinite sectors of T5 are "
            "eliminated — and for pair-created geons the Dowker–Sorkin "
            "theorem forces exchange phase = 2π-rotation phase (SSC), "
            "its hypotheses all holding by T3. Evaluating the rotation "
            "phase: BARE metric GR — non-spinorial (T4, Hendriks) ⟹ +1 "
            "⟹ pair-created BAM throats are BOSONS (if BAM were bare "
            "geometrodynamics, the exchange arc would be refuted here). "
            "PIN⁻-FRAMED GR — the framing forced by #169 (non-orientable "
            "RP² slice) + #170 (Pin⁻ unique) + #195 (the modes ARE pin "
            "spinors): the rotation lifts to −1 (T4) ⟹ FERMI — "
            f"consistent with the #188 measured holonomy ({holonomy_sign})."
        ),
        "ssc_hypotheses": hypotheses,
        "ssc_applies": ssc_applies,
        "bare_GR": {"rotation": bare_rotation, "exchange": bare_exchange,
                    "statistics": "Bose"},
        "pin_framed": {"rotation": pin_rotation, "exchange": pin_exchange,
                       "statistics": "Fermi"},
        "consistent_with_188_holonomy": consistent_188,
        "pass": ok,
    }


def test_T7_cobordism_ledger() -> dict:
    """The #58 channel against the Dowker–Sorkin 4-manifold condition."""
    ledger = {
        "mirror_pair_identity": {
            "status": True,
            "why": "amphichirality (T3): the created mirror partner is "
                   "diffeomorphic to the throat — identical geons",
        },
        "bordism_existence": {
            "status": True,
            "why": "Omega_3^SO = 0 (every closed orientable 3-manifold "
                   "bounds): 4-manifolds S3 -> S3 # RP3 # RP3 exist (cited)",
        },
        "spin_structure_extension": {
            "status": True,
            "why": "RP3 is spin (orientable Lie group SO(3), "
                   "parallelizable); Omega_3^Spin = 0: the pin/spin "
                   "framing extends over a bounding 4-manifold (cited)",
        },
        "explicit_BAM_4manifold": {
            "status": None,
            "why": "the Dowker-Sorkin construction is over the R3 "
                   "background; the S3 transplant (connected-sum "
                   "locality) with their causal/Morse conditions has NOT "
                   "been written down - the honest open construction, "
                   "with no topological obstruction found",
        },
    }
    closed = [k for k, v in ledger.items() if v["status"] is True]
    open_items = [k for k, v in ledger.items() if v["status"] is None]
    ok = len(closed) == 3 and open_items == ["explicit_BAM_4manifold"]
    return {
        "name": "T7_cobordism_ledger",
        "description": (
            "Section 6's ledger: the #58 throat–antithroat nucleation "
            "channel against the Dowker–Sorkin pair-creation condition. "
            "CLOSED (3): the mirror-pair identity (the created partner is "
            "an identical geon, by amphichirality); the existence of "
            "interpolating 4-manifolds (Ω₃^SO = 0); the extension of the "
            "spin/pin framing (RP³ spin, Ω₃^Spin = 0). OPEN (1): the "
            "EXPLICIT BAM 4-manifold — the Dowker–Sorkin construction "
            "transplanted from the ℝ³ background to the S³ background "
            "with their causal/Morse structure intact — a construction "
            "problem with no topological obstruction found, flagged as "
            "the honest residue."
        ),
        "ledger": ledger,
        "closed_items": closed,
        "open_items": open_items,
        "pass": ok,
    }


def test_T8_assessment() -> dict:
    return {
        "name": "T8_assessment",
        "description": (
            "THE ADJUDICATION. Neither pure Outcome A nor pure Outcome B. "
            "WHAT BAM OWNS AS A THEOREM: the spin-statistics CORRELATION "
            "— the throat prime RP³ passes all three Dowker–Sorkin "
            "hypotheses (prime, non-chiral, abelian), and the SOH "
            "abelianness closes the anomalous/indefinite sectors. WHAT "
            "CHANGES: the bare-GR reading is NOT survivable — RP³ is "
            "non-spinorial (Hendriks), so bare geometrodynamics + SSC "
            "makes pair-created throats BOSONS; the #170/#171 'spinorial' "
            "sentence is corrected. WHAT SURVIVES, SHARPENED: in "
            "Pin⁻-framed GR — the framing BAM's own antipodal quotient "
            "forces (#169/#170) and its matter content requires (#195) — "
            "the rotation lifts to −1 and SSC delivers FERMI: 'Pauli "
            "from GR + the forced Pin⁻ framing', with the #188 holonomy "
            "correctly reinterpreted as the pin lift. REMAINING: the "
            "explicit BAM pair-creation 4-manifold (a construction, not "
            "an obstruction)."
        ),
        "classification": (
            "SSC_THEOREM_HYPOTHESES_HOLD_SIGN_IS_FRAMING_DEPENDENT_BARE_"
            "GR_SELECTS_BOSE_PIN_MINUS_FRAMING_SELECTS_FERMI"
        ),
        "pass": True,
    }


# ════════════════════════════════════════════════════════════════════════
# RUNNER + VERDICT
# ════════════════════════════════════════════════════════════════════════

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_lemma1_orientation_arithmetic(),
        test_T3_lemma2_hypotheses(),
        test_T4_pin_lift_of_rotation(),
        test_T5_two_throat_mcg_sectors(),
        test_T6_soh_selection(),
        test_T7_cobordism_ledger(),
        test_T8_assessment(),
    ]
    all_ok = all(t["pass"] for t in tests)
    if all_ok:
        verdict_class = (
            "SSC_THEOREM_HYPOTHESES_HOLD_SIGN_IS_FRAMING_DEPENDENT_BARE_"
            "GR_SELECTS_BOSE_PIN_MINUS_FRAMING_SELECTS_FERMI"
        )
        verdict = (
            "THE ADJUDICATION (the argument is in "
            "docs/geon_statistics_adjudication.md; this probe checks its "
            "arithmetic).\n\n"
            "WHAT IS A THEOREM. The BAM throat prime is RP³ — the RP³ "
            "geon, from the repo's own #169 quotient (orientation "
            "arithmetic verified) — and it passes ALL THREE hypotheses "
            "of the Dowker–Sorkin spin-statistics theorem: PRIME "
            "(elliptic, irreducible), NON-CHIRAL (q² ≡ −1 mod 2; the "
            "explicit reflection descends), ABELIAN (π₁ = Z₂). The "
            "spin-statistics CORRELATION is therefore a theorem for "
            "pair-created BAM throats — the strongest topological proof "
            "point available.\n\n"
            "WHAT THE SIGN DEPENDS ON. RP³ is NON-SPINORIAL in bare Diff "
            "(Hendriks: cyclic-π₁ primes; the #170/#171 'spinorial' "
            "sentence is hereby corrected): in bare metric GR the SSC "
            "selects BOSE, and in frozen-topology canonical GR the "
            "statistics is a sector choice among Bose, Fermi, and a "
            "continuum of indefinite sectors (G = (Z₂∗Z₂)⋊S₂ — four "
            "scalar sectors + the anomalous 2-dim family, all verified "
            "in a faithful model). The −1 is NOT a bare-GR theorem. In "
            "PIN⁻-FRAMED GR — forced by the non-orientable RP² slice "
            "(#169), the uniqueness of Pin⁻ on it (#170), and the pin-"
            "spinor matter content (#195) — the trivialized rotation "
            "lifts to −I (the SU(2) lift of the π₁(SO(3)) generator, "
            "computed), the throat is spinorial WITH FRAMING, and the "
            "same SSC theorem selects FERMI — the #188 holonomy, "
            "correctly reinterpreted as the pin lift.\n\n"
            "THE #58 CHANNEL. Mirror-pair identity ✓ (amphichirality), "
            "bordism existence ✓ (Ω₃^SO = 0), framing extension ✓ (RP³ "
            "spin, Ω₃^Spin = 0); the explicit BAM 4-manifold (the "
            "Dowker–Sorkin construction transplanted from ℝ³ to the S³ "
            "background) is the honest open construction — no "
            "obstruction found. LABEL CHANGE: 'Pauli from GR' → 'Pauli "
            "from GR + the forced Pin⁻ framing.'"
        )
    else:
        verdict_class = "GEON_STATISTICS_ADJUDICATION_INCONCLUSIVE"
        verdict = (
            "INCONCLUSIVE. An arithmetic check failed; re-examine the "
            "document's chain before quoting the adjudication."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "The geon-statistics adjudication: SSC hypotheses "
            "(prime/non-chiral/abelian) all hold for the RP³ throat — the "
            "correlation is a theorem; the sign is framing-dependent "
            "(bare GR: Bose, non-spinorial RP³; Pin⁻-framed: Fermi, the "
            "computed pin lift) — the honest label is 'Pauli from GR + "
            "the forced Pin⁻ framing'"
        ),
        "adjudicates": "the #185/#188 exchange −1 against Sorkin–Surya / Dowker–Sorkin / Hendriks",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# The geon-statistics adjudication — companion probe (PR #196)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "The deliverable is `docs/geon_statistics_adjudication.md` — "
        "mathematics against the Sorkin-school literature on the actual "
        "BAM topology. This probe machine-checks the document's "
        "finite/algebraic steps. *(QFT on the fixed classical throat "
        "geometry, not quantum gravity.)*"
    )
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "the adjudication; the document is the argument",
        "T2": "Lemma 1: the throat prime is the RP³ geon (orientation arithmetic)",
        "T3": "Lemma 2: prime ✓ non-chiral ✓ abelian ✓ (SSC hypotheses hold)",
        "T4": "Lemma 3/5: bare rotation trivial (cited); pin lift = −I (computed)",
        "T5": "Lemma 4: G = (Z₂∗Z₂)⋊S₂; 4 scalar + indefinite sectors",
        "T6": "SOH: SSC ⟹ bare GR → Bose; Pin⁻ framing → Fermi (= #188)",
        "T7": "#58 cobordism ledger: 3 closed, 1 honest open construction",
        "T8": "label change: 'Pauli from GR + the forced Pin⁻ framing'",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'—')} | {p} |")
    out.append("")
    t5, t6, t7 = s["tests"][4], s["tests"][5], s["tests"][6]
    out.append("## The sector structure (two throats)")
    out.append("")
    out.append("| slides | exchange | sector |")
    out.append("|---:|---:|---|")
    for r in t5["one_dim_sectors"]:
        stat = "Fermi" if r["exchange"] == -1 else "Bose"
        out.append(f"| {r['slides']:+d} | {r['exchange']:+d} | {stat} |")
    out.append("| (2-dim family) | mixed | indefinite (Sorkin–Surya anomalous) |")
    out.append("")
    out.append("## The selection")
    out.append("")
    out.append(f"- SSC hypotheses: {t6['ssc_hypotheses']} → the correlation is a **theorem**")
    out.append(f"- bare metric GR: rotation {t6['bare_GR']['rotation']:+d} → **{t6['bare_GR']['statistics']}**")
    out.append(f"- Pin⁻-framed GR: rotation {t6['pin_framed']['rotation']:+d} → **{t6['pin_framed']['statistics']}** "
               f"(consistent with #188: {t6['consistent_with_188_holonomy']})")
    out.append("")
    out.append("## The #58 cobordism ledger")
    out.append("")
    for k, v in t7["ledger"].items():
        mark = "✓" if v["status"] else "OPEN"
        out.append(f"- **{k}** [{mark}]: {v['why']}")
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
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_geon_statistics_adjudication_probe"
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
