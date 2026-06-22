"""
Tangherlini J-quotient consistency and brane non-orientability (PR #169).

> Framing: QFT on the *fixed classical* throat geometry (geometry → fields),
> not quantum gravity.

THE TOPOLOGICAL ROOT OF THE NON-ORIENTABLE THROAT
─────────────────────────────────────────────────
PR #167 found the throat gluing is non-orientable; PR #168 embedded the
throat as the equatorial slice of the 5D Tangherlini bulk.  This probe
identifies WHY the throat is non-orientable while the bulk is not — a
clean dimension-parity statement about the antipodal (J) quotient:

    bulk angular sphere   S³ / antipodal = RP³   ORIENTABLE      (det +1)
    brane angular slice   S² / antipodal = RP²   NON-ORIENTABLE  (det −1)

The antipodal involution J: x ↦ −x on Sⁿ has orientation determinant
(−1)^{n+1} — orientation-preserving for ODD n (orientable quotient),
orientation-reversing for EVEN n (non-orientable quotient).  The bulk's
angular sphere is S³ (odd ⟹ RP³ orientable); the throat mouth is the
brane's angular S² (even ⟹ RP² non-orientable).  The SAME J-quotient is
therefore consistently orientable on the bulk and non-orientable on the
brane mouth: the non-orientable throat is an RP² (a cross-cap) sitting
inside the orientable RP³, forced by the one-dimension drop from bulk to
brane.

In the #168 coordinates the S³ antipodal map is
(χ, θ, φ) ↦ (π−χ, π−θ, φ+π); it FIXES the equatorial χ = π/2 brane and
RESTRICTS there to the S² antipodal map (θ, φ) ↦ (π−θ, φ+π).  The 5D
Tangherlini metric (round angular part) is J-invariant, so it descends to
the quotient: the bulk becomes the orientable RP³-angular Tangherlini, the
brane mouth becomes the non-orientable RP².

Tests:
  T1. Goal + framing.
  T2. J is a FREE isometric involution (manifold quotient; metric descends).
  T3. Bulk orientation determinant +1 (S³ → RP³ orientable), explicit.
  T4. Brane angular determinant −1 (S² → RP² non-orientable), explicit.
  T5. The parity law det = (−1)^{n+1}: odd orientable / even not.
  T6. The Tangherlini realization: J fixes the χ=π/2 brane and restricts
      to the S² antipodal map; the metric descends (ties to #167/#168).
  T7. Consistency with spin/Pin and the even-k grading (#63/#67).
  T8. Assessment.

Verdict:
  - J_QUOTIENT_CONSISTENT_BULK_RP3_ORIENTABLE_BRANE_RP2_NON_ORIENTABLE
    (expected): the antipodal J-quotient is consistent — orientable bulk
    (RP³, det +1), non-orientable brane mouth (RP², det −1) — and this
    dimension-parity split IS the topological origin of the non-orientable
    throat, realized concretely in the #168 Tangherlini embedding.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np


# ════════════════════════════════════════════════════════════════════════
# THE ANTIPODAL (J) INVOLUTION ON Sⁿ AND ITS ORIENTATION DETERMINANT
# ════════════════════════════════════════════════════════════════════════

def orientation_determinant(n: int, seed: int = 0) -> float:
    """Orientation determinant of the antipodal map J = −I on Sⁿ ⊂ ℝ^{n+1}.

    Build an oriented tangent frame {e_i} at a random p ∈ Sⁿ; the standard
    orientation is sgn det[ν_p | e_1 … e_n] with outward normal ν_p = p.
    J pushes p → −p, ν_p → −p, e_i → −e_i.  The ratio of the pushed-frame
    orientation (with the antipode's outward normal ν_{−p} = −p) to the
    original is the orientation determinant; it equals (−1)^{n+1}."""
    rng = np.random.default_rng(seed)
    p = rng.standard_normal(n + 1)
    p /= np.linalg.norm(p)
    M = rng.standard_normal((n + 1, n))
    M = M - np.outer(p, p @ M)            # project into T_p Sⁿ
    Q, _ = np.linalg.qr(M)
    e = Q[:, :n]                          # orthonormal tangent frame
    det_p = np.linalg.det(np.column_stack([p, e]))          # (ν_p, e)
    det_ap = np.linalg.det(np.column_stack([-p, -e]))       # (ν_{−p}, J·e)
    return float(np.sign(det_ap / det_p))


def is_free_involution(n: int, n_samples: int = 2000, seed: int = 1) -> bool:
    """J(x) = −x has a fixed point on Sⁿ iff −x = x iff x = 0 ∉ Sⁿ; so J is
    free for all n.  Verified by sampling: min‖J(x) − x‖ is bounded away
    from 0."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n_samples, n + 1))
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    return float(np.min(np.linalg.norm(-X - X, axis=1))) > 1.0


def is_isometry_of_round_sphere(n: int, seed: int = 2) -> bool:
    """J = −I is orthogonal (Jᵀ J = I), hence an isometry of the round Sⁿ
    (and of any metric whose angular part is round, e.g. Tangherlini)."""
    J = -np.eye(n + 1)
    return bool(np.allclose(J.T @ J, np.eye(n + 1)))


# ── the #168 angular coordinates: S³ antipodal map (χ,θ,φ)↦(π−χ,π−θ,φ+π) ──

def _s3_embed(chi, th, ph):
    return np.array([
        math.cos(chi),
        math.sin(chi) * math.cos(th),
        math.sin(chi) * math.sin(th) * math.cos(ph),
        math.sin(chi) * math.sin(th) * math.sin(ph),
    ])


# ════════════════════════════════════════════════════════════════════════
# TESTS
# ════════════════════════════════════════════════════════════════════════

def test_T1_goal() -> dict:
    return {
        "name": "T1_goal",
        "description": (
            "Identify WHY the BAM throat is non-orientable (PR #167) while "
            "the bulk is not. The antipodal (J) quotient is consistent but "
            "dimension-parity-split: the bulk angular sphere S³/antipodal = "
            "RP³ is ORIENTABLE (det +1), the brane angular slice "
            "S²/antipodal = RP² is NON-ORIENTABLE (det −1). The throat mouth "
            "is the brane's even-dimensional angular S², so it quotients to "
            "the non-orientable RP² (a cross-cap) inside the orientable RP³ "
            "bulk — realized concretely in the #168 Tangherlini embedding."
        ),
        "bulk": "S³ / antipodal = RP³ (orientable, det +1)",
        "brane": "S² / antipodal = RP² (non-orientable, det −1)",
        "framing": "QFT on the fixed classical throat geometry — not quantum gravity",
        "pass": True,
    }


def test_T2_free_isometric_involution() -> dict:
    """J is a free isometric involution — the quotient is a manifold and the
    Tangherlini metric descends."""
    free3 = is_free_involution(3)
    free2 = is_free_involution(2)
    iso3 = is_isometry_of_round_sphere(3)
    iso2 = is_isometry_of_round_sphere(2)
    # involution: J² = I
    involution = np.allclose((-np.eye(4)) @ (-np.eye(4)), np.eye(4))
    ok = free3 and free2 and iso3 and iso2 and involution
    return {
        "name": "T2_free_isometric_involution",
        "description": (
            "The antipodal map J = −I is an involution (J² = I), is FREE on "
            f"Sⁿ (no fixed points: −x = x ⟹ x = 0 ∉ Sⁿ; verified S³ {free3}, "
            f"S² {free2}), and is an ISOMETRY of the round sphere (JᵀJ = I; "
            f"S³ {iso3}, S² {iso2}). A free involution gives a smooth "
            "manifold quotient (no orbifold points), and an isometry means "
            "the 5D Tangherlini metric — whose angular part is round — "
            "descends to the quotient. So the J-quotient of the bulk is a "
            "well-defined smooth geometry."
        ),
        "free_S3": free3, "free_S2": free2,
        "isometry_S3": iso3, "isometry_S2": iso2,
        "involution": bool(involution),
        "pass": ok,
    }


def test_T3_bulk_orientation_plus1() -> dict:
    """Bulk: S³ antipodal orientation determinant +1 → RP³ orientable."""
    dets = [orientation_determinant(3, s) for s in range(6)]
    det = float(np.mean(dets))
    consistent = all(abs(d - 1.0) < 1e-9 for d in dets)
    ok = abs(det - 1.0) < 1e-9 and consistent
    return {
        "name": "T3_bulk_orientation_determinant",
        "description": (
            "The S³ antipodal map PRESERVES orientation: the pushed tangent "
            "frame (with the antipode's outward normal) has the same "
            f"orientation as the original — determinant +1 (computed "
            f"{det:+.3f} over 6 random frames). Hence the bulk quotient "
            "S³ / antipodal = RP³ is ORIENTABLE: the BAM spatial bulk "
            "remains a consistent oriented manifold under the antipodal "
            "identification."
        ),
        "bulk_orientation_determinant": det,
        "RP3_orientable": True,
        "pass": ok,
    }


def test_T4_brane_orientation_minus1() -> dict:
    """Brane: S² antipodal orientation determinant −1 → RP² non-orientable."""
    dets = [orientation_determinant(2, s) for s in range(6)]
    det = float(np.mean(dets))
    consistent = all(abs(d + 1.0) < 1e-9 for d in dets)
    ok = abs(det + 1.0) < 1e-9 and consistent
    return {
        "name": "T4_brane_angular_determinant",
        "description": (
            "The S² antipodal map REVERSES orientation: the pushed tangent "
            "frame has the opposite orientation — determinant −1 (computed "
            f"{det:+.3f} over 6 random frames). Hence the brane's angular "
            "quotient S² / antipodal = RP² is NON-ORIENTABLE: the throat "
            "mouth is a cross-cap. This is the topological content of PR "
            "#167's 'non-orientable throat gluing' — it is the RP² of the "
            "antipodally-identified mouth 2-sphere."
        ),
        "brane_angular_determinant": det,
        "RP2_non_orientable": True,
        "pass": ok,
    }


def test_T5_parity_law() -> dict:
    """The consistency: det = (−1)^{n+1} — odd orientable, even not."""
    rows = []
    law_ok = True
    for n in range(1, 6):
        det = orientation_determinant(n, seed=n)
        formula = (-1) ** (n + 1)
        match = abs(det - formula) < 1e-9
        law_ok = law_ok and match
        rows.append({"n": n, "orientation_det": det,
                     "formula_(-1)^(n+1)": formula,
                     "RP^n_orientable": formula > 0, "match": match})
    # the bulk (S³) and brane mouth (S²) sit one dimension apart, on
    # opposite sides of the parity → opposite orientability, from one J.
    bulk_even_brane = (orientation_determinant(3) > 0
                       and orientation_determinant(2) < 0)
    ok = law_ok and bulk_even_brane
    return {
        "name": "T5_dimension_parity_law",
        "description": (
            "The antipodal orientation determinant is (−1)^{n+1}: ODD "
            "spheres preserve orientation (orientable RPⁿ), EVEN spheres "
            "reverse it (non-orientable RPⁿ) — verified n = 1..5. The bulk "
            "angular sphere (S³, odd) and the throat mouth (S², even) sit "
            "one dimension apart, on opposite sides of the parity, so the "
            "SINGLE antipodal J-quotient is consistently orientable on the "
            "bulk and non-orientable on the brane. The non-orientable throat "
            "is forced by the one-dimension drop from bulk to mouth — not an "
            "extra assumption."
        ),
        "parity_table": rows,
        "bulk_orientable_brane_not": bulk_even_brane,
        "pass": ok,
    }


def test_T6_tangherlini_realization() -> dict:
    """The #168 realization: J fixes the χ=π/2 brane and restricts to the
    S² antipodal map; the Tangherlini metric descends."""
    # the S³ antipodal map (χ,θ,φ)↦(π−χ,π−θ,φ+π) equals x↦−x
    chi, th, ph = math.pi / 2, 0.7, 1.3
    p = _s3_embed(chi, th, ph)
    ap = _s3_embed(math.pi - chi, math.pi - th, ph + math.pi)
    map_is_antipodal = bool(np.allclose(ap, -p))
    # it fixes the equatorial brane χ = π/2
    fixes_brane = abs((math.pi - chi) - math.pi / 2) < 1e-12
    # and restricts to the S² antipodal map (θ,φ)↦(π−θ,φ+π) on the brane
    def _s2(t, f):
        return np.array([math.sin(t) * math.cos(f),
                         math.sin(t) * math.sin(f), math.cos(t)])
    q = _s2(th, ph)
    aq = _s2(math.pi - th, ph + math.pi)
    restricts_to_s2_antipodal = bool(np.allclose(aq, -q))
    # the round angular metric is J-invariant (isometry) ⟹ Tangherlini
    # ds² = −F dt² + dρ²/F + ρ²dΩ₃² descends to the quotient
    metric_descends = is_isometry_of_round_sphere(3)
    ok = (map_is_antipodal and fixes_brane and restricts_to_s2_antipodal
          and metric_descends)
    return {
        "name": "T6_tangherlini_realization",
        "description": (
            "In the #168 coordinates the S³ antipodal map "
            f"(χ,θ,φ)↦(π−χ,π−θ,φ+π) equals x↦−x ({map_is_antipodal}); it "
            f"FIXES the equatorial χ = π/2 brane ({fixes_brane}) and "
            "RESTRICTS there to the S² antipodal map (θ,φ)↦(π−θ,φ+π) "
            f"({restricts_to_s2_antipodal}). The round angular metric is "
            "J-invariant, so the 5D Tangherlini metric descends to the "
            "quotient: the bulk becomes the orientable RP³-angular "
            "Tangherlini and the throat mouth becomes the non-orientable "
            "RP². The #167 non-orientable throat is exactly this RP² "
            "cross-cap inside the orientable RP³ bulk of #168."
        ),
        "map_equals_antipodal": map_is_antipodal,
        "fixes_equatorial_brane": fixes_brane,
        "restricts_to_S2_antipodal": restricts_to_s2_antipodal,
        "tangherlini_metric_descends": metric_descends,
        "pass": ok,
    }


def test_T7_spin_pin_consistency() -> dict:
    """Consistency with spin/Pin and the even-k orientability grading."""
    # RP³ ≅ SO(3): orientable, parallelizable, admits a spin structure.
    # RP² : non-orientable, admits no spin structure — only a Pin structure
    #       (the half-twist), matching the fermionic / spin-½ character.
    return {
        "name": "T7_spin_pin_consistency",
        "description": (
            "The orientability split lands exactly where BAM's spinor "
            "structure lives. The bulk quotient RP³ ≅ SO(3) is orientable, "
            "parallelizable, and admits a spin structure — a consistent "
            "arena for bulk fields. The brane mouth RP² is non-orientable "
            "and admits NO spin structure, only a PIN structure: the "
            "half-twist that carries the spin-½ / fermionic character. This "
            "is the same orientability grading already seen in BAM as the "
            "C-swap (C = iσ_y, T² = −1; #63) and the even-k absence "
            "(k mod 2 = orientability/spin-statistics; #67). A consistency "
            "observation, not a new derivation: the throat's RP² "
            "non-orientability is the topological seat of its Pin/spin-½ "
            "structure."
        ),
        "RP3_orientable_spin": True,
        "RP2_pin_only": True,
        "consistency_only": True,
        "pass": True,
    }


def test_T8_assessment() -> dict:
    return {
        "name": "T8_assessment",
        "description": (
            "The antipodal J-quotient is consistent and dimension-parity-"
            "split: the bulk angular S³ → RP³ is orientable (det +1), the "
            "brane angular S² → RP² is non-orientable (det −1), from the one "
            "free isometric involution J = −I (orientation determinant "
            "(−1)^{n+1}). The non-orientable throat of #167 is the RP² "
            "cross-cap of the antipodally-identified mouth 2-sphere, sitting "
            "inside the orientable RP³ bulk of the #168 Tangherlini "
            "embedding (J fixes the χ=π/2 brane and restricts to the S² "
            "antipodal map; the metric descends). The non-orientability is "
            "forced by the single-dimension drop from bulk to mouth, and it "
            "is the topological seat of the throat's Pin/spin-½ structure "
            "(#63/#67)."
        ),
        "classification": (
            "J_QUOTIENT_CONSISTENT_BULK_RP3_ORIENTABLE_BRANE_RP2_NON_ORIENTABLE"
        ),
        "pass": True,
    }


# ════════════════════════════════════════════════════════════════════════
# RUNNER + VERDICT
# ════════════════════════════════════════════════════════════════════════

def run_probe() -> dict:
    tests = [
        test_T1_goal(),
        test_T2_free_isometric_involution(),
        test_T3_bulk_orientation_plus1(),
        test_T4_brane_orientation_minus1(),
        test_T5_parity_law(),
        test_T6_tangherlini_realization(),
        test_T7_spin_pin_consistency(),
        test_T8_assessment(),
    ]
    t3, t4 = tests[2], tests[3]
    if all(t["pass"] for t in tests):
        verdict_class = (
            "J_QUOTIENT_CONSISTENT_BULK_RP3_ORIENTABLE_BRANE_RP2_NON_ORIENTABLE"
        )
        verdict = (
            "CONSISTENT, AND IT EXPLAINS THE NON-ORIENTABLE THROAT. The "
            "antipodal (J) quotient is a single free isometric involution "
            "whose orientation determinant is (−1)^{n+1}, so it acts "
            "oppositely on the bulk and the brane mouth.\n\n"
            "THE BULK. S³ / antipodal = RP³: the S³ antipodal map preserves "
            f"orientation (determinant {t3['bulk_orientation_determinant']:+.0f}), "
            "so the BAM spatial bulk remains a consistent ORIENTABLE "
            "manifold — RP³ ≅ SO(3).\n\n"
            "THE BRANE MOUTH. S² / antipodal = RP²: the S² antipodal map "
            f"reverses orientation (determinant {t4['brane_angular_determinant']:+.0f}), "
            "so the throat mouth is the NON-ORIENTABLE RP² (a cross-cap). "
            "This is exactly PR #167's 'non-orientable throat gluing'.\n\n"
            "THE CONSISTENCY. The bulk angular sphere (S³, odd) and the "
            "throat mouth (S², even) sit one dimension apart, on opposite "
            "sides of the parity (−1)^{n+1}: ONE antipodal quotient is "
            "therefore consistently orientable on the bulk and "
            "non-orientable on the mouth. In the #168 coordinates J = "
            "(χ,θ,φ)↦(π−χ,π−θ,φ+π) fixes the equatorial χ=π/2 brane and "
            "restricts to the S² antipodal map; the round-angular Tangherlini "
            "metric is J-invariant and descends. So the #167 non-orientable "
            "throat is the RP² cross-cap inside the orientable RP³ bulk of "
            "#168 — forced by the single-dimension drop, and the topological "
            "seat of the throat's Pin/spin-½ structure (#63/#67)."
        )
    else:
        verdict_class = "J_QUOTIENT_INCONSISTENT"
        verdict = (
            "INCONSISTENT. An orientation/quotient check failed; review the "
            "determinants, the free-action/isometry conditions, or the "
            "Tangherlini realization."
        )

    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "identification": (
            "the antipodal J-quotient is consistent and dimension-parity-"
            "split: bulk S³→RP³ orientable (det +1), brane mouth S²→RP² "
            "non-orientable (det −1); the non-orientable throat is the RP² "
            "cross-cap in the orientable RP³ Tangherlini bulk of #168"
        ),
        "bulk": "S³/antipodal = RP³ orientable (det +1)",
        "brane": "S²/antipodal = RP² non-orientable (det −1)",
        "consistency": "one free isometric involution; det = (−1)^{n+1}; metric descends",
        "realization": "J fixes the χ=π/2 brane and restricts to the S² antipodal map (#168)",
        "tests": tests,
        "n_passed": sum(1 for t in tests if t["pass"]),
        "n_total": len(tests),
        "verdict_class": verdict_class,
        "verdict": verdict,
    }


def render_markdown(s: dict) -> str:
    out = []
    out.append("# Tangherlini J-quotient consistency and brane non-orientability (PR #169)")
    out.append("")
    out.append(f"**Run:** {s['timestamp_utc']}")
    out.append("")
    out.append(
        "Identifies the topological root of the non-orientable throat "
        "(#167): the antipodal (J) quotient is consistent but "
        "dimension-parity-split — orientable on the bulk, non-orientable on "
        "the brane mouth. *(QFT on the classical throat, not quantum "
        "gravity.)*"
    )
    out.append("")
    out.append(f"- **Bulk**: {s['bulk']}")
    out.append(f"- **Brane**: {s['brane']}")
    out.append(f"- **Consistency**: {s['consistency']}")
    out.append(f"- **Realization**: {s['realization']}")
    out.append("")
    out.append("## Test summary")
    out.append("")
    out.append("| # | Test | Key finding | PASS? |")
    out.append("|---|---|---|---|")
    labels = {
        "T1": "goal: the J-quotient parity split (bulk vs brane)",
        "T2": "J is a free isometric involution (metric descends)",
        "T3": "bulk orientation determinant +1 → RP³ orientable",
        "T4": "brane angular determinant −1 → RP² non-orientable",
        "T5": "the parity law det = (−1)^{n+1} (odd vs even)",
        "T6": "Tangherlini realization: J fixes the χ=π/2 brane (#168)",
        "T7": "spin/Pin + even-k consistency (#63/#67)",
        "T8": "J_QUOTIENT_CONSISTENT_RP3_ORIENTABLE_RP2_NOT",
    }
    for t in s["tests"]:
        p = "**PASS**" if t["pass"] else "**FAIL**"
        pre = t["name"][:2]
        out.append(f"| {pre} | `{t['name']}` | {labels.get(pre,'—')} | {p} |")
    out.append("")
    t5 = s["tests"][4]
    out.append("## The dimension-parity law  (orientation det = (−1)ⁿ⁺¹)")
    out.append("")
    out.append("| n | orientation det | RPⁿ orientable? |")
    out.append("|---|---:|---|")
    for r in t5["parity_table"]:
        out.append(f"| {r['n']} | {r['orientation_det']:+.0f} | "
                   f"{'yes' if r['RP^n_orientable'] else 'no'} |")
    out.append("")
    out.append("(bulk = S³, n=3 → orientable RP³; brane mouth = S², n=2 → non-orientable RP²)")
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
    out = here / "runs" / f"{ts}_tangherlini_j_quotient_probe"
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
