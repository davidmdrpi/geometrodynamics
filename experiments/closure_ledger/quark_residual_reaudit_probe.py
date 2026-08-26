"""The downstream re-derivation PR #271 deferred: the quark residual sector.

THE RESULT, FIRST
─────────────────
PR #271 corrected the 5D radial scalar operator and reopened the LEPTON
sector's geometric closure. The three residual knobs of the frozen v3 QUARK
lock are derived from the same eigensolver, so the expectation was that they
would break the same way.

They do not. ALL THREE MOVE TOWARD THEIR LOCKED VALUES.

    residual                          locked    legacy            corrected
    pinhole    Sum V_max[1..5]        22.25     22.008 (-1.09%)   22.331 (+0.36%)
    transport  mean <u1|dV|u2>         0.54      0.5447 (+0.88%)   0.5438 (+0.70%)
    resistance transport*ln(a5/a1)     0.14      0.1407 (+0.49%)   0.1400 (-0.02%)

WHY THE TWO SECTORS PART COMPANY
────────────────────────────────
One barrier feeds both sectors, and the correction moves it the same way for
both. What differs is how hard each Hamiltonian leans on it:

    lepton   d ln m_mu / d ln gamma    = -17.471
    quark    d ln m_s  / d ln pinhole  =  +4.798

a factor of 3.64. The lepton's corrected residual (-0.75%) is a MEASURED 15.2%
muon error; the quark's (+0.36%) is a measured 1.79% strange error.

    A PERCENTAGE AGREEMENT BETWEEN A GEOMETRIC QUANTITY AND A FITTED KNOB
    MEANS NOTHING UNTIL MULTIPLIED BY THE ELASTICITY OF WHAT IT FEEDS.

That is the transferable lesson, and it is why #271's headline and this one
point in opposite directions without contradicting each other.

A NUMBER THIS ROUND CORRECTS
────────────────────────────
The README row written in #271 reads "the residual improves to -0.75%, but
d ln m_mu/d ln gamma = -17.5 at the lock makes that a 9% muon error". The
elasticity and the residual are both right; the 9% is not. It came from a
generic illustration in the #271 probe ("a half-percent error in gamma is a
nine-percent error in the muon"), applied to a residual that is not a half
percent. Linearising -0.75% gives +14.0%; the locked block actually returns
+15.2%, and #271's own A/B/C table carried that value the whole time (T7).

THE COMPOSITE IS NOT THE SUM OF ITS PARTS
─────────────────────────────────────────
Substituting all three derived residuals at once, retuning NOTHING, the
corrected set scores WORSE than the legacy set (3.78% vs 3.44%) even though
each corrected residual is individually closer to its lock. The legacy triple
enjoys a partial cancellation its pinhole error alone would not earn (3.61%
alone -> 3.44% together); the corrected triple has none (3.40% -> 3.78%).

Neither composite reaches the fitted lock's 1.61%. "Each knob is derived to
within 1%" and "the derived knobs reproduce the ladder" are DIFFERENT CLAIMS
and only the first was ever established -- under either operator. The
correction did not create that gap, it exposed it.

THE ONE THING THAT GENUINELY IMPROVES
─────────────────────────────────────
Read as demands on R_OUTER, the legacy operator has both sectors asking for
MORE than 1.26 (1.2723 and 1.2874), putting the canonical value outside the
bracket they define. The corrected operator has them STRADDLE it -- quark
1.2564 below, lepton 1.2679 above, 1.26 at 31% of the way across.

That is evidence for R_OUTER = 1.26, and it is not the evidence #271 reopened:
that was a single-sector fixed point, this is a two-sector bracket. It is also
weak -- two numbers 0.91% apart bracket anything between them. Recorded as
suggestive, NOT as a restored derivation.

NOT RE-RUN, AND WHY
───────────────────
The same exclusion as #271: the four quark shell-index axioms are expressible
in k_5 = 5 alone, the Hopf transport phase is topological, and the flavor-CP
layer inherits the v3 eigenvalues by construction. None reads the radial
scalar operator. Proximity is not dependence.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from geometrodynamics.qcd.residual_audit import (  # noqa: E402
    LOCKED_PINHOLE, LOCKED_RESISTANCE, LOCKED_TRANSPORT,
    measure_the_cross_sector_r_outer_bracket,
    measure_the_ladder_without_retuning,
    measure_the_lepton_resistance_selector, measure_the_n_stability,
    measure_the_pinhole_elasticity, measure_the_quark_ledger,
    measure_the_three_residuals_under_both_operators,
    measure_the_two_sectors_side_by_side)


def run_probe() -> dict:
    checks: List[dict] = []

    res = measure_the_three_residuals_under_both_operators()
    checks.append({
        "id": "Q1",
        "name": "*** all three quark residuals move TOWARD their locked values ***",
        "detail": res,
        "pass": bool(res["all_three_move_toward_the_lock"])})

    lad = measure_the_ladder_without_retuning()
    checks.append({
        "id": "Q2",
        "name": "the composite reverses the per-knob ordering, and neither reaches the lock",
        "detail": lad,
        "pass": bool(lad["corrected_pinhole_alone_is_better"]
                     and lad["corrected_composite_is_worse"]
                     and lad["legacy_composite_beats_its_own_pinhole"]
                     and lad["neither_composite_reaches_the_lock"])})

    ela = measure_the_pinhole_elasticity()
    checks.append({
        "id": "Q3",
        "name": "the quark ladder is 3.6x less stiff than the lepton chain",
        "detail": ela,
        "pass": bool(ela["stiffest_species"] == "s"
                     and ela["correction_halves_the_miss"]
                     and ela["ladder_error_at_what_it_wants"]
                     < ela["ladder_error_at_the_v3_lock"])})

    brk = measure_the_cross_sector_r_outer_bracket()
    checks.append({
        "id": "Q4",
        "name": "only the corrected operator brackets the canonical R_OUTER",
        "detail": brk,
        "pass": bool(brk["only_the_corrected_operator_brackets_1_26"]
                     and brk["correction_narrows_the_split"])})

    sel = measure_the_lepton_resistance_selector()
    checks.append({
        "id": "Q5",
        "name": "7pi/100 survives, but not for the reason the README gives",
        "detail": sel,
        "pass": bool(sel["competitor_beat_the_closed_form_under_legacy"]
                     and sel["closed_form_wins_under_the_correction"])})

    nst = measure_the_n_stability()
    checks.append({
        "id": "Q6",
        "name": "N still drifts: beta remains the model's last fit knob",
        "detail": nst,
        "pass": bool(nst["N_drifts_under_every_mode"])})

    two = measure_the_two_sectors_side_by_side()
    checks.append({
        "id": "Q7",
        "name": "the two sectors side by side, with MEASURED errors",
        "detail": two,
        "pass": bool(two["elasticity_ratio"] > 3.0
                     and two["lepton"]["measured_error_pct"] > 10.0
                     and two["quark"]["measured_error_pct"] < 3.0)})

    led = measure_the_quark_ledger()
    checks.append({
        "id": "Q8", "name": "the ledger", "detail": led,
        "pass": bool(len(led["entries"]) >= 10)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    L: List[str] = [
        "# The quark residual sector, re-derived from the corrected operator",
        "",
        f"**{summary['passed']}/{summary['total']} checks pass.**",
        "",
        "PR #271 corrected the 5D radial scalar operator and reopened the "
        "**lepton** sector's geometric closure. The frozen v3 **quark** lock "
        "derives three residual knobs from the same eigensolver, so the "
        "expectation was that they would break the same way. They do not.",
        "",
    ]

    res = detail("Q1")
    L += ["## Q1 — all three residuals move toward their locked values", "",
          "| residual | locked | legacy | corrected | legacy off | corrected off | drift |",
          "|--|--|--|--|--|--|--|"]
    for r in res["rows"]:
        L.append(f"| `{r['residual']}` | `{r['locked']}` | "
                 f"`{r['legacy']:.6f}` | `{r['corrected']:.6f}` | "
                 f"`{r['legacy_rel_pct']:+.2f}%` | "
                 f"**`{r['corrected_rel_pct']:+.2f}%`** | "
                 f"`{r['drift_pct']:+.3f}%` |")
    L += ["", f"*{res['note'].capitalize()}.*", ""]

    two = detail("Q7")
    L += ["## Q7 — why the two sectors part company", "",
          "One barrier feeds both, and the correction moves it the same way "
          "for both. What differs is how hard each Hamiltonian leans on it.",
          "",
          "| sector | scalar | corrected residual | elasticity | linearised | **measured** |",
          "|--|--|--|--|--|--|"]
    for key in ("lepton", "quark"):
        s = two[key]
        L.append(f"| {key} (`m_{s['species']}`) | {s['scalar']} | "
                 f"`{s['residual_pct']:+.2f}%` | `{s['elasticity']:+.3f}` | "
                 f"`{s['linearised_error_pct']:+.2f}%` | "
                 f"**`{s['measured_error_pct']:+.2f}%`** |")
    L += ["",
          f"The elasticity ratio is `{two['elasticity_ratio']:.2f}`. "
          f"{two['why_the_verdicts_differ'].capitalize()}.",
          "",
          "> **A percentage agreement between a geometric quantity and a "
          "fitted knob means nothing until multiplied by the elasticity of "
          "what it feeds.**",
          "",
          f"**A number this round corrects.** The README row written in #271 "
          f"quoted `{two['readme_said_pct']:.0f}%` for the muon error at the "
          f"corrected `γ`. The elasticity and the residual in that row are "
          f"both right; the `9%` is not — it came from a generic illustration "
          f"in the #271 probe (*a half-percent error in γ is a nine-percent "
          f"error in the muon*) applied to a residual that is not a half "
          f"percent. Linearising gives "
          f"`{two['lepton']['linearised_error_pct']:+.2f}%`; the locked block "
          f"returns `{two['corrected_to_pct']:+.2f}%`, and #271's own A/B/C "
          f"table carried that value the whole time.", ""]

    lad = detail("Q2")
    L += ["## Q2 — the composite is not the sum of its parts", "",
          "Every substitution below goes into the **locked** quark "
          "Hamiltonian with **nothing retuned**, `d`-anchored, max relative "
          "error over `{s, c, b, t}`.", "",
          "| substitution | legacy-derived | corrected-derived |",
          "|--|--|--|"]
    for knob in ("pinhole", "transport", "resistance"):
        L.append(f"| `{knob}` alone | "
                 f"`{100*lad['one_at_a_time_max_rel_err']['legacy'][knob]:.2f}%` | "
                 f"`{100*lad['one_at_a_time_max_rel_err']['scalar_correct'][knob]:.2f}%` |")
    L.append(f"| **all three** | "
             f"`{100*lad['composite_max_rel_err']['legacy']:.2f}%` | "
             f"**`{100*lad['composite_max_rel_err']['scalar_correct']:.2f}%`** |")
    L.append(f"| *the fitted v3 lock* | "
             f"`{100*lad['locked_max_rel_err']:.2f}%` | "
             f"`{100*lad['locked_max_rel_err']:.2f}%` |")
    L += ["",
          "The ordering **reverses**: each corrected residual is individually "
          "closer to its lock, yet the corrected triple scores worse. The "
          "legacy triple enjoys a partial cancellation its pinhole error alone "
          "would not earn; the corrected triple has none.",
          "",
          "> **\"Each knob is derived to within 1 %\" and \"the derived knobs "
          "reproduce the ladder\" are different claims, and only the first was "
          "ever established** — under either operator. The correction did not "
          "create that gap, it exposed it.", ""]

    ela = detail("Q3")
    L += ["## Q3 — the elasticity, and the pinhole the ladder actually wants",
          "", "| species | local `d ln m/d ln pinhole` at the lock | secant "
          "across the two derived pinholes |", "|--|--|--|"]
    for s in ("s", "c", "b", "t"):
        L.append(f"| `{s}` | "
                 f"`{ela['local_d_ln_m_d_ln_pinhole_at_the_lock'][s]:+.4f}` | "
                 f"`{ela['secant_across_the_two_derived_pinholes'][s]:+.4f}` |")
    d = ela["distance_from_what_it_wants_pct"]
    L += ["",
          f"Scanning the locked ladder over the pinhole alone, it is minimised "
          f"at `{ela['pinhole_the_ladder_wants']:.6f}` "
          f"(`{100*ela['ladder_error_at_what_it_wants']:.3f}%`), not at the "
          f"fitted `22.25` (`{100*ela['ladder_error_at_the_v3_lock']:.3f}%`). "
          f"Measured against the ladder's **own** optimum rather than a fitted "
          f"knob, the correction roughly halves the miss:", "",
          "| pinhole | distance from what the ladder wants |", "|--|--|",
          f"| v3 lock `22.25` | `{d['v3_lock']:+.3f}%` |",
          f"| legacy `22.008` | `{d['legacy']:+.3f}%` |",
          f"| corrected `22.331` | **`{d['corrected']:+.3f}%`** |", ""]

    brk = detail("Q4")
    L += ["## Q4 — the one thing that genuinely improves: cross-sector `R_OUTER`",
          "",
          "Each sector, read as a demand on `R_OUTER` through its own locked "
          "scalar:", "",
          "| operator | lepton (`γ=22.5`) | quark (`pinhole=22.25`) | split | "
          "brackets `1.26`? |", "|--|--|--|--|--|"]
    for name, label in (("legacy", "legacy"), ("scalar_correct", "corrected")):
        b = brk["per_operator"][name]
        L.append(f"| {label} | `{b['lepton_r_outer']:.5f}` "
                 f"(`{b['lepton_vs_canonical_pct']:+.2f}%`) | "
                 f"`{b['quark_r_outer']:.5f}` "
                 f"(`{b['quark_vs_canonical_pct']:+.2f}%`) | "
                 f"`{b['split_pct_of_mean']:.3f}%` | "
                 f"{'**yes**' if b['brackets_canonical'] else 'no'} |")
    corr = brk["per_operator"]["scalar_correct"]
    L += ["",
          f"Under the legacy operator both sectors demand **more** than `1.26`, "
          f"putting the canonical value outside the bracket they define. Under "
          f"the corrected operator they **straddle** it — the quark sector "
          f"wants less, the lepton sector more, and `1.26` sits at "
          f"`{corr['canonical_position_in_bracket']:.2f}` of the way across.",
          "",
          f"This is evidence *for* `R_OUTER = 1.26`, and it is **not** the "
          f"evidence #271 reopened: that was a single-sector fixed point, this "
          f"is a two-sector bracket. It is also weak — {brk['caveat']}. "
          f"Recorded as suggestive, not as a restored derivation.", ""]

    sel = detail("Q5")
    L += ["## Q5 — the last lepton claim that reads a radial eigenvalue", "",
          "README: *\"Resistance = 7π/100 — selected over `4·(ω−1)` by "
          "`R_OUTER` bisection.\"*", "",
          "| candidate | value | off the locked `0.217869` |", "|--|--|--|",
          f"| `7π/100` (operator-independent) | "
          f"`{sel['seven_pi_over_100']:.6f}` | "
          f"`{sel['closed_form_rel_pct']:+.3f}%` |"]
    for name, label in (("legacy", "legacy"), ("scalar_correct", "corrected")):
        r = sel["per_operator"][name]
        L.append(f"| `4(ω−1)`, {label} | `{r['four_omega_minus_one']:.6f}` | "
                 f"`{r['competitor_rel_pct']:+.3f}%` |")
    L += ["",
          f"Under the legacy operator the **rejected** candidate fitted "
          f"*better* than the selected one, so the choice rested entirely on "
          f"the `R_OUTER` bisection — which #271 reopened. Under the "
          f"correction the competitor degrades to "
          f"`{sel['per_operator']['scalar_correct']['competitor_rel_pct']:+.2f}%` "
          f"and `7π/100` wins on proximity too.", "",
          f"> **{sel['verdict'].capitalize()}.**", ""]

    nst = detail("Q6")
    L += ["## Q6 — `N` still drifts", "",
          "| residuals pinned to | `N` range | width |", "|--|--|--|"]
    for mode, r in nst["per_mode"].items():
        L.append(f"| {mode} | `[{r['range'][0]}, {r['range'][1]}]` | "
                 f"`{r['width']}` |")
    L += ["",
          f"{nst['verdict'].capitalize()}. Pinning the residuals to *corrected* "
          f"geometry does not stabilise `N` any more than pinning them to "
          f"legacy geometry did. The v3 lock's own conclusion stands unchanged; "
          f"only the digits move, `N` at PDG masses included.", ""]

    led = detail("Q8")
    L += ["## Q8 — the ledger", "",
          "| claim | verdict | evidence |", "|--|--|--|"]
    for e in led["entries"]:
        L.append(f"| {e['claim']} | **{e['verdict']}** | {e['evidence']} |")
    L += ["", f"**Headline.** {led['headline'].capitalize()}.", "",
          "**Not re-run, and why.** The four quark shell-index axioms are "
          "expressible in `k₅ = 5` alone, the Hopf transport phase is "
          "topological, and the flavor-CP layer inherits the v3 eigenvalues by "
          "construction. None of them reads the radial scalar operator. "
          "Proximity is not dependence.",
          "",
          "**Still open.** Nothing here derives `22.25`, `0.54` or `0.14` — "
          "each remains a fitted knob that a geometric quantity lands near. "
          "The composite gap (Q2) is the sharp successor target: a residual "
          "set that is individually right and jointly wrong points at a "
          "missing correlation between the three, not at three independent "
          "near-misses.", ""]
    return "\n".join(L)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, tuple):
        return list(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = (Path(__file__).resolve().parent / "runs"
           / f"{stamp}_quark_residual_reaudit_probe")
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2,
                                               default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
