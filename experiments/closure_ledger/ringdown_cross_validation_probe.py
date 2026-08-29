"""Settling PR #270's ringdown cross-validation.

THE VERDICT
───────────
PR #270 built two horizon-penetrating time-domain codes for a test scalar on a
fixed D=5 Tangherlini background. Both stable, both converged, and they
disagreed: real parts within 0.3% at l=1, DAMPING RATES DIFFERING BY 37%. #270
refused to quote a frequency -- correctly, "a converged number is not a correct
number" -- and named its own prime suspect: the Kerr-Schild inner cut.

THE KERR-SCHILD CODE WAS RIGHT. THE TORTOISE CODE'S DAMPING WAS WRONG.

An independent Gundlach-Price-Pullin characteristic evolution, written from
scratch and sharing no code with either, gives at l=1

    this round (h -> 0)   1.01612 - 0.36244i
    #270 Kerr-Schild      1.01622 - 0.36231i     agrees to ~1e-4
    #270 tortoise         1.01876 - 0.26404i     damping 37% away

Agreement to 1e-4 between independent implementations, against a 37% miss, is
not ambiguous. #270's suspicion pointed at the wrong code.

WHY D=5 MAKES THE ANSWER CHECKABLE
──────────────────────────────────
Two exact facts do most of the work.

  * r* = r + (1/2)ln((r-1)/(r+1)) -> r, with NO logarithmic tail (unlike 4D),
    so the far field is a pure Hankel wave with no Coulomb phase.
  * V -> [(l+1)^2 - 1/4]/r^2 exactly -- the same flat-limit Bessel identity
    PR #271 used to settle which radial operator was correct.

Together they give an EXACT eikonal asymptote to judge every solver against:
r_ph = sqrt(2), Omega_c = 1/2, lambda = 1/sqrt(2), so

    omega -> 0.5(l+1) - 0.353553i     as l -> infinity

and the characteristic evolution tracks it:

    l=0   0.5354 - 0.3842i        eikonal  0.5 - 0.35355i
    l=1   1.0161 - 0.3624i        eikonal  1.0 - 0.35355i
    l=2   1.5106 - 0.3575i        eikonal  1.5 - 0.35355i

Real parts approach 0.5(l+1) from above, damping approaches -0.35355. That is
what a correct solver must do. THE TORTOISE DAMPING -0.264 IS NOT NEAR THAT
ASYMPTOTE AND NEVER APPROACHES IT.

FOUR LINES OF EVIDENCE, ALL AGAINST THE TORTOISE DAMPING
────────────────────────────────────────────────────────
    characteristic evolution   -0.36244   (this round, converged in h)
    #270 Kerr-Schild           -0.36231
    first-order WKB            -0.36095
    exact eikonal asymptote    -0.35355
    ---------------------------------------
    #270 tortoise              -0.26404   <- excluded

WHAT DID NOT WORK, REPORTED AS SUCH
───────────────────────────────────
FREQUENCY-DOMAIN SHOOTING FAILS, and this round REPRODUCED the failure rather
than fixing it. Matching an ingoing Frobenius series at the horizon to
sqrt(r) H^(1) at large r gives a root that moves with every knob: 1.229-0.152i
-> 1.173-0.102i -> 1.155-0.214i across eps. The QNM boundary-value problem is
exponentially ill-conditioned in real r -- for Im w < 0 the outgoing piece grows
like e^{|Im w| R} and swamps the incoming coefficient one is trying to zero.
#270's diagnosis stands.

SIXTH-ORDER WKB IS UNUSABLE BY FINITE DIFFERENCES. It needs V^(6) at the peak,
and central differences amplify roundoff as h^-6: refining the grid makes the
answer DIVERGE (9.0 -> 18.6 -> 623). First-order WKB is well conditioned and is
what is used, with its accuracy STATED: 0.4% on the damping, 13% on the real
part at l=1, improving to 6.8% at l=2. Low-l WKB is simply inaccurate on the
real part -- a known limitation, not a discrepancy between solvers.

SCOPE
─────
A test scalar on a FIXED exact Tangherlini background -- no backreaction.
Fundamental (n=0) mode only. l=0 carries a wider uncertainty: its barrier is
weakest and the power-law tail contaminates the fit earliest.

The autopsy is NOT available: neither #270 code was landed in the tree, only
their reported numbers. This round establishes WHICH number is right and by how
much the other is excluded; it cannot say which line of the unlanded code did it.
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

from geometrodynamics.tangherlini.ringdown import (
    measure_against_the_published_spectrum,
    measure_the_extraction_systematics,  # noqa: E402
    measure_the_background_asymptotics_are_exact,
    measure_the_characteristic_scheme_converges,
    measure_the_cross_validation_verdict, measure_the_eikonal_invariants,
    measure_the_fundamental_modes, measure_the_ringdown_ledger,
    measure_what_did_not_work)


def run_probe() -> dict:
    checks: List[dict] = []

    asy = measure_the_background_asymptotics_are_exact()
    checks.append({
        "id": "R0", "name": "the D=5 asymptotics are exact: no log tail, Bessel V",
        "detail": asy,
        "pass": bool(asy["no_logarithmic_tail"]
                     and asy["the_tail_is_asymptotically_bessel"])})

    eik = measure_the_eikonal_invariants()
    checks.append({
        "id": "R1", "name": "the eikonal invariants are exact",
        "detail": eik,
        "pass": bool(eik["r_photon_squared_is_exactly_two"]
                     and eik["omega_c_is_exactly_one_half"]
                     and eik["lyapunov_is_one_over_sqrt_two"])})

    cnv = measure_the_characteristic_scheme_converges()
    checks.append({
        "id": "R2", "name": "the independent solver converges in step size",
        "detail": cnv, "pass": bool(cnv["all_converging"])})

    fun = measure_the_fundamental_modes()
    checks.append({
        "id": "R3", "name": "the frequencies track the exact eikonal asymptote",
        "detail": fun,
        "pass": bool(fun["every_real_part_sits_above_the_eikonal"]
                     and fun["every_damping_within_10_percent_of_the_asymptote"])})

    ver = measure_the_cross_validation_verdict()
    checks.append({
        "id": "R4",
        "name": "*** Kerr-Schild confirmed; the tortoise damping excluded ***",
        "detail": ver,
        "pass": bool(ver["kerr_schild_is_confirmed"]
                     and ver["tortoise_damping_is_excluded"])})

    neg = measure_what_did_not_work()
    checks.append({
        "id": "R5", "name": "two honest negatives, reproduced not papered over",
        "detail": neg,
        "pass": bool(neg["frequency_domain_shooting"]["so_pr_270s_diagnosis_stands"])})

    pub = measure_against_the_published_spectrum()
    checks.append({
        "id": "R7",
        "name": "*** the spectrum matches a published high-precision "
                "calculation ***",
        "detail": pub,
        "pass": bool(pub["every_mode_within_0p3_percent"]
                     and pub["ell_1_and_2_within_0p05_percent"])})

    sysm = measure_the_extraction_systematics()
    checks.append({
        "id": "R8",
        "name": "where the error floor actually lives, by varying knobs",
        "detail": sysm,
        "pass": bool(sysm["the_window_dominates"]
                     and sysm["extraction_band_exceeds_step_refinement"])})

    led = measure_the_ringdown_ledger()
    checks.append({"id": "R6", "name": "the ledger", "detail": led,
                   "pass": bool(len(led["entries"]) >= 7)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def _gap(ver: dict, key: str) -> float:
    """Damping error of `ver[key]` against the published reference, in percent."""
    reference = ver["published_reference"][1]
    return 100.0 * abs(ver[key][1] - reference) / abs(reference)


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    ver = detail("R4")
    L: List[str] = [
        "# Settling PR #270's ringdown cross-validation", "",
        f"**{summary['passed']}/{summary['total']} checks pass.**", "",
        "PR #270 built two time-domain codes for a test scalar on a fixed "
        "`D = 5` Tangherlini background. Both stable, both converged, and they "
        "disagreed — real parts within `0.3%`, and damping rates apart by **37% of "
        "the smaller value** (equivalently, the wrong one was **27% low**). "
        "#270 refused to quote a frequency, correctly, and named its "
        "own prime suspect: the Kerr–Schild inner cut.", "",
        "> ## The Kerr–Schild code was right. The tortoise code's damping was wrong.",
        "",
        "All damping errors below are relative to the **published** value.", "",
        "| source | `ω` at `ℓ = 1` | damping error |", "|--|--|--|",
        f"| **published** (continued fractions + Hill) | "
        f"`{ver['published_reference'][0]:.8f} "
        f"{ver['published_reference'][1]:+.8f}i` | — |",
        f"| #270 Kerr–Schild | "
        f"`{ver['pr_270_kerr_schild'][0]:.5f} "
        f"{ver['pr_270_kerr_schild'][1]:+.5f}i` | "
        f"**`{_gap(ver, 'pr_270_kerr_schild'):.3f}%`** |",
        f"| **this round** (characteristic, `h → 0`) | "
        f"`{ver['this_round_characteristic'][0]:.5f} "
        f"{ver['this_round_characteristic'][1]:+.5f}i` | "
        f"**`{_gap(ver, 'this_round_characteristic'):.3f}%`** |",
        f"| #270 tortoise | `{ver['pr_270_tortoise'][0]:.5f} "
        f"{ver['pr_270_tortoise'][1]:+.5f}i` | "
        f"**`{_gap(ver, 'pr_270_tortoise'):.1f}%`** ← excluded |", "",
        "An independent Gundlach–Price–Pullin evolution, written from scratch "
        "and sharing no code with either #270 code, lands on a spectrum "
        "computed elsewhere by continued fractions. **#270's suspicion pointed "
        "at the wrong code.**", "",
        "> **Note the ordering.** #270's Kerr–Schild code is ~6× *more* accurate "
        "than this round's solver. The characteristic scheme arbitrated "
        "correctly; it did not out-resolve.", "",
        "**Naming the denominator.** The tortoise damping is "
        f"`{ver['the_denominator_is_named']['tortoise_relative_error_against_published']:.1f}%` "
        "off in the conventional sense (divided by the published value), and "
        "equivalently the correct damping is "
        f"`{ver['the_denominator_is_named']['correct_damping_is_larger_than_tortoise_by']:.1f}%` "
        "*larger* than the tortoise value (divided by the tortoise value). "
        "PR #270 measured the latter because it had two codes and no reference. "
        "Both are true; neither should be quoted bare.", "",
    ]

    asy = detail("R0")
    L += ["## R0 — why `D = 5` makes the answer checkable", "",
          "Two exact facts do most of the work.", "",
          "**No logarithmic tail.** The correction has the exact closed form "
          "`r* − r = −artanh(1/r) = −1/r − 1/(3r³) − ⋯`, so `−1/r` is the "
          "*leading behaviour*, not an equality. Every term decays — unlike "
          "4D's growing `2M ln r` — so the far field is a pure Hankel wave with "
          "no Coulomb phase. The series predicts the next coefficient, and it "
          "is checked rather than assumed:", "",
          "| `r` | `r(r*−r) + 1` | predicted `−1/(3r²)` | rel. error |",
          "|--|--|--|--|"]
    for radius, scaled, pred, err in zip(
            [50.0, 200.0, 1000.0], asy["tortoise_minus_r_times_r"],
            asy["predicted_next_term_minus_one_over_three_r_squared"],
            asy["next_coefficient_relative_errors"]):
        L.append(f"| `{radius:.0f}` | `{scaled + 1.0:.4e}` | `{pred:.4e}` | "
                 f"`{err:.1e}` |")
    L += ["",
          "**The tail is *asymptotically* Bessel.** Exactly, "
          "`V = L/r² + (9/4−L)/r⁴ − (9/4)/r⁶` with `L = (ℓ+1)² − ¼`, so "
          "`√r H⁽¹⁾` is the *leading* outgoing solution, not the exact one at "
          "finite `r`. `V·r²` against `L`:", "",
          "| `ℓ` | limit | `V·r²` at `r = 1000` | rel. error |", "|--|--|--|--|"]
    for b in asy["bessel_tail"]:
        L.append(f"| {b['ell']} | `{b['limit']}` | "
                 f"`{b['V_times_r2_at_far_radii'][-1]:.6f}` | "
                 f"`{b['relative_error_at_1000']:.1e}` |")
    L += ["", f"*{asy['the_same_identity_that_fixed_the_operator']}*",
          ""]

    eik = detail("R1")
    L += ["## R1 — an exact asymptote to judge every solver against", "",
          "| quantity | value | exact? |", "|--|--|--|",
          f"| `r_ph²` | `{eik['r_photon_squared']:.15f}` | **2 exactly** |",
          f"| `Ω_c` | `{eik['omega_c']:.15f}` | **1/2 exactly** |",
          f"| `λ` | `{eik['lyapunov']:.6f}` | **1/√2** |", "",
          f"> `{eik['the_law']}`", ""]

    cnv = detail("R2")
    L += ["## R2 — the independent solver converges", "",
          "Gundlach–Price–Pullin on the null grid `u = t − r*`, `v = t + r*`. "
          "**No spatial boundary conditions are applied or needed** — the "
          "domain of dependence is the null diamond, so the horizon and "
          "infinity are limits rather than boundaries. That is exactly why this "
          "construction is immune to the excision question that both the "
          "frequency-domain shoot and the Kerr–Schild inner cut raise.", "",
          "| `ℓ` | " + " | ".join(f"`h = {s}`" for s in cnv["steps"])
          + " | successive `\\|Δ\\|` |", "|--|" + "--|" * (len(cnv["steps"]) + 1)]
    for row in cnv["rows"]:
        cells = []
        for s in cnv["steps"]:
            v = row["by_step"][str(s)]
            cells.append("—" if v is None else f"`{v[0]:.7f} {v[1]:+.7f}i`")
        L.append(f"| {row['ell']} | " + " | ".join(cells) + " | "
                 + ", ".join(f"`{d:.1e}`" for d in row["successive_differences"])
                 + " |")
    L.append("")

    fun = detail("R3")
    L += ["## R3 — the frequencies, against the exact asymptote", "",
          "| `ℓ` | characteristic | eikonal | 1st-order WKB |", "|--|--|--|--|"]
    for row in fun["rows"]:
        w = row["omega"]
        s = "—" if w is None else f"`{w[0]:.6f} {w[1]:+.6f}i`"
        L.append(f"| {row['ell']} | {s} | "
                 f"`{row['eikonal'][0]:.6f} {row['eikonal'][1]:+.6f}i` | "
                 f"`{row['wkb_first_order'][0]:.6f} "
                 f"{row['wkb_first_order'][1]:+.6f}i` |")
    L += ["", f"**{fun['the_pattern_a_correct_solver_must_show']}** "
          f"*{fun['ell_zero_is_least_certain']}*", ""]

    L += ["## R4 — four lines of evidence, all against the tortoise damping", "",
          "| source | `Im ω` at `ℓ = 1` |", "|--|--|"]
    for name, val in ver["damping_lines_of_evidence"].items():
        mark = "  ← **excluded**" if "tortoise" in name else ""
        L.append(f"| {name} | `{val:+.5f}`{mark} |")
    L += ["", f"> **{ver['verdict']}**", "",
          f"**What this round cannot do.** {ver['what_this_round_cannot_do']}",
          ""]

    neg = detail("R5")
    L += ["## R5 — what did not work", "", "**Frequency-domain shooting.** "
          f"{neg['frequency_domain_shooting']['status']}. Roots across `ε`: "
          + ", ".join(f"`{r}`" for r in
                      neg["frequency_domain_shooting"]["roots_across_epsilon"])
          + f". {neg['frequency_domain_shooting']['why']}", "",
          "**Sixth-order WKB.** "
          f"{neg['sixth_order_wkb']['status']}. Under refinement: "
          + " → ".join(f"`{v}`" for v in
                       neg["sixth_order_wkb"]["values_under_refinement"])
          + f". {neg['sixth_order_wkb']['why']}", "",
          "**First-order WKB accuracy**, stated rather than assumed: damping "
          f"{neg['first_order_wkb_accuracy']['damping_at_ell_1']}, real part "
          f"{neg['first_order_wkb_accuracy']['real_part_at_ell_1']} at `ℓ = 1`, "
          f"{neg['first_order_wkb_accuracy']['real_part_at_ell_2']} at `ℓ = 2`. "
          f"{neg['first_order_wkb_accuracy']['reading']}", ""]

    pub = detail("R7")
    L += ["## R7 — against a published high-precision spectrum", "",
          "The strongest check available, and external to this repository "
          "entirely. Source: " + pub["source"] + ".", "",
          f"| `ℓ` | characteristic (`h = {pub['step_size_of_the_rows']}`) | "
          "published | real err | damping err |",
          "|--|--|--|--|--|"]
    for row in pub["rows"]:
        if row["characteristic"] is None:
            continue
        L.append(
            f"| {row['ell']} | `{row['characteristic'][0]:.6f} "
            f"{row['characteristic'][1]:+.6f}i` | "
            f"`{row['published'][0]:.8f} {row['published'][1]:+.8f}i` | "
            f"`{100*row['real_relative_error']:.4f}%` | "
            f"`{100*row['damping_relative_error']:.4f}%` |")
    L += ["", "**" + pub["the_reframing"] + "**", "",
          "*" + pub["who_is_most_accurate"] + "*", ""]

    ref = pub["refinement_versus_truth"]
    L += ["### What the reference exposed about this solver", "",
          "| `h` | `Im ω` | distance to published |", "|--|--|--|"]
    for step, value, distance in zip(ref["steps"], ref["damping_by_step"],
                                     ref["distance_to_published_by_step"]):
        if value is None:
            continue
        L.append(f"| `{step}` | `{value:.7f}` | `{distance:.2e}` |")
    L += ["",
          f"Last successive difference `{ref['last_successive_difference']:.2e}`, "
          f"actual distance from the finest value to the published one "
          f"`{ref['distance_from_finest_to_published']:.2e}` — "
          f"**understated by `{ref['understatement_factor']:.1f}×`**, and the "
          "`h = 0.05` value lands *closer* than the `h = 0.025` one.", "",
          "> " + pub["the_lesson"], ""]

    sysm = detail("R8")
    L += ["## R8 — where the error floor actually lives", "",
          "An earlier draft *asserted* the residual was extraction "
          "systematics. Nothing was varied behind that. Varying them:", "",
          "| extraction window | damping rel. error |", "|--|--|"]
    for row in sysm["by_extraction_window"]:
        err = row["damping_relative_error"]
        L.append(f"| `{tuple(row['window'])}` | "
                 + (f"`{100*err:.4f}%`" if err is not None else "—") + " |")
    L += ["", "| observer `r*` | damping rel. error |", "|--|--|"]
    for row in sysm["by_observer_radius"]:
        err = row["damping_relative_error"]
        L.append(f"| `{row['observer_rstar']}` | "
                 + (f"`{100*err:.4f}%`" if err is not None else "—") + " |")
    band = sysm["band_over_reasonable_choices"]
    L += ["",
          f"**The window dominates** by orders of magnitude — late windows "
          f"admit the power-law tail. **`t_max` is bit-irrelevant** "
          f"({'confirmed' if sysm['t_max_is_irrelevant'] else 'NOT confirmed'}), "
          "because the extraction window sits well inside it.", "",
          f"Over reasonable choices the band is "
          f"`{100*band['min']:.4f}%`–`{100*band['max']:.4f}%`, against a "
          f"step-refinement difference of "
          f"`{100*sysm['step_refinement_last_relative_difference']:.4f}%` — "
          f"**{sysm['how_many_times_larger']:.1f}× larger**. That is why step "
          "refinement alone was the wrong error bar.", "",
          "> " + sysm["what_this_corrects"], ""]

    led = detail("R6")
    L += ["## R6 — the ledger", "", "| claim | verdict | evidence |", "|--|--|--|"]
    for e in led["entries"]:
        L.append(f"| {e['claim']} | **{e['verdict']}** | {e['evidence']} |")
    L += ["", f"**The standing lesson held.** {led['the_standing_lesson_held']}",
          "", f"**The lesson this round adds.** {led['the_lesson_this_round_adds']}",
          "", f"**Still open.** {led['still_open']}",
          "", f"**The next object.** {led['the_next_object']}", ""]
    return "\n".join(L)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, tuple):
        return list(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = (Path(__file__).resolve().parent / "runs"
           / f"{stamp}_ringdown_cross_validation_probe")
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(json.dumps(summary, indent=2,
                                               default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0 if summary["passed"] == summary["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
