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

from geometrodynamics.tangherlini.ringdown import (  # noqa: E402
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
                     and asy["the_tail_is_exactly_bessel"])})

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

    led = measure_the_ringdown_ledger()
    checks.append({"id": "R6", "name": "the ledger", "detail": led,
                   "pass": bool(len(led["entries"]) >= 7)})

    return {"checks": checks,
            "passed": sum(1 for c in checks if c["pass"]),
            "total": len(checks)}


def render_markdown(summary: dict) -> str:
    def detail(cid: str) -> dict:
        return next(c for c in summary["checks"] if c["id"] == cid)["detail"]

    ver = detail("R4")
    L: List[str] = [
        "# Settling PR #270's ringdown cross-validation", "",
        f"**{summary['passed']}/{summary['total']} checks pass.**", "",
        "PR #270 built two time-domain codes for a test scalar on a fixed "
        "`D = 5` Tangherlini background. Both stable, both converged, and they "
        "disagreed — real parts within `0.3%`, **damping rates differing by "
        "37%**. #270 refused to quote a frequency, correctly, and named its "
        "own prime suspect: the Kerr–Schild inner cut.", "",
        "> ## The Kerr–Schild code was right. The tortoise code's damping was wrong.",
        "",
        "| source | `ω` at `ℓ = 1` | gap to this round |", "|--|--|--|",
        f"| **this round** (characteristic, `h → 0`) | "
        f"`{ver['this_round_characteristic'][0]:.5f} "
        f"{ver['this_round_characteristic'][1]:+.5f}i` | — |",
        f"| #270 Kerr–Schild | "
        f"`{ver['pr_270_kerr_schild'][0]:.5f} "
        f"{ver['pr_270_kerr_schild'][1]:+.5f}i` | "
        f"**`{ver['gap_to_kerr_schild']['imag_percent']:.3f}%`** damping |",
        f"| #270 tortoise | `{ver['pr_270_tortoise'][0]:.5f} "
        f"{ver['pr_270_tortoise'][1]:+.5f}i` | "
        f"**`{ver['gap_to_tortoise']['imag_percent']:.1f}%`** damping |", "",
        "An independent Gundlach–Price–Pullin evolution, written from scratch "
        "and sharing no code with either, agrees with Kerr–Schild to `~1e-4` "
        "and misses the tortoise damping by 37%. **#270's suspicion pointed at "
        "the wrong code.**", "",
    ]

    asy = detail("R0")
    L += ["## R0 — why `D = 5` makes the answer checkable", "",
          "Two exact facts do most of the work.", "",
          f"**No logarithmic tail.** `r* − r` at `r = 50, 200, 1000` is "
          f"`{', '.join('%.1e' % v for v in asy['tortoise_minus_r_at_far_radii'])}` "
          f"— unlike 4D, `r* → r`, so the far field is a pure Hankel wave with "
          f"no Coulomb phase.", "",
          "**The tail is exactly Bessel.** `V·r²` against `(ℓ+1)² − ¼`:", "",
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

    led = detail("R6")
    L += ["## R6 — the ledger", "", "| claim | verdict | evidence |", "|--|--|--|"]
    for e in led["entries"]:
        L.append(f"| {e['claim']} | **{e['verdict']}** | {e['evidence']} |")
    L += ["", f"**The standing lesson held.** {led['the_standing_lesson_held']}",
          "", f"**Still open.** {led['still_open']}", ""]
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
