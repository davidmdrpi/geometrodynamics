"""One scalar field on one surface: the construction v66 should have used.

> Scope: a LINEAR scalar on a FIXED round background.  The crossing rule that
> glues R_outer to R_inner is a REPRESENTATION choice, not a derived boundary
> condition.  The gain is a DISPLAY amplitude -- and the fixed-energy check is
> the one place that distinction is made to do work.

THE CORRECTION
──────────────
v66 drew TWO curves and asked whether their images meet through the seam.  Two
curves in one frame are two surfaces, and reading their overlap as a connection
is a statement about a picture, not about a field.  If the two contributions are
two pieces of ONE scalar deformation of ONE surface there is only ever one curve

    u = s_A u_A + s_B u_B          r = R_mid + eps u

and the question is whether that single curve reaches R_outer at one theta and
R_inner at another.

THE REPAIR COSTS NOTHING, WHICH IS THE SURPRISE
────────────────────────────────────────────────
delta = r_A - r_B = eps(s_A u_A - s_B u_B) -- the quantity v66 plotted as a
"separation between two curves" -- IS the one-surface deformation with the
second sign flipped.  The same array, to R_mid's own rounding.  So every number
v66 reported survives, with the two configurations SWAPPING NAMES:

    v66 "like-signed"  ==  one surface OPPOSED
    v66 "opposed"      ==  one surface LIKE-signed

which inverts v66's headline: its cheapest-when-co-located result belongs to the
LIKE pair, and its identically-zero bisector is the node of the OPPOSED field.

WHAT IS CHECKED
───────────────
T1  *** THE REPAIR IS FREE. *** v66's separation and the one-surface field are
    the same array to 2 ulps of R_mid, at five offsets over a whole run.  The
    residue is the mid-radius, not the fields, and is reported as such.

T2  *** COINCIDENCE CANCELS, EXACTLY. *** alpha = 0 with opposite orientation
    gives u == 0 at every time -- zero, not small -- so no gain connects it.
    The like pair is unaffected.

T3  THE MONOCHROMATIC LAW.  B = 2A|sin(m alpha/2)|, measured on a grid and
    compared with the closed form; the first optimum is alpha* = pi/m, HALF A
    WAVELENGTH.  The antipode is simply the m = 1 member of that family.

T4  *** THE ANTIPODE IS PARITY-DEPENDENT. *** sin(m pi/2) is +-1 for odd m and
    0 for even m, so opposite foci are maximal for odd modes and CANCEL for
    even ones.  Checked twice: on the plane-wave factor, and on the exact S^3
    zonal harmonics through Z_n(pi) = (-1)^n.

T5  *** WHERE THE TWO MODELS PART COMPANY. *** A zonal harmonic is CENTRED
    (Z_n(0) = 1, |Z_n| <= 1), so the opposed pair obeys |Z_A - Z_B| <= 2 with
    equality only where one focus sees +1 and the other -1 -- exactly the
    antipode with odd n.  For zonal modes alpha* = pi for EVERY odd n and it
    SATURATES the bound; half a wavelength reaches only a fraction of it.  For
    even n the antipode cancels and nothing reaches the bound at all.  The
    parity carries across the two models; the LOCATION of the optimum does not.

T6  A PULSE IS NOT A MODE.  v46's field is a launched pulse, so the 1/|sin|
    divergence is confined to about its own width: past that the two pulses
    stop overlapping, nothing cancels, and the threshold SATURATES instead of
    continuing down.  The coincident cancellation is real; the law governing
    its approach is not.

T7  THE CHORD SHRINKS WITH FREQUENCY.  At the optimum the two extrema are pi/m
    apart, so L = sqrt(D^2 + 4 R_out R_in sin^2(pi/2m)) falls from 2.000 to the
    purely radial gap 0.520.  Same deformation, shorter connection.

T8  FIXED ENERGY IS NOT FIXED AMPLITUDE.  E ~ omega^2 A^2, so at fixed energy
    A ~ 1/omega and the attainable span falls as fast as the chord.  A
    frequency slider cannot hold displacement fixed and be read as physics.
    No favourable frequency is claimed -- that needs an energy normalisation
    and a packet focusing law this model does not contain.

    python -m experiments.closure_ledger.one_surface_probe
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional

import numpy as np

from geometrodynamics.viz.one_surface import (
    measure_a_pulse_saturates_where_a_mode_diverges,
    measure_how_one_surface_answers_two_fronts,
    measure_coincident_foci_cancel_exactly,
    measure_fixed_energy_reverses_the_preference,
    measure_the_antipode_is_parity_dependent,
    measure_the_chord_shrinks_with_frequency,
    measure_the_optimum_is_half_a_wavelength,
    measure_the_two_curve_picture_was_one_field_all_along,
    measure_the_zonal_optimum_is_the_antipode,
)


def run_probe() -> dict:
    checks: List[dict] = []

    same = measure_the_two_curve_picture_was_one_field_all_along()
    checks.append({"id": "T1", "name": "*** the repair costs nothing ***",
                   "detail": same,
                   "pass": bool(same["they_are_the_same_array"])})

    zero = measure_coincident_foci_cancel_exactly()
    checks.append({"id": "T2", "name": "*** coincidence cancels exactly ***",
                   "detail": zero,
                   "pass": bool(zero["the_opposed_field_is_identically_zero"]
                                and zero["no_gain_connects_it"]
                                and zero["the_like_pair_is_unaffected"])})

    half = measure_the_optimum_is_half_a_wavelength()
    checks.append({"id": "T3", "name": "the optimum is half a wavelength",
                   "detail": half,
                   "pass": bool(half["the_closed_form_is_the_measured_amplitude"]
                                and half["the_first_optimum_is_half_a_wavelength"])})

    par = measure_the_antipode_is_parity_dependent()
    checks.append({"id": "T4",
                   "name": "*** the antipode is parity-dependent ***",
                   "detail": par,
                   "pass": bool(par["odd_modes_are_maximal_at_the_antipode"]
                                and par["even_modes_cancel_at_the_antipode"]
                                and par["the_zonal_spectrum_agrees"]
                                and par["zonal_odd_doubles"]
                                and par["zonal_even_cancels"])})

    zon = measure_the_zonal_optimum_is_the_antipode()
    checks.append({"id": "T5",
                   "name": "*** the zonal optimum IS the antipode ***",
                   "detail": zon,
                   "pass": bool(zon["the_bound_is_two"]
                                and zon["odd_orders_peak_at_the_antipode"]
                                and zon["and_they_saturate_the_bound"]
                                and zon["half_a_wavelength_does_not_reproduce_it"]
                                and zon["even_orders_cancel_at_the_antipode"])})

    pul = measure_a_pulse_saturates_where_a_mode_diverges()
    checks.append({"id": "T6", "name": "a pulse saturates where a mode diverges",
                   "detail": pul,
                   "pass": bool(pul["the_pulse_threshold_saturates"]
                                and pul["it_still_cancels_at_coincidence"])})

    ch = measure_the_chord_shrinks_with_frequency()
    checks.append({"id": "T7", "name": "the chord shrinks with frequency",
                   "detail": ch,
                   "pass": bool(ch["the_closed_form_is_the_law_of_cosines"]
                                and ch["the_chord_shrinks_monotonically"]
                                and ch["the_span_is_the_same_at_every_mode"])})

    en = measure_fixed_energy_reverses_the_preference()
    checks.append({"id": "T8", "name": "fixed energy is not fixed amplitude",
                   "detail": en,
                   "pass": bool(en["span_is_flat_at_fixed_amplitude"]
                                and en["span_falls_like_one_over_omega"]
                                and en["the_chord_and_the_span_both_shrink"])})

    two = measure_how_one_surface_answers_two_fronts()
    checks.append({"id": "T9",
                   "name": "*** the offset turns cancellation OFF, not "
                           "interference on ***",
                   "detail": two,
                   "pass": bool(two["the_contributions_barely_overlap"]
                                and two["the_field_peak_sits_on_a_contribution_peak"]
                                and two["coincidence_is_the_only_total_overlap"])})

    return {
        "probe": "one_surface",
        "question": "read as ONE scalar field on ONE surface rather than two "
                    "overlaid curves, where do two oppositely-oriented "
                    "refocusing contributions best span the bulk gap?",
        "answer": "not at 'the antipode' in general -- at half a wavelength "
                  "for a travelling mode, and at the antipode for every ODD "
                  "zonal mode, which is where the two models disagree; "
                  "coincident foci cancel exactly in both",
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "total": len(checks),
    }


def render_markdown(summary: dict) -> str:
    d = {c["id"]: c["detail"] for c in summary["checks"]}
    lines = [
        "# One-surface probe — one field, one curve, one question",
        "",
        f"**Question.** {summary['question']}",
        "",
        f"**Answer.** {summary['answer']}",
        "",
        f"**{summary['passed']}/{summary['total']} checks pass.**",
        "",
        "| id | check | result |",
        "|----|-------|--------|",
    ]
    for c in summary["checks"]:
        lines.append(f"| {c['id']} | {c['name']} | "
                     f"{'PASS' if c['pass'] else 'FAIL'} |")

    lines += [
        "",
        "## The repair costs nothing",
        "",
        "| | |",
        "|--|--|",
        f"| v66 separation vs one-surface field | "
        f"`{d['T1']['residue_in_ulps']:.1f}` ulp of `R_mid` |",
        f"| the opposed field at `α = 0` | "
        f"`{d['T2']['worst_absolute_field']:.1e}` — exactly zero |",
        f"| its span threshold | `{d['T2']['opposed_span_gain']}` |",
        f"| the like pair's, unaffected | "
        f"`{d['T2']['like_signed_cheapest_span_gain']:.4f}` |",
        "",
        f"> {d['T1']['the_labels_swap']}",
        "",
        "## The monochromatic law, and the parity",
        "",
        "| mode `m` | `α*` = half a wavelength | at the antipode | verdict |",
        "|--|--|--|--|",
    ]
    by_mode = {r["mode"]: r for r in d["T3"]["rows"]}
    for r in d["T4"]["rows"]:
        opt = by_mode.get(r["mode"])
        oc = (f"`{opt['first_optimum_over_pi']:.4f}π`" if opt else "—")
        lines.append(f"| `{r['mode']}` | {oc} | "
                     f"`{r['plane_wave_factor_at_pi']:.4f}` | "
                     f"**{r['verdict']}** |")
    lines += [
        "",
        f"> {d['T4']['why']}",
        "",
        "## Where the two models part company",
        "",
        "| zonal order `n` | `ω = n+1` | measured `α*` | peak strength |"
        " at the antipode | at half a wavelength |",
        "|--|--|--|--|--|--|",
    ]
    for r in d["T5"]["rows"]:
        lines.append(
            f"| `{r['order']}` | `{r['omega']}` | "
            f"`{r['measured_optimum_over_pi']:.4f}π` | "
            f"`{r['peak_strength']:.4f}` | "
            f"`{r['strength_at_the_antipode']:.2e}` | "
            f"`{r['strength_at_half_wavelength']:.4f}` |")
    lines += [
        "",
        f"> {d['T5']['why']}",
        "",
        f"> {d['T5']['the_kernel_is_odd']}",
        "",
        "## A pulse is not a mode",
        "",
        f"Pulse width `{d['T6']['pulse_width']}` rad. The threshold plateaus at "
        f"`{d['T6']['plateau_threshold']:.4f}` "
        f"(spread `{d['T6']['plateau_spread']:.1e}`).",
        "",
        "| offset `α/π` | in pulse widths | monochromatic `A_req` |"
        " pulse threshold | ratio |",
        "|--|--|--|--|--|",
    ]
    for r in d["T6"]["rows"]:
        lines.append(f"| `{r['offset_over_pi']:.2f}` | "
                     f"`{r['offset_in_pulse_widths']:.2f}` | "
                     f"`{r['monochromatic_required_amplitude']:.4f}` | "
                     f"`{r['pulse_threshold']:.4f}` | "
                     f"`{r['ratio']:.2f}×` |")
    lines += [
        "",
        f"> {d['T6']['why']}",
        "",
        "## The chord, and what it costs",
        "",
        "| mode `m` | separation | span / `A` | bulk chord | `L/D` |"
        " span at fixed energy | spans the gap? |",
        "|--|--|--|--|--|--|--|",
    ]
    en = {r["mode"]: r for r in d["T8"]["rows"]}
    for r in d["T7"]["rows"]:
        e = en.get(r["mode"])
        lines.append(
            f"| `{r['mode']}` | `{r['optimal_separation_over_pi']:.4f}π` | "
            f"`{r['radial_span_over_A']:.4f}` | `{r['bulk_chord']:.4f}` | "
            f"`{r['chord_over_gap']:.3f}` | "
            + (f"`{e['span_at_fixed_energy']:.4f}` | "
               f"{'yes' if e['spans_the_gap_at_fixed_energy'] else '**no**'} |"
               if e else "— | — |"))
    lines += [
        "",
        f"The chord falls from `{d['T7']['chord_at_the_fundamental']:.4f}` to "
        f"`{d['T7']['chord_at_the_highest_mode']:.4f}`, with the limit the "
        f"purely radial gap `{d['T7']['limit_is_the_radial_gap']}`.",
        "",
        f"> {d['T8']['the_tradeoff']}",
        "",
        f"At fixed energy the highest mode that still spans the gap is "
        f"`m = {d['T8']['highest_mode_that_still_spans_at_fixed_energy']}`.",
        "",
        f"**Not claimed:** {d['T8']['what_is_not_claimed']}",
        "",
        "## Where each front sits on the one surface",
        "",
        "| offset `α/π` | at `t/π` | peak `c_A` | peak `u` | amplification |"
        " overlap arc | `σ` of `c_A` | `σ` of `c_B` | `σ` of `u` |",
        "|--|--|--|--|--|--|--|--|--|",
    ]
    for r in d["T9"]["rows"]:
        lines.append(
            f"| `{r['offset_over_pi']:.2f}` | `{r['at_time_over_pi']:.2f}` | "
            f"`{r['peak_contribution']:.4f}` | `{r['peak_field']:.4f}` | "
            f"`{r['amplification']:.4f}` | `{r['overlap_arc']:.3f}` | "
            f"`{r['sigma_of_peak_a_over_pi']:+.2f}π` | "
            f"`{r['sigma_of_peak_b_over_pi']:+.2f}π` | "
            f"`{r['sigma_of_peak_field_over_pi']:+.2f}π` |")
    lo, hi = d["T9"]["amplification_when_apart"]
    lines += [
        "",
        f"Apart, the total is `{lo:.3f}–{hi:.3f}×` **one** contribution, and it "
        "peaks exactly where a single contribution does.",
        "",
        f"> {d['T9']['what_the_offset_buys']}",
        "",
    ]
    return "\n".join(lines)


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    return str(o)


def main(argv: Optional[list] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = here / "runs" / f"{ts}_one_surface_probe"
    out.mkdir(parents=True, exist_ok=True)
    (out / "probe.json").write_text(
        json.dumps(summary, indent=2, default=_json_default))
    (out / "probe.md").write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main(sys.argv[1:]))
