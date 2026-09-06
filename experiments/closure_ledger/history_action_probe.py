"""Probe: does a classical BAM action select the oriented branch?

Against `docs/history_action_prereg.md` (commit `a33a901`), written before
`geometrodynamics/bulk/history_action.py`.
Run:  python -m experiments.closure_ledger.history_action_probe
"""

from __future__ import annotations

import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from geometrodynamics.bulk import history_action as ha   # noqa: E402


def run_probe() -> dict:
    oracle = ha.morse_bott_oracle()
    degen = ha.class_function_degeneracy()
    addit = ha.additive_functionals_have_no_critical_points()
    kappas = {k: ha.saddle_branch_ratio(v) for k, v in
              {"pi/4": math.pi / 4, "3pi/4": 3 * math.pi / 4,
               "1 (hbar=1)": 1.0, "pi/2": math.pi / 2}.items()}
    masses = ha.morse_bott_component_masses(
        [0.0, 0.0, 1.0], [math.sin(1.0), 0.0, math.cos(1.0)])
    offcrit = ha.no_off_closure_critical_points(trials=6)
    amp = ha.amplitude_dependence()
    excis = ha.excision_estimate()
    discrete = ha.discrete_symmetry_extension()
    orbits = ha.sector_orbits()
    fibre = ha.fibre_action_is_weight_blind()
    homog = ha.detector_response_homogeneity()
    quad = ha.quadratic_readouts_disagree()
    locsq = ha.local_square_mean_is_closed_form()
    radial = ha.radial_action_compatibility()
    gate = ha.source_observable_signalling()
    v = ha.verdicts()

    checks = [
        ("O1-O4  the frozen oracle reproduces", oracle["passes"]),
        ("O4     saddle magnitude IS the positive coarea density",
         oracle["O4_saddle_magnitude_is_coarea_rel"] < 1e-3),
        ("A-ctl  every F(cos theta) shares the critical set",
         degen["all_share_the_critical_set"]),
        ("A1     theta additive, S_H not",
         addit["theta_additive"] and addit["S_H_not_additive"]),
        ("A2     the additive functional is nowhere stationary",
         addit["theta_has_no_critical_point"]),
        ("A3     Crit(S_H) = Gamma: no critical points OFF closure either",
         offcrit["no_off_closure_critical_points"]),
        ("F1     component masses are unequal; only the PHASE is undetermined",
         masses["masses_are_unequal"]
         and all(abs(abs(r["phase_factor"]) - 1.0) < 1e-12
                 for r in kappas.values())),
        ("F1b    the masses reproduce BOTH candidate aggregations exactly",
         masses["oriented_identity_residual"] < 1e-9
         and masses["positive_count_identity_residual"] < 1e-9),
        ("F1c    ...but CONDITIONAL on the round-5 Haar amplitude a = 1",
         amp["masses_are_conditional_on_the_measure"]),
        ("F1d    the singular punctures carry no mass: excised disc is O(eps^2)",
         excis["excision_is_safe"]),
        ("F2     phase factor real iff 4kappa/pi odd",
         kappas["pi/4"]["phase_factor_is_real"]
         and kappas["3pi/4"]["phase_factor_is_real"]
         and not kappas["1 (hbar=1)"]["phase_factor_is_real"]),
        ("B1     two sector orbits at every angle but pi/2",
         orbits["forced_only_at_right_angle"]),
        ("B2     r not forced at any CHSH angle",
         not orbits["forced_at_any_chsh_angle"]),
        ("B3     fibre symmetries cannot change sector weights",
         fibre["fibre_blind"]),
        ("B4     no IDENTIFIED discrete operation mixes like and unlike",
         discrete["no_identified_operation_mixes_like_and_unlike"]),
        ("C1     the six audited stress/flux quantities are degree 2",
         homog["all_quadratic"]),
        ("C2     but two ordinary quadratics disagree: no readout is named",
         quad["the_two_quadratics_disagree"]),
        ("C3     <D_s^2> = (1+c)(2+c) closed form", locsq["closed_form_holds"]),
        ("D1     the radial action is a genuine one-form integral",
         radial["radial_action_is_a_one_form_integral"]
         and radial["both_legs_nonzero"]),
        ("E1     source density is exactly antipodally even",
         gate["density_is_antipodally_even_residual"] < 1e-15),
        ("E2     odd means/signs blind, but full readout laws separate",
         gate["odd_means_and_signs_are_blind"]
         and not gate["odd_full_distributions_are_blind"]
         and gate["some_even_functions_separate_the_ensembles"]),
        ("E2b    but NOT every even observable separates (constants, x.x)",
         gate["not_every_even_observable_separates"]),
        ("E4     no operational readout is claimed",
         not gate["operational_readout_constructed"]
         and not gate["bam_couplings_shown_even_in_x"]),
        ("E3     non-coplanar supports are mutually singular",
         gate["non_coplanar_supports_mutually_singular"]),
    ]
    return {
        "oracle": oracle, "degeneracy": degen, "additivity": addit,
        "kappa": kappas, "masses": masses, "off_closure": offcrit,
        "amplitude": amp, "excision": excis, "discrete": discrete,
        "orbits": orbits, "fibre": fibre, "local_square": locsq,
        "homogeneity": homog, "quadratic": quad, "radial": radial,
        "gate": gate, "verdicts": v,
        "ledger": ha.dependency_ledger(),
        "checks": checks,
        "passed": sum(1 for _, ok in checks if ok), "total": len(checks),
    }


def render(s: dict) -> str:
    L = ["# Round 7 — does a classical BAM action select the oriented branch?", "",
         f"**{s['passed']}/{s['total']} checks pass.**", "",
         "## Five independent verdicts", ""]
    for k, val in s["verdicts"].items():
        L.append(f"* **{k}** — `{val}`")
    L += ["", "## Checks", ""]
    for name, ok in s["checks"]:
        L.append(f"* {'PASS' if ok else 'FAIL'}  {name}")

    o, a, q, g = s["oracle"], s["additivity"], s["quadratic"], s["gate"]
    L += ["", "## The candidate, and why it is not an action", "",
          f"* `-1/2 Tr G = -cos(Omega/2) = -D/sqrt(D^2+N^2)` to `{o['O1_closed_form']:.1e}`",
          f"* closure set is the critical manifold: `|grad S_H|` <= `{o['O2_grad_on_closure']:.1e}`",
          f"* transverse curvature `sgn(D)|uxv|^2/D^2`, rel `{o['O3_transverse_curvature_rel']:.1e}`;"
          f" index is `sgn D`: {o['O3_index_is_sign_D']}",
          f"* **saddle magnitude = positive coarea density**, rel "
          f"`{o['O4_saddle_magnitude_is_coarea_rel']:.1e}`",
          "",
          f"* `theta` additive to `{a['theta_is_additive_residual']:.1e}`;"
          f" `S_H` additivity defect `{a['S_H_additivity_defect']:.3f}`",
          f"* `min |grad theta|` on closure = `{a['min_|grad theta|_on_closure']:.4f}` > 0",
          ""]
    m = s["masses"]
    L += ["", f"* component masses `M_0 = {m['M_0']:.6f}`, `M_pi = {m['M_pi']:.6f}`,"
          f" ratio `{m['mass_ratio']:.6f}` — **not** 1",
          f"* `(M_0 - M_pi)|uxv| = int D` (residual `{m['oriented_identity_residual']:.1e}`)"
          f" and `(M_0 + M_pi)|uxv| = int |D|` (residual "
          f"`{m['positive_count_identity_residual']:.1e}`): stationary phase supplies"
          " both candidate magnitudes; only the relative phase is open",
          f"* conditional on the round-5 Haar amplitude `a = 1`: `1 + 0.5 x.u`"
          f" moves the oriented sum to `{s['amplitude']['rows'][1]['oriented']:.4f}`"
          f" from `{s['amplitude']['rows'][0]['oriented']:.4f}`",
          f"* the singular punctures carry no mass: excised disc `O(eps^2)`,"
          f" worst relative error `{s['excision']['worst_relative_error']:.2e}`",
          f"* no critical points off closure: `min |grad theta| = "
          f"{s['off_closure']['min_grad_theta_off_closure']:.4f}`,"
          " and the only candidate needs `x_p = -sec(gamma/2)`, outside the sphere",
          ""]
    L += ["| kappa | 4k/pi | phase real? | selects |", "|---|---|---|---|"]
    for name, r in s["kappa"].items():
        L.append(f"| {name} | {r['4kappa_over_pi']:.4f} | {r['phase_factor_is_real']} "
                 f"| {r['selects']} |")

    L += ["", "## B — sector orbits", "",
          "| gamma | group order | orbits | r forced |", "|---|---|---|---|"]
    for r in s["orbits"]["rows"]:
        L.append(f"| {r['gamma']:.4f} | {r['group_order']} | {r['n_orbits']} "
                 f"| {r['r_forced']} |")

    L += ["", "## C — measured homogeneity of the six audited quantities", "",
          "| coupling | where | degree |", "|---|---|---|"]
    for r in s["homogeneity"]["rows"]:
        L.append(f"| `{r['coupling']}` | `{r['where']}` | {r['measured_degree']:.6f} |")
    L += ["", f"* linear readout: `S_max = {q['linear']['S_max']:.6f}` (Tsirelson)",
          f"* quadratic, square-of-integral: `S_max = {q['square_of_integral']['S_max']:.6f}`"
          " = `8 sqrt2/3`",
          f"* quadratic, integral-of-square: `S_max = {q['integral_of_square']['S_max']:.6f}`",
          "* both keep marginals at exactly `1/2`; **they disagree**, so degree-2"
          " homogeneity does not name a readout"]

    L += ["", "## E — causality gate", "",
          f"* conditioned density antipodally even to `{g['density_is_antipodally_even_residual']:.1e}`",
          f"* odd means/signs blind (spread `{g['odd_observable_spread']:.1e}`); "
          f"projection variances separate (spread `{g['odd_projection_variance_spread']:.4f}`)",
          f"* *some* even functions separate (spread "
          f"`{g['some_even_observables_separate']:.4f}`); constants and `x.x` do not"
          f" (spread `{g['blind_even_observable_spread']:.1e}`)",
          "* no map from `x` to field configurations, and no two-boundary-compatible"
          " source readout, is constructed — a hazard, not a channel",
          f"* non-coplanar total variation `{g['non_coplanar_total_variation']}`"
          " — mutually singular, a one-shot signal",
          "", "## Dependency ledger", ""]
    for row in s["ledger"]:
        L.append(f"* `{row['input']}` — **{row['status']}** ({row['where']})")
    return "\n".join(L) + "\n"


def main() -> int:
    s = run_probe()
    text = render(s)
    print(text)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    outdir = os.path.join(os.path.dirname(__file__), "runs", f"{stamp}_history_action_probe")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "probe.json"), "w") as h:
        json.dump(s, h, indent=2, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else str(o))
    with open(os.path.join(outdir, "probe.md"), "w") as h:
        h.write(text)
    return 0 if s["passed"] == s["total"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
