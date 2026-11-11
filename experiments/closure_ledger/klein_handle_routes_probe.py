"""Two routes to the sheet pairing, on equal terms. Freeze ``fd78693``.

Reports separate verdict fields. Exits nonzero when a required numerical
check fails. A frozen prediction contradicted by a resolved measurement is
reported as such and is not a numerical failure.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import klein_handle_electric as ke
from geometrodynamics.bulk import klein_handle_response as kh
from geometrodynamics.bulk import klein_handle_routes as kr

PREREG = "fd78693"
SCOPE = ("Round 14's lattice electrodynamics on the doubly reversed seam. "
         "Trapped flux and a matter current compared at matched total "
         "energy, over 12 draws at each of two sizes, with gauge-invariant "
         "observables only. Classical throughout; lattice units; no norm, "
         "detector law or statement about which route BAM uses.")
ZERO, EXACT = 1e-9, 1e-12
SIZES = {"primary": kr.PRIMARY, "replicate": kr.REPLICATE}
CONTRADICTED = {
    "screening_rule":
        "Frozen P8 required the early flux change to match "
        "-1/2 g^2 Phi (B_0/|C|) t^2 within 5 % at t = 0.1 for every resolved "
        "draw. One replicate draw misses by 6.4 %. The sign rule holds in "
        "every resolved draw. The formula is the leading order: its error "
        "falls as t^2 and is below 1 % at t = 0.04. The tolerance was "
        "mis-set, and so was my first explanation of the miss.",
}


def measure(config: dict, long_steps: int, numerics: bool) -> dict:
    net = kh.Network(**config)
    model = ke.Electric(net)
    rng = np.random.default_rng(kr.SEED)
    ds = kr.draws(model, rng)
    out = {
        "config": config, "sites": net.size, "draws": len(ds),
        "loop_edges": kr.loop_edges(model),
        "P1": kr.rest_flux_report(model, ds),
        "P2_P3": kr.scaling_report(model, ds),
        "P4": kr.flux_fraction_report(model, ds),
        "P5": kr.coupling_report(net, ds),
        "P6": kr.sector_report(model, ds),
        "P7": kr.interference_report(model, ds),
        "P8": kr.screening_report(model, ds),
        "P8_leading_order": kr.screening_convergence_report(model, ds),
        "P9": kr.reversal_report(model, ds, rng),
        "long_window": kr.long_window_report(model, ds, long_steps),
    }
    if numerics:
        out["numerics"] = kr.numerical_report(model, ds[0], long_steps)
    return out


def evaluate(m: dict) -> dict:
    p1, sc, p4, p5, p6, p7, p8, lo, p9, lw = (
        m[k] for k in ("P1", "P2_P3", "P4", "P5", "P6", "P7", "P8",
                       "P8_leading_order", "P9", "long_window"))
    ratio = sc["ratio"]

    def band(values, centre, width):
        return all(abs(v - centre) < width for v in values)

    out = {
        "flux_route_requirement": {
            "P1 flux with the field at rest gives no sheet charge": (
                p1["max_sheet_charge_at_rest"] < EXACT
                and p1["deck_fixed_error"] < EXACT
                and p1["ct0_fixed_error"] < EXACT),
            "P1 the flux still evolves, and a moving field does give one": (
                p1["min_flux_change_at_rest"] > 1e-12
                and p1["min_sheet_charge_moving"] > 1e-10),
        },
        "energy_scaling": {
            "P2 matter-current sheet charge is linear in energy": (
                abs(sc["M"]["slope_of_median"] - 1) < 0.1
                and abs(sc["M0"]["slope_of_median"] - 1) < 0.1),
            "P2 flux-route sheet charge goes as energy to the three halves":
                abs(sc["F"]["slope_of_median"] - 1.5) < 0.15,
        },
        "short_time_ranking": {
            "P3 the ratio grows as the square root of energy": (
                abs(ratio["slope_of_median"] - 0.5) < 0.2),
            "P3 the median ratio is below one at all four energies": (
                len(ratio["median"]) == 4
                and all(r < 1 for r in ratio["median"])),
        },
        "flux_fraction": {
            "P4 ratios to the one-third value are within 10 %": all(
                abs(a / b - 1) < 0.1
                for a, b in zip(p4["median_ratio"], p4["predicted_ratio"])),
        },
        "coupling_dependence": {
            "P5 matter route does not depend on the coupling": band(
                p5["M0_slopes"], 0, 0.1) and len(p5["M0_slopes"]) == 2,
            "P5 flux route is first order in the coupling": band(
                p5["F_slopes"], 1, 0.1) and len(p5["F_slopes"]) == 2,
            "P5 seam flux is g times the sheet charge": (
                p5["gauss_identity_error"] < ZERO),
        },
        "electric_sectors": {
            "P6 matter current at rest has no loop flux, holonomy, loop "
            "current or swap-even field": all(
                p6["M0"][k] < EXACT for k in (
                    "handle_flux", "loop_holonomy", "loop_current",
                    "swap_even_field")),
            "P6 it does have sheet charge, and the flux route has loop flux":
                (p6["M0"]["sheet_charge"] > 1e-9
                 and p6["F"]["handle_flux"] > 1e-4),
        },
        "interference": {
            "P7 additivity defect is second order in the flux": (
                band(p7["defect_slopes"], 2, 0.3)
                and len(p7["defect_slopes"]) == 2),
            "P7 reversing the flux alone changes the history": all(
                v > 1e-9 for v in p7["median_flux_reversal_difference"]),
            "P7 reversing both orienting data is a symmetry": (
                p7["joint_reversal_error"] < EXACT),
        },
        "screening_rule": {
            "P8 sign of the early flux change is minus the sign of the bond "
            "order": (p8["resolved_draws"] > 0
                      and p8["sign_agreement"] == p8["resolved_draws"]),
            "P8 early flux change within 5 % of the formula (frozen "
            "tolerance)": p8["worst_ratio_error"] < 0.05,
        },
        "screening_leading_order": {
            "post-freeze: the formula's error falls as t squared": (
                abs(lo["error_slope"] - 2) < 0.3),
            "post-freeze: the error is below 1 % at t = 0.04": (
                lo["max_error"][0] < 0.01),
            "post-freeze: the rate of the bond order vanishes by parity": (
                lo["bond_rate_over_bond"] < EXACT),
            "post-freeze: the sign rule holds at all three times": (
                lo["sign_agreement_all_times"] == lo["resolved_draws"]),
        },
        "reversal_requirements": {
            "P9 matter current at rest is fixed by the time-preserving "
            "fibre-reversing deck": (
                p9["rest_data_fixed_by_alice_deck"] < EXACT
                and p9["rest_data_fixed_by_klein_deck"] < EXACT),
            "P9 its history stays in that fixed set and separates charge": (
                p9["alice_fixed_set_invariance"] < ZERO
                and p9["rest_sheet_charge"] > 1e-9),
            "P9 ordinary handle keeps the sheet charge at zero": (
                p9["ordinary_fixed_error"] < EXACT
                and p9["ordinary_gauss"] < EXACT
                and p9["ordinary_sheet_charge"] < EXACT
                and p9["ordinary_local_charge"] > 1e-7),
        },
        "long_time_ranking": {
            "P10 the long-window ratio exceeds the short-window ratio at "
            "every energy": (len(lw["rows"]) == 4 and all(
                r["ratio_long_rms"] > r["ratio_short_peak"]
                for r in lw["rows"])),
        },
    }
    numerics = {"long-window Gauss law holds": all(
        r["final_gauss"] < ZERO for r in lw["rows"])}
    if "numerics" in m:
        num = m["numerics"]
        numerics["rms sheet charge is stable under a halved step"] = all(
            num[k]["rms_change_halved_step"] < 1e-2 for k in num)
        numerics["rms sheet charge agrees with DOP853"] = all(
            num[k]["rms_change_dop853"] < 1e-2
            and num[k]["dop853_gauss"] < ZERO for k in num)
    out["long_window_numerics"] = numerics
    return out


def run_probe(progress=lambda s: None, sizes=None, long_steps=kr.LONG_STEPS,
              numerics=True):
    sizes = sizes or SIZES
    measured, per_size = {}, {}
    for name, config in sizes.items():
        progress(f"{name}: {config}")
        measured[name] = measure(config, long_steps, numerics)
        per_size[name] = evaluate(measured[name])
    first = next(iter(sizes))
    groups = {f: list(g) for f, g in per_size[first].items()}
    checks = {c: all(bool(per_size[s][f][c]) for s in sizes)
              for f, names in groups.items() for c in names}
    counter = {k: v for k, v in CONTRADICTED.items()
               if not all(checks[c] for c in groups[k])}
    v = kr.verdict(checks, groups, counter)
    return {
        "public_preregistration": PREREG, "scope": SCOPE, "seed": kr.SEED,
        "long_window": long_steps * kr.STEP,
        "measurements": measured,
        "checks_by_size": {s: {c: bool(val) for g in per_size[s].values()
                               for c, val in g.items()} for s in sizes},
        "checks": checks, "check_groups": groups,
        "contradicted_predictions": counter,
        "checks_passed": v["numerical_validity"] == "PASS",
        "verdict": v,
    }


def render(report):
    v = report["verdict"]
    sizes = list(report["measurements"])
    lines = ["# Two routes to the sheet pairing, on equal terms", "",
             report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. Seed "
             f"`{report['seed']}`. Long window `[0, "
             f"{report['long_window']:g}]`.", "",
             "| verdict field | value |", "|---|---|"]
    lines += [f"| {k} | **{v[k]}** |" for k in kr.VERDICT_FIELDS]
    for field, text in report["contradicted_predictions"].items():
        lines += ["", f"**Frozen prediction contradicted ({field}).** {text}"]
    if v["failed_checks"] or v["missing_checks"]:
        lines += ["", "**Outside tolerance:** "
                  + ", ".join(f"`{c}`" for c in v["failed_checks"]),
                  "", "**Not run:** "
                  + ", ".join(f"`{c}`" for c in v["missing_checks"])]
    lines += ["", "## Short window: sheet charge at matched energy", "",
              "| size | energy | matter current | at rest | flux route | "
              "ratio F/M (median, range) |", "|---|---:|---:|---:|---:|---|"]
    for size in sizes:
        sc = report["measurements"][size]["P2_P3"]
        for j, e in enumerate(sc["energies"]):
            r = sc["ratio"]
            lines.append(
                f"| {size} | {e:g} | {sc['M']['median_peak'][j]:.3e} | "
                f"{sc['M0']['median_peak'][j]:.3e} | "
                f"{sc['F']['median_peak'][j]:.3e} | {r['median'][j]:.4f} "
                f"({r['min'][j]:.4f} to {r['max'][j]:.4f}) |")
    lines += ["", "| size | slope M | slope M0 | slope F | slope of ratio |",
              "|---|---:|---:|---:|---:|"]
    for size in sizes:
        sc = report["measurements"][size]["P2_P3"]
        lines.append(f"| {size} | {sc['M']['slope_of_median']:.4f} | "
                     f"{sc['M0']['slope_of_median']:.4f} | "
                     f"{sc['F']['slope_of_median']:.4f} | "
                     f"{sc['ratio']['slope_of_median']:.4f} |")
    lines += ["", "## Long window (measurements, except the ranking)", "",
              "| size | energy | rms M | rms F | rms B | ratio long | ratio "
              "short | F larger | flux retained (min) | current half-life |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for size in sizes:
        m = report["measurements"][size]
        for r in m["long_window"]["rows"]:
            lines.append(
                f"| {size} | {r['energy']:g} | {r['median_rms']['M']:.3e} | "
                f"{r['median_rms']['F']:.3e} | {r['median_rms']['B']:.3e} | "
                f"{r['ratio_long_rms']:.3f} | {r['ratio_short_peak']:.4f} | "
                f"{r['draws_flux_route_larger_long']} of {m['draws']} | "
                f"{r['flux_retention_min']:.4f} | "
                f"{r['median_current_half_life']:.2f} |")
    lines += ["", "| Check | Pass |", "|---|---|"]
    lines += [f"| {k} | {val} |" for k, val in report["checks"].items()]
    lines += ["", f"{sum(report['checks'].values())} of "
              f"{len(report['checks'])} checks hold; "
              f"{len(v['contradicted_checks'])} frozen prediction check "
              "contradicted, reported above.", "",
              "Matched total energy is one definition of equal terms. "
              "Neither orienting datum is derived, and nothing here is a "
              "probability.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    report = run_probe(progress=lambda s: print(s, flush=True))
    summary = render(report)
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "probe.json").write_text(
            json.dumps(report, indent=2, allow_nan=False, default=float) + "\n")
        (args.output_dir / "probe.md").write_text(summary)
    print(summary)
    return 0 if report["checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
