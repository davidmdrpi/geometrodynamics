"""Temporal symmetry of exchange across the seam. Freeze ``95d7a80``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.bulk import klein_handle_parity as kp
from geometrodynamics.bulk import klein_handle_response as kh

PREREG = "95d7a80"
SCOPE = ("Round 15's ungauged spinor doublet on four seam types, responding "
         "mouths, generic deck-fixed data. Which parity in time the sheet "
         "value of each conserved quantity has. Classical; the reading for "
         "BAM rests on round 15's unreviewed matching.")
SIZES = {"primary": kp.PRIMARY, "replicate": kp.REPLICATE}


def run_probe(progress=lambda s: None, sizes=None, draws=kp.DRAWS):
    sizes = sizes or SIZES
    measured = {}
    for name, config in sizes.items():
        progress(f"{name}: {config}")
        measured[name] = kp.parity_report(
            kh.Network(**config), np.random.default_rng(kp.SEED), draws)
    rule, live = {}, {}
    for rows in measured.values():
        for r in rows:
            key = (f"{r['internal']} seam, time "
                   f"{'reversed' if r['s_t'] == -1 else 'kept'}, "
                   f"{r['generator']}: {r['predicted']}")
            ok = kp.satisfies(r)
            sound = (r["symmetry_defect"] < 1e-9 and r["max_drift"] < 1e-9)
            rule[f"S1 {key}"] = rule.get(f"S1 {key}", True) and ok["rule"] \
                and sound
            live[f"S2 {key}"] = live.get(f"S2 {key}", True) and ok["live"]
    first = next(iter(measured.values()))
    pick = {(r["internal"], r["s_t"], r["generator"]): r["predicted"]
            for r in first}
    opposite = {"S3 charge is odd and all three spin components even on the "
                "transport seam with time reversal": (
                    pick[("J", -1, "i")] == "odd" and all(
                        pick[("J", -1, g)] == "even"
                        for g in kp.GENERATORS[1:])
                    and all(v for k, v in rule.items() if "J seam, time rev"
                            in k) and all(v for k, v in live.items()
                                          if "J seam, time rev" in k))}
    checks = {**rule, **live, **opposite}
    groups = {"parity_rule": list(rule), "non_vacuous": list(live),
              "charge_against_spin": list(opposite)}
    v = kp.verdict(checks, groups)
    return {"public_preregistration": PREREG, "scope": SCOPE, "seed": kp.SEED,
            "measurements": measured, "checks": checks,
            "check_groups": groups,
            "checks_passed": v["numerical_validity"] == "PASS", "verdict": v}


def render(report):
    v = report["verdict"]
    lines = ["# Temporal symmetry of exchange across the seam", "",
             report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. Seed "
             f"`{report['seed']}`.", "", "| verdict field | value |",
             "|---|---|"]
    lines += [f"| {k} | **{v[k]}** |" for k in kp.VERDICT_FIELDS]
    for size, rows in report["measurements"].items():
        lines += ["", f"## {size}", "",
                  "| seam | time | quantity | predicted | even part | odd "
                  "part | total |", "|---|---|---|---|---|---|---:|"]
        for r in rows:
            if r["generator"] not in ("i", "sigma_z"):
                continue
            lines.append(
                f"| {r['internal']} | "
                f"{'reversed' if r['s_t'] == -1 else 'kept'} | {r['kind']} | "
                f"{r['predicted']} | {r['min_even_part']:.1e} to "
                f"{r['max_even_part']:.1e} | {r['min_odd_part']:.1e} to "
                f"{r['max_odd_part']:.1e} | {r['max_abs_total']:.1e} |")
    lines += ["", f"Passed {sum(report['checks'].values())}/"
              f"{len(report['checks'])} required checks.", "",
              "The field type and the time extension select the parity; the "
              "geometry alone does not.", ""]
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
