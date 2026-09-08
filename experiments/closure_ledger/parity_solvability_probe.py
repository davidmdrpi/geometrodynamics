"""Does constraint solvability force antipodal parity? Freeze ``495f1f1``.

Reports separate verdict fields and exits nonzero when any required check
fails.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves import parity_solvability as ps


PREDICTED_FREE = ([1, 4], [2, 5], [3, 6])
PREDICTED_OBSTRUCTED = ([1, 2], [2, 3], [3, 4])
PARITY_PURE = ([1, 3], [2, 4], [1, 3, 5], [2, 4, 6])


def run_probe(progress=lambda s: None):
    progress("selection rule table")
    table = ps.selection_rule_table(ps.MAX_DEGREE)

    progress("degree-set scan")
    scan = ps.scan_degree_sets(list(PREDICTED_FREE) + list(PREDICTED_OBSTRUCTED)
                               + list(PARITY_PURE), samples=200)
    by_set = {tuple(r["degrees"]): r for r in scan}

    progress("independent improved-stress cross-check, full frozen coverage")
    # Correction C6. The first version sampled seven hand-picked degree sets
    # with field data only. Frozen check 2 demands every pair through degree 6
    # AND nonzero momentum. That is 21 unordered pairs times three data kinds.
    rng = np.random.default_rng(ps.SEED + 3)
    agreement = []
    for n in range(1, ps.MAX_DEGREE + 1):
        for m in range(n, ps.MAX_DEGREE + 1):
            degrees = [n] if n == m else [n, m]
            for kind in ("field", "momentum", "both"):
                base = ps.random_state(degrees, rng)
                field = base if kind in ("field", "both") else {}
                momentum = (ps.random_state(degrees, rng)
                            if kind in ("momentum", "both") else {})
                agreement.append({"degrees": degrees, "kind": kind,
                                  **ps.route_agreement(field, momentum)})

    progress("bilinearity through the independent route")
    halving = {str(d): ps.bilinearity_in_amplitude(d, samples=200)
               for d in PREDICTED_OBSTRUCTED}

    progress("quadrature exactness")
    exactness = ps.quadrature_exactness()

    progress("small-projection counterexample")
    small = ps.small_projection_counterexample()

    progress("complete momentum audit")
    audit = ps.complete_momentum_audit()

    progress("maximality")
    maximal = [ps.subspace_maximality(n, n + 1) for n in (1, 2, 3)]
    graphs = [ps.graph_subspace_search(n, n + 1) for n in (1, 2)]

    progress("momentum sector")
    momentum = ps.momentum_sector_report()

    checks = {
        "selection rule holds for every pair up to degree 6":
            table["selection_rule_holds"],
        "predicted free mixed-parity sets have no obstruction": all(
            by_set[tuple(d)]["max_abs_dipole"] < 1e-12 for d in PREDICTED_FREE),
        "predicted obstructed sets do obstruct": all(
            by_set[tuple(d)]["max_abs_dipole"] > 1e-3 for d in PREDICTED_OBSTRUCTED),
        "parity-pure sets have no obstruction": all(
            by_set[tuple(d)]["max_abs_dipole"] < 1e-12 for d in PARITY_PURE),
        "overlap reduction matches the inherited improved stress": all(
            (a["relative"] < 1e-10 if not a["below_floor"] else a["absolute"] < 1e-12)
            for a in agreement),
        "the obstruction is bilinear, measured independently": all(
            v["max_halving_error"] < 1e-8 for v in halving.values())
            and exactness["exact"],
        "no nonzero subspace evades an adjacent partner": all(
            m["kernel_is_trivial"] for m in maximal),
        "no nonzero graph subspace evades an adjacent pair": all(
            g["only_trivial_graph"] for g in graphs),
        "a small-projection mixed-parity subspace evades an adjacent pair":
            small["evades"],
        "the complete six Killing and four gradient charges are audited":
            audit["incomplete_audit_would_report_zero"]
            and audit["gradient_sector_is_independent"]
            and audit["gradient_obeys_adjacency"]
            and audit["killing_is_diagonal_in_degree"],
        "momentum charges are not a parity condition":
            momentum["parity_pure_kills_dipole"]
            and momentum["parity_pure_can_carry_charge"],
    }
    checks = {k: bool(v) for k, v in checks.items()}
    assert set(checks) == set(ps.REQUIRED_CHECKS), (
        "probe checks must match the frozen required set exactly")
    note = ("KILLING_DIAGONAL_IN_DEGREE_AND_GRADIENT_CKV_ADJACENT; "
            + momentum["note"])
    return {
        "public_preregistration": ps.PUBLIC_PREREG,
        "seed": ps.SEED,
        "selection_rule": table, "degree_set_scan": scan,
        "route_agreement": agreement, "halving": halving,
        "maximality": maximal, "graph_subspaces": graphs,
        "momentum_sector": momentum, "complete_momentum_audit": audit,
        "quadrature_exactness": exactness,
        "small_projection_counterexample": small,
        "checks": checks, "checks_passed": all(checks.values()),
        "verdict": ps.verdict(checks, note),
        "implementation_corrections": [
            "C1: the overlap route summed ordered pairs (n,m) and (m,n), each "
            "with the full frozen P2 coefficient, double counting the total "
            "cross term. Found by the frozen improved-stress cross-check, "
            "which disagreed by exactly 2.000. One overall factor of 1/2 fixes "
            "it. No frozen prediction changed: a global factor cannot move a "
            "zero, so the selection rule and both counterexample classes are "
            "unaffected; only reported magnitudes halve."],
        "scope": "Linearized Hamiltonian constraint at order phi^2 on the round "
                 "S^3, zero supporting-matter perturbation, linear subspaces of "
                 "pure-multiplet data. Not the nonlinear constraint, not "
                 "evolution, not a variety-level characterization.",
        "triangle_map": "NOT_DERIVED", "readout": "NOT_DERIVED",
    }


def render(report):
    v = report["verdict"]
    t = report["selection_rule"]
    lines = ["# Does constraint solvability force antipodal parity?", "",
             report["scope"], "",
             f"Public freeze: `{report['public_preregistration']}`. "
             f"Seed `{report['seed']}`.", "",
             "| verdict field | value |", "|---|---|"]
    for key in ps.VERDICT_FIELDS:
        lines.append(f"| {key} | **{v[key]}** |")
    if v.get("failed_checks"):
        lines += ["", "**Required checks failed; every field is `UNRESOLVED`.** "
                  "Failing: " + ", ".join(f"`{c}`" for c in v["failed_checks"]), ""]
    lines += ["", "## Selection rule", "",
              f"Minimum adjacent overlap norm `{t['min_adjacent_norm']:.6g}`; "
              f"maximum non-adjacent `{t['max_non_adjacent_norm']:.3e}`.", "",
              "| degrees | mixed parity | adjacent pair | max \\|P^A\\| |",
              "|---|---|---|---:|"]
    for r in report["degree_set_scan"]:
        lines.append(f"| `{r['degrees']}` | {r['mixed_parity']} | "
                     f"{r['has_adjacent_pair']} | {r['max_abs_dipole']:.3e} |")
    lines += ["", "## Momentum sector (no prediction was frozen)", "",
              "| degree | parity | Hamiltonian dipole | Killing charge |",
              "|---:|---|---:|---:|"]
    for r in report["momentum_sector"]["rows"]:
        lines.append(f"| {r['degree']} | {r['parity']} | "
                     f"{r['hamiltonian_dipole']:.3e} | "
                     f"{r['killing_charge_aligned']:.6f} |")
    lines += ["", "| Required check | Pass |", "|---|---|"]
    lines += [f"| {k} | {val} |" for k, val in report["checks"].items()]
    lines += ["", f"Passed {sum(report['checks'].values())}/"
              f"{len(report['checks'])} required checks.", "",
              "### Implementation corrections", ""]
    lines += [f"- {c}" for c in report["implementation_corrections"]]
    lines += ["", "Sufficiency of antipodal parity is not a derivation of it. "
              "No counting function, Born law, triangle map or readout follows.", ""]
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
