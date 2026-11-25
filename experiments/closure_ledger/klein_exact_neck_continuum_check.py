"""Independent re-measurement of the exact-neck held-out cells with PR #360's continuum solver.

Not part of the freeze and not a gate. The frozen verdict is the finite-element
one in runs/20261124_klein_exact_neck/probe.json. This check re-measures every
gated and diagnostic cell with #360's Wronskian/shooting solver
(alice_ground_split.measure), which shares no code with #355's finite-element
solver. It then reports:
* the agreement of the two solvers;
* the frozen gates scored on the continuum splittings;
* the rim-kink ratio on the continuum splittings.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "experiments/closure_ledger/runs/20261124_klein_exact_neck"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRED = _load("klein_exact_neck_predictions")
CONT = _load("alice_ground_split")


def continuum_delta(a, m):
    """|omega_+ - omega_-| of the lowest Klein pair from #360's solver."""
    c = CONT.measure(a, m)
    mean = (np.sqrt(c["lambda_N"]) + np.sqrt(c["lambda_D"])) / 2
    return float(c["split"] * mean), float(mean)


def check(archive):
    rows = []
    for key in ("gated", "diagnostic"):
        for g in archive[key]:
            a, m = g["a"], g["m"]
            pred = PRED.prediction(a, m)
            delta, mean = continuum_delta(a, m)
            r, r_h = delta / pred["predicted_delta"] - 1, delta / pred["herring_delta"] - 1
            kink = abs(pred["neck_log_derivative"] / pred["band_log_derivative"] - 1)
            rows.append(dict(a=a, m=m, gated=key == "gated", continuum_delta=delta, mean_omega=mean,
                             fem_over_continuum_minus_1=g["sizes"]["4000"]["delta"] / delta - 1,
                             r_minus_1=r, r_h_minus_1=r_h, neck_fraction=abs(r) / abs(r_h),
                             bound=pred["bound"], kink=kink, r_minus_1_over_kink=r / kink))
    gated = [x for x in rows if x["gated"]]
    summary = dict(
        worst_solver_difference=max(abs(x["fem_over_continuum_minus_1"]) for x in gated),
        worst_solver_difference_all=max(abs(x["fem_over_continuum_minus_1"]) for x in rows),
        fem_below_continuum_everywhere=all(x["fem_over_continuum_minus_1"] < 0 for x in rows),
        G1_on_continuum=all(abs(x["r_minus_1"]) <= x["bound"] for x in gated),
        G2_on_continuum=all(x["neck_fraction"] <= 0.2 for x in gated),
        worst_bound_use=max(abs(x["r_minus_1"]) / x["bound"] for x in gated),
        worst_neck_fraction=max(x["neck_fraction"] for x in gated),
        kink_ratio_gated_m_ge_2=[min(x["r_minus_1_over_kink"] for x in gated if x["m"] >= 2),
                                 max(x["r_minus_1_over_kink"] for x in gated if x["m"] >= 2)])
    return dict(solver="experiments/closure_ledger/alice_ground_split.py (PR #360), measure()",
                freeze=archive["freeze"], cells=rows, summary=summary)


def markdown(rep):
    s = rep["summary"]
    lines = ["# Exact-neck held-out cells: continuum cross-check (not a gate)", "",
             f"Solver: {rep['solver']}. Finite-element values from `probe.json`, N = 4000.", "",
             "| a | m | continuum Δ | FEM/continuum − 1 | r − 1 | \\|r − 1\\|/\\|r_H − 1\\| | (r − 1)/kink |",
             "|---|---|---|---|---|---|---|"]
    for x in rep["cells"]:
        tag = "" if x["gated"] else " (diagnostic)"
        lines.append(f"| {x['a']}{tag} | {x['m']} | {x['continuum_delta']:.5e} | "
                     f"{x['fem_over_continuum_minus_1']:+.1e} | {x['r_minus_1']:+.5f} | "
                     f"{x['neck_fraction']:.3f} | {x['r_minus_1_over_kink']:.3f} |")
    lines += ["",
              f"Worst solver difference on gated cells: {s['worst_solver_difference']:.1e} "
              f"(all cells {s['worst_solver_difference_all']:.1e}); finite element below the continuum "
              f"everywhere: {s['fem_below_continuum_everywhere']}.",
              f"Frozen gates scored on the continuum: G1 {'pass' if s['G1_on_continuum'] else 'FAIL'} "
              f"(worst use of bound {s['worst_bound_use']:.3f}), G2 {'pass' if s['G2_on_continuum'] else 'FAIL'} "
              f"(worst {s['worst_neck_fraction']:.3f}).",
              f"(r − 1)/kink, gated m ≥ 2: {s['kink_ratio_gated_m_ge_2'][0]:.3f}–{s['kink_ratio_gated_m_ge_2'][1]:.3f}.",
              ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=RUN)
    args = parser.parse_args(argv)
    rep = check(json.loads((RUN / "probe.json").read_text()))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "continuum_check.json").write_text(json.dumps(rep, indent=1) + "\n")
    (args.out / "continuum_check.md").write_text(markdown(rep))
    print(markdown(rep))


if __name__ == "__main__":
    main()
