"""Held-out test of the exact-neck Klein splitting law, frozen in docs/klein_exact_neck_prereg.md.

Uses PR #355's solver unchanged, through the measurement of
klein_splitting_law_probe.py. First reproduces the 17 in-sample splittings of
#355's and #362's archives, then measures the gated and diagnostic cells at
N = 2000 and 4000 and applies the frozen gates and verdict rules.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FREEZE = "834ac70"
SIZES = (2000, 4000)
ARCHIVE_355 = ROOT / "experiments/closure_ledger/runs/20261118_alice_handle_ports/probe.json"
ARCHIVE_362 = ROOT / "experiments/closure_ledger/runs/20261123_klein_splitting_law/probe.json"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRED = _load("klein_exact_neck_predictions")
SPLIT = _load("klein_splitting_law_probe")


def archived_deltas():
    """The 17 in-sample splittings, by cell and size, from their own archives."""
    out = {}
    a355 = json.loads(ARCHIVE_355.read_text())
    for n in map(str, SIZES):
        for row in a355["per_size"][n]["measured"]["P4"]:
            pair = row["pairs"][0]
            out.setdefault((row["a"], row["m"]), {})[n] = abs(pair["omega_cos"] - pair["omega_sin"])
    a362 = json.loads(ARCHIVE_362.read_text())
    for g in a362["gated"] + a362["diagnostic"]:
        out[(g["a"], g["m"])] = {n: g["sizes"][n]["delta"] for n in map(str, SIZES)}
    return out


def cell_report(a, m, measured):
    """measured: {size: SPLIT.measure output}. Ratios for n = 0 and the excited pairs."""
    pred = PRED.prediction(a, m)
    row = dict(a=a, m=m, predicted_delta=pred["predicted_delta"], herring_delta=pred["herring_delta"],
               bound=pred["bound"], reversal_time=float(np.pi / (2 * pred["predicted_delta"])),
               neck_log_derivative=pred["neck_log_derivative"],
               band_log_derivative=pred["band_log_derivative"], sizes={})
    for n, got in measured.items():
        delta = got["pairs"][0]["delta"]
        r, r_h = delta / pred["predicted_delta"], delta / pred["herring_delta"]
        excited = {str(p["n"]): p["delta"] / PRED.prediction(a, m, p["n"])["predicted_delta"]
                   for p in got["pairs"][1:]}
        row["sizes"][str(n)] = dict(delta=delta, r_minus_1=float(r - 1), r_h_minus_1=float(r_h - 1),
                                    neck_fraction=float(abs(r - 1) / abs(r_h - 1)),
                                    purity=got["purity"], excited_ratio=excited)
    d0, d1 = (row["sizes"][str(n)]["delta"] for n in SIZES)
    row["resolution"] = float(abs(d0 / d1 - 1))
    row["resolution_signed"] = float(d0 / d1 - 1)
    return row


def verdict(gated):
    """The frozen gates and verdict rules, applied to the gated cell reports."""
    per_size = {}
    for n in map(str, SIZES):
        g1 = [(g["a"], g["m"]) for g in gated if not abs(g["sizes"][n]["r_minus_1"]) <= g["bound"]]
        g2 = [(g["a"], g["m"]) for g in gated if not g["sizes"][n]["neck_fraction"] <= 0.2]
        per_size[n] = dict(G1=not g1, G2=not g2, G1_misses=g1, G2_misses=g2)
    g3_misses = [(g["a"], g["m"]) for g in gated if not g["resolution"] <= 1e-3]
    gates = dict(per_size=per_size, G3=not g3_misses, G3_misses=g3_misses,
                 worst_resolution=float(max(g["resolution"] for g in gated)))
    both_fail = [k for k in ("G1", "G2") if not any(per_size[n][k] for n in per_size)]
    if both_fail:
        return gates, "COUNTEREXAMPLE", both_fail
    if all(per_size[n]["G1"] and per_size[n]["G2"] for n in per_size) and not g3_misses:
        return gates, "CONFIRMED", []
    if g3_misses:
        return gates, "UNRESOLVED_NUMERICALLY", ["G3"]
    # A G1/G2 miss at one size only: it shrank if only the coarse size missed.
    missed = [k for k in ("G1", "G2") if not all(per_size[n][k] for n in per_size)]
    shrank = all(per_size[str(SIZES[-1])][k] for k in missed)
    return gates, ("UNRESOLVED_NUMERICALLY" if shrank else "NO_VERDICT"), missed


def run(cache=None, progress=print):
    started = time.time()
    report = dict(freeze=FREEZE, solver="geometrodynamics/bulk/alice_handle_ports.py at 07d58de",
                  sizes=list(SIZES))

    def measure(a, m, n):
        progress(f"a={a} m={m} N={n}")
        return SPLIT.cached(cache, f"cell_{a}_{m}_{n}", lambda: SPLIT.measure(a, m, n))

    # 1. Reproduction check on the 17 in-sample cells.
    archive, repro, worst = archived_deltas(), [], 0.
    for a, m in PRED.IN_SAMPLE:
        for n in SIZES:
            got = measure(a, m, n)["pairs"][0]["delta"]
            ref = archive[(a, m)][str(n)]
            rel = abs(got / ref - 1)
            worst = max(worst, rel)
            repro.append(dict(a=a, m=m, N=n, delta=got, archive_delta=ref, relative=rel))
    report["reproduction"] = dict(cells=repro, worst_relative=worst, passed=worst <= 1e-6)
    if not report["reproduction"]["passed"]:
        report.update(verdict="INCOMPLETE", reasons=["reproduction check"], seconds=time.time() - started)
        return report

    # 2. Gated and diagnostic cells.
    gated = [cell_report(a, m, {n: measure(a, m, n) for n in SIZES}) for a, m in PRED.GATED]
    diagnostic = [cell_report(a, m, {n: measure(a, m, n) for n in SIZES}) for a, m in PRED.DIAGNOSTIC]
    report.update(gated=gated, diagnostic=diagnostic)

    # 3. Frozen gates.
    gates, v, reasons = verdict(gated)
    gates["max_abs_r_minus_1_over_a4"] = float(max(abs(g["sizes"]["4000"]["r_minus_1"]) / g["a"] ** 4
                                                   for g in gated))
    report.update(gates=gates, verdict=v, reasons=reasons, seconds=time.time() - started)
    return report


def markdown(rep):
    lines = [f"# Exact-neck Klein splitting law: held-out test against freeze `{rep['freeze']}`", "",
             f"**Verdict: {rep['verdict']}**" + (f" ({', '.join(rep['reasons'])})" if rep["reasons"] else ""),
             "", f"Reproduction of the 17 in-sample splittings: worst relative difference "
             f"{rep['reproduction']['worst_relative']:.2e} (gate 1e-6).", ""]
    if "gated" not in rep:
        return "\n".join(lines)
    lines += ["| a | m | predicted Δ = E | measured Δ (N=4000) | r − 1, N=2000 | r − 1, N=4000 | bound | "
              "\\|r − 1\\|/\\|r_H − 1\\|, worse size | \\|Δ₂₀₀₀/Δ₄₀₀₀ − 1\\| |",
              "|---|---|---|---|---|---|---|---|---|"]
    for g in rep["gated"] + rep["diagnostic"]:
        s2, s4 = g["sizes"]["2000"], g["sizes"]["4000"]
        tag = "" if g in rep["gated"] else " (diagnostic)"
        lines.append(f"| {g['a']}{tag} | {g['m']} | {g['predicted_delta']:.4e} | {s4['delta']:.4e} | "
                     f"{s2['r_minus_1']:+.5f} | {s4['r_minus_1']:+.5f} | {g['bound']:.2e} | "
                     f"{max(s2['neck_fraction'], s4['neck_fraction']):.3f} | {g['resolution']:.1e} |")
    per = rep["gates"]["per_size"]
    lines += ["", "| gate | N=2000 | N=4000 |", "|---|---|---|"]
    for k in ("G1", "G2"):
        lines.append(f"| {k} | {'pass' if per['2000'][k] else 'FAIL ' + str(per['2000'][k + '_misses'])} | "
                     f"{'pass' if per['4000'][k] else 'FAIL ' + str(per['4000'][k + '_misses'])} |")
    lines += [f"| G3 (resolution) | {'pass' if rep['gates']['G3'] else 'FAIL'} "
              f"(worst {rep['gates']['worst_resolution']:.1e}) | |", "",
              f"max \\|r − 1\\|/a⁴ over gated cells at N=4000: {rep['gates']['max_abs_r_minus_1_over_a4']:.3f}.",
              "", "Excited pairs, measured/E at N=4000 (diagnostic):", "",
              "| a | m | n = 1 | n = 2 |", "|---|---|---|---|"]
    for g in rep["gated"] + rep["diagnostic"]:
        ex = g["sizes"]["4000"]["excited_ratio"]
        lines.append(f"| {g['a']} | {g['m']} | {ex['1']:.4f} | {ex['2']:.4f} |")
    lines += ["", f"Run time {rep['seconds']:.0f} s.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "experiments/closure_ledger/runs/20261124_klein_exact_neck")
    parser.add_argument("--cache", type=Path, default=None)
    args = parser.parse_args(argv)
    rep = run(args.cache, lambda s: print(s, flush=True))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "probe.json").write_text(json.dumps(rep, indent=1) + "\n")
    (args.out / "probe.md").write_text(markdown(rep))
    print(markdown(rep))


if __name__ == "__main__":
    main()
