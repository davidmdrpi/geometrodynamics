"""Held-out test of the Klein splitting law, frozen in docs/klein_splitting_law_prereg.md.

Uses PR #355's solver unchanged. First reproduces the six in-sample splittings
of #355's archive, then measures the gated cells at N = 2000 and 4000 and
applies the frozen gates and verdict rules.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.special import gammaln, lpmv

from geometrodynamics.bulk import alice_handle_ports as ah

ROOT = Path(__file__).resolve().parents[2]
FREEZE = "56055ab"
SIZES = (2000, 4000)
ARCHIVE_355 = ROOT / "experiments/closure_ledger/runs/20261118_alice_handle_ports/probe.json"
_spec = importlib.util.spec_from_file_location(
    "klein_splitting_law_predictions", Path(__file__).with_name("klein_splitting_law_predictions.py"))
PRED = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(PRED)


def cached(cache, name, compute):
    if cache is None:
        return compute()
    path = Path(cache) / f"{name}.json"
    if path.exists():
        return json.loads(path.read_text())
    out = compute()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out))
    return out


def measure(a, m, n, pairs=3):
    """Lowest Klein pairs, paired as in #355's p4_p5."""
    lam, _, cos, purity = ah.klein_sectors(ah.Sector(ah.Loop(a, n), m, ah.KLEIN))
    ic, is_ = np.flatnonzero(cos), np.flatnonzero(~cos)
    rows = []
    for k in range(pairs):
        wp, wm = np.sqrt(lam[ic[k]]), np.sqrt(lam[is_[k]])
        rows.append(dict(n=k, omega_cos=float(wp), omega_sin=float(wm), delta=float(abs(wp - wm))))
    return dict(a=a, m=m, N=n, purity=float(purity), pairs=rows)


def excited_prediction(a, m, k):
    """Diagnostic only: the same Herring construction with the exact band state P_{m+k}^m."""
    ell = m + k
    omega = np.sqrt(ell * (ell + 1))
    if omega >= m / np.sin(a):
        return None
    kappa = np.sqrt(m * m - omega * omega * np.sin(a) ** 2)
    norm = 2 * np.exp(gammaln(ell + m + 1) - gammaln(ell - m + 1)) / (2 * ell + 1)
    w = PRED.neck_action(a, m, omega)
    return float(8 * kappa * lpmv(m, ell, np.cos(a)) ** 2 * np.exp(-w) / norm / (2 * omega))


def run(cache=None, progress=print):
    started = time.time()
    archive = json.loads(ARCHIVE_355.read_text())
    report = dict(freeze=FREEZE, solver="geometrodynamics/bulk/alice_handle_ports.py at 07d58de",
                  sizes=list(SIZES))

    # 1. Reproduction check on #355's in-sample cells.
    repro, worst = [], 0.
    for a, m in PRED.IN_SAMPLE:
        for n in SIZES:
            progress(f"in-sample a={a} m={m} N={n}")
            got = cached(cache, f"cell_{a}_{m}_{n}", lambda: measure(a, m, n))
            ref = next(r for r in archive["per_size"][str(n)]["measured"]["P4"]
                       if r["a"] == a and r["m"] == m)["pairs"][0]
            ref_delta = abs(ref["omega_cos"] - ref["omega_sin"])
            rel = abs(got["pairs"][0]["delta"] / ref_delta - 1)
            worst = max(worst, rel)
            repro.append(dict(a=a, m=m, N=n, delta=got["pairs"][0]["delta"],
                              archive_delta=ref_delta, relative=rel))
    report["reproduction"] = dict(cells=repro, worst_relative=worst, passed=worst <= 1e-6)
    if not report["reproduction"]["passed"]:
        report.update(verdict="INCOMPLETE", reasons=["reproduction check"])
        return report

    # 2. Gated and diagnostic cells.
    def cell(a, m):
        pred = PRED.prediction(a, m)
        row = dict(a=a, m=m, predicted_delta=pred["predicted_delta"],
                   reversal_time=float(np.pi / (2 * pred["predicted_delta"])), sizes={})
        for n in SIZES:
            progress(f"a={a} m={m} N={n}")
            got = cached(cache, f"cell_{a}_{m}_{n}", lambda: measure(a, m, n))
            r = got["pairs"][0]["delta"] / pred["predicted_delta"]
            excited = {}
            for p in got["pairs"][1:]:
                h = excited_prediction(a, m, p["n"])
                excited[str(p["n"])] = None if h is None else p["delta"] / h
            row["sizes"][str(n)] = dict(delta=got["pairs"][0]["delta"], ratio=float(r),
                                        m_r_minus_1=float(m * (r - 1)), purity=got["purity"],
                                        excited_ratio=excited)
        d0, d1 = (row["sizes"][str(n)]["delta"] for n in SIZES)
        row["resolution"] = float(abs(d0 / d1 - 1))
        return row

    gated = [cell(a, m) for a, m in PRED.GATED]
    diagnostic = [cell(a, m) for a, m in PRED.DIAGNOSTIC]
    report.update(gated=gated, diagnostic=diagnostic)

    # 3. Frozen gates, per resolution.
    gates = {}
    for n in map(str, SIZES):
        c = {a: np.mean([g["sizes"][n]["m_r_minus_1"] for g in gated if g["a"] == a])
             for a in (0.2, 0.4)}
        g1 = all(0 < g["sizes"][n]["m_r_minus_1"] <= 0.25 for g in gated)
        g2 = (c[0.2] < c[0.4]) and all(0.048 <= g["sizes"][n]["m_r_minus_1"] <= 0.182
                                       for g in gated if g["a"] == 0.4)
        gates[n] = dict(G1=bool(g1), G2=bool(g2), c_02=float(c[0.2]), c_04=float(c[0.4]))
    g3 = all(g["resolution"] <= 1e-3 for g in gated)
    report["gates"] = dict(per_size=gates, G3=bool(g3),
                           worst_resolution=float(max(g["resolution"] for g in gated)))

    both_fail = [k for k in ("G1", "G2") if not any(gates[n][k] for n in gates)]
    if both_fail:
        verdict, reasons = "COUNTEREXAMPLE", both_fail
    elif all(gates[n]["G1"] and gates[n]["G2"] for n in gates) and g3:
        verdict, reasons = "CONFIRMED", []
    elif not g3:
        verdict, reasons = "UNRESOLVED_NUMERICALLY", ["G3"]
    else:
        # A G1/G2 miss at one size only. It shrank (to nothing) if only the coarse
        # size missed; that is UNRESOLVED_NUMERICALLY. A miss that appears only at
        # the finer size grew, which the freeze does not name.
        missed = [k for k in ("G1", "G2") if not all(gates[n][k] for n in gates)]
        shrank = all(gates[str(SIZES[-1])][k] for k in missed)
        verdict, reasons = ("UNRESOLVED_NUMERICALLY" if shrank else "NO_VERDICT"), missed
    report.update(verdict=verdict, reasons=reasons, seconds=time.time() - started)
    return report


def markdown(rep):
    lines = [f"# Klein splitting law: held-out test against freeze `{rep['freeze']}`", "",
             f"**Verdict: {rep['verdict']}**" + (f" ({', '.join(rep['reasons'])})" if rep["reasons"] else ""),
             "", f"Reproduction of #355's six in-sample splittings: worst relative difference "
             f"{rep['reproduction']['worst_relative']:.2e} (gate 1e-6).", ""]
    if "gated" not in rep:
        return "\n".join(lines)
    lines += ["| a | m | predicted Δ | measured Δ (N=4000) | r = Δ/H | m(r−1) N=2000 | m(r−1) N=4000 | \\|Δ₂₀₀₀/Δ₄₀₀₀ − 1\\| |",
              "|---|---|---|---|---|---|---|---|"]
    for g in rep["gated"] + rep["diagnostic"]:
        s2, s4 = g["sizes"]["2000"], g["sizes"]["4000"]
        tag = "" if g in rep["gated"] else " (diagnostic)"
        lines.append(f"| {g['a']}{tag} | {g['m']} | {g['predicted_delta']:.4e} | {s4['delta']:.4e} | "
                     f"{s4['ratio']:.4f} | {s2['m_r_minus_1']:.3f} | {s4['m_r_minus_1']:.3f} | "
                     f"{g['resolution']:.1e} |")
    lines += ["", "| gate | N=2000 | N=4000 |", "|---|---|---|"]
    for k in ("G1", "G2"):
        lines.append(f"| {k} | {'pass' if rep['gates']['per_size']['2000'][k] else 'FAIL'} | "
                     f"{'pass' if rep['gates']['per_size']['4000'][k] else 'FAIL'} |")
    lines += [f"| G3 (resolution) | {'pass' if rep['gates']['G3'] else 'FAIL'} "
              f"(worst {rep['gates']['worst_resolution']:.1e}) | |", "",
              "c(a) = mean m(r−1) over gated m: " + ", ".join(
                  f"N={n}: c(0.2) = {v['c_02']:.3f}, c(0.4) = {v['c_04']:.3f}"
                  for n, v in rep["gates"]["per_size"].items()), "",
              f"Run time {rep['seconds']:.0f} s.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "experiments/closure_ledger/runs/20261123_klein_splitting_law")
    parser.add_argument("--cache", type=Path, default=None)
    args = parser.parse_args(argv)
    rep = run(args.cache, lambda s: print(s, flush=True))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "probe.json").write_text(json.dumps(rep, indent=1) + "\n")
    (args.out / "probe.md").write_text(markdown(rep))
    print(markdown(rep))


if __name__ == "__main__":
    main()
