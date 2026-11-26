"""Held-out run of docs/alice_vcurve_prereg.md (freeze 291ab5f) on (0.5, 2) and (0.35, 2)."""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FREEZE = "291ab5f"
OUT = ROOT / "experiments/closure_ledger/runs/20261125_alice_vcurve"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V = _load("alice_vcurve_probe")
Q = _load("alice_vcurve_predictions")
SPECTRA = {}


def spectrum(a, n):
    if (a, n) not in SPECTRA:
        SPECTRA.clear()
        SPECTRA[(a, n)] = V.BASE.solve(a, Q.M, n)
    return SPECTRA[(a, n)]


def measure(cache, a, tau, kind, eta, n=2000, keep=150, progress=print):
    key = f"{a}_{n}_{keep}_{kind}_{eta}_{tau:.6f}"
    path = Path(cache) / f"{key}.json" if cache else None
    if path is not None and path.exists():
        return json.loads(path.read_text())
    progress(key)
    try:
        r = V.point(spectrum(a, n), a, Q.M, n, keep, tau, kind, eta)
        r.pop("U"), r.pop("V")
    except RuntimeError as err:
        # The requested depletion is not reachable on the weak-coupling branch
        # at this window: record the largest depletion that is.
        geo = V.Sweep(a, Q.M, n, keep, spectrum(a, n), tau, kind)
        strengths = 0.05 * 8.0 / tau * 1.5 ** np.arange(40)
        reached = [geo.absorbed_plus(g, geo.mean) for g in strengths]
        r = dict(a=a, m=Q.M, N=n, K=keep, tau=tau, profile=kind, eta=eta, unreachable=True,
                 max_eta=float(max(reached)), at_strength=float(strengths[int(np.argmax(reached))]), error=str(err))
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(r))
    return r


def run(cache=None, progress=print):
    started = time.time()
    frozen = json.loads((OUT / "frozen_predictions.json").read_text())
    pred = {(r["a"], r["eta"], round(r["tau"], 6)): r for r in frozen["vcurve"]}
    rep = dict(freeze=FREEZE, cells=[Q.A1, Q.A2], m=Q.M, grid=list(Q.GRID))

    # V-curves, sorted by cell so each spectrum is solved once
    vcurve = []
    for a in (Q.A1, Q.A2):
        for eta in Q.ETAS:
            for tau in Q.GRID:
                r = measure(cache, a, tau, "sin2", eta, progress=progress)
                p = pred[(a, eta, round(tau, 6))]
                r.update(predicted=p["total"], predicted_switching=p["switching"], predicted_exchange=p["exchange"])
                if not r.get("unreachable"):
                    r["ratio"] = r["max_abs_defect"] / p["total"]
                vcurve.append(r)
        if a == Q.A1:
            scan = [measure(cache, a, float(t), "sin2", 0.5, progress=progress) for t in Q.NULL_SCAN]
        else:
            envelope = [measure(cache, a, tau, kind, 0.5, progress=progress)
                        for tau in Q.ENVELOPE_TAUS for kind in Q.PROFILES]
        geo = V.Sweep(a, Q.M, 2000, 150, spectrum(a, 2000), Q.GRID[0], "sin2")
        rep.setdefault("undetected", {})[str(a)] = geo.run(0.0, geo.mean)["max_abs_defect"]
    rep.update(vcurve=vcurve, null_scan=scan, envelope=envelope)

    # Numerics at the measured minima: K = 80 and N = 4000
    minima = {}
    for a in (Q.A1, Q.A2):
        for eta in Q.ETAS:
            rows = [r for r in vcurve if r["a"] == a and r["eta"] == eta and not r.get("unreachable")]
            minima[(a, eta)] = min(rows, key=lambda r: r["max_abs_defect"])
    checks = []
    for a in (Q.A1, Q.A2):
        for eta in Q.ETAS:
            base = minima[(a, eta)]
            k80 = measure(cache, a, base["tau"], "sin2", eta, keep=80, progress=progress)
            checks.append(dict(a=a, eta=eta, tau=base["tau"], K150=base["max_abs_defect"], K80=k80["max_abs_defect"]))
    for a in (Q.A1, Q.A2):
        for eta in Q.ETAS:
            base = minima[(a, eta)]
            n4 = measure(cache, a, base["tau"], "sin2", eta, n=4000, progress=progress)
            next(c for c in checks if c["a"] == a and c["eta"] == eta)["N4000"] = n4["max_abs_defect"]
    rep["numerics"] = checks
    rep.update(gates(rep, minima, frozen))
    rep["seconds"] = time.time() - started
    return rep


def gates(rep, minima, frozen):
    g = {}
    unreachable = [(r["a"], r["eta"], r["tau"], r["max_eta"]) for r in rep["vcurve"] + rep["null_scan"] + rep["envelope"]
                   if r.get("unreachable")]
    vc = [r for r in rep["vcurve"] if not r.get("unreachable")]
    g["G1"] = dict(passed=all(1 / 3 <= r["ratio"] <= 3 for r in vc),
                   worst=[min(r["ratio"] for r in vc), max(r["ratio"] for r in vc)],
                   misses=[(r["a"], r["eta"], r["tau"], r["ratio"]) for r in vc if not 1 / 3 <= r["ratio"] <= 3])
    g2 = []
    for (a, eta), best in minima.items():
        taus = list(Q.GRID)
        i = taus.index(best["tau"])
        prows = [r for r in frozen["vcurve"] if r["a"] == a and r["eta"] == eta]
        j = taus.index(min(prows, key=lambda r: r["total"])["tau"])
        g2.append(dict(a=a, eta=eta, measured_index=i, predicted_index=j, interior=0 < i < len(taus) - 1, near=abs(i - j) <= 1))
    g["G2"] = dict(passed=all(x["interior"] and x["near"] for x in g2), cells=g2)
    g3 = []
    for eta in Q.ETAS:
        e1, e2 = minima[(Q.A1, eta)]["max_abs_defect"], minima[(Q.A2, eta)]["max_abs_defect"]
        g3.append(dict(eta=eta, E_min_a1=e1, E_min_a2=e2, ratio=e1 / e2, below_1e4=e2 <= 1e-4, ok=(e2 <= e1 / 3) and e2 <= 1e-4))
    g["G3"] = dict(passed=all(x["ok"] for x in g3), cells=g3)
    taus = np.array([r["tau"] for r in rep["null_scan"]])
    vals = np.array([r["max_abs_defect"] for r in rep["null_scan"]])
    g4 = []
    for tk in frozen["null_taus"]:
        win = np.abs(taus - tk) <= 0.64
        i = np.argmin(np.where(win, vals, np.inf))
        g4.append(dict(tau_k=tk, tau_min=float(taus[i]), offset=float(taus[i] - tk), depth=float(vals[i] / vals[win].max()),
                       ok=bool(abs(taus[i] - tk) <= 0.15 and vals[i] <= 0.4 * vals[win].max())))
    g["G4"] = dict(passed=all(x["ok"] for x in g4), nulls=g4)
    env = {(r["tau"], r["profile"]): r["max_abs_defect"] for r in rep["envelope"]}
    penv = {(r["tau"], r["profile"]): r["total"] for r in frozen["envelope"]}
    g5 = []
    for tau in Q.ENVELOPE_TAUS:
        for kind in ("sin4", "bump"):
            m_ratio, p_ratio = env[(tau, kind)] / env[(tau, "sin2")], penv[(tau, kind)] / penv[(tau, "sin2")]
            g5.append(dict(tau=tau, profile=kind, measured_ratio=m_ratio, predicted_ratio=p_ratio,
                           ok=bool(1 / 2.5 <= m_ratio / p_ratio <= 2.5)))
    ranking = env[(16.0, "sin4")] <= 0.2 * env[(16.0, "sin2")] and env[(16.0, "bump")] > env[(16.0, "sin2")]
    g["G5"] = dict(passed=all(x["ok"] for x in g5) and bool(ranking), ratios=g5, ranking_at_16=bool(ranking))
    num = rep["numerics"]
    k_ok = all(abs(c["K80"] - c["K150"]) <= max(1e-6, 0.1 * c["K150"]) for c in num)
    n_ok = all(abs(c["N4000"] - c["K150"]) <= max(1e-6, 0.2 * c["N4000"]) for c in num)
    balance = max(r["balance"] for r in vc + rep["null_scan"] + rep["envelope"] if not r.get("unreachable"))
    und = max(rep["undetected"].values())
    g["G6"] = dict(passed=bool(k_ok and n_ok and balance <= 1e-8 and und <= 1e-8),
                   modal=bool(k_ok), size=bool(n_ok), balance=balance, undetected=und)
    g["unreachable_points"] = unreachable
    failed = [k for k in ("G1", "G2", "G3", "G4", "G5") if not g[k]["passed"]]
    if failed:
        verdict, reasons = "COUNTEREXAMPLE", failed
    elif unreachable:
        # A frozen point could not be measured at its depletion: a case the
        # prereg does not name.
        verdict, reasons = "NO_VERDICT", ["unreachable depletion at " + ", ".join(
            f"a={u[0]} eta={u[1]} tau={u[2]:.2f} (max {u[3]:.3f})" for u in unreachable)]
    elif g["G6"]["passed"]:
        verdict, reasons = "CONFIRMED", []
    else:
        verdict, reasons = "UNRESOLVED_NUMERICALLY", ["G6"]
    return dict(gates=g, verdict=verdict, reasons=reasons)


def markdown(rep):
    L = [f"# Alice-handle exchange: V-curve, aperture scaling and envelope control (freeze `{rep['freeze']}`)", "",
         f"**Verdict: {rep['verdict']}**" + (f" ({', '.join(rep['reasons'])})" if rep["reasons"] else ""), "",
         "Largest |u + v − 1|, sin², N = 2000, K = 150. Each cell is measured (predicted).", "",
         "| cell | η | " + " | ".join(f"{t:.2f}" for t in rep["grid"]) + " |",
         "|---|---|" + "---|" * len(rep["grid"])]
    for a in rep["cells"]:
        for eta in (0.5, 0.9):
            rows = sorted([r for r in rep["vcurve"] if r["a"] == a and r["eta"] == eta], key=lambda r: r["tau"])
            L.append(f"| {a} | {eta} | " + " | ".join(
                (f"η unreachable (max {r['max_eta']:.3f})" if r.get("unreachable") else f"{r['max_abs_defect']:.2e} ({r['predicted']:.1e})")
                for r in rows) + " |")
    g = rep["gates"]
    L += ["", "| gate | result |", "|---|---|",
          f"| G1 predictor within ×3 | {'pass' if g['G1']['passed'] else 'FAIL'} (measured/predicted {g['G1']['worst'][0]:.2f}–{g['G1']['worst'][1]:.2f}) |",
          f"| G2 V shape | {'pass' if g['G2']['passed'] else 'FAIL'} |",
          "| G3 aperture scaling | " + ("pass" if g["G3"]["passed"] else "FAIL") + " (" + "; ".join(
              f"η={c['eta']}: {c['E_min_a1']:.2e} → {c['E_min_a2']:.2e}, ×{c['ratio']:.1f}" for c in g["G3"]["cells"]) + ") |",
          "| G4 spectral nulls | " + ("pass" if g["G4"]["passed"] else "FAIL") + " (" + "; ".join(
              f"{n['tau_k']:.3f}: min at {n['tau_min']:.1f}, depth {n['depth']:.2f}" for n in g["G4"]["nulls"]) + ") |",
          "| G5 envelope | " + ("pass" if g["G5"]["passed"] else "FAIL") + " (" + "; ".join(
              f"τ={r['tau']:.0f} {r['profile']}/sin² {r['measured_ratio']:.3f} vs {r['predicted_ratio']:.3f}" for r in g["G5"]["ratios"]) + ") |",
          f"| G6 numerics | {'pass' if g['G6']['passed'] else 'FAIL'} (K {g['G6']['modal']}, N {g['G6']['size']}, balance {g['G6']['balance']:.1e}, undetected {g['G6']['undetected']:.1e}) |",
          "", f"Run time {rep['seconds']:.0f} s.", ""]
    return "\n".join(L)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, default=None)
    args = parser.parse_args(argv)
    rep = run(args.cache, lambda s: print(s, flush=True))
    (OUT / "probe.json").write_text(json.dumps(rep, indent=1) + "\n")
    (OUT / "probe.md").write_text(markdown(rep))
    print(markdown(rep))


if __name__ == "__main__":
    main()
