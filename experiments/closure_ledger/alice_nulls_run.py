"""Held-out run of docs/alice_nulls_prereg.md (freeze f03b034) on (0.45, 2) and (0.55, 2)."""
import argparse
import copy
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FREEZE = "f03b034"
OUT = ROOT / "experiments/closure_ledger/runs/20261126_alice_nulls"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V = _load("alice_vcurve_probe")
P = _load("alice_nulls_predictions")


def doublet_only(geo):
    """The same geometry with the modal basis cut to the doublet."""
    g = copy.copy(geo)
    idx = np.array([geo.iA, geo.iB])
    g.keep, g.omega_all = 2, geo.omega_all[idx]
    g.W1, g.W0 = geo.W1[np.ix_(idx, idx)].copy(), geo.W0[np.ix_(idx, idx)].copy()
    g.iA, g.iB = 0, 1
    return g


def identical_profiles(d):
    """Doublet-only geometry with the two mode profiles made identical inside the strip."""
    h = copy.copy(d)
    w = 0.5 * (d.W1[0, 0] + d.W1[1, 1])
    h.W1 = np.diag([w, w])
    h.W0 = np.array([[0.0, np.sign(d.W0[0, 1]) * w], [np.sign(d.W0[1, 0]) * w, 0.0]])
    return h


def evaluate(geo, eta):
    g0 = geo.calibrate(eta, geo.mean)
    r = geo.run(g0, geo.mean)
    return dict(defect=r["max_abs_defect"], g0=g0, rate_over_wbar=g0 / geo.mean, balance=r["balance"])


def cached(cache, key, fn, progress):
    path = Path(cache) / f"{key}.json" if cache else None
    if path is not None and path.exists():
        return json.loads(path.read_text())
    progress(key)
    out = fn()
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(out))
    return out


def run(cache=None, progress=print):
    started = time.time()
    frozen = json.loads((OUT / "frozen_predictions.json").read_text())
    pred = {(r["a"], r["eta"], round(r["tau"], 6)): r for r in frozen["nulls"]}
    rep = dict(freeze=FREEZE, cells=list(P.CELLS), m=P.M, null_ks=list(P.NULL_KS), points=[], boundary=[], cell_info=[])
    omega = np.sqrt(P.M * (P.M + 1))
    for a in P.CELLS:
        spec = V.BASE.solve(a, P.M, 2000)
        base = V.Sweep(a, P.M, 2000, 150, spec, 8.0, "sin2")
        wbar = base.mean
        eps = float(base.W1[base.iA, base.iA] / base.W1[base.iB, base.iB] - 1)
        und = V.Sweep(a, P.M, 2000, 150, spec, 8.0, "sin2").run(0.0, wbar)["max_abs_defect"]
        rep["cell_info"].append(dict(a=a, wbar=float(wbar), split=float(base.split), eps=eps,
                                     eps_over_split=float(abs(eps) / (base.split / wbar)), undetected=float(und)))
        for eta in P.ETAS:
            for k in P.NULL_KS:
                tau = float(k * np.pi / wbar)
                key = f"{a}_{eta}_{k}"
                geo = V.Sweep(a, P.M, 2000, 150, spec, tau, "sin2")
                row = cached(cache, key, lambda: dict(
                    full=evaluate(geo, eta), doublet=evaluate(doublet_only(geo), eta),
                    identical=evaluate(identical_profiles(doublet_only(geo)), eta),
                    K80=evaluate(V.Sweep(a, P.M, 2000, 80, spec, tau, "sin2"), eta)), progress)
                p = pred[(a, eta, round(float(k * np.pi / omega), 6))]
                row.update(a=a, eta=eta, k=k, tau=tau, predicted=p)
                rep["points"].append(row)
        for tau in P.BOUNDARY_TAUS:
            def boundary():
                g = V.Sweep(a, P.M, 2000, 150, spec, tau, "sin2")
                strengths = 0.05 * 8.0 / tau * 1.25 ** np.arange(50)
                reached = [g.absorbed_plus(x, g.mean) for x in strengths]
                i = int(np.argmax(reached))
                return dict(eta_max=float(reached[i]), rate_over_wbar_at_max=float(strengths[i] / g.mean))
            b = cached(cache, f"{a}_boundary_{tau}", boundary, progress)
            b.update(a=a, tau=tau, in_sample=P.BOUNDARY_IN_SAMPLE[tau])
            rep["boundary"].append(b)
        del spec
    # N = 4000 at the k = 8 points
    for a in P.CELLS:
        spec4 = None
        for eta in P.ETAS:
            def fine():
                nonlocal spec4
                spec4 = spec4 or V.BASE.solve(a, P.M, 4000)
                g = V.Sweep(a, P.M, 4000, 150, spec4, 8.0, "sin2")
                return evaluate(V.Sweep(a, P.M, 4000, 150, spec4, float(8 * np.pi / g.mean), "sin2"), eta)
            r = cached(cache, f"{a}_{eta}_8_N4000", fine, progress)
            next(x for x in rep["points"] if x["a"] == a and x["eta"] == eta and x["k"] == 8)["N4000"] = r
    rep.update(gates(rep))
    rep["seconds"] = time.time() - started
    return rep


def gates(rep):
    pts = rep["points"]
    outside = [(x["a"], x["eta"], x["k"], x["full"]["rate_over_wbar"]) for x in pts if x["full"]["rate_over_wbar"] > P.DOMAIN_RATE]
    g = {}
    n1 = [(x["a"], x["eta"], x["k"], x["full"]["defect"] / x["predicted"]["total"]) for x in pts]
    g["N1"] = dict(passed=all(0.5 <= r <= 2 for *_, r in n1), ratios=n1)
    n2 = [(x["a"], x["eta"], x["k"], abs(x["doublet"]["defect"] / x["full"]["defect"] - 1)) for x in pts]
    g["N2"] = dict(passed=all(d <= 0.03 for *_, d in n2), deviations=n2)
    eps_ok = all(3 <= c["eps_over_split"] <= 5 for c in rep["cell_info"])
    n3 = [(x["a"], x["k"], (x["doublet"]["defect"] - x["identical"]["defect"]) / x["predicted"]["mismatch"])
          for x in pts if x["eta"] == 0.5]
    g["N3"] = dict(passed=eps_ok and all(0.5 <= r <= 2 for *_, r in n3), eps_ok=eps_ok, ratios=n3)
    n4 = []
    for x in pts:
        if x["a"] == 0.45 and x["k"] in (5, 8):
            ex = P.terms(x["a"], P.M, x["tau"], "sin2", x["eta"])["exchange"]
            n4.append((x["eta"], x["k"], (x["identical"]["defect"] - ex) * (np.sqrt(6) * x["tau"]) ** 3 / P.C[x["eta"]]))
    g["N4"] = dict(passed=all(0.5 <= r <= 2 for *_, r in n4), ratios=n4)
    n5 = [(b["a"], b["tau"], b["eta_max"] - b["in_sample"]) for b in rep["boundary"]]
    g["N5"] = dict(passed=all(abs(d) <= 0.005 for *_, d in n5), deviations=n5)
    k_ok = all(abs(x["K80"]["defect"] - x["full"]["defect"]) <= max(1e-6, 0.05 * x["full"]["defect"]) for x in pts)
    n_ok = all(abs(x["N4000"]["defect"] - x["full"]["defect"]) <= max(1e-6, 0.1 * x["N4000"]["defect"]) for x in pts if "N4000" in x)
    balance = max(max(x[k]["balance"] for k in ("full", "doublet", "identical", "K80")) for x in pts)
    und = max(c["undetected"] for c in rep["cell_info"])
    g["N6"] = dict(passed=bool(k_ok and n_ok and balance <= 1e-8 and und <= 1e-8), modal=bool(k_ok), size=bool(n_ok),
                   balance=float(balance), undetected=float(und))
    g["outside_domain"] = outside
    failed = [k for k in ("N1", "N2", "N3", "N4", "N5") if not g[k]["passed"]]
    if outside:
        verdict, reasons = "NO_VERDICT", ["primary point outside the domain"]
    elif failed:
        verdict, reasons = "COUNTEREXAMPLE", failed
    elif g["N6"]["passed"]:
        verdict, reasons = "CONFIRMED", []
    else:
        verdict, reasons = "UNRESOLVED_NUMERICALLY", ["N6"]
    return dict(gates=g, verdict=verdict, reasons=reasons)


def markdown(rep):
    g = rep["gates"]
    L = [f"# Null residual of the Alice-handle exchange: held-out test against freeze `{rep['freeze']}`", "",
         f"**Verdict: {rep['verdict']}**" + (f" ({', '.join(rep['reasons'])})" if rep["reasons"] else ""), "",
         "Defect at the exact sin² nulls (full field, N = 2000, K = 150), measured against predicted.", "",
         "| cell | η | k | τ | measured | predicted | ratio | doublet only | identical profiles | peak rate / ω̄ |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for x in rep["points"]:
        L.append(f"| {x['a']} | {x['eta']} | {x['k']} | {x['tau']:.3f} | {x['full']['defect']:.3e} | {x['predicted']['total']:.3e} | "
                 f"{x['full']['defect'] / x['predicted']['total']:.2f} | {x['doublet']['defect']:.3e} | {x['identical']['defect']:.3e} | "
                 f"{x['full']['rate_over_wbar']:.3f} |")
    L += ["", "| gate | result |", "|---|---|"]
    r1 = [r for *_, r in g["N1"]["ratios"]]
    L.append(f"| N1 null values within ×2 | {'pass' if g['N1']['passed'] else 'FAIL'} ({min(r1):.2f}–{max(r1):.2f}) |")
    d2 = [d for *_, d in g["N2"]["deviations"]]
    L.append(f"| N2 doublet only = full within 3% | {'pass' if g['N2']['passed'] else 'FAIL'} (max {max(d2):.2%}) |")
    r3 = [r for *_, r in g["N3"]["ratios"]]
    L.append(f"| N3 mismatch | {'pass' if g['N3']['passed'] else 'FAIL'} (ε/(Δ/ω̄) " + ", ".join(
        f"{c['eps_over_split']:.2f}" for c in rep["cell_info"]) + f"; ratios {min(r3):.2f}–{max(r3):.2f}) |")
    r4 = [r for *_, r in g["N4"]["ratios"]]
    L.append(f"| N4 intrinsic (ωτ)⁻³ law | {'pass' if g['N4']['passed'] else 'FAIL'} ({min(r4):.2f}–{max(r4):.2f}) |")
    d5 = [d for *_, d in g["N5"]["deviations"]]
    L.append(f"| N5 depletion boundary | {'pass' if g['N5']['passed'] else 'FAIL'} (max |Δη_max| {max(abs(x) for x in d5):.4f}) |")
    L.append(f"| N6 numerics | {'pass' if g['N6']['passed'] else 'FAIL'} (K {g['N6']['modal']}, N {g['N6']['size']}, "
             f"balance {g['N6']['balance']:.1e}, undetected {g['N6']['undetected']:.1e}) |")
    L += ["", f"Points outside the domain: {g['outside_domain'] or 'none'}.", "", f"Run time {rep['seconds']:.0f} s.", ""]
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
