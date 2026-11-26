"""Held-out test of docs/alice_exchange_prereg.md (freeze 0fb3b6e).

The lowest Klein doublet of #355's handle exchanges its +m and -m circulating
channels over pi/Delta. A detector strip on the band absorbs during 0<t<tau,
either from the +m channel only (counter-rotating dashpot) or from both
(static dashpot). The +m doublet energy is read at the window centre on the
undamped reference (u) and pi/Delta later on the detected evolution (v). Every
reading is a quadratic form in the initial doublet state, so the largest
|u + v - 1| over all states is an eigenvalue.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy import linalg
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from geometrodynamics.bulk import alice_handle_ports as ah

ROOT = Path(__file__).resolve().parents[2]
FREEZE = "0fb3b6e"
SIZES = (2000, 4000)
KEEPS = (150, 300)
_spec = importlib.util.spec_from_file_location(
    "alice_exchange_predictions", Path(__file__).with_name("alice_exchange_predictions.py"))
PRED = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(PRED)
TAU, ETAS, HALF = PRED.TAU, PRED.ETAS, PRED.STRIP_HALF_WIDTH


def rotate_pairs(X):
    """(A, B) -> (-B, A) at every node, in each node's own frame."""
    Y = np.empty_like(X)
    Y[0::2], Y[1::2] = -X[1::2], X[0::2]
    return Y


def solve(a, m, n):
    loop = ah.Loop(a, n)
    sec = ah.Sector(loop, m, ah.KLEIN)
    lam, _, V = sec.modes
    return loop, sec, lam, V


class Geometry:
    """One cell at one size: modes, the detector strip and the doublet forms."""

    def __init__(self, a, m, n, keep, spectrum=None):
        loop, sec, lam, V = spectrum if spectrum is not None else solve(a, m, n)
        cos = sec.sector_of(V) > 0.5
        self.iA, self.iB = int(np.flatnonzero(cos)[0]), int(np.flatnonzero(~cos)[0])
        if max(self.iA, self.iB) >= keep:
            raise ValueError("doublet outside the kept modes")
        self.a, self.m, self.n, self.keep = a, m, n, keep
        self.omega_all = np.sqrt(lam[:keep])
        Vk = V[:, :keep]
        theta = a + loop.ell
        w = np.zeros(n)
        band = np.flatnonzero(loop.on_band)
        x = (theta[band] - np.pi / 2) / HALF
        w[band] = np.where(np.abs(x) < 1, np.cos(0.5 * np.pi * x) ** 2, 0.0)
        weight = np.repeat(loop.mass * w, 2)
        self.W1 = Vk.T @ (weight[:, None] * Vk)
        self.W0 = Vk.T @ (weight[:, None] * rotate_pairs(Vk))
        G = Vk.T @ (np.repeat(loop.mass, 2)[:, None] * rotate_pairs(Vk))
        self.wA, self.wB = self.omega_all[self.iA], self.omega_all[self.iB]
        self.split = abs(self.wA - self.wB)
        self.mean = 0.5 * (self.wA + self.wB)
        # doublet phase coordinates s = (qA, qB, pA, pB)
        d = [self.iA, self.iB]
        self.E_d = np.diag([self.wA ** 2 / 2, self.wB ** 2 / 2, 0.5, 0.5])
        H = np.pi * m * G[np.ix_(d, d)]               # L = p^T H q
        L_sym = np.zeros((4, 4))
        L_sym[2:, :2] = H / 2
        L_sym[:2, 2:] = H.T / 2
        mu, xi = linalg.eigh(L_sym, self.E_d)
        plus = xi[:, mu > 0]
        self.channel_eigen = mu.tolist()
        self.P_plus = (self.E_d @ plus) @ (self.E_d @ plus).T
        # E-orthonormal basis of the doublet: s = B c
        self.Bmat = np.diag([np.sqrt(2) / self.wA, np.sqrt(2) / self.wB, np.sqrt(2), np.sqrt(2)])
        plus_c = np.linalg.solve(self.Bmat, plus)
        self.plus_c = plus_c / np.sqrt(np.sum(plus_c ** 2, axis=0))
        self.band_overlap = float(abs(H[0, 1]) / (np.pi * m))

    def basis(self):
        """Modal (q, p) of the four E-orthonormal doublet basis states, shape (2K, 4)."""
        Y = np.zeros((2 * self.keep, 4))
        for j, (row, idx) in enumerate([(0, self.iA), (0, self.iB), (1, self.iA), (1, self.iB)]):
            Y[row * self.keep + idx, j] = self.Bmat[j, j]
        return Y

    def free(self, Y, t, omega=None):
        omega = self.omega_all if omega is None else omega
        ph = ah.phase(omega, np.array([float(t)]))[:, 0]
        c, s = np.cos(ph), np.sin(ph)
        q, p = Y[:self.keep], Y[self.keep:]
        return np.vstack([c[:, None] * q + (s / omega)[:, None] * p,
                          -(omega * s)[:, None] * q + c[:, None] * p])

    def doublet(self, Y):
        k = self.keep
        return np.vstack([Y[self.iA], Y[self.iB], Y[k + self.iA], Y[k + self.iB]])

    def energy_form(self, Y):
        q, p = Y[:self.keep], Y[self.keep:]
        return 0.5 * (p.T @ p + q.T @ (self.omega_all[:, None] ** 2 * q))

    def window(self, g0, rotation):
        """Integrate the four basis states through the detector window."""
        k = self.keep
        om2 = self.omega_all ** 2

        def rhs(t, y):
            Y = y[:8 * k].reshape(2 * k, 4)
            q, p = Y[:k], Y[k:]
            g = g0 * float(PRED.window(t))
            force = -g * (self.W1 @ p + rotation * (self.W0 @ q))
            dY = np.vstack([p, -om2[:, None] * q + force])
            dD = -0.5 * (p.T @ force + force.T @ p)        # absorbed-energy form
            return np.concatenate([dY.ravel(), dD.ravel()])
        y0 = np.concatenate([self.basis().ravel(), np.zeros(16)])
        sol = solve_ivp(rhs, (0.0, TAU), y0, method="DOP853", rtol=1e-11, atol=1e-13)
        if not sol.success:
            raise RuntimeError(sol.message)
        y = sol.y[:, -1]
        return y[:8 * k].reshape(2 * k, 4), y[8 * k:].reshape(4, 4)

    def absorbed_plus(self, g0, rotation):
        _, D = self.window(g0, rotation)
        return float(np.trace(self.plus_c.T @ D @ self.plus_c) / 2)

    def calibrate(self, eta, rotation):
        """The smallest g0 at which the +m channel loses eta (the weak-coupling branch).

        Absorption is not monotonic in g0: a very stiff medium drags the field
        with it and stops absorbing. Scan upward from zero, then bisect.
        """
        f = lambda g: self.absorbed_plus(g, rotation) - eta
        lo, g = 0.0, 1e-3
        while f(g) < 0:
            lo, g = g, 1.5 * g
            if g > 1e3:
                raise RuntimeError("depletion eta not reached")
        return brentq(f, lo, g, xtol=1e-13, rtol=1e-12)

    def run(self, g0, rotation, split_for_timing=None):
        split = self.split if split_for_timing is None else split_for_timing
        Y, D = self.window(g0, rotation)
        balance = float(np.max(np.abs(self.energy_form(Y) + D - np.eye(4))))
        centre = TAU / 2
        S_u = self.doublet(self.free(self.basis(), centre))
        U = S_u.T @ self.P_plus @ S_u
        t2 = centre + np.pi / split
        Y2 = self.free(Y, t2 - TAU)
        S_v = self.doublet(Y2)
        Vf = S_v.T @ self.P_plus @ S_v
        sym = lambda X: 0.5 * (X + X.T)
        defect = np.linalg.eigvalsh(sym(U + Vf - np.eye(4)))
        total = self.energy_form(Y2)
        dE = S_v.T @ self.E_d @ S_v
        leak = np.linalg.eigvalsh(sym(total - dE))
        return dict(max_abs_defect=float(np.max(np.abs(defect))), defect_range=[float(defect.min()), float(defect.max())],
                    leakage_max=float(leak.max()), balance=balance, t2=float(t2), U=U.tolist(), V=Vf.tolist(), D=D.tolist())


def cell_report(a, m, n, keep, progress=print, spectrum=None):
    geo = Geometry(a, m, n, keep, spectrum)
    out = dict(a=a, m=m, N=n, K=keep, omega_A=float(geo.wA), omega_B=float(geo.wB), split=float(geo.split),
               mean=float(geo.mean), band_overlap=geo.band_overlap, channel_eigenvalues=geo.channel_eigen,
               predicted_split=PRED.NECK.prediction(a, m)["predicted_delta"])
    out["undetected"] = geo.run(0.0, geo.mean)
    out["selective"] = {}
    for eta in ETAS:
        progress(f"a={a} m={m} N={n} K={keep} selective eta={eta}")
        g0 = geo.calibrate(eta, geo.mean)
        res = geo.run(g0, geo.mean)
        res.update(g0=g0, eta_check=geo.absorbed_plus(g0, geo.mean),
                   exact_neck_timing=geo.run(g0, geo.mean, out["predicted_split"])["max_abs_defect"])
        U, D = np.array(res["U"]), np.array(res["D"])
        res["record_vs_eta_u"] = float(np.max(np.abs(np.linalg.eigvalsh(0.5 * (D + D.T) - eta * U))))
        out["selective"][str(eta)] = res
    progress(f"a={a} m={m} N={n} K={keep} static control")
    g0 = geo.calibrate(0.5, 0.0)
    res = geo.run(g0, 0.0)
    res.update(g0=g0, eta_check=geo.absorbed_plus(g0, 0.0))
    out["control"] = res
    return out


def gates(cells):
    by = {(c["a"], c["m"], c["N"], c["K"]): c for c in cells}
    primary = KEEPS[0]
    per_size = {}
    for n in SIZES:
        rows = [c for c in cells if c["N"] == n and c["K"] == primary]
        g1 = all(c["selective"][str(e)]["max_abs_defect"] <= 0.05 * e for c in rows for e in ETAS)
        small = next(c for c in rows if (c["a"], c["m"]) == (0.3, 2))
        g2 = all(small["selective"][str(e)]["max_abs_defect"] <= 1e-4 for e in ETAS)
        g3 = all(c["control"]["max_abs_defect"] >= 0.25 for c in rows)
        per_size[str(n)] = dict(G1=g1, G2=g2, G3=g3)
    checks = dict(undetected=max(c["undetected"]["max_abs_defect"] for c in cells),
                  balance=max(max([c["control"]["balance"]] + [r["balance"] for r in c["selective"].values()])
                              for c in cells))
    k_ok, n_ok = True, True
    for (a, m) in PRED.CELLS:
        for n in SIZES:
            lo, hi = by[(a, m, n, KEEPS[0])], by[(a, m, n, KEEPS[1])]
            pairs = [(lo["selective"][str(e)]["max_abs_defect"], hi["selective"][str(e)]["max_abs_defect"]) for e in ETAS]
            pairs.append((lo["control"]["max_abs_defect"], hi["control"]["max_abs_defect"]))
            k_ok &= all(abs(x - y) <= max(1e-5, 0.1 * y) for x, y in pairs)
        c2, c4 = by[(a, m, SIZES[0], primary)], by[(a, m, SIZES[1], primary)]
        for e in ETAS:
            x, y = c2["selective"][str(e)]["max_abs_defect"], c4["selective"][str(e)]["max_abs_defect"]
            n_ok &= abs(x - y) <= max(1e-4, 0.2 * y)
    g4 = checks["undetected"] <= 1e-8 and checks["balance"] <= 1e-8 and k_ok and n_ok
    checks.update(modal_agreement=bool(k_ok), size_agreement=bool(n_ok), G4=bool(g4))
    both = [k for k in ("G1", "G2", "G3") if not any(per_size[n][k] for n in per_size)]
    if both:
        verdict, reasons = "COUNTEREXAMPLE", both
    elif all(all(v.values()) for v in per_size.values()) and g4:
        verdict, reasons = "CONFIRMED", []
    elif not g4 or any(not all(v.values()) for v in per_size.values()):
        verdict = "UNRESOLVED_NUMERICALLY"
        reasons = (["G4"] if not g4 else []) + sorted({k for v in per_size.values() for k, ok in v.items() if not ok})
    else:
        verdict, reasons = "NO_VERDICT", []
    return dict(per_size=per_size, checks=checks), verdict, reasons


def run(progress=print, cells=PRED.CELLS, sizes=SIZES, keeps=KEEPS, cache=None):
    started = time.time()
    results = []
    for a, m in cells:
        for n in sizes:
            spectrum = None
            for keep in keeps:
                path = None if cache is None else Path(cache) / f"cell_{a}_{m}_{n}_{keep}.json"
                if path is not None and path.exists():
                    results.append(json.loads(path.read_text()))
                    continue
                if spectrum is None:
                    progress(f"a={a} m={m} N={n}: modes")
                    spectrum = solve(a, m, n)
                row = cell_report(a, m, n, keep, progress, spectrum)
                if path is not None:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(json.dumps(row))
                results.append(row)
            del spectrum
    report = dict(freeze=FREEZE, solver="geometrodynamics/bulk/alice_handle_ports.py at 07d58de",
                  tau=TAU, etas=list(ETAS), strip_half_width=HALF, cells=results)
    if tuple(cells) == tuple(PRED.CELLS) and tuple(sizes) == SIZES and tuple(keeps) == KEEPS:
        g, verdict, reasons = gates(results)
        report.update(gates=g, verdict=verdict, reasons=reasons)
    report["seconds"] = time.time() - started
    return report


def markdown(rep):
    lines = [f"# Alice-handle exchange under detector depletion: test against freeze `{rep['freeze']}`", "",
             f"**Verdict: {rep.get('verdict', 'n/a')}**" + (f" ({', '.join(rep['reasons'])})" if rep.get("reasons") else ""), "",
             "Largest |u + v − 1| over all doublet states (K = 150).", "",
             "| a | m | N | Δ (FEM) | Δτ | selective η=0.1 | η=0.5 | η=0.9 | static η=0.5 | leakage (η=0.5) |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for c in rep["cells"]:
        if c["K"] != KEEPS[0]:
            continue
        s = c["selective"]
        lines.append(f"| {c['a']} | {c['m']} | {c['N']} | {c['split']:.4e} | {c['split'] * rep['tau']:.2e} | "
                     f"{s['0.1']['max_abs_defect']:.2e} | {s['0.5']['max_abs_defect']:.2e} | {s['0.9']['max_abs_defect']:.2e} | "
                     f"{c['control']['max_abs_defect']:.3f} | {s['0.5']['leakage_max']:.1e} |")
    if "gates" in rep:
        lines += ["", "| gate | N=2000 | N=4000 |", "|---|---|---|"]
        for k in ("G1", "G2", "G3"):
            lines.append(f"| {k} | {'pass' if rep['gates']['per_size']['2000'][k] else 'FAIL'} | "
                         f"{'pass' if rep['gates']['per_size']['4000'][k] else 'FAIL'} |")
        ch = rep["gates"]["checks"]
        lines += [f"| G4 | {'pass' if ch['G4'] else 'FAIL'} (undetected {ch['undetected']:.1e}, balance {ch['balance']:.1e}, "
                  f"K agreement {ch['modal_agreement']}, N agreement {ch['size_agreement']}) | |"]
    lines += ["", f"Run time {rep['seconds']:.0f} s.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "experiments/closure_ledger/runs/20261125_alice_exchange")
    parser.add_argument("--cache", type=Path, default=None)
    args = parser.parse_args(argv)
    rep = run(lambda s: print(s, flush=True), cache=args.cache)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "probe.json").write_text(json.dumps(rep, indent=1) + "\n")
    (args.out / "probe.md").write_text(markdown(rep))
    print(markdown(rep))


if __name__ == "__main__":
    main()
