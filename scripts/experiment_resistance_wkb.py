#!/usr/bin/env python3
"""
scripts/experiment_resistance_wkb.py
=====================================

User-named milestone (axioms §8): derive resistance from a WKB
tunneling integral, since resistance is literally the exponential
suppression coefficient in the same-partition off-diagonal:

    H_off(k1, k2) = -transport · exp(-resistance · dk) · cos(phase · dk)

The WKB factor for tunneling under a barrier V_barrier above an
energy ω² is exp(-∫ √(V_barrier(r*) - ω²) dr*), with the integral
running between the two classical turning points (where V = ω²).

If resistance = α controls exp(-α·dk), it must equal κ / dk, where

    κ_{l1,l2} = ∫ √(V_barrier(r*) - ω²) dr*    (over the forbidden region)

For a single resistance to apply across pairs (1,3), (3,5), (1,5),
the ratio κ_{l1,l2} / dk must be approximately constant.  The
model uses dk = max(l1, l2) at the lock (winding_mode = "max"), so
we test that convention plus dk = |l2 - l1| as a control.

Conventions probed:
  ω²       ∈ {ω(l1)², ω(l2)², (ω(l1)² + ω(l2)²)/2, ω(l1)·ω(l2)}
  V_barrier ∈ {V_{l1}, V_{l2}, (V_{l1}+V_{l2})/2, max_{r}(V_{l1},V_{l2})}
  dk       ∈ {max(l1,l2), l2-l1}
  domain   ∈ {classical forbidden region, full mode domain}

For each combination we report κ per pair and the resulting
resistance candidate (mean of κ/dk).  Joint-fit at the candidate.
"""

from __future__ import annotations

import os
import sys
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import argparse
import json
import math
from dataclasses import asdict, replace
from typing import Callable

import numpy as np

from geometrodynamics.qcd.quark_spectrum import (
    OBSERVED_MASSES_MEV as _PDG_MASSES,
    QUARK_ACTION_BASE,
    QUARK_SPECIES,
    QuarkParams,
    extract_physical_spectrum,
)
from geometrodynamics.tangherlini.radial import (
    solve_radial_modes,
    r_to_rstar,
    rstar_to_r,
    V_tangherlini,
)
from geometrodynamics.constants import R_MID, R_OUTER, R_INNER

trap = np.trapezoid

# ── Constraint-reduced template ───────────────────────────────────
_TEMPLATE = QuarkParams(
    action_base=QUARK_ACTION_BASE,
    beta=466 * math.pi / 2.0,
    gamma_q=0.10,
    u_q_form="k_minus_2",
    phase=0.0,
    transport=0.54,
    pinhole=22.25,
    resistance=0.14,
    partition_mixing=0.0,
    winding_mode="max",
    resistance_model="exponential",
    depth_cost_mode="tunnel_only",
    uplift_mode="partition_asymmetric",
    uplift_asymmetry=24.0 / 25.0,
    spectrum_zero_mode="min_eigenvalue",
    chi_q_k3=20.0,
    eta_k3k5_minus=5.0,
)


def _err(observed, params):
    target = np.array([observed[s] for s in QUARK_SPECIES], dtype=float)
    skip = QUARK_SPECIES.index("u")
    try:
        sm = extract_physical_spectrum(params)
    except Exception:
        return float("inf")
    if sm["d"] <= 1e-6:
        return float("inf")
    scale = observed["d"] / sm["d"]
    pred = np.array([sm[s] * scale for s in QUARK_SPECIES], dtype=float)
    rel = np.abs(pred - target) / target
    return float(np.max([rel[i] for i in range(6) if i != skip]))


def joint_fit_resistance(resistance_value: float) -> tuple[int, float]:
    base = replace(_TEMPLATE, resistance=resistance_value)
    cur = base
    cur_err = _err(_PDG_MASSES, cur)
    free = [
        ("beta", np.arange(380, 561, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("transport", np.arange(0.40, 0.80, 0.01),
            lambda p, v: replace(p, transport=v)),
        ("pinhole", np.arange(18.0, 27.0, 0.25),
            lambda p, v: replace(p, pinhole=v)),
    ]
    for _ in range(6):
        improved = False
        for name, vals, apply in free:
            best_v, best_e = None, cur_err
            for v in vals:
                e = _err(_PDG_MASSES, apply(cur, float(v)))
                if e < best_e:
                    best_e, best_v = e, float(v)
            if best_v is not None:
                cur = apply(cur, best_v)
                cur_err = best_e
                improved = True
        if not improved:
            break
    return int(round(cur.beta * 2.0 / math.pi)), cur_err


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    args = p.parse_args()

    # ── Solve modes ──────────────────────────────────────────────
    omegas, funcs = {}, {}
    for l in [1, 3, 5]:
        oms, fns, _ = solve_radial_modes(l=l, N=80, n_modes=4)
        omegas[l] = float(oms[0])
        funcs[l] = fns[0]
    r_phys = funcs[1]['r_phys']
    r_star = np.array([r_to_rstar(rv, R_MID) for rv in r_phys])

    # Use a denser uniform tortoise grid for accurate WKB integration
    N_DENSE = 1024
    rs_dense = np.linspace(r_star[0], r_star[-1], N_DENSE)
    r_dense = np.array([rstar_to_r(s, R_MID) for s in rs_dense])

    print("Mode frequencies:")
    for l, w in omegas.items():
        print(f"  ω(l={l}, n=0) = {w:.4f}  (ω² = {w**2:.4f})")

    pairs = [(1, 3), (3, 5), (1, 5)]

    def V_l(rv, l): return V_tangherlini(rv, l, R_MID)

    def kappa_integral(omega2: float, V_arr: np.ndarray, x: np.ndarray) -> tuple[float, float, float]:
        """
        ∫ √(V - ω²) dx over the region V > ω².
        Returns (κ, r*_inner, r*_outer) of the forbidden region.
        """
        diff = V_arr - omega2
        forbidden = diff > 0
        if not np.any(forbidden):
            return 0.0, float("nan"), float("nan")
        # find contiguous forbidden region (assume one)
        idx = np.where(forbidden)[0]
        i_lo, i_hi = idx[0], idx[-1]
        x_seg = x[i_lo:i_hi + 1]
        d_seg = diff[i_lo:i_hi + 1]
        d_seg = np.maximum(d_seg, 0.0)  # numerical safety at boundaries
        kappa = float(trap(np.sqrt(d_seg), x_seg))
        return kappa, float(x[i_lo]), float(x[i_hi])

    # ── Run all conventions ──────────────────────────────────────
    omega_choices = {
        "omega²(l1)":      lambda l1, l2: omegas[l1] ** 2,
        "omega²(l2)":      lambda l1, l2: omegas[l2] ** 2,
        "mean ω²":         lambda l1, l2: 0.5 * (omegas[l1]**2 + omegas[l2]**2),
        "ω(l1)·ω(l2)":     lambda l1, l2: omegas[l1] * omegas[l2],
    }
    barrier_choices = {
        "V_l1":         lambda l1, l2: V_l(r_dense, l1),
        "V_l2":         lambda l1, l2: V_l(r_dense, l2),
        "mean V":       lambda l1, l2: 0.5 * (V_l(r_dense, l1) + V_l(r_dense, l2)),
        "max V":        lambda l1, l2: np.maximum(V_l(r_dense, l1), V_l(r_dense, l2)),
    }
    dk_choices = {
        "dk=max(l1,l2)": lambda l1, l2: max(l1, l2),
        "dk=l2-l1":      lambda l1, l2: l2 - l1,
    }

    print()
    print("=" * 100)
    print("WKB κ_{l1,l2} = ∫ √(V_barrier - ω²) dr* over forbidden region, on tortoise grid")
    print("=" * 100)
    print(f"{'omega2':<14}  {'V_barrier':<10}  {'dk':<14}  "
          f"{'(1,3)':>8}  {'(3,5)':>8}  {'(1,5)':>8}  {'mean':>8}  "
          f"{'spread':>8}  {'|mean-0.14|':>12}")
    print("-" * 100)
    all_results: dict[str, dict] = {}
    for o_name, o_fn in omega_choices.items():
        for v_name, v_fn in barrier_choices.items():
            for dk_name, dk_fn in dk_choices.items():
                per_pair = {}
                for (l1, l2) in pairs:
                    omega2 = o_fn(l1, l2)
                    V_arr = v_fn(l1, l2)
                    dk = dk_fn(l1, l2)
                    kappa, _, _ = kappa_integral(omega2, V_arr, rs_dense)
                    per_pair[f"({l1},{l2})"] = {
                        "kappa": kappa,
                        "dk": dk,
                        "alpha": kappa / dk if dk > 0 else float("nan"),
                    }
                alphas = [per_pair[f"({l1},{l2})"]["alpha"] for (l1, l2) in pairs]
                mean = float(np.mean(alphas))
                spread = float(np.max(alphas) - np.min(alphas))
                d = abs(mean - 0.14)
                key = f"{o_name} | {v_name} | {dk_name}"
                all_results[key] = {
                    "per_pair": per_pair,
                    "alpha_mean": mean,
                    "alpha_spread": spread,
                }
                marker = " ←" if d < 0.02 and spread < 0.05 else ""
                print(f"{o_name:<14}  {v_name:<10}  {dk_name:<14}  "
                      f"{alphas[0]:>+8.4f}  {alphas[1]:>+8.4f}  {alphas[2]:>+8.4f}  "
                      f"{mean:>+8.4f}  {spread:>+8.4f}  {d:>12.4f}{marker}")

    # ── Joint-fit the best candidates ─────────────────────────────
    print()
    print("=" * 80)
    print("Joint-fit at most-promising candidates (alpha_mean closest to 0.14)")
    print("=" * 80)
    sorted_candidates = sorted(
        all_results.items(),
        key=lambda kv: (abs(kv[1]["alpha_mean"] - 0.14)
                        + 0.5 * kv[1]["alpha_spread"]),
    )[:8]
    print(f"{'convention':<55}  {'alpha':>8}  {'spread':>8}  "
          f"{'N':>5}  {'err':>10}")
    print("-" * 90)
    for key, res in sorted_candidates:
        a = res["alpha_mean"]
        if 0.06 <= a <= 0.30:
            N, e = joint_fit_resistance(a)
        else:
            N, e = None, None
        N_s = f"{N}" if N else "—"
        e_s = f"{e:.4e}" if e is not None else "—"
        res["joint_fit_N"] = N
        res["joint_fit_err"] = e
        print(f"{key:<55}  {a:>+8.4f}  {res['alpha_spread']:>8.4f}  "
              f"{N_s:>5}  {e_s:>10}")

    out = {
        "omegas": {str(l): w for l, w in omegas.items()},
        "all_results": all_results,
    }
    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nFull results written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
