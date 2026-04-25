#!/usr/bin/env python3
"""
scripts/experiment_transport_overlap.py
========================================

User-named milestone (axioms §8): drop the broad Hopf-scalar
search and compute the overlap integral

    T_{l,l+2} = ∫ u_{l,0}(r*) u_{l+2,0}(r*) w(r*) dr*

on the tortoise grid that just gave pinhole = ∑V_max within 1%.
This is naturally an off-diagonal matrix element, structurally
exactly the role "transport" plays in the QCD Hamiltonian's
same-partition coupling.

The eigensolver returns ground-state radial wave functions
u_{l,0}(r*) for l = 1, 2, ..., 5.  The model's same-partition
off-diagonal couples basis states (k=1,p) ↔ (k=3,p) ↔ (k=5,p),
which we identify with l = 1, 3, 5 (odd partial waves only — by
the same parity choice that gives the closure pass-counts).
The three relevant pairs are therefore (1,3), (3,5), (1,5).

Weighting choices probed:
  w = 1               — raw normalised overlap
  w = V_{l+2} − V_l   — exact QM perturbation operator (the
                        Hamiltonian shift between two l shells;
                        this IS the off-diagonal matrix element
                        of the cross-l Hamiltonian variation)
  w = V(r, l_avg)     — average effective potential
  w = (l₂² − l₁²)/r²  — pure centrifugal difference

Each candidate is also joint-fit-tested: pin transport to that
value, let N + remaining residuals descend, report (N, max_rel_err).
A candidate is "viable" if joint-fit err stays close to the
unpinned baseline (0.0161) and within 2× of the pinhole-pinned
joint-fit err (0.0295) for fairness.
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
    V_tangherlini,
)
from geometrodynamics.constants import R_MID

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


def joint_fit_transport(transport_value: float) -> tuple[int, float]:
    base = replace(_TEMPLATE, transport=transport_value)
    cur = base
    cur_err = _err(_PDG_MASSES, cur)
    free = [
        ("beta", np.arange(380, 561, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("pinhole", np.arange(18.0, 27.0, 0.25),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.06, 0.22, 0.005),
            lambda p, v: replace(p, resistance=v)),
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

    # ── Solve modes for l = 1, 2, 3, 4, 5 ─────────────────────────
    funcs = {}
    for l in [1, 2, 3, 4, 5]:
        oms, fns, _ = solve_radial_modes(l=l, N=80, n_modes=4)
        funcs[l] = (oms[0], fns[0])
    r_phys = funcs[1][1]['r_phys']
    r_star = np.array([r_to_rstar(rv, R_MID) for rv in r_phys])

    def f_metric(rv): return 1.0 - (R_MID / rv) ** 2

    # Pairs that appear as same-partition off-diagonals in the model
    same_partition_pairs = [(1, 3), (3, 5), (1, 5)]

    # Weighting candidates
    def make_weights():
        return {
            "w=1 (raw)":           np.ones_like(r_phys),
            "w=1/r^2":             1.0 / r_phys ** 2,
            "w=f(r)/r^2":          f_metric(r_phys) / r_phys ** 2,
            "w=V(r,l_avg)":        None,            # filled per-pair
            "w=V(r,l_max)":        None,            # filled per-pair
            "w=V_l2-V_l1":         None,            # filled per-pair
            "w=(l2-l1)*f/r^2":     None,            # filled per-pair
            "w=cos(pi r/2)":       np.cos(math.pi * r_phys / 2.0),
            "w=sin(pi r/2)":       np.sin(math.pi * r_phys / 2.0),
        }

    def per_pair_weight(name: str, l1: int, l2: int):
        if name == "w=V(r,l_avg)":
            return V_tangherlini(r_phys, (l1 + l2) // 2, R_MID)
        if name == "w=V(r,l_max)":
            return V_tangherlini(r_phys, l2, R_MID)
        if name == "w=V_l2-V_l1":
            return f_metric(r_phys) * (l2*(l2+2) - l1*(l1+2)) / r_phys**2
        if name == "w=(l2-l1)*f/r^2":
            return (l2 - l1) * f_metric(r_phys) / r_phys ** 2
        return None

    def normed_overlap(u1, u2, w, x):
        n1 = np.sqrt(trap(u1*u1, x))
        n2 = np.sqrt(trap(u2*u2, x))
        return float(trap(u1*u2*w, x) / (n1 * n2))

    # ── Compute per-pair overlaps for each weighting on tortoise ──
    print("\n" + "=" * 82)
    print("Pairwise weighted overlaps T_{l1,l2} on tortoise grid (r*)")
    print("=" * 82)
    weight_names = list(make_weights().keys())
    pair_table: dict[str, dict[str, float]] = {}
    for name in weight_names:
        per_pair = {}
        for (l1, l2) in same_partition_pairs:
            u1 = funcs[l1][1]['u_half']
            u2 = funcs[l2][1]['u_half']
            base_w = make_weights()[name]
            w = base_w if base_w is not None else per_pair_weight(name, l1, l2)
            T = normed_overlap(u1, u2, w, r_star)
            per_pair[f"({l1},{l2})"] = T
        per_pair["mean"] = float(np.mean(list(per_pair.values())))
        per_pair["geomean_signed"] = float(
            np.sign(np.prod(list(per_pair.values()))) *
            abs(np.prod(list(per_pair.values()))) ** (1.0 / 3.0)
        )
        pair_table[name] = per_pair

    print(f"{'weight':<25}  {'(1,3)':>8}  {'(3,5)':>8}  {'(1,5)':>8}  "
          f"{'mean':>8}  {'|mean-0.54|':>12}")
    print("-" * 82)
    for name, per_pair in pair_table.items():
        m = per_pair['mean']
        d = abs(m - 0.54)
        marker = " ← target" if d < 0.05 else ""
        print(f"{name:<25}  {per_pair['(1,3)']:>+8.4f}  "
              f"{per_pair['(3,5)']:>+8.4f}  {per_pair['(1,5)']:>+8.4f}  "
              f"{m:>+8.4f}  {d:>12.4f}{marker}")

    # ── Joint-fit each (mean) value ───────────────────────────────
    print("\n" + "=" * 82)
    print("Joint-fit at transport = mean(weighted overlap)")
    print("=" * 82)
    print(f"{'weight':<25}  {'mean overlap':>12}  {'N':>5}  {'max rel err':>12}")
    print("-" * 82)
    joint_results: dict[str, dict] = {}
    for name, per_pair in pair_table.items():
        v = per_pair['mean']
        if 0.30 <= v <= 0.80:
            N, e = joint_fit_transport(v)
        else:
            N, e = None, None
        N_s = f"{N}" if N else "—"
        e_s = f"{e:.4e}" if e is not None else "—"
        joint_results[name] = {"mean_overlap": v, "N": N, "err": e}
        print(f"{name:<25}  {v:>+12.4f}  {N_s:>5}  {e_s:>12}")

    # ── Identify the strongest candidate ───────────────────────────
    print("\n" + "=" * 82)
    print("Best candidates (joint-fit err < 0.04)")
    print("=" * 82)
    best = sorted(
        ((name, r) for name, r in joint_results.items()
         if r["err"] is not None and r["err"] < 0.04),
        key=lambda kv: kv[1]["err"],
    )
    for name, r in best:
        diff = abs(r["mean_overlap"] - 0.54)
        print(f"  {name}: T={r['mean_overlap']:.4f} (Δ={diff:+.4f}), "
              f"N={r['N']}, err={r['err']:.4e}")

    out = {
        "pair_overlaps": pair_table,
        "joint_fit": joint_results,
        "fitted_transport": 0.54,
    }
    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nFull results written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
