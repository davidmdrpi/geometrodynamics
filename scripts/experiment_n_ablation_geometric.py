#!/usr/bin/env python3
"""
scripts/experiment_n_ablation_geometric.py
============================================

User-named decisive test: re-run the N-stability ablation, but
this time with all three residuals pinned to their geometric
derivations (transport from V-difference overlap, pinhole from
∑V_max, resistance from transport·ln(α_q ratio)).

If N=466 remains stable under anchor changes and per-species
mass perturbations, it becomes a real interpretive target.
If it still drifts, then N remains phenomenological while the
residual sector is largely geometrized.

Comparison points:
  - Original ablation (residuals FREE):    [432, 510], width 78
  - All-broad-pinned ablation:              [538, 540], width  2 (err 11.7%)
  - This run (all-geometric-pinned):        ?

Free knob: N only.  All other knobs at the constraint-reduced
defaults (k_5 = 5 readings + γ_q = 1/10 + phase = 0).
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
from geometrodynamics.tangherlini.alpha_q import derive_alpha_q
from geometrodynamics.constants import R_MID, R_OUTER, R_INNER

trap = np.trapezoid


def derive_residuals() -> dict:
    """
    Compute the three derived geometric residuals from
    tangherlini.radial + tangherlini.alpha_q on the eigensolver
    tortoise grid.
    """
    # Modes for l = 1..5 needed for both pinhole and transport
    funcs = {}
    omegas = {}
    modes_for_alpha = {}
    for l in range(1, 6):
        oms, fns, _ = solve_radial_modes(l=l, N=80, n_modes=4)
        funcs[l] = fns[0]
        omegas[l] = float(oms[0])
        modes_for_alpha[l] = {"omega": oms, "funcs": fns}

    r_phys = funcs[1]['r_phys']
    r_star = np.array([r_to_rstar(rv, R_MID) for rv in r_phys])

    # ── Transport: mean ⟨u_l | V_{l+2} − V_l | u_{l+2}⟩ over the
    # three same-partition pairs on tortoise grid ───────────────
    def f_metric(rv): return 1.0 - (R_MID / rv) ** 2
    pair_overlaps = []
    for (l1, l2) in [(1, 3), (3, 5), (1, 5)]:
        u1 = funcs[l1]['u_half']
        u2 = funcs[l2]['u_half']
        w = f_metric(r_phys) * (l2*(l2+2) - l1*(l1+2)) / r_phys**2
        n1 = np.sqrt(trap(u1*u1, r_star))
        n2 = np.sqrt(trap(u2*u2, r_star))
        T = float(trap(u1*u2*w, r_star) / (n1*n2))
        pair_overlaps.append(T)
    transport = float(np.mean(pair_overlaps))

    # ── Pinhole: Σ V_max(l=1..5) on tortoise grid ────────────────
    rs_min = r_to_rstar(R_MID + 5e-4, R_MID)
    rs_max = r_to_rstar(R_OUTER - 5e-4, R_MID)
    N = 80
    x = np.cos(np.pi * np.arange(N + 1) / N)
    L = (rs_max - rs_min) / 2.0
    rsg = rs_min + L * (1.0 - x)
    from geometrodynamics.tangherlini.radial import rstar_to_r
    rg = np.array([rstar_to_r(s, R_MID) for s in rsg])
    pinhole = float(sum(np.max(V_tangherlini(rg, l, R_MID))
                        for l in range(1, 6)))

    # ── Resistance: transport · ln(α_q(k_5) / α_q(k_1)) ──────────
    aq = derive_alpha_q(modes_for_alpha)
    ln_ratio = math.log(aq[(5, 0)] / aq[(1, 0)])
    resistance = transport * ln_ratio

    return {
        "transport": transport,
        "pinhole": pinhole,
        "resistance": resistance,
        "ln_alpha_q_ratio": float(ln_ratio),
        "pair_overlaps": pair_overlaps,
        "alpha_q_5_0": float(aq[(5, 0)]),
    }


def err_factory(observed: dict[str, float], anchor_species: str):
    target = np.array([observed[s] for s in QUARK_SPECIES], dtype=float)
    skip = [QUARK_SPECIES.index("u")]

    def err(params: QuarkParams) -> float:
        try:
            sm = extract_physical_spectrum(params)
        except Exception:
            return float("inf")
        anchor = sm[anchor_species]
        if anchor <= 1e-6:
            return float("inf")
        scale = observed[anchor_species] / anchor
        pred = np.array([sm[s] * scale for s in QUARK_SPECIES], dtype=float)
        rel = np.abs(pred - target) / target
        keep = [i for i in range(len(QUARK_SPECIES)) if i not in skip]
        return float(np.max(rel[keep]))
    return err


def best_N(template: QuarkParams,
           err: Callable[[QuarkParams], float]) -> tuple[int, float]:
    """Coordinate descent over N only — all other knobs are pinned."""
    cur = template
    cur_err = err(cur)
    grid = np.arange(380, 561, 2)
    for _ in range(3):
        improved = False
        best_v, best_e = None, cur_err
        for v in grid:
            e = err(replace(cur, beta=float(v) * math.pi / 2.0))
            if e < best_e:
                best_e, best_v = e, float(v)
        if best_v is not None:
            cur = replace(cur, beta=best_v * math.pi / 2.0)
            cur_err = best_e
            improved = True
        if not improved:
            break
    return int(round(cur.beta * 2.0 / math.pi)), cur_err


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    args = p.parse_args()

    derived = derive_residuals()
    print("Derived residuals (from tangherlini machinery on tortoise grid):")
    for k in ("transport", "pinhole", "resistance", "ln_alpha_q_ratio"):
        print(f"  {k:>20}: {derived[k]:.6f}")

    template = QuarkParams(
        action_base=QUARK_ACTION_BASE,
        beta=466 * math.pi / 2.0,
        gamma_q=0.10,
        u_q_form="k_minus_2",
        phase=0.0,
        transport=derived["transport"],
        pinhole=derived["pinhole"],
        resistance=derived["resistance"],
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

    perturbations = {
        "PDG":           _PDG_MASSES,
        "PDG_x_1.10":    {s: m*1.10 for s, m in _PDG_MASSES.items()},
        "PDG_x_0.90":    {s: m*0.90 for s, m in _PDG_MASSES.items()},
        "c_x_1.10":      {**_PDG_MASSES, "c": _PDG_MASSES["c"]*1.10},
        "b_x_1.10":      {**_PDG_MASSES, "b": _PDG_MASSES["b"]*1.10},
        "t_x_1.10":      {**_PDG_MASSES, "t": _PDG_MASSES["t"]*1.10},
        "t_x_0.90":      {**_PDG_MASSES, "t": _PDG_MASSES["t"]*0.90},
        "all_perturb_5pct": {
            "u": _PDG_MASSES["u"]*1.05, "d": _PDG_MASSES["d"]*0.95,
            "s": _PDG_MASSES["s"]*1.05, "c": _PDG_MASSES["c"]*0.95,
            "b": _PDG_MASSES["b"]*1.05, "t": _PDG_MASSES["t"]*0.95,
        },
    }

    print("\n" + "=" * 60)
    print("N-ablation (all three residuals fixed at derived values)")
    print("=" * 60)
    results: dict = {"derived": derived}

    print("\n── (1) Mass perturbations (anchor=d) ──")
    print(f"{'scenario':<22}  {'best N':>7}  {'err':>10}")
    n_values_pert = []
    for label, masses in perturbations.items():
        e = err_factory(masses, "d")
        N, err_val = best_N(template, e)
        results[label] = {"N": N, "err": err_val}
        if err_val < 0.5:
            n_values_pert.append(N)
        print(f"{label:<22}  {N:>7d}  {err_val:>10.4e}")
    if n_values_pert:
        print(f"  N range = [{min(n_values_pert)}, {max(n_values_pert)}], "
              f"width = {max(n_values_pert) - min(n_values_pert)}")

    print("\n── (2) Anchor-species ablation (PDG masses) ──")
    print(f"{'anchor':<22}  {'best N':>7}  {'err':>10}")
    n_values_anchor = []
    for sp in ("d", "s", "c", "b", "t"):
        e = err_factory(_PDG_MASSES, sp)
        N, err_val = best_N(template, e)
        key = f"anchor={sp}"
        results[key] = {"N": N, "err": err_val}
        if err_val < 0.5:
            n_values_anchor.append(N)
        print(f"{key:<22}  {N:>7d}  {err_val:>10.4e}")
    if n_values_anchor:
        print(f"  N range = [{min(n_values_anchor)}, {max(n_values_anchor)}], "
              f"width = {max(n_values_anchor) - min(n_values_anchor)}")

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    all_n = [r["N"] for k, r in results.items()
             if k != "derived" and r["err"] < 0.5]
    if all_n:
        print(f"All ablations (well-fit only):  "
              f"N range [{min(all_n)}, {max(all_n)}], "
              f"width = {max(all_n) - min(all_n)}, "
              f"unique values = {sorted(set(all_n))}")

    print("\nReference comparison:")
    print(f"  Residuals FREE (original ablation):       N range [432, 510], width 78")
    print(f"  All broad scalars (Hopf 0.5, etc):        N range [538, 540], width 2 (err 11.7%)")
    if all_n:
        print(f"  All geometric derivations (this run):     N range [{min(all_n)}, {max(all_n)}], "
              f"width {max(all_n) - min(all_n)}")

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(results, fh, indent=2)
        print(f"\nFull results written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
