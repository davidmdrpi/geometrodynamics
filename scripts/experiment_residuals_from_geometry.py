#!/usr/bin/env python3
"""
scripts/experiment_residuals_from_geometry.py
==============================================

User-named milestone: replace one residual knob at a time with
geometry-derived quantities from the existing codebase, then
re-run the N-ablation.  If N stops drifting after those are
derived, it may become meaningful again; if it still drifts,
β should stay explicitly phenomenological.

Mapping (this script's choices, documented inline):

  transport  ←  hopf.hopf_connection(0) = 1/2
                (canonical Hopf connection at χ=0; the
                equator-to-pole holonomy coefficient on S³)

  resistance ←  (α_q(5,0) − α_q(1,0)) / 2 ≈ 0.1473
                (half-range of throat-flux ratios across the
                closure pass-count band; α_q from
                tangherlini.alpha_q.derive_alpha_q)

  pinhole    ←  Σ_{l=1}^{5} V_max(l) ≈ 21.80
                (cumulative Tangherlini centrifugal-barrier
                height summed over partial waves up to the
                heaviest pass count k_5 = 5)

For each substitution mode we hold the four shell-index
constraints (ε=24/25, η=5, χ=20, phase=0) and γ_q = 1/10 fixed,
plus pin the substituted residual to its derived value.  Then we
let only N (and any non-substituted residuals) vary, and check
how N moves under per-species mass perturbations.

Stability metric: range of N across the four mass scenarios
(PDG, c×1.10, b×1.10, t×1.10).  Baseline (no substitution) is
the §8 N-ablation result.
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
from geometrodynamics.hopf.connection import hopf_connection
from geometrodynamics.tangherlini.radial import (
    solve_radial_modes,
    V_tangherlini,
)
from geometrodynamics.tangherlini.alpha_q import derive_alpha_q
from geometrodynamics.constants import R_MID, R_OUTER, R_INNER


# ─── Derive the substitution scalars ───────────────────────────────
def derive_scalars() -> dict:
    # Hopf scalar
    transport_geom = float(hopf_connection(0.0))  # = 1/2

    # alpha_q range over k=1..5 (= partial waves l=1..5)
    modes = {l: None for l in range(1, 6)}
    for l in modes:
        oms, fns, _ = solve_radial_modes(l, N=80, n_modes=4)
        modes[l] = {"omega": oms, "funcs": fns}
    aq = derive_alpha_q(modes)
    # range from k=1 to k=5, halved (per-pair "step size")
    resistance_geom = (aq[(5, 0)] - aq[(1, 0)]) / 2.0

    # Sum of centrifugal-barrier maxima over l=1..5
    rs = np.linspace(R_INNER + 0.01, R_OUTER - 0.01, 400)
    pinhole_geom = float(sum(np.max(V_tangherlini(rs, l, R_MID))
                             for l in range(1, 6)))

    return {
        "transport_geom": transport_geom,
        "resistance_geom": float(resistance_geom),
        "pinhole_geom":   pinhole_geom,
        "alpha_q_table":  {f"{k}": float(v) for k, v in aq.items()},
    }


# ─── Constraint-reduced template ───────────────────────────────────
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


def _err_for_masses(observed: dict[str, float]) -> Callable[[QuarkParams], float]:
    target = np.array([observed[s] for s in QUARK_SPECIES], dtype=float)
    skip = [QUARK_SPECIES.index("u")]

    def err(params: QuarkParams) -> float:
        try:
            sm = extract_physical_spectrum(params)
        except Exception:
            return float("inf")
        anchor_val = sm["d"]
        if anchor_val <= 1e-6:
            return float("inf")
        scale = observed["d"] / anchor_val
        pred = np.array([sm[s] * scale for s in QUARK_SPECIES], dtype=float)
        rel = np.abs(pred - target) / target
        keep = [i for i in range(len(QUARK_SPECIES)) if i not in skip]
        return float(np.max(rel[keep]))
    return err


def coord_descent(
    base: QuarkParams,
    err: Callable[[QuarkParams], float],
    free_axes: list[tuple[str, np.ndarray, Callable]],
    max_rounds: int = 6,
) -> tuple[QuarkParams, float]:
    cur = base
    cur_err = err(cur)
    for _ in range(max_rounds):
        improved = False
        for name, vals, apply in free_axes:
            best_v, best_e = None, cur_err
            for v in vals:
                e = err(apply(cur, float(v)))
                if e < best_e:
                    best_e = e
                    best_v = float(v)
            if best_v is not None:
                cur = apply(cur, best_v)
                cur_err = best_e
                improved = True
        if not improved:
            break
    return cur, cur_err


_BETA_GRID = np.arange(380, 561, 2)
_TRANSPORT_GRID = np.arange(0.40, 0.80, 0.01)
_PINHOLE_GRID = np.arange(18.0, 27.0, 0.25)
_RESISTANCE_GRID = np.arange(0.06, 0.22, 0.005)


def _free_axes_for_mode(mode: str) -> list:
    """
    mode = 'baseline' (all residuals free)
         | 'transport' (transport pinned, others free)
         | 'resistance' (resistance pinned, others free)
         | 'pinhole' (pinhole pinned, others free)
         | 'all' (all three pinned, only N free)
    """
    axes = [
        ("beta", _BETA_GRID,
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
    ]
    if mode in ("baseline", "resistance", "pinhole"):
        axes.append(("transport", _TRANSPORT_GRID,
                     lambda p, v: replace(p, transport=v)))
    if mode in ("baseline", "transport", "pinhole"):
        axes.append(("resistance", _RESISTANCE_GRID,
                     lambda p, v: replace(p, resistance=v)))
    if mode in ("baseline", "transport", "resistance"):
        axes.append(("pinhole", _PINHOLE_GRID,
                     lambda p, v: replace(p, pinhole=v)))
    return axes


def best_N(lock: QuarkParams) -> int:
    return int(round(lock.beta * 2.0 / math.pi))


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    args = p.parse_args()

    scalars = derive_scalars()
    print("Derived scalars:")
    for k, v in scalars.items():
        if k != "alpha_q_table":
            print(f"  {k:>20}: {v:>10.6f}")
    print(f"  alpha_q table: {scalars['alpha_q_table']}")
    print()

    perturbations = {
        "PDG":      _PDG_MASSES,
        "c_x_1.10": {**_PDG_MASSES, "c": _PDG_MASSES["c"] * 1.10},
        "b_x_1.10": {**_PDG_MASSES, "b": _PDG_MASSES["b"] * 1.10},
        "t_x_1.10": {**_PDG_MASSES, "t": _PDG_MASSES["t"] * 1.10},
    }

    def base_for_mode(mode: str) -> QuarkParams:
        b = _TEMPLATE
        if mode in ("transport", "all"):
            b = replace(b, transport=scalars["transport_geom"])
        if mode in ("resistance", "all"):
            b = replace(b, resistance=scalars["resistance_geom"])
        if mode in ("pinhole", "all"):
            b = replace(b, pinhole=scalars["pinhole_geom"])
        return b

    results: dict[str, dict] = {}
    print(f"{'mode':<14}  {'mass set':<10}  {'N':>5}  {'err':>10}")
    print("─" * 50)
    for mode in ("baseline", "transport", "resistance", "pinhole", "all"):
        results[mode] = {}
        for pert_name, masses in perturbations.items():
            err = _err_for_masses(masses)
            base = base_for_mode(mode)
            free = _free_axes_for_mode(mode)
            lock, e = coord_descent(base, err, free)
            results[mode][pert_name] = {
                "N": best_N(lock), "err": e,
                "transport": lock.transport,
                "resistance": lock.resistance,
                "pinhole": lock.pinhole,
            }
            print(f"{mode:<14}  {pert_name:<10}  "
                  f"{best_N(lock):>5d}  {e:>10.4e}")
        ns = [results[mode][p]["N"] for p in perturbations]
        print(f"  -> N range = [{min(ns)}, {max(ns)}], "
              f"width = {max(ns) - min(ns)}")
        print()

    print("=" * 60)
    print("N-stability comparison")
    print("=" * 60)
    print(f"{'substitution mode':<20}  {'N range':>14}  {'width':>6}")
    print("─" * 60)
    for mode in ("baseline", "transport", "resistance", "pinhole", "all"):
        ns = [results[mode][p]["N"] for p in perturbations]
        print(f"{mode:<20}  [{min(ns):>4}, {max(ns):>4}]  {max(ns)-min(ns):>6}")

    out = {"scalars": scalars, "results": results,
           "perturbations": list(perturbations.keys())}
    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nFull results written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
