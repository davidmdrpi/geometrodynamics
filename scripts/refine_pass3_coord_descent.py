#!/usr/bin/env python3
"""
scripts/refine_pass3_coord_descent.py
======================================

Pass 3: coordinate-descent refinement starting from the pass-2 best
point.  Each axis is swept on a fine 1D grid with the others held
fixed at the running best; an axis update fires only if it improves
the max-rel-err.  The cycle repeats until no axis improves or the
configured iteration cap is hit.

Faster than a full 8-D fine grid — needs only O(N_axes × resolution)
evaluations per round, vs O(resolution^N_axes) for a grid scan.
Trades exhaustive coverage for local convergence; the basin probes
already established that the pass-2 point sits in a real basin, so
local descent is now the right tool.

Goal (per user milestone): max_rel_err (excl. u) below 0.10.
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
from typing import Callable, Optional

import numpy as np

from geometrodynamics.qcd.quark_spectrum import (
    OBSERVED_MASSES_MEV,
    QUARK_ACTION_BASE,
    QUARK_SPECIES,
    QuarkParams,
    extract_physical_spectrum,
)


_ANCHOR_SPECIES = "d"
_ANCHOR_MASS_MEV = OBSERVED_MASSES_MEV[_ANCHOR_SPECIES]
_TARGET = np.array(
    [OBSERVED_MASSES_MEV[s] for s in QUARK_SPECIES], dtype=float,
)


def _max_rel_err_excl_u(params: QuarkParams) -> float:
    try:
        species = extract_physical_spectrum(params)
    except Exception:
        return float("inf")
    anchor = species[_ANCHOR_SPECIES]
    if anchor <= 1e-6:
        return float("inf")
    scale = _ANCHOR_MASS_MEV / anchor
    predicted = np.array([species[s] * scale for s in QUARK_SPECIES], dtype=float)
    non_anchor = np.array(
        [i for i, s in enumerate(QUARK_SPECIES) if s != "u"], dtype=int,
    )
    rel = np.abs(predicted - _TARGET) / _TARGET
    return float(np.max(rel[non_anchor]))


def _predict(params: QuarkParams) -> np.ndarray:
    species = extract_physical_spectrum(params)
    scale = _ANCHOR_MASS_MEV / species[_ANCHOR_SPECIES]
    return np.array([species[s] * scale for s in QUARK_SPECIES], dtype=float)


def coord_descent(
    base: QuarkParams,
    axes: list[tuple[str, np.ndarray, Callable[[QuarkParams, float], QuarkParams]]],
    max_rounds: int = 6,
    verbose: bool = False,
) -> tuple[QuarkParams, float, list[dict]]:
    current = base
    current_err = _max_rel_err_excl_u(current)
    history: list[dict] = []
    if verbose:
        print(f"start err = {current_err:.4e}")

    for round_idx in range(max_rounds):
        improved = False
        for axis_name, values, apply in axes:
            best_v = None
            best_e = current_err
            for v in values:
                trial = apply(current, float(v))
                e = _max_rel_err_excl_u(trial)
                if e < best_e:
                    best_e = e
                    best_v = float(v)
            if best_v is not None:
                current = apply(current, best_v)
                history.append({
                    "round": round_idx,
                    "axis": axis_name,
                    "new_value": best_v,
                    "new_err": best_e,
                })
                if verbose:
                    print(f"  round {round_idx}  {axis_name:>20} -> "
                          f"{best_v:.4g}, err = {best_e:.4e}")
                current_err = best_e
                improved = True
        if not improved:
            if verbose:
                print(f"converged after {round_idx} round(s)")
            break

    return current, current_err, history


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--max-rounds", type=int, default=6)
    args = p.parse_args()

    base = QuarkParams(
        action_base=QUARK_ACTION_BASE,
        beta=400 * math.pi / 2.0,
        gamma_q=0.10,
        u_q_form="k_minus_2",
        phase=0.001,
        transport=0.6,
        pinhole=22.0,
        resistance=0.15,
        partition_mixing=0.0,
        winding_mode="max",
        resistance_model="exponential",
        depth_cost_mode="tunnel_only",
        uplift_mode="partition_asymmetric",
        uplift_asymmetry=0.95,
        spectrum_zero_mode="min_eigenvalue",
        chi_q_k3=20.0,
        eta_k3k5_minus=5.0,
    )

    # 1D grids tuned to the basin widths from basin_probe_v2:
    #   N width ~100, chi width ~1.6, eta width ~9.5
    axes = [
        # Integer N: descend with step 5 first, refined to step 1 in
        # later rounds; we just give a single fine grid and let
        # descent visit it once per round.
        ("integer_winding", np.arange(360, 461, 4),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("uplift_asymmetry", np.arange(0.80, 1.20, 0.02),
            lambda p, v: replace(p, uplift_asymmetry=v)),
        ("chi_q_k3", np.arange(18.0, 22.0, 0.1),
            lambda p, v: replace(p, chi_q_k3=v)),
        ("eta_k3k5_minus", np.arange(0.0, 12.0, 0.25),
            lambda p, v: replace(p, eta_k3k5_minus=v)),
        ("gamma_q", np.arange(0.04, 0.20, 0.01),
            lambda p, v: replace(p, gamma_q=v)),
        ("transport", np.arange(0.30, 1.20, 0.05),
            lambda p, v: replace(p, transport=v)),
        ("pinhole", np.arange(15.0, 30.0, 0.5),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.08, 0.30, 0.01),
            lambda p, v: replace(p, resistance=v)),
        ("phase", np.arange(0.0001, 0.005, 0.0002),
            lambda p, v: replace(p, phase=v)),
    ]

    final, final_err, history = coord_descent(
        base, axes, max_rounds=args.max_rounds, verbose=args.verbose,
    )

    pred = _predict(final)
    rel_errs = np.abs(pred - _TARGET) / _TARGET

    print()
    print("=" * 72)
    print(f"final max_rel_err (excl. u): {final_err:.4e}")
    print("=" * 72)
    print(f"N (= 2*beta/pi):     {final.beta * 2.0 / math.pi:.2f}")
    print(f"uplift_asymmetry:    {final.uplift_asymmetry:.4f}")
    print(f"chi_q_k3:            {final.chi_q_k3:.4f}")
    print(f"eta_k3k5_minus:      {final.eta_k3k5_minus:.4f}")
    print(f"gamma_q:             {final.gamma_q:.4f}")
    print(f"transport:           {final.transport:.4f}")
    print(f"pinhole:             {final.pinhole:.4f}")
    print(f"resistance:          {final.resistance:.4f}")
    print(f"phase:               {final.phase:.6f}")
    print()
    print("Predicted vs observed (MeV):")
    for sp, p_, o, e in zip(QUARK_SPECIES, pred, _TARGET, rel_errs):
        print(f"  {sp}: pred={p_:>14.4f}  obs={o:>10.4f}  rel_err={e:.4e}")

    if args.output_json:
        out = {
            "final_max_rel_err_excl_u": final_err,
            "final_params": asdict(final),
            "predicted_masses_mev": pred.tolist(),
            "target_masses_mev": _TARGET.tolist(),
            "rel_err_per_species": rel_errs.tolist(),
            "species_order": list(QUARK_SPECIES),
            "history": history,
        }
        with open(args.output_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nFull result written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
