#!/usr/bin/env python3
"""
scripts/basin_probe_topological_locks.py
==========================================

Credibility test for the pass-2 best point: are N=400, χ=20, η=5
basin features or grid coincidences?

For each of the three topological-lock candidates (N, χ_q_k3,
η_k3k5_minus), this script holds the other seven axes at the pass-2
best point and sweeps the target axis across a fine grid wider than
the pass-2 grid spacing.  A clean basin shows:

  * Error rises smoothly on both sides of the central value.
  * The 2×-best-error window has nontrivial width (not a needle).
  * The minimum within the swept window matches (or improves on)
    the pass-2 grid value.

A grid coincidence would show:

  * The pass-2 value sits on a slope (better off-axis points exist
    very close by).
  * The basin is asymmetric / pathological, or the 2×-best window
    collapses to a single sample.

Output: per-axis CSV-style tables to stdout plus a JSON summary.
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
from typing import Optional

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


# Pass-2 best point (from docs/calibration_runs/experiment_refined_k3k5_pass2.json)
_BEST = QuarkParams(
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


def _max_rel_err_excl_u(params: QuarkParams) -> float:
    try:
        species = extract_physical_spectrum(params)
    except Exception:
        return float("nan")
    anchor = species[_ANCHOR_SPECIES]
    if anchor <= 1e-6:
        return float("nan")
    scale = _ANCHOR_MASS_MEV / anchor
    predicted = np.array([species[s] * scale for s in QUARK_SPECIES], dtype=float)
    non_anchor = np.array(
        [i for i, s in enumerate(QUARK_SPECIES) if s != "u"], dtype=int,
    )
    rel = np.abs(predicted - _TARGET) / _TARGET
    return float(np.max(rel[non_anchor]))


def probe_axis(
    axis_name: str,
    values: np.ndarray,
    apply: callable,
    label_fmt: str = "{:.4g}",
    verbose: bool = True,
) -> dict[str, object]:
    errs = []
    for v in values:
        params = apply(_BEST, float(v))
        errs.append(_max_rel_err_excl_u(params))
    errs_arr = np.asarray(errs, dtype=float)
    valid = ~np.isnan(errs_arr)
    if not np.any(valid):
        raise RuntimeError(f"no valid points along {axis_name}")
    valid_vals = values[valid]
    valid_errs = errs_arr[valid]
    best_idx = int(np.argmin(valid_errs))
    best_v = float(valid_vals[best_idx])
    best_err = float(valid_errs[best_idx])
    threshold = 2.0 * best_err

    # Find contiguous window around best below threshold
    window_lo = window_hi = best_v
    sorted_pairs = sorted(zip(valid_vals.tolist(), valid_errs.tolist()))
    sorted_vs = [v for v, _ in sorted_pairs]
    sorted_es = [e for _, e in sorted_pairs]
    bi = sorted_vs.index(best_v)
    i = bi
    while i > 0 and sorted_es[i - 1] <= threshold:
        i -= 1
        window_lo = sorted_vs[i]
    i = bi
    while i < len(sorted_vs) - 1 and sorted_es[i + 1] <= threshold:
        i += 1
        window_hi = sorted_vs[i]

    if verbose:
        print(f"\n── basin probe along {axis_name} ──")
        print(f"{'value':>14}  {'max_rel_err':>12}")
        for v, e in zip(valid_vals.tolist(), valid_errs.tolist()):
            mark = "  <- best" if v == best_v else ""
            print(f"{label_fmt.format(v):>14}  {e:>12.4e}{mark}")
        print(f"\n  best at {axis_name} = {label_fmt.format(best_v)}, "
              f"err = {best_err:.4e}")
        print(f"  2×best window: [{label_fmt.format(window_lo)}, "
              f"{label_fmt.format(window_hi)}]")

    return {
        "axis": axis_name,
        "values": [float(v) for v in valid_vals],
        "errors": [float(e) for e in valid_errs],
        "best_value": best_v,
        "best_error": best_err,
        "two_x_best_window": [float(window_lo), float(window_hi)],
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    args = p.parse_args()

    print("Pass-2 best point baseline:")
    base_err = _max_rel_err_excl_u(_BEST)
    print(f"  max_rel_err (excl. u) = {base_err:.4e}")

    # ── N (integer winding) ─────────────────────────────────────────
    # Sweep N ∈ [200, 800] in steps of 25; pass-2 was N=400.
    n_values = np.arange(200, 825, 25)
    n_result = probe_axis(
        "integer_winding",
        n_values,
        apply=lambda p, v: replace(p, beta=v * math.pi / 2.0),
        label_fmt="{:.0f}",
    )

    # ── χ (chi_q_k3) ────────────────────────────────────────────────
    # Sweep χ ∈ [10, 22.4] in steps of 0.2.  At χ ≈ 22.5 the d
    # eigenvalue crosses zero and the d-anchor scaling breaks down;
    # the smooth basin lives below that level-crossing boundary.
    chi_values = np.arange(10.0, 22.45, 0.2)
    chi_result = probe_axis(
        "chi_q_k3",
        chi_values,
        apply=lambda p, v: replace(p, chi_q_k3=v),
        label_fmt="{:.2f}",
    )

    # ── η (eta_k3k5_minus) ─────────────────────────────────────────
    # Sweep η ∈ [0, 18] in steps of 0.5; pass-2 was η=5.  η ≥ 18
    # crosses a similar level-crossing boundary.
    eta_values = np.arange(0.0, 18.05, 0.5)
    eta_result = probe_axis(
        "eta_k3k5_minus",
        eta_values,
        apply=lambda p, v: replace(p, eta_k3k5_minus=v),
        label_fmt="{:.2f}",
    )

    summary = {
        "pass2_best_params": asdict(_BEST),
        "pass2_baseline_err": base_err,
        "probes": {
            "integer_winding": n_result,
            "chi_q_k3":        chi_result,
            "eta_k3k5_minus":  eta_result,
        },
    }

    print("\n" + "=" * 72)
    print("BASIN PROBE SUMMARY")
    print("=" * 72)
    print(f"{'axis':>20}  {'best':>10}  {'best_err':>10}  {'2x window':>14}  {'width':>8}")
    for axis_key in ("integer_winding", "chi_q_k3", "eta_k3k5_minus"):
        r = summary["probes"][axis_key]
        lo, hi = r["two_x_best_window"]
        print(f"{axis_key:>20}  {r['best_value']:>10.4g}  "
              f"{r['best_error']:>10.4e}  "
              f"[{lo:>5.4g}, {hi:>5.4g}]  {hi - lo:>8.4g}")

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(summary, fh, indent=2)
        print(f"\nFull summary written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
