#!/usr/bin/env python3
"""
scripts/lock_quark_beta_probe.py
==================================

Step 4 of the v3 calibration pipeline: with the integer-winding lock
and action_base fixed, optimize the residual continuous knobs against
the six-species spectrum and emit the final locked parameter point.

Mirrors ``scripts/lock_beta_50pi_probe.py`` in the lepton sector.

Discipline (v3 §5):
    β and action_base are HARD-LOCKED here; only the residual knobs
    (phase, transport, pinhole, resistance, partition_mixing, γ_q)
    are optimized.  If the fit is bad with β locked, do NOT unlock β
    — revisit the spec or report that the integer-winding discovery
    was a false positive.
"""

from __future__ import annotations

# ── Repo-root path shim (r2) ────────────────────────────────────────────
import os
import sys
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
# ────────────────────────────────────────────────────────────────────────

import argparse
import itertools
import json
import math
from dataclasses import asdict
from typing import Optional

import numpy as np

from geometrodynamics.qcd.quark_spectrum import (
    OBSERVED_MASSES_MEV,
    QUARK_ACTION_BASE,
    QUARK_ANCHOR_MASS_MEV,
    QUARK_ANCHOR_SPECIES,
    QUARK_SPECIES,
    QuarkParams,
    extract_physical_spectrum,
)


_TARGET_MASSES_MEV = np.array(
    [OBSERVED_MASSES_MEV[s] for s in QUARK_SPECIES], dtype=float,
)


def _predict(params: QuarkParams) -> Optional[np.ndarray]:
    try:
        species_map = extract_physical_spectrum(params)
    except Exception:
        return None
    anchor_val = species_map[QUARK_ANCHOR_SPECIES]
    if anchor_val <= 1e-6:
        return None
    scale = QUARK_ANCHOR_MASS_MEV / anchor_val
    return np.array([species_map[s] * scale for s in QUARK_SPECIES], dtype=float)


def _residual(params: QuarkParams) -> float:
    """RMS log residual across all six species.  Returns inf if unphysical."""
    predicted = _predict(params)
    if predicted is None:
        return float("inf")
    log_residuals = np.log(predicted / _TARGET_MASSES_MEV)
    return float(np.sqrt(np.mean(log_residuals ** 2)))


def action_base_from_label(label: str) -> float:
    label = label.lower().strip()
    if label == "pi":
        return math.pi
    if label == "2pi":
        return 2.0 * math.pi
    if label == "pi_over_2":
        return math.pi / 2.0
    try:
        return float(label)
    except ValueError as exc:
        raise ValueError(
            f"Unknown action-base label: {label!r}.  Use 'pi', '2pi', "
            f"'pi_over_2', or a numeric value."
        ) from exc


def optimize_residual_knobs(
    integer_winding: int,
    action_base: float,
    n_points: int = 12,
    verbose: bool = False,
) -> dict[str, object]:
    beta = integer_winding * math.pi / 2.0

    grid = {
        "phase":             np.linspace(0.0005, 0.003,  n_points),
        "transport":         np.linspace(0.3,    5.0,    n_points),
        "pinhole":           np.linspace(10.0,   40.0,   n_points),
        "resistance":        np.linspace(0.10,   0.40,   n_points),
        "partition_mixing":  np.linspace(0.0,    1.0,    n_points),
        "gamma_q":           np.linspace(0.01,   0.5,    n_points),
    }
    axis_names = list(grid.keys())
    axis_values = [grid[name] for name in axis_names]
    total = int(np.prod([len(v) for v in axis_values]))

    if verbose:
        print(f"Optimizing {total} grid points with β = {integer_winding}·π/2 "
              f"= {beta:.4f}, action_base = {action_base:.4f}")

    best = {"residual": float("inf"), "params": None, "predicted": None}
    rejected = 0
    counter = 0
    for combo in itertools.product(*axis_values):
        counter += 1
        kwargs = {k: float(v) for k, v in zip(axis_names, combo)}
        params = QuarkParams(
            action_base=action_base,
            beta=beta,
            u_q_form="k_minus_2",
            winding_mode="max",
            resistance_model="exponential",
            depth_cost_mode="tunnel_only",
            **kwargs,
        )
        r = _residual(params)
        if not math.isfinite(r):
            rejected += 1
            continue
        if r < best["residual"]:
            best = {
                "residual": r,
                "params": params,
                "predicted": _predict(params),
            }
            if verbose and counter % max(1, total // 40) == 0:
                print(f"  [{counter}/{total}] residual={r:.4e}")

    if best["params"] is None:
        raise RuntimeError(
            f"optimization found no valid point "
            f"(rejected {rejected} / {total})"
        )

    rel_errs = np.abs(best["predicted"] - _TARGET_MASSES_MEV) / _TARGET_MASSES_MEV

    return {
        "integer_winding": integer_winding,
        "action_base": action_base,
        "locked_params": asdict(best["params"]),
        "residual_log_rms": best["residual"],
        "predicted_masses_mev": best["predicted"].tolist(),
        "target_masses_mev": _TARGET_MASSES_MEV.tolist(),
        "rel_err_per_species": rel_errs.tolist(),
        "max_rel_err": float(np.max(rel_errs)),
        "species_order": list(QUARK_SPECIES),
        "rejected_points": rejected,
        "total_points": total,
    }


def emit_module_lock_snippet(result: dict) -> str:
    lp = result["locked_params"]
    lines = [
        "# Populated by scripts/lock_quark_beta_probe.py",
        f"#   integer_winding = {result['integer_winding']}",
        f"#   action_base     = {result['action_base']}",
        f"#   max_rel_err     = {result['max_rel_err']:.4e}",
        "",
        "LOCKED_QUARK_PARAMS = QuarkParams(",
    ]
    for key, val in lp.items():
        lines.append(f"    {key}={val!r},")
    lines.append(")")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Hard-lock β and action_base, optimize residual knobs.",
    )
    parser.add_argument("--integer-winding", type=int, required=True)
    parser.add_argument("--action-base-label", type=str, default="pi",
                        help="'pi', '2pi', 'pi_over_2', or numeric value")
    parser.add_argument("--n-points", type=int, default=12)
    parser.add_argument("--output-json", type=str, default=None)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    action_base = action_base_from_label(args.action_base_label)

    result = optimize_residual_knobs(
        integer_winding=args.integer_winding,
        action_base=action_base,
        n_points=args.n_points,
        verbose=args.verbose,
    )

    print()
    print("=" * 72)
    print("FINAL LOCK")
    print("=" * 72)
    print(f"integer_winding  = {result['integer_winding']}")
    print(f"action_base      = {result['action_base']:.6f}")
    print(f"residual (log)   = {result['residual_log_rms']:.4e}")
    print(f"max rel error    = {result['max_rel_err']:.4e}")
    print(f"rejected points  = {result['rejected_points']} / {result['total_points']}")
    print()
    print("Predicted (MeV):")
    for species, mass in zip(result["species_order"], result["predicted_masses_mev"]):
        print(f"  {species}: {mass:>14.4f}")
    print("Relative error per species:")
    for species, err in zip(result["species_order"], result["rel_err_per_species"]):
        print(f"  {species}: {err:.4e}")
    print()
    print("=" * 72)
    print("Paste into quark_spectrum.py:")
    print("=" * 72)
    print(emit_module_lock_snippet(result))

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"\nFull lock result written to {args.output_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
