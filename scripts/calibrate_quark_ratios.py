#!/usr/bin/env python3
"""
scripts/calibrate_quark_ratios.py
==================================

Step 1 of the v3 calibration pipeline: coarse grid scan on the residual
continuous knobs of the shelled-closure Hamiltonian.

For each point on the grid, builds the 6×6 Hamiltonian, extracts the
physical spectrum (with adiabatic species labeling and spectrum-zero
shift), and reports the relative error against the observed masses
after anchoring the lightest species to its observed mass.

Parallel to ``scripts/calibrate_muon_ratio.py`` in the lepton sector.

Usage
─────
    python scripts/calibrate_quark_ratios.py --n-points 16
    python scripts/calibrate_quark_ratios.py --n-points 8 --verbose

Output
──────
Writes best-fit parameter point to stdout, optionally to JSON.  Does
NOT write LOCKED_QUARK_PARAMS back into the module — that happens only
after ``sweep_quark_beta.py`` and ``lock_quark_beta_probe.py``.
"""

from __future__ import annotations

# ── Repo-root path shim (r2) ────────────────────────────────────────────
# Mirrors the convention used by the lepton calibration CLIs so this
# script can be run without PYTHONPATH wrangling.
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


# ════════════════════════════════════════════════════════════════════════
# SCAN CONFIGURATION
# ════════════════════════════════════════════════════════════════════════

_TARGET_MASSES_MEV: np.ndarray = np.array(
    [OBSERVED_MASSES_MEV[s] for s in QUARK_SPECIES], dtype=float,
)
_TARGET_RATIOS_TO_ANCHOR: np.ndarray = (
    _TARGET_MASSES_MEV / OBSERVED_MASSES_MEV[QUARK_ANCHOR_SPECIES]
)


def _ratio_residual(params: QuarkParams) -> tuple[float, Optional[np.ndarray], Optional[str]]:
    """
    Extract physical spectrum, compute relative error on mass ratios
    against observed values.

    Returns (max_rel_err, predicted_masses_mev, reject_reason).
    If the parameter point is unphysical (anchor ≤ 0 after shift, or
    adiabatic tracking fails), returns (inf, None, reason).
    """
    try:
        species_map = extract_physical_spectrum(params)
    except Exception as exc:
        return float("inf"), None, f"extract failed: {exc}"

    anchor_val = species_map[QUARK_ANCHOR_SPECIES]
    if anchor_val <= 1e-6:
        return float("inf"), None, f"anchor non-positive: {anchor_val:.4e}"

    scale = QUARK_ANCHOR_MASS_MEV / anchor_val
    predicted = np.array(
        [species_map[s] * scale for s in QUARK_SPECIES], dtype=float,
    )
    predicted_ratios = predicted / QUARK_ANCHOR_MASS_MEV
    rel_err = np.abs(predicted_ratios - _TARGET_RATIOS_TO_ANCHOR) / _TARGET_RATIOS_TO_ANCHOR
    return float(np.max(rel_err)), predicted, None


def _default_grid(n_points: int) -> dict[str, np.ndarray]:
    """
    5-axis coarse grid.  Ranges seeded from the lepton locked baseline
    and widened modestly.  If the best point lands on an edge, widen
    that axis in a follow-up run.
    """
    return {
        "phase":             np.linspace(0.0005, 0.002, n_points),
        "transport":         np.linspace(0.5,    5.0,   n_points),
        "pinhole":           np.linspace(15.0,   30.0,  n_points),
        "resistance":        np.linspace(0.15,   0.30,  n_points),
        "partition_mixing":  np.linspace(0.0,    0.5,   n_points),
    }


_STRUCTURAL_LOCKS = dict(
    action_base=QUARK_ACTION_BASE,
    beta=0.0,
    u_q_form="k_minus_2",
    winding_mode="max",
    resistance_model="exponential",
    depth_cost_mode="tunnel_only",
)


def run_scan(
    n_points: int = 8,
    gamma_q_values: Optional[list[float]] = None,
    verbose: bool = False,
) -> dict[str, object]:
    if gamma_q_values is None:
        gamma_q_values = [0.05, 0.1, 0.2, 0.5]

    grid = _default_grid(n_points)
    axis_names = list(grid.keys())
    axis_values = [grid[name] for name in axis_names]

    total = int(np.prod([len(v) for v in axis_values])) * len(gamma_q_values)
    if verbose:
        print(f"Scanning {total} points ...")

    best_err = float("inf")
    best_params: Optional[QuarkParams] = None
    best_predicted: Optional[np.ndarray] = None
    rejected = 0

    counter = 0
    for gamma_q in gamma_q_values:
        for combo in itertools.product(*axis_values):
            counter += 1
            kwargs = dict(zip(axis_names, combo))
            params = QuarkParams(
                action_base=_STRUCTURAL_LOCKS["action_base"],
                beta=_STRUCTURAL_LOCKS["beta"],
                gamma_q=float(gamma_q),
                u_q_form=_STRUCTURAL_LOCKS["u_q_form"],
                winding_mode=_STRUCTURAL_LOCKS["winding_mode"],
                resistance_model=_STRUCTURAL_LOCKS["resistance_model"],
                depth_cost_mode=_STRUCTURAL_LOCKS["depth_cost_mode"],
                **{k: float(v) for k, v in kwargs.items()},
            )
            rel_err, predicted, reject_reason = _ratio_residual(params)
            if reject_reason is not None:
                rejected += 1
                continue

            if rel_err < best_err:
                best_err = rel_err
                best_params = params
                best_predicted = predicted
                if verbose:
                    print(
                        f"  [{counter}/{total}] new best: "
                        f"max_rel_err={rel_err:.4f}, γ_q={gamma_q:.3f}"
                    )

    if best_params is None:
        raise RuntimeError(
            f"scan found no valid points out of {total} "
            f"(rejected {rejected}).  Try widening the grid ranges "
            f"or reducing γ_q to keep the anchor above action_base."
        )

    return {
        "best_max_rel_err": best_err,
        "best_params": asdict(best_params),
        "best_predicted_masses_mev": (
            best_predicted.tolist() if best_predicted is not None else None
        ),
        "target_masses_mev": _TARGET_MASSES_MEV.tolist(),
        "species_order": list(QUARK_SPECIES),
        "scan_points_total": counter,
        "rejected_points": rejected,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Coarse grid scan on shelled-closure residual knobs.",
    )
    parser.add_argument("--n-points", type=int, default=8,
                        help="grid density per axis (default: 8)")
    parser.add_argument("--output-json", type=str, default=None,
                        help="optional path to write the full scan summary")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    result = run_scan(n_points=args.n_points, verbose=args.verbose)

    print()
    print("=" * 72)
    print("Best point")
    print("=" * 72)
    print(f"max rel err vs observed ratios: {result['best_max_rel_err']:.4e}")
    print(f"scan points:                    {result['scan_points_total']}")
    print(f"rejected (unphysical):          {result['rejected_points']}")
    print()
    print("Predicted (MeV):")
    for species, mass in zip(result["species_order"], result["best_predicted_masses_mev"]):
        print(f"  {species}: {mass:>12.4f}")
    print("Observed (MeV):")
    for species, mass in zip(result["species_order"], result["target_masses_mev"]):
        print(f"  {species}: {mass:>12.4f}")
    print()
    print("Best parameters:")
    for k, v in result["best_params"].items():
        if isinstance(v, float):
            print(f"  {k:>20s}: {v:.6f}")
        else:
            print(f"  {k:>20s}: {v}")

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"\nFull result written to {args.output_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
