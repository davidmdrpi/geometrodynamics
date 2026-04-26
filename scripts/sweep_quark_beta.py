#!/usr/bin/env python3
"""
scripts/sweep_quark_beta.py
============================

Step 2 of the v3 calibration pipeline: hunt for the integer-winding
β-uplift lock, parallel to the lepton sector's `4β = 200π = 100·(2π)`
τ-discovery.

Governing principle (v3 spec §0, §5 step 3):
    The heaviest-species uplift must be an integer multiple of 2π.
    For the lepton τ, that integer was 100.  For the heaviest shelled-
    closure species, expect a much larger integer; this script sweeps
    candidate integer windings.

Parallel to ``scripts/sweep_k_uplift_beta.py`` in the lepton sector.
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
import json
import math
from dataclasses import asdict, replace
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


def _build_and_solve(params: QuarkParams) -> Optional[np.ndarray]:
    """Physical spectrum anchored to m_u.  Returns None if unphysical."""
    try:
        species_map = extract_physical_spectrum(params)
    except Exception:
        return None
    anchor_val = species_map[QUARK_ANCHOR_SPECIES]
    if anchor_val <= 1e-6:
        return None
    scale = QUARK_ANCHOR_MASS_MEV / anchor_val
    return np.array([species_map[s] * scale for s in QUARK_SPECIES], dtype=float)


def _relative_error_vector(predicted: np.ndarray) -> np.ndarray:
    return np.abs(predicted - _TARGET_MASSES_MEV) / _TARGET_MASSES_MEV


def integer_to_beta(N: int) -> float:
    """
    β = N · π / 2, so that `4β = N · (2π)`.  The factor 4 is the
    `max(0, k-3)²` coefficient at k=5.  Mirrors the lepton sector's
    `4β_lepton = 100·(2π)` when N=100 there (β_lepton = 50π).
    """
    return N * math.pi / 2.0


def sweep_integers(
    base_params: QuarkParams,
    integer_min: int,
    integer_max: int,
    log_spaced: bool = True,
    n_samples: int = 200,
    verbose: bool = False,
) -> dict[str, object]:
    if log_spaced:
        log_min = math.log(max(integer_min, 1))
        log_max = math.log(integer_max)
        integers = sorted(set(
            int(round(math.exp(x)))
            for x in np.linspace(log_min, log_max, n_samples)
        ))
    else:
        stride = max(1, (integer_max - integer_min) // n_samples)
        integers = list(range(integer_min, integer_max + 1, stride))

    best = {
        "N": None,
        "beta": None,
        "max_rel_err": float("inf"),
        "predicted": None,
        "rel_err_vector": None,
    }
    rejected = 0

    for N in integers:
        beta = integer_to_beta(N)
        params = replace(base_params, beta=beta)
        predicted = _build_and_solve(params)
        if predicted is None:
            rejected += 1
            continue
        rel_err = _relative_error_vector(predicted)
        max_err = float(np.max(rel_err))

        if max_err < best["max_rel_err"]:
            best = {
                "N": N,
                "beta": beta,
                "max_rel_err": max_err,
                "predicted": predicted.tolist(),
                "rel_err_vector": rel_err.tolist(),
            }
            if verbose:
                heavy_err = rel_err[-1]
                print(
                    f"  N={N:>6d}  β={beta:>12.4f}  "
                    f"max_err={max_err:.4e}  heavy_err={heavy_err:.4e}"
                )

    if best["N"] is None:
        raise RuntimeError(
            f"sweep found no valid integer winding "
            f"(rejected {rejected} / {len(integers)})"
        )

    return {
        "best_integer_winding": best["N"],
        "best_beta": best["beta"],
        "best_max_rel_err": best["max_rel_err"],
        "predicted_masses_mev": best["predicted"],
        "rel_err_per_species": best["rel_err_vector"],
        "species_order": list(QUARK_SPECIES),
        "target_masses_mev": _TARGET_MASSES_MEV.tolist(),
        "integer_range": [integer_min, integer_max],
        "samples_tried": len(integers),
        "rejected_points": rejected,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Sweep candidate integer-winding β-uplift values.",
    )
    parser.add_argument("--integer-min", type=int, default=10)
    parser.add_argument("--integer-max", type=int, default=100000)
    parser.add_argument("--n-samples", type=int, default=400)
    parser.add_argument("--linear", action="store_true")
    parser.add_argument("--input-json", type=str, default=None,
                        help="output from calibrate_quark_ratios.py")
    parser.add_argument("--output-json", type=str, default=None)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.input_json:
        with open(args.input_json, "r") as fh:
            step1 = json.load(fh)
        allowed = set(QuarkParams.__dataclass_fields__.keys())
        base_params = QuarkParams(
            **{k: v for k, v in step1["best_params"].items() if k in allowed}
        )
    else:
        base_params = QuarkParams(
            action_base=QUARK_ACTION_BASE,
            gamma_q=0.1,
            transport=1.0,
            partition_mixing=0.1,
        )

    result = sweep_integers(
        base_params=base_params,
        integer_min=args.integer_min,
        integer_max=args.integer_max,
        log_spaced=not args.linear,
        n_samples=args.n_samples,
        verbose=args.verbose,
    )

    print()
    print("=" * 72)
    print(f"Best integer winding: N = {result['best_integer_winding']}")
    print(f"  β = N · π/2 = {result['best_beta']:.6f}")
    print(f"  max rel err: {result['best_max_rel_err']:.4e}")
    print(f"  rejected:    {result['rejected_points']} / {result['samples_tried']}")
    print("=" * 72)
    print("Predicted (MeV):")
    for species, mass in zip(result["species_order"], result["predicted_masses_mev"]):
        print(f"  {species}: {mass:>12.4f}")
    print("Relative error per species:")
    for species, err in zip(result["species_order"], result["rel_err_per_species"]):
        print(f"  {species}: {err:.4e}")

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"\nSweep summary written to {args.output_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
