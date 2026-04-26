#!/usr/bin/env python3
"""
scripts/map_basin_quark_uplift.py
===================================

Step 3 of the v3 calibration pipeline: probe the basin around the
integer-winding β-lock from step 2 to verify that the lock sits on a
genuine attractor, not a narrow needle.

Mirrors ``scripts/map_basin_k_uplift.py`` in the lepton sector.

v3 spec §5 step 5 criterion: topological locks should be robust —
small deviations from the integer winding should degrade the fit
smoothly.  If the basin is a needle, the integer-winding claim is
suspect and the spec needs revisiting.
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
from dataclasses import replace
from typing import Optional

import numpy as np

from geometrodynamics.qcd.quark_spectrum import (
    OBSERVED_MASSES_MEV,
    QUARK_ANCHOR_MASS_MEV,
    QUARK_ANCHOR_SPECIES,
    QUARK_SPECIES,
    QuarkParams,
    extract_physical_spectrum,
)


_TARGET_MASSES_MEV = np.array(
    [OBSERVED_MASSES_MEV[s] for s in QUARK_SPECIES], dtype=float,
)


def _max_rel_err(params: QuarkParams) -> float:
    try:
        species_map = extract_physical_spectrum(params)
    except Exception:
        return float("nan")
    anchor_val = species_map[QUARK_ANCHOR_SPECIES]
    if anchor_val <= 1e-6:
        return float("nan")
    scale = QUARK_ANCHOR_MASS_MEV / anchor_val
    predicted = np.array([species_map[s] * scale for s in QUARK_SPECIES], dtype=float)
    rel_err = np.abs(predicted - _TARGET_MASSES_MEV) / _TARGET_MASSES_MEV
    return float(np.max(rel_err))


def probe_basin(
    base_params: QuarkParams,
    center_integer: int,
    half_width: int,
    n_points: int,
) -> dict[str, object]:
    integers_raw = np.linspace(
        center_integer - half_width,
        center_integer + half_width,
        n_points,
    )
    integers = sorted({int(round(x)) for x in integers_raw if round(x) >= 1})

    errors = []
    for N in integers:
        beta = N * math.pi / 2.0
        params = replace(base_params, beta=beta)
        errors.append(_max_rel_err(params))

    errors_arr = np.asarray(errors, dtype=float)
    valid = ~np.isnan(errors_arr)
    if not np.any(valid):
        raise RuntimeError("no valid points in basin scan")

    valid_integers = [n for n, v in zip(integers, valid) if v]
    valid_errors = errors_arr[valid]
    best_idx = int(np.argmin(valid_errors))
    best_N = valid_integers[best_idx]
    best_err = float(valid_errors[best_idx])

    # Basin width within 2× of best error
    threshold = 2.0 * best_err
    width_lo = width_hi = center_integer
    if center_integer in integers:
        center_idx = integers.index(center_integer)
        if not math.isnan(errors_arr[center_idx]):
            i = center_idx
            while i > 0 and not math.isnan(errors_arr[i - 1]) and errors_arr[i - 1] <= threshold:
                i -= 1
                width_lo = integers[i]
            i = center_idx
            while i < len(integers) - 1 and not math.isnan(errors_arr[i + 1]) and errors_arr[i + 1] <= threshold:
                i += 1
                width_hi = integers[i]

    return {
        "center_integer": center_integer,
        "best_integer": best_N,
        "best_max_rel_err": best_err,
        "basin_within_2x_of_best": [width_lo, width_hi],
        "integers_probed": integers,
        "errors_per_integer": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Basin-width probe around the shelled-closure β-lock.",
    )
    parser.add_argument("--integer-winding", type=int, required=True)
    parser.add_argument("--half-width", type=int, default=20)
    parser.add_argument("--n-points", type=int, default=41)
    parser.add_argument("--input-json", type=str, default=None)
    parser.add_argument("--output-json", type=str, default=None)
    args = parser.parse_args()

    # If a step-1 JSON is given, inherit its parameter point; otherwise
    # use defensive defaults that keep the anchor physical.
    if args.input_json:
        with open(args.input_json, "r") as fh:
            data = json.load(fh)
        allowed = set(QuarkParams.__dataclass_fields__.keys())
        if "best_params" in data:
            base_params = QuarkParams(
                **{k: v for k, v in data["best_params"].items() if k in allowed}
            )
        else:
            base_params = QuarkParams(gamma_q=0.1, transport=1.0, partition_mixing=0.1)
    else:
        base_params = QuarkParams(gamma_q=0.1, transport=1.0, partition_mixing=0.1)

    result = probe_basin(
        base_params=base_params,
        center_integer=args.integer_winding,
        half_width=args.half_width,
        n_points=args.n_points,
    )

    print()
    print("=" * 72)
    print(f"Basin probe around N = {result['center_integer']}")
    print("=" * 72)
    print(f"best N in window:       {result['best_integer']}")
    print(f"best max rel err:       {result['best_max_rel_err']:.4e}")
    lo, hi = result["basin_within_2x_of_best"]
    print(f"basin (≤ 2× best err):  [{lo}, {hi}]  (width = {hi - lo})")

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"\nBasin probe written to {args.output_json}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
