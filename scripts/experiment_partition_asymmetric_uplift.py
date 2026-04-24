#!/usr/bin/env python3
"""
scripts/experiment_partition_asymmetric_uplift.py
==================================================

Follow-up experiment for the negative calibration result logged in
docs/quark_axioms.md §8.  Tests candidate #2 from the §8 "Verdict":

    Does a partition-dependent k=5 uplift break the b/t degeneracy?

Extends the minimal v3 ansatz with
    uplift = beta * (k-3)^2 * (1 + epsilon * sigma(p))
so that at k=5 the partition-"+" state (t) receives 4*beta*(1+eps)
and partition-"-" (b) receives 4*beta*(1-eps).  The minimal ansatz is
recovered at epsilon = 0.

Runs a coarse grid over (epsilon, N) plus the residual knobs and
reports whether the best max_rel_err materially improves over the
§8 baseline.

Not part of the locked pipeline.  Output goes to stdout + JSON so
next session can decide whether this becomes a spec revision or
gets abandoned.
"""

from __future__ import annotations

import os
import sys
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

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


_TARGET = np.array(
    [OBSERVED_MASSES_MEV[s] for s in QUARK_SPECIES], dtype=float,
)


def _predict(params: QuarkParams) -> Optional[np.ndarray]:
    try:
        species = extract_physical_spectrum(params)
    except Exception:
        return None
    anchor = species[QUARK_ANCHOR_SPECIES]
    if anchor <= 1e-6:
        return None
    scale = QUARK_ANCHOR_MASS_MEV / anchor
    return np.array([species[s] * scale for s in QUARK_SPECIES], dtype=float)


def _max_rel_err(predicted: np.ndarray) -> float:
    return float(np.max(np.abs(predicted - _TARGET) / _TARGET))


def run(
    epsilons: np.ndarray,
    integer_windings: np.ndarray,
    gamma_qs: np.ndarray,
    transports: np.ndarray,
    resistances: np.ndarray,
    pinholes: np.ndarray,
    verbose: bool = False,
) -> dict[str, object]:
    axes = [
        ("uplift_asymmetry", epsilons),
        ("integer_winding",  integer_windings),
        ("gamma_q",          gamma_qs),
        ("transport",        transports),
        ("resistance",       resistances),
        ("pinhole",          pinholes),
    ]
    total = int(np.prod([len(v) for _, v in axes]))
    if verbose:
        print(f"scanning {total} points "
              f"(eps x N x gamma_q x transport x resistance x pinhole)")

    best_err = float("inf")
    best_params: Optional[QuarkParams] = None
    best_predicted: Optional[np.ndarray] = None
    rejected = 0
    for combo in itertools.product(*[v for _, v in axes]):
        eps, N, gq, tr, rs, pn = (float(x) for x in combo)
        params = QuarkParams(
            action_base=QUARK_ACTION_BASE,
            beta=N * math.pi / 2.0,
            gamma_q=gq,
            u_q_form="k_minus_2",
            phase=0.001,
            transport=tr,
            pinhole=pn,
            resistance=rs,
            partition_mixing=0.0,
            winding_mode="max",
            resistance_model="exponential",
            depth_cost_mode="tunnel_only",
            uplift_mode="partition_asymmetric",
            uplift_asymmetry=eps,
        )
        predicted = _predict(params)
        if predicted is None:
            rejected += 1
            continue
        err = _max_rel_err(predicted)
        if err < best_err:
            best_err = err
            best_params = params
            best_predicted = predicted
            if verbose:
                print(f"  eps={eps:+.2f}  N={int(N):>5d}  "
                      f"gamma_q={gq:.3f}  err={err:.4e}")
    if best_params is None:
        raise RuntimeError("no valid point found")

    rel_errs = np.abs(best_predicted - _TARGET) / _TARGET
    return {
        "best_max_rel_err": best_err,
        "best_params": asdict(best_params),
        "predicted_masses_mev": best_predicted.tolist(),
        "target_masses_mev": _TARGET.tolist(),
        "rel_err_per_species": rel_errs.tolist(),
        "species_order": list(QUARK_SPECIES),
        "scan_points_total": total,
        "rejected_points": rejected,
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    # Deliberately compact grid for tractability; the scan is
    # exploratory, not final.  ~1.3k points per epsilon slice.
    result = run(
        epsilons=np.linspace(-0.9, 0.9, 10),
        integer_windings=np.array([10, 30, 60, 100, 150, 250, 500, 1000]),
        gamma_qs=np.linspace(0.01, 0.5, 6),
        transports=np.array([0.5, 1.5, 3.0]),
        resistances=np.array([0.15, 0.22, 0.30]),
        pinholes=np.array([15.0, 22.5, 30.0]),
        verbose=args.verbose,
    )

    print()
    print("=" * 72)
    print(f"best max_rel_err:   {result['best_max_rel_err']:.4e}")
    print(f"rejected:           {result['rejected_points']} / "
          f"{result['scan_points_total']}")
    print(f"best uplift_asymm:  "
          f"{result['best_params']['uplift_asymmetry']:+.4f}")
    print(f"best beta:          {result['best_params']['beta']:.4f}  "
          f"(N = {result['best_params']['beta'] * 2.0 / math.pi:.1f})")
    print(f"best gamma_q:       {result['best_params']['gamma_q']:.4f}")
    print()
    print("Predicted vs observed (MeV):")
    for species, pred, obs, err in zip(
        result["species_order"],
        result["predicted_masses_mev"],
        result["target_masses_mev"],
        result["rel_err_per_species"],
    ):
        print(f"  {species}: pred={pred:>14.4f}  obs={obs:>10.4f}  "
              f"rel_err={err:.4e}")

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"\nResult written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
