#!/usr/bin/env python3
"""
scripts/experiment_refined_k3k5.py
===================================

Follow-up experiment for docs/quark_axioms.md §8 next-session milestone:

    Refine the edge-limited optimum, add one targeted off-diagonal
    coupling, see whether max_rel_err drops below about 0.3.

Inherits from experiment 3's best point (N=150, ε=+0.9, χ=15,
γ_q=0.1, transport=0.5, pinhole=15, resistance=0.15) and extends:

  1. ε past the 0.9 grid edge (into ε > 1, which pushes (5,−)
     uplift toward zero under min_eigenvalue zero).
  2. χ past the 15 grid edge (up to 140, since χ was also edge-limited).
  3. A new knob ``eta_k3k5_minus``: a targeted off-diagonal amplitude
     that lives only on the (3,−)↔(5,−) element of the Hamiltonian.
     Physical motivation: s and t were the outliers in experiment 3,
     and both sit in the partition-"−" block at k=3 and k=5. Level
     repulsion between them is the minimal structural fix, and a
     single targeted coupling is the minimal way to implement it.
     Default 0.0 recovers experiment 3.

Anchor: d quark (u = 0 by construction under min_eigenvalue zero).
Goal: max_rel_err (excluding u) below 0.3 without broadening
structural assumptions.
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
    QUARK_SPECIES,
    QuarkParams,
    extract_physical_spectrum,
)


_ANCHOR_SPECIES = "d"
_ANCHOR_MASS_MEV = OBSERVED_MASSES_MEV[_ANCHOR_SPECIES]
_TARGET = np.array(
    [OBSERVED_MASSES_MEV[s] for s in QUARK_SPECIES], dtype=float,
)


def _predict(params: QuarkParams) -> Optional[np.ndarray]:
    try:
        species = extract_physical_spectrum(params)
    except Exception:
        return None
    anchor = species[_ANCHOR_SPECIES]
    if anchor <= 1e-6:
        return None
    scale = _ANCHOR_MASS_MEV / anchor
    return np.array([species[s] * scale for s in QUARK_SPECIES], dtype=float)


def _max_rel_err_excl_u(predicted: np.ndarray) -> float:
    non_anchor = np.array(
        [i for i, s in enumerate(QUARK_SPECIES) if s != "u"], dtype=int,
    )
    rel = np.abs(predicted - _TARGET) / _TARGET
    return float(np.max(rel[non_anchor]))


def run(axes: list[tuple[str, np.ndarray]], verbose: bool = False) -> dict[str, object]:
    total = int(np.prod([len(v) for _, v in axes]))
    if verbose:
        print(f"scanning {total} points (anchor = d = {_ANCHOR_MASS_MEV} MeV)")

    best_err = float("inf")
    best_params: Optional[QuarkParams] = None
    best_predicted: Optional[np.ndarray] = None
    rejected = 0
    for combo in itertools.product(*[v for _, v in axes]):
        kwargs = dict(zip([n for n, _ in axes], (float(x) for x in combo)))
        N = kwargs.pop("integer_winding")
        params = QuarkParams(
            action_base=QUARK_ACTION_BASE,
            beta=N * math.pi / 2.0,
            u_q_form="k_minus_2",
            phase=0.001,
            partition_mixing=0.0,
            winding_mode="max",
            resistance_model="exponential",
            depth_cost_mode="tunnel_only",
            uplift_mode="partition_asymmetric",
            spectrum_zero_mode="min_eigenvalue",
            **kwargs,
        )
        predicted = _predict(params)
        if predicted is None:
            rejected += 1
            continue
        err = _max_rel_err_excl_u(predicted)
        if err < best_err:
            best_err = err
            best_params = params
            best_predicted = predicted
            if verbose:
                kv = ", ".join(
                    f"{k}={v:g}" for k, v in kwargs.items()
                    if k in ("uplift_asymmetry", "chi_q_k3",
                             "eta_k3k5_minus", "gamma_q")
                )
                print(f"  N={int(N):>4d}  {kv}  err={err:.4e}")
    if best_params is None:
        raise RuntimeError("no valid point found")

    rel_errs = np.abs(best_predicted - _TARGET) / _TARGET
    return {
        "anchor_species": _ANCHOR_SPECIES,
        "anchor_mass_mev": _ANCHOR_MASS_MEV,
        "best_max_rel_err_excl_u": best_err,
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

    # Pass 1: hold transport/pinhole/resistance at exp-3 best and
    # sweep the primary four axes (N, ε, χ, η) plus γ_q.  The
    # residual knobs can be re-optimized in a follow-up once the
    # primary basin is located.
    axes = [
        ("integer_winding",   np.array([80, 120, 150, 180, 220, 280])),
        ("uplift_asymmetry",  np.array([0.85, 0.95, 1.05, 1.15, 1.30, 1.50, 1.70, 1.90])),
        ("chi_q_k3",          np.array([10.0, 18.0, 28.0, 40.0, 55.0, 75.0, 100.0, 140.0])),
        ("eta_k3k5_minus",    np.array([0.0, 0.5, 1.0, 2.0, 3.5, 5.0, 7.5, 10.0])),
        ("gamma_q",           np.array([0.05, 0.10, 0.20, 0.40])),
        ("transport",         np.array([0.5])),
        ("pinhole",           np.array([15.0])),
        ("resistance",        np.array([0.15])),
    ]

    result = run(axes, verbose=args.verbose)

    print()
    print("=" * 72)
    print(f"best max_rel_err (excl. u): "
          f"{result['best_max_rel_err_excl_u']:.4e}")
    print(f"rejected:           {result['rejected_points']} / "
          f"{result['scan_points_total']}")
    bp = result["best_params"]
    print(f"best N:             {bp['beta'] * 2.0 / math.pi:.1f}")
    print(f"best eps:           {bp['uplift_asymmetry']:+.4f}")
    print(f"best chi:           {bp['chi_q_k3']:.4f}")
    print(f"best eta(3,5)-:     {bp['eta_k3k5_minus']:.4f}")
    print(f"best gamma_q:       {bp['gamma_q']:.4f}")
    print(f"best transport:     {bp['transport']:.4f}")
    print(f"best pinhole:       {bp['pinhole']:.4f}")
    print(f"best resistance:    {bp['resistance']:.4f}")
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
