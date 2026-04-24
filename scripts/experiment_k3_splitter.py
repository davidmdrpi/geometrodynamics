#!/usr/bin/env python3
"""
scripts/experiment_k3_splitter.py
==================================

Plan step 4 from the user's 2026-04-24 session recommendations.
Takes the experiment_min_eigenvalue_zero.py setup (min_eigenvalue
zero + partition_asymmetric uplift + d-anchor) and adds one more
knob: chi_q_k3, a k=3-specific partition split that targets the
c/s degeneracy surfaced by the previous scan.

Structure of the new term (see _diagonal_entry):
    k3_split = chi_q_k3 * sigma(p) if k == 3 else 0

So at k=3, "+" (c) gets +chi_q_k3 and "-" (s) gets -chi_q_k3,
leaving k=1 and k=5 alone.  Zero chi_q_k3 recovers the previous
experiment exactly.

Runs a 7-axis coarse grid.  Writes JSON for next-session review.
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
_NON_ANCHOR_IDX = np.array(
    [i for i, s in enumerate(QUARK_SPECIES) if s != "u"], dtype=int,
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


def _max_rel_err(predicted: np.ndarray) -> float:
    rel = np.abs(predicted - _TARGET) / _TARGET
    return float(np.max(rel[_NON_ANCHOR_IDX]))


def run(
    integer_windings: np.ndarray,
    uplift_asymmetries: np.ndarray,
    gamma_qs: np.ndarray,
    chi_k3s: np.ndarray,
    transports: np.ndarray,
    pinholes: np.ndarray,
    resistances: np.ndarray,
    verbose: bool = False,
) -> dict[str, object]:
    axes = [
        ("integer_winding",   integer_windings),
        ("uplift_asymmetry",  uplift_asymmetries),
        ("gamma_q",           gamma_qs),
        ("chi_q_k3",          chi_k3s),
        ("transport",         transports),
        ("pinhole",           pinholes),
        ("resistance",        resistances),
    ]
    total = int(np.prod([len(v) for _, v in axes]))
    if verbose:
        print(f"scanning {total} points (anchor = d = "
              f"{_ANCHOR_MASS_MEV} MeV)")

    best_err = float("inf")
    best_params: Optional[QuarkParams] = None
    best_predicted: Optional[np.ndarray] = None
    rejected = 0
    for combo in itertools.product(*[v for _, v in axes]):
        N, eps, gq, chi, tr, pn, rs = (float(x) for x in combo)
        params = QuarkParams(
            action_base=QUARK_ACTION_BASE,
            beta=N * math.pi / 2.0,
            gamma_q=gq,
            chi_q_k3=chi,
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
            spectrum_zero_mode="min_eigenvalue",
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
                print(f"  N={int(N):>5d}  eps={eps:+.2f}  "
                      f"gq={gq:.2f}  chi={chi:.2f}  err={err:.4e}")
    if best_params is None:
        raise RuntimeError("no valid point found")

    rel_errs = np.abs(best_predicted - _TARGET) / _TARGET
    return {
        "anchor_species": _ANCHOR_SPECIES,
        "anchor_mass_mev": _ANCHOR_MASS_MEV,
        "best_max_rel_err_excluding_u": best_err,
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

    # Compact 7-axis grid.  chi_q_k3 explores 0..30 since the target
    # c/s = 13.6 requires the k=3 split to be comparable to the
    # k_cost + pinhole at k=3.
    result = run(
        integer_windings=np.array([60, 150, 500, 1500, 5000]),
        uplift_asymmetries=np.linspace(-0.8, 0.9, 6),
        gamma_qs=np.array([0.1, 0.5, 1.5, 3.0]),
        chi_k3s=np.array([0.0, 1.0, 3.0, 8.0, 15.0, 30.0]),
        transports=np.array([0.5, 2.0]),
        pinholes=np.array([15.0, 25.0]),
        resistances=np.array([0.15, 0.25]),
        verbose=args.verbose,
    )

    print()
    print("=" * 72)
    print(f"best max_rel_err (excl. u): "
          f"{result['best_max_rel_err_excluding_u']:.4e}")
    print(f"rejected:           {result['rejected_points']} / "
          f"{result['scan_points_total']}")
    bp = result["best_params"]
    print(f"best N:             {bp['beta'] * 2.0 / math.pi:.1f}")
    print(f"best eps:           {bp['uplift_asymmetry']:+.4f}")
    print(f"best gamma_q:       {bp['gamma_q']:.4f}")
    print(f"best chi_q_k3:      {bp['chi_q_k3']:.4f}")
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
