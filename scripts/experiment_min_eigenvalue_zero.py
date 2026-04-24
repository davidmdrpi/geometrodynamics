#!/usr/bin/env python3
"""
scripts/experiment_min_eigenvalue_zero.py
==========================================

Follow-up experiment to the partition-asymmetric uplift scan
(see docs/quark_axioms.md §8 "Follow-up experiment").  Combines:

  1. spectrum_zero_mode = "min_eigenvalue"  (candidate #3 from §3.5)
       Sets the zero to min(eigenvalues) after diagonalization, so u
       becomes massless by construction.  This removes the positivity
       constraint that kept gamma_q pinned small.
  2. uplift_mode = "partition_asymmetric"   (candidate #2, retained)
       Keeps the b/t-splitting knob from the previous experiment.
  3. Anchor species = d                     (required by #1)
       With u at exactly 0 under min_eigenvalue zero, the lightest
       non-anchor mass is d.  Scale such that predicted_d matches
       observed_d = 4.67 MeV.

Scan: (integer_winding, uplift_asymmetry, gamma_q, transport, pinhole)
with gamma_q pushed up to 5.0 now that positivity is no longer a cap.

Non-breaking: uses the opt-in modes; defaults untouched.  Reports a
JSON summary so next session can decide if this becomes a spec
revision.
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


def _max_rel_err(predicted: np.ndarray) -> float:
    # Skip u in the max (it's 0 by construction under min_eigenvalue
    # zero, so dividing by observed u is the "full error")
    non_anchor = np.array(
        [i for i, s in enumerate(QUARK_SPECIES) if s != "u"], dtype=int,
    )
    rel = np.abs(predicted - _TARGET) / _TARGET
    return float(np.max(rel[non_anchor]))


def run(
    integer_windings: np.ndarray,
    uplift_asymmetries: np.ndarray,
    gamma_qs: np.ndarray,
    transports: np.ndarray,
    pinholes: np.ndarray,
    resistances: np.ndarray,
    verbose: bool = False,
) -> dict[str, object]:
    axes = [
        ("integer_winding",   integer_windings),
        ("uplift_asymmetry",  uplift_asymmetries),
        ("gamma_q",           gamma_qs),
        ("transport",         transports),
        ("pinhole",           pinholes),
        ("resistance",        resistances),
    ]
    total = int(np.prod([len(v) for _, v in axes]))
    if verbose:
        print(f"scanning {total} points ... "
              f"(anchor = d = {_ANCHOR_MASS_MEV} MeV)")

    best_err = float("inf")
    best_params: Optional[QuarkParams] = None
    best_predicted: Optional[np.ndarray] = None
    rejected = 0
    for combo in itertools.product(*[v for _, v in axes]):
        N, eps, gq, tr, pn, rs = (float(x) for x in combo)
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
                print(f"  N={int(N):>5d}  eps={eps:+.3f}  "
                      f"gamma_q={gq:.3f}  transport={tr:.2f}  "
                      f"max_err={err:.4e}")
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

    # Compact but expressive grid.  gamma_q range goes well past the
    # r_q positivity bound that blocked the previous experiments.
    result = run(
        integer_windings=np.array([30, 60, 120, 250, 500, 1000, 2500, 5000]),
        uplift_asymmetries=np.linspace(-0.9, 0.9, 10),
        gamma_qs=np.array([0.05, 0.3, 1.0, 2.0, 3.5, 5.0]),
        transports=np.array([0.5, 1.5, 3.0]),
        pinholes=np.array([15.0, 22.5, 30.0]),
        resistances=np.array([0.15, 0.22, 0.30]),
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
