#!/usr/bin/env python3
"""
scripts/experiment_constraint_search.py
========================================

Constraint-reduction pass on the pass-3 lock
(see docs/quark_axioms.md §8 "Pass 3").

The pass-3 best point sits at:
    N=460, eps=0.96, chi=19.8, eta=5.0,
    gamma_q=0.10, transport=0.55, pinhole=22.0, resistance=0.14,
    phase=0.0049    (max_rel_err = 0.0165, anchor=d).

Several of these numbers are suspiciously clean:
    chi · eta = 19.8 · 5.0 = 99    (~100, lepton tau winding)
    chi / eta = 19.8 / 5.0 = 3.96  (~4)
    eps = 0.96 = 24/25 exactly
    phase = 0.0049 ≈ 1/200 = 0.005

This script tests each candidate relation by fixing it and
re-running coordinate descent over the remaining axes.  If the
relation is structurally meaningful, the constrained best should
match the unconstrained 1.6% within scan resolution.  If not, the
residual error will rise — quantifying how "free" each candidate is.

A relation that holds the lock open at 1.6% with one fewer free
knob is a real reduction.  A relation that costs 5x in error is
a coincidence.
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

_PASS3_BEST = QuarkParams(
    action_base=QUARK_ACTION_BASE,
    beta=460 * math.pi / 2.0,
    gamma_q=0.10,
    u_q_form="k_minus_2",
    phase=0.0049,
    transport=0.55,
    pinhole=22.0,
    resistance=0.14,
    partition_mixing=0.0,
    winding_mode="max",
    resistance_model="exponential",
    depth_cost_mode="tunnel_only",
    uplift_mode="partition_asymmetric",
    uplift_asymmetry=0.96,
    spectrum_zero_mode="min_eigenvalue",
    chi_q_k3=19.8,
    eta_k3k5_minus=5.0,
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


def coord_descent(
    base: QuarkParams,
    free_axes: list[tuple[str, np.ndarray, Callable[[QuarkParams, float], QuarkParams]]],
    max_rounds: int = 6,
) -> tuple[QuarkParams, float]:
    """Descent over free_axes only; constrained axes stay locked."""
    current = base
    current_err = _max_rel_err_excl_u(current)
    for _ in range(max_rounds):
        improved = False
        for axis_name, values, apply in free_axes:
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
                current_err = best_e
                improved = True
        if not improved:
            break
    return current, current_err


# ── candidate constraints ───────────────────────────────────────────
# Each constraint takes the unconstrained pass-3 best, applies a
# constraint surface, then optimizes the remaining axes.

def _free_axes_full() -> list:
    """All 9 axes from pass 3, used as the 'no constraint' baseline."""
    return [
        ("integer_winding", np.arange(440, 481, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("uplift_asymmetry", np.arange(0.92, 1.00, 0.01),
            lambda p, v: replace(p, uplift_asymmetry=v)),
        ("chi_q_k3", np.arange(19.5, 20.1, 0.05),
            lambda p, v: replace(p, chi_q_k3=v)),
        ("eta_k3k5_minus", np.arange(4.0, 6.0, 0.1),
            lambda p, v: replace(p, eta_k3k5_minus=v)),
        ("gamma_q", np.arange(0.06, 0.16, 0.01),
            lambda p, v: replace(p, gamma_q=v)),
        ("transport", np.arange(0.45, 0.70, 0.02),
            lambda p, v: replace(p, transport=v)),
        ("pinhole", np.arange(20.0, 25.0, 0.5),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.10, 0.18, 0.01),
            lambda p, v: replace(p, resistance=v)),
        ("phase", np.arange(0.001, 0.008, 0.0005),
            lambda p, v: replace(p, phase=v)),
    ]


def constraint_chi_eta_product(constant: float) -> tuple[QuarkParams, list]:
    """Enforce chi · eta = constant.  Free knob: chi (eta = const/chi)."""
    base = replace(_PASS3_BEST, chi_q_k3=19.8, eta_k3k5_minus=constant / 19.8)
    def apply_chi(p, v):
        return replace(p, chi_q_k3=v, eta_k3k5_minus=constant / v)
    free = [
        ("integer_winding", np.arange(440, 481, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("uplift_asymmetry", np.arange(0.92, 1.00, 0.01),
            lambda p, v: replace(p, uplift_asymmetry=v)),
        ("chi (with eta=K/chi)", np.arange(18.0, 22.0, 0.1), apply_chi),
        ("gamma_q", np.arange(0.06, 0.16, 0.01),
            lambda p, v: replace(p, gamma_q=v)),
        ("transport", np.arange(0.45, 0.70, 0.02),
            lambda p, v: replace(p, transport=v)),
        ("pinhole", np.arange(20.0, 25.0, 0.5),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.10, 0.18, 0.01),
            lambda p, v: replace(p, resistance=v)),
        ("phase", np.arange(0.001, 0.008, 0.0005),
            lambda p, v: replace(p, phase=v)),
    ]
    return base, free


def constraint_chi_eta_ratio(ratio: float) -> tuple[QuarkParams, list]:
    """Enforce chi / eta = ratio.  Free knob: eta (chi = ratio*eta)."""
    base = replace(_PASS3_BEST, eta_k3k5_minus=5.0, chi_q_k3=ratio * 5.0)
    def apply_eta(p, v):
        return replace(p, eta_k3k5_minus=v, chi_q_k3=ratio * v)
    free = [
        ("integer_winding", np.arange(440, 481, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("uplift_asymmetry", np.arange(0.92, 1.00, 0.01),
            lambda p, v: replace(p, uplift_asymmetry=v)),
        ("eta (with chi=R*eta)", np.arange(4.0, 6.0, 0.1), apply_eta),
        ("gamma_q", np.arange(0.06, 0.16, 0.01),
            lambda p, v: replace(p, gamma_q=v)),
        ("transport", np.arange(0.45, 0.70, 0.02),
            lambda p, v: replace(p, transport=v)),
        ("pinhole", np.arange(20.0, 25.0, 0.5),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.10, 0.18, 0.01),
            lambda p, v: replace(p, resistance=v)),
        ("phase", np.arange(0.001, 0.008, 0.0005),
            lambda p, v: replace(p, phase=v)),
    ]
    return base, free


def constraint_eps_fixed(value: float) -> tuple[QuarkParams, list]:
    """Pin uplift_asymmetry to a specified clean value."""
    base = replace(_PASS3_BEST, uplift_asymmetry=value)
    free = [
        ("integer_winding", np.arange(440, 481, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        # eps NOT in free list
        ("chi_q_k3", np.arange(19.5, 20.1, 0.05),
            lambda p, v: replace(p, chi_q_k3=v)),
        ("eta_k3k5_minus", np.arange(4.0, 6.0, 0.1),
            lambda p, v: replace(p, eta_k3k5_minus=v)),
        ("gamma_q", np.arange(0.06, 0.16, 0.01),
            lambda p, v: replace(p, gamma_q=v)),
        ("transport", np.arange(0.45, 0.70, 0.02),
            lambda p, v: replace(p, transport=v)),
        ("pinhole", np.arange(20.0, 25.0, 0.5),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.10, 0.18, 0.01),
            lambda p, v: replace(p, resistance=v)),
        ("phase", np.arange(0.001, 0.008, 0.0005),
            lambda p, v: replace(p, phase=v)),
    ]
    return base, free


def constraint_phase_fixed(value: float) -> tuple[QuarkParams, list]:
    """Pin phase to a specified clean value."""
    base = replace(_PASS3_BEST, phase=value)
    free = [
        ("integer_winding", np.arange(440, 481, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("uplift_asymmetry", np.arange(0.92, 1.00, 0.01),
            lambda p, v: replace(p, uplift_asymmetry=v)),
        ("chi_q_k3", np.arange(19.5, 20.1, 0.05),
            lambda p, v: replace(p, chi_q_k3=v)),
        ("eta_k3k5_minus", np.arange(4.0, 6.0, 0.1),
            lambda p, v: replace(p, eta_k3k5_minus=v)),
        ("gamma_q", np.arange(0.06, 0.16, 0.01),
            lambda p, v: replace(p, gamma_q=v)),
        ("transport", np.arange(0.45, 0.70, 0.02),
            lambda p, v: replace(p, transport=v)),
        ("pinhole", np.arange(20.0, 25.0, 0.5),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.10, 0.18, 0.01),
            lambda p, v: replace(p, resistance=v)),
    ]
    return base, free


def constraint_chi_eta_int_ratio() -> tuple[QuarkParams, list]:
    """Enforce both chi=20 and eta=5 (clean integers, ratio 4)."""
    base = replace(_PASS3_BEST, chi_q_k3=20.0, eta_k3k5_minus=5.0)
    free = [
        ("integer_winding", np.arange(440, 481, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("uplift_asymmetry", np.arange(0.92, 1.00, 0.01),
            lambda p, v: replace(p, uplift_asymmetry=v)),
        # chi and eta both pinned
        ("gamma_q", np.arange(0.06, 0.16, 0.01),
            lambda p, v: replace(p, gamma_q=v)),
        ("transport", np.arange(0.45, 0.70, 0.02),
            lambda p, v: replace(p, transport=v)),
        ("pinhole", np.arange(20.0, 25.0, 0.5),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.10, 0.18, 0.01),
            lambda p, v: replace(p, resistance=v)),
        ("phase", np.arange(0.001, 0.008, 0.0005),
            lambda p, v: replace(p, phase=v)),
    ]
    return base, free


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    args = p.parse_args()

    # Unconstrained baseline (re-runs descent with the same finer
    # grids used inside the constraints, so all comparisons are
    # apples-to-apples)
    base = _PASS3_BEST
    free = _free_axes_full()
    final, err = coord_descent(base, free, max_rounds=4)
    print(f"\n[baseline] no constraint: max_rel_err = {err:.4e}")

    results = {"baseline": {"err": err, "params": asdict(final)}}

    # Constraints to test
    tests = [
        ("chi_eta_product=99",  constraint_chi_eta_product(99.0)),
        ("chi_eta_product=100", constraint_chi_eta_product(100.0)),
        ("chi_eta_product=99.16", constraint_chi_eta_product(19.8 * 5.008080)),
        ("chi_eta_ratio=4",     constraint_chi_eta_ratio(4.0)),
        ("chi_eta_ratio=3.96",  constraint_chi_eta_ratio(3.96)),
        ("chi=20, eta=5",       constraint_chi_eta_int_ratio()),
        ("eps=0.96 (=24/25)",   constraint_eps_fixed(24.0 / 25.0)),
        ("eps=1.0",             constraint_eps_fixed(1.0)),
        ("eps=0.95",            constraint_eps_fixed(0.95)),
        ("phase=0.005 (=1/200)", constraint_phase_fixed(0.005)),
        ("phase=0.0",            constraint_phase_fixed(0.0)),
    ]

    for label, (base_c, free_c) in tests:
        final_c, err_c = coord_descent(base_c, free_c, max_rounds=4)
        ratio = err_c / err if err > 0 else float("inf")
        print(f"[{label:<22}] err = {err_c:.4e}  "
              f"({ratio:5.2f}× baseline)")
        results[label] = {
            "err": err_c,
            "ratio_to_baseline": ratio,
            "params": asdict(final_c),
        }

    print("\n" + "=" * 72)
    print("Constraint-reduction summary")
    print("=" * 72)
    print(f"{'constraint':>26}  {'err':>10}  {'ratio':>8}")
    for label in ["baseline"] + [t[0] for t in tests]:
        e = results[label]["err"]
        r = results[label].get("ratio_to_baseline", 1.0)
        print(f"{label:>26}  {e:>10.4e}  {r:>8.2f}")

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(results, fh, indent=2)
        print(f"\nFull results written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
