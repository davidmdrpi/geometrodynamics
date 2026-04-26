#!/usr/bin/env python3
"""
scripts/experiment_n_ablation.py
=================================

N-stability ablation, addressing the user's milestone:

    Run ablations that hold the four shell-index constraints fixed
    and test whether N shifts under alternate anchor conventions,
    alternate observed-mass sets, and alternate spectrum-zero
    definitions.  If N=466 is stable across those choices, it
    becomes a real target for topological interpretation.  If it
    moves, then N is an effective compensator for the residual
    transport/pinhole/resistance sector.

The four shell-index constraints (from §8 constraint-reduction):
    uplift_asymmetry  = 24/25     (= 1 − 1/k_5²)
    eta_k3k5_minus    = 5         (= k_5)
    chi_q_k3          = 20        (= (k_5−1)·k_5)
    phase             = 0
    gamma_q           = 1/10      (held empirical-clean)

Free knobs in every ablation: N, transport, pinhole, resistance.

Three ablation families:

  (1) Anchor species: d (baseline), s, c, b, t.
      For each, scale the predicted spectrum to that species'
      observed mass instead of d.  If N is set by topology, the
      best N should not depend on which species is the anchor.

  (2) Observed-mass perturbations: PDG (baseline), all-up-10%,
      all-down-10%, c+10%, b+10%, t+10%, t-10%.  A uniform
      multiplicative shift should leave optimal N fixed (it just
      rescales the MeV anchor).  Per-species perturbations test
      whether N is moving to compensate for individual masses.

  (3) Spectrum-zero: min_eigenvalue (baseline) vs second_min.
      The "second_min" mode shifts to the second-lightest
      eigenvalue, treating d as zero by construction (instead of
      u).  Implemented inline; reuses the extraction machinery.

Output: per-ablation best N + max_rel_err, plus a summary table.
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
    OBSERVED_MASSES_MEV as _PDG_MASSES,
    QUARK_ACTION_BASE,
    QUARK_SPECIES,
    QuarkParams,
    extract_physical_spectrum,
    build_quark_hamiltonian,
)


# ── Constraint-reduced template (everything except the 4 free knobs) ──
_TEMPLATE = QuarkParams(
    action_base=QUARK_ACTION_BASE,
    beta=466 * math.pi / 2.0,
    gamma_q=0.10,
    u_q_form="k_minus_2",
    phase=0.0,
    transport=0.54,
    pinhole=22.25,
    resistance=0.14,
    partition_mixing=0.0,
    winding_mode="max",
    resistance_model="exponential",
    depth_cost_mode="tunnel_only",
    uplift_mode="partition_asymmetric",
    uplift_asymmetry=24.0 / 25.0,
    spectrum_zero_mode="min_eigenvalue",
    chi_q_k3=20.0,
    eta_k3k5_minus=5.0,
)


def _err_factory(
    observed: dict[str, float],
    anchor_species: str,
    spectrum_zero_strategy: str = "min_eigenvalue",
) -> Callable[[QuarkParams], float]:
    """
    Returns an error function for the given (observed-masses, anchor,
    spectrum-zero) ablation choice.  Skip u from the max if u is at 0
    by construction (anchor==d under min_eigenvalue zero).
    """
    target = np.array([observed[s] for s in QUARK_SPECIES], dtype=float)
    anchor_idx = QUARK_SPECIES.index(anchor_species)
    anchor_obs = observed[anchor_species]

    if spectrum_zero_strategy == "min_eigenvalue":
        # Anchor sits at non-zero only if anchor != min_species (= u).
        # Skip u from max-error in min_eigenvalue mode.
        skip_indices = [QUARK_SPECIES.index("u")]
        params_override = lambda p: replace(p, spectrum_zero_mode="min_eigenvalue")
    elif spectrum_zero_strategy == "second_min":
        # Use the second-smallest eigenvalue as the zero — implemented
        # by computing the spectrum, sorting, and shifting by sorted[1].
        # Requires bypassing the lib's spectrum_zero_mode dispatch.
        skip_indices = [QUARK_SPECIES.index("u"), QUARK_SPECIES.index("d")]
        params_override = lambda p: replace(p, spectrum_zero_mode="min_eigenvalue")
    elif spectrum_zero_strategy == "action_base":
        skip_indices: list[int] = []
        params_override = lambda p: replace(p, spectrum_zero_mode="action_base")
    else:
        raise ValueError(f"Unknown spectrum_zero_strategy: {spectrum_zero_strategy}")

    def err(params: QuarkParams) -> float:
        try:
            species_map = extract_physical_spectrum(params_override(params))
        except Exception:
            return float("inf")
        if spectrum_zero_strategy == "second_min":
            # Re-shift by second-min eigenvalue
            vals = sorted(species_map.values())
            second_min = vals[1]
            species_map = {s: v - second_min for s, v in species_map.items()}
        anchor_val = species_map[anchor_species]
        if anchor_val <= 1e-6:
            return float("inf")
        scale = anchor_obs / anchor_val
        predicted = np.array(
            [species_map[s] * scale for s in QUARK_SPECIES], dtype=float,
        )
        rel = np.abs(predicted - target) / target
        # Mask the species that are anchored at 0 by construction
        keep = [i for i in range(len(QUARK_SPECIES)) if i not in skip_indices]
        return float(np.max(rel[keep]))
    return err


def coord_descent(
    base: QuarkParams,
    err: Callable[[QuarkParams], float],
    max_rounds: int = 6,
) -> tuple[QuarkParams, float]:
    """4 free axes: N, transport, pinhole, resistance.  Returns (lock, err)."""
    free = [
        ("beta",       np.arange(420, 521, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
        ("transport",  np.arange(0.40, 0.80, 0.01),
            lambda p, v: replace(p, transport=v)),
        ("pinhole",    np.arange(18.0, 27.0, 0.25),
            lambda p, v: replace(p, pinhole=v)),
        ("resistance", np.arange(0.08, 0.22, 0.005),
            lambda p, v: replace(p, resistance=v)),
    ]
    cur = base
    cur_err = err(cur)
    for _ in range(max_rounds):
        improved = False
        for name, vals, apply in free:
            best_v, best_e = None, cur_err
            for v in vals:
                e = err(apply(cur, float(v)))
                if e < best_e:
                    best_e = e
                    best_v = float(v)
            if best_v is not None:
                cur = apply(cur, best_v)
                cur_err = best_e
                improved = True
        if not improved:
            break
    return cur, cur_err


def best_N(lock: QuarkParams) -> int:
    return int(round(lock.beta * 2.0 / math.pi))


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    args = p.parse_args()

    # ── Baseline ─────────────────────────────────────────────────────
    err = _err_factory(_PDG_MASSES, "d", "min_eigenvalue")
    lock, e = coord_descent(_TEMPLATE, err)
    print(f"\n[baseline] anchor=d, PDG masses, min_eig:")
    print(f"  best N = {best_N(lock):>4d}  err = {e:.4e}")
    print(f"  (transport={lock.transport:.3f}, pinhole={lock.pinhole:.2f}, "
          f"resistance={lock.resistance:.4f})")
    results = {"baseline": {"N": best_N(lock), "err": e,
                             "params": asdict(lock)}}

    # ── (1) Anchor-species ablation ─────────────────────────────────
    print("\n── (1) Anchor-species ablation (PDG masses, min_eig) ──")
    for sp in ("s", "c", "b", "t"):
        err = _err_factory(_PDG_MASSES, sp, "min_eigenvalue")
        lock, e = coord_descent(_TEMPLATE, err)
        print(f"  anchor={sp}:  best N = {best_N(lock):>4d}  err = {e:.4e}")
        results[f"anchor={sp}"] = {"N": best_N(lock), "err": e}

    # ── (2) Observed-mass perturbations ─────────────────────────────
    print("\n── (2) Observed-mass perturbations (anchor=d, min_eig) ──")
    perturbations = {
        "PDG_x_1.10":  {s: m * 1.10 for s, m in _PDG_MASSES.items()},
        "PDG_x_0.90":  {s: m * 0.90 for s, m in _PDG_MASSES.items()},
        "c_x_1.10":    {**_PDG_MASSES, "c": _PDG_MASSES["c"] * 1.10},
        "b_x_1.10":    {**_PDG_MASSES, "b": _PDG_MASSES["b"] * 1.10},
        "t_x_1.10":    {**_PDG_MASSES, "t": _PDG_MASSES["t"] * 1.10},
        "t_x_0.90":    {**_PDG_MASSES, "t": _PDG_MASSES["t"] * 0.90},
        "all_perturb_5pct": {  # rng-free, deterministic ±5% per species
            "u": _PDG_MASSES["u"] * 1.05,
            "d": _PDG_MASSES["d"] * 0.95,
            "s": _PDG_MASSES["s"] * 1.05,
            "c": _PDG_MASSES["c"] * 0.95,
            "b": _PDG_MASSES["b"] * 1.05,
            "t": _PDG_MASSES["t"] * 0.95,
        },
    }
    for label, masses in perturbations.items():
        err = _err_factory(masses, "d", "min_eigenvalue")
        lock, e = coord_descent(_TEMPLATE, err)
        print(f"  {label:<22}: best N = {best_N(lock):>4d}  err = {e:.4e}")
        results[label] = {"N": best_N(lock), "err": e}

    # ── (3) Spectrum-zero ablation ──────────────────────────────────
    print("\n── (3) Spectrum-zero ablation (anchor=d, PDG masses) ──")
    for strat in ("second_min",):  # action_base will fail for d-anchor
        err = _err_factory(_PDG_MASSES, "d", strat)
        lock, e = coord_descent(_TEMPLATE, err)
        print(f"  zero={strat}:  best N = {best_N(lock):>4d}  err = {e:.4e}")
        results[f"zero={strat}"] = {"N": best_N(lock), "err": e}

    # ── Summary ─────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("N-ablation summary")
    print("=" * 72)
    print(f"{'ablation':>26}  {'best N':>8}  {'err':>12}")
    for k, v in results.items():
        print(f"{k:>26}  {v['N']:>8d}  {v['err']:>12.4e}")
    print()

    n_values = [v["N"] for k, v in results.items() if v["err"] < 0.5]
    if n_values:
        print(f"N range across well-fit ablations: "
              f"[{min(n_values)}, {max(n_values)}]  "
              f"(median {sorted(n_values)[len(n_values)//2]})")

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(results, fh, indent=2)
        print(f"\nFull results written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
