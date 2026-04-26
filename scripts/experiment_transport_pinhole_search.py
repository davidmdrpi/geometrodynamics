#!/usr/bin/env python3
"""
scripts/experiment_transport_pinhole_search.py
================================================

User-named milestone (axioms §8 next-session priorities, refined):

  1. Transport: a 1D geometric search.  transport = ½ cos(χ) is the
     Hopf-connection ansatz.  Since the fitted value 0.54 exceeds
     the connection's max of ½, normalisation may be wrong, not just
     χ.  Test connection / holonomy / curvature / overlap forms and
     find which (if any) reproduces 0.54 at a clean χ.

  2. Pinhole: refine ∑V_max.  The minimal candidate is 21.80
     (sum over l=1..5).  Variants to try:
       - only odd l (l=1, 3, 5)
       - degeneracy-weighted (2l+1)
       - evaluated at classical turning points instead of V_max
       - on the same tortoise-coordinate grid as the eigensolver

For each candidate, also do a joint-fit check: pin the residual
to that value (others free), let N drift, report (N, max_rel_err)
under PDG masses.  A candidate is "viable" if the joint-fit
error stays within 2× of the unpinned baseline (≈ 0.032).
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
from typing import Callable

import numpy as np

from geometrodynamics.qcd.quark_spectrum import (
    OBSERVED_MASSES_MEV as _PDG_MASSES,
    QUARK_ACTION_BASE,
    QUARK_SPECIES,
    QuarkParams,
    extract_physical_spectrum,
)
from geometrodynamics.hopf.connection import (
    hopf_connection,
    hopf_holonomy,
    hopf_curvature,
)
from geometrodynamics.tangherlini.radial import (
    solve_radial_modes,
    V_tangherlini,
    r_to_rstar,
    rstar_to_r,
)
from geometrodynamics.constants import R_MID, R_OUTER, R_INNER

# ── Constraint-reduced template ───────────────────────────────────
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


def _err_for_masses(observed: dict[str, float]) -> Callable[[QuarkParams], float]:
    target = np.array([observed[s] for s in QUARK_SPECIES], dtype=float)
    skip = [QUARK_SPECIES.index("u")]

    def err(params: QuarkParams) -> float:
        try:
            sm = extract_physical_spectrum(params)
        except Exception:
            return float("inf")
        if sm["d"] <= 1e-6:
            return float("inf")
        scale = observed["d"] / sm["d"]
        pred = np.array([sm[s] * scale for s in QUARK_SPECIES], dtype=float)
        rel = np.abs(pred - target) / target
        keep = [i for i in range(len(QUARK_SPECIES)) if i not in skip]
        return float(np.max(rel[keep]))
    return err


def _coord_descent(base, err, free, max_rounds=6):
    cur, cur_err = base, err(base)
    for _ in range(max_rounds):
        improved = False
        for name, vals, apply in free:
            best_v, best_e = None, cur_err
            for v in vals:
                e = err(apply(cur, float(v)))
                if e < best_e:
                    best_e, best_v = e, float(v)
            if best_v is not None:
                cur = apply(cur, best_v)
                cur_err = best_e
                improved = True
        if not improved:
            break
    return cur, cur_err


# ── Joint-fit helper: pin one residual, let N + others vary ───────
def joint_fit(
    template: QuarkParams,
    pinned_field: str,
    pinned_value: float,
    masses: dict[str, float],
) -> tuple[int, float]:
    base = replace(template, **{pinned_field: pinned_value})
    err = _err_for_masses(masses)
    free_names = ["transport", "pinhole", "resistance"]
    free_names.remove(pinned_field)
    free = [
        ("beta", np.arange(380, 561, 2),
            lambda p, v: replace(p, beta=v * math.pi / 2.0)),
    ]
    grids = {
        "transport": np.arange(0.40, 0.80, 0.01),
        "pinhole": np.arange(18.0, 27.0, 0.25),
        "resistance": np.arange(0.06, 0.22, 0.005),
    }
    for f in free_names:
        free.append((f, grids[f],
                     lambda p, v, _f=f: replace(p, **{_f: v})))
    lock, err_val = _coord_descent(base, err, free)
    return int(round(lock.beta * 2.0 / math.pi)), err_val


# ── (1) Transport-form sweep ──────────────────────────────────────
def transport_search() -> dict:
    """
    For each candidate functional form f(χ), find:
      - The χ that solves f(χ) = 0.54 (if reachable)
      - The numerical f at a few canonical χ values
      - Joint-fit (N, err) when transport is pinned to f(χ_canonical)
    """
    candidates = [
        ("(1/2)*cos(chi)",     lambda chi: 0.5 * math.cos(chi),       0.5),
        ("cos(chi)",           lambda chi: math.cos(chi),             1.0),
        ("(1/2)*sin(chi)",     lambda chi: 0.5 * math.sin(chi),       0.5),
        ("sin(chi)",           lambda chi: math.sin(chi),             1.0),
        ("pi*cos(chi)",        lambda chi: math.pi * math.cos(chi),   math.pi),
        ("(1/pi)*pi*cos(chi)", lambda chi: math.cos(chi),             1.0),  # same as cos
        ("cos²(chi)",          lambda chi: math.cos(chi) ** 2,        1.0),
        ("sin²(chi)",          lambda chi: math.sin(chi) ** 2,        1.0),
    ]
    results = {}
    for label, f, max_amp in candidates:
        # Solve f(chi) = 0.54
        if max_amp >= 0.54:
            # Brute-force scan to find chi
            chis = np.linspace(0.0, math.pi, 1024)
            values = np.array([f(c) for c in chis])
            idx = int(np.argmin(np.abs(values - 0.54)))
            chi_match = float(chis[idx])
            value_at_match = float(values[idx])
        else:
            chi_match = None
            value_at_match = None

        # Sample at canonical χ values and pick closest to 0.54
        canonical_chis = {
            "0":     0.0,
            "pi/6":  math.pi / 6,
            "pi/4":  math.pi / 4,
            "pi/3":  math.pi / 3,
            "pi/2":  math.pi / 2,
        }
        canonical_values = {k: float(f(v)) for k, v in canonical_chis.items()}

        # Pick the canonical χ that gives the value closest to 0.54
        best_canon_label = min(canonical_values,
                               key=lambda k: abs(canonical_values[k] - 0.54))
        best_canon_value = canonical_values[best_canon_label]

        # Joint fit at the closest canonical
        if 0.30 <= best_canon_value <= 0.78:
            N, err = joint_fit(_TEMPLATE, "transport", best_canon_value, _PDG_MASSES)
        else:
            N, err = None, None

        results[label] = {
            "max_amplitude": max_amp,
            "chi_solving_0.54": chi_match,
            "value_at_chi_solve": value_at_match,
            "canonical_values": canonical_values,
            "best_canonical_chi": best_canon_label,
            "best_canonical_value": best_canon_value,
            "joint_fit_N": N,
            "joint_fit_err": err,
        }
    return results


# ── (2) Pinhole refinement search ─────────────────────────────────
def pinhole_search() -> dict:
    """
    Compute V_max(l) for l=1..5 and combine in several refined ways.
    For each combination, also run a joint fit to see which gives
    the cleanest match to the fitted 22.25.
    """
    # Solve modes
    modes = {}
    rs = np.linspace(R_INNER + 0.005, R_OUTER - 0.005, 600)
    V_max_per_l = {}
    for l in range(1, 7):
        V_max_per_l[l] = float(np.max(V_tangherlini(rs, l, R_MID)))
        oms, _, _ = solve_radial_modes(l, N=80, n_modes=4)
        modes[l] = oms

    # Turning point: V(r*) = ω². For ground state ω(l, 0).
    def turning_point_V(l: int) -> float:
        """V at the outer classical turning point r* of the l ground state."""
        omega_sq = float(modes[l][0]) ** 2
        # Find r* > sqrt(rs²)·1.05 where V(r*) = ω². V is monotonic
        # decreasing past V_max so the outer turning point is the
        # outer of the two roots of V - ω² = 0.
        Vs = V_tangherlini(rs, l, R_MID)
        diff = Vs - omega_sq
        # Find sign changes
        sign_changes = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]
        if len(sign_changes) == 0:
            return float("nan")
        # Outer turning point: last sign change (largest r)
        i = sign_changes[-1]
        r_tp = float(rs[i])
        return float(V_tangherlini(r_tp, l, R_MID))

    # Sum on tortoise-coordinate grid: same as the eigensolver's grid
    def tortoise_max_V(l: int) -> float:
        rs_min = r_to_rstar(R_MID + 5e-4, R_MID)
        rs_max = r_to_rstar(R_OUTER - 5e-4, R_MID)
        N = 80
        x, _ = (np.cos(np.pi * np.arange(N + 1) / N), None)
        L = (rs_max - rs_min) / 2.0
        rsg = rs_min + L * (1.0 - x)
        rg = np.array([rstar_to_r(s, R_MID) for s in rsg])
        Vg = V_tangherlini(rg, l, R_MID)
        return float(np.max(Vg))

    candidates = {
        "Sum_l=1..5_V_max":
            sum(V_max_per_l[l] for l in range(1, 6)),
        "Sum_l=1..5_V_max_extended_l=6":
            sum(V_max_per_l[l] for l in range(1, 7)),
        "Sum_odd_l_in_{1,3,5}_V_max":
            V_max_per_l[1] + V_max_per_l[3] + V_max_per_l[5],
        "Sum_odd_l_with_factor_2":
            2.0 * (V_max_per_l[1] + V_max_per_l[3] + V_max_per_l[5]),
        "Sum_l=1..5_(2l+1)_V_max":
            sum((2 * l + 1) * V_max_per_l[l] for l in range(1, 6)),
        "Sum_l=1..5_l_V_max":
            sum(l * V_max_per_l[l] for l in range(1, 6)),
        "Sum_l=1..5_V_at_turning":
            sum(turning_point_V(l) for l in range(1, 6)),
        "Sum_l=1..5_V_on_tortoise_grid":
            sum(tortoise_max_V(l) for l in range(1, 6)),
        "V_max_l=5_only":
            V_max_per_l[5],
        "V_max_l=5_x_pi":
            math.pi * V_max_per_l[5],
        "V_max_l=5_x_e":
            math.e * V_max_per_l[5],
    }

    aux = {
        "V_max_per_l": V_max_per_l,
        "omega_per_l": {l: float(modes[l][0]) for l in range(1, 7)},
    }

    # Joint fit for each candidate within plausible range
    results = {}
    for label, value in candidates.items():
        if 5.0 <= value <= 50.0:
            N, err = joint_fit(_TEMPLATE, "pinhole", value, _PDG_MASSES)
        else:
            N, err = None, None
        results[label] = {
            "value": float(value),
            "diff_from_fitted": float(value - 22.25),
            "rel_diff_from_fitted": float((value - 22.25) / 22.25),
            "joint_fit_N": N,
            "joint_fit_err": err,
        }

    return {"candidates": results, "aux": aux}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output-json", default=None)
    args = p.parse_args()

    print("\n" + "=" * 72)
    print("(1) Transport — geometric form search")
    print("=" * 72)
    print(f"Fitted value: 0.54.  hopf_connection(0) = 0.5 (max) — too small.")
    print()
    transport_results = transport_search()
    print(f"{'form':<22}  {'max_amp':>8}  "
          f"{'chi=0.54':>10}  {'best_can':>10}  {'val':>6}  "
          f"{'N':>4}  {'err':>10}")
    print("-" * 90)
    for label, r in transport_results.items():
        chi_str = f"{r['chi_solving_0.54']:.4f}" if r['chi_solving_0.54'] else "N/A"
        N_str = f"{r['joint_fit_N']}" if r['joint_fit_N'] else "—"
        err_str = f"{r['joint_fit_err']:.4e}" if r['joint_fit_err'] is not None else "—"
        print(f"{label:<22}  {r['max_amplitude']:>8.4f}  "
              f"{chi_str:>10}  {r['best_canonical_chi']:>10}  "
              f"{r['best_canonical_value']:>6.3f}  "
              f"{N_str:>4}  {err_str:>10}")

    print("\n" + "=" * 72)
    print("(2) Pinhole — refined ∑V_max construction search")
    print("=" * 72)
    print(f"Fitted value: 22.25.")
    print()
    pinhole_results = pinhole_search()
    aux = pinhole_results["aux"]
    print("V_max(l) per l:")
    for l, v in aux["V_max_per_l"].items():
        print(f"  l={l}: V_max = {v:>8.4f}, omega(l,n=0) = {aux['omega_per_l'][l]:>6.4f}")
    print()
    print(f"{'candidate':<40}  {'value':>10}  {'rel_diff':>10}  "
          f"{'N':>4}  {'err':>10}")
    print("-" * 90)
    for label, r in pinhole_results["candidates"].items():
        N_str = f"{r['joint_fit_N']}" if r['joint_fit_N'] else "—"
        err_str = f"{r['joint_fit_err']:.4e}" if r['joint_fit_err'] is not None else "—"
        rel = r['rel_diff_from_fitted']
        print(f"{label:<40}  {r['value']:>10.4f}  {rel:>+10.2%}  "
              f"{N_str:>4}  {err_str:>10}")

    out = {
        "transport_search": transport_results,
        "pinhole_search": pinhole_results,
    }
    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"\nFull results written to {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
