#!/usr/bin/env python3
"""Dense local scan for locked-action tau refinement near a beta anchor."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geometrodynamics.tangherlini.lepton_spectrum import (  # noqa: E402
    S3_ACTION_BASE,
    derive_geometric_beta,
)
from scripts.calibrate_muon_ratio import TAU_MEV, calibrate_grid  # noqa: E402
from scripts.map_basin_k_uplift import evaluate  # noqa: E402


@dataclass(frozen=True)
class LocalBest:
    beta: float
    phase_per_pass: float
    transport_strength: float
    hard_pinhole_gamma: float
    resistance_scale: float
    mu_rel_error: float
    tau_pred_mev: float
    tau_rel_error: float
    basin_mu_err_lo: float
    basin_mu_err_hi: float


def _parse_int_list(raw: str) -> list[int]:
    vals = [int(v.strip()) for v in raw.split(",") if v.strip()]
    if not vals:
        raise ValueError("beta-family list is empty")
    return vals


def refine_locked_scan(
    *,
    beta_min: float,
    beta_max: float,
    beta_step: float,
    phase_min: float,
    phase_max: float,
    phase_steps: int,
    transport_min: float,
    transport_max: float,
    transport_steps: int,
    pinhole_min: float,
    pinhole_max: float,
    pinhole_steps: int,
    resistance_min: float,
    resistance_max: float,
    n_points: int,
    action_base: float,
) -> list[LocalBest]:
    if beta_step <= 0.0:
        raise ValueError("beta_step must be positive")
    betas = np.arange(beta_min, beta_max + 0.5 * beta_step, beta_step)
    out: list[LocalBest] = []
    for beta in betas:
        candidates = calibrate_grid(
            phase_steps=phase_steps,
            transport_min=transport_min,
            transport_max=transport_max,
            transport_steps=transport_steps,
            resistance_min=resistance_min,
            resistance_max=resistance_max,
            pinhole_min=pinhole_min,
            pinhole_max=pinhole_max,
            pinhole_steps=pinhole_steps,
            n_points=n_points,
            depth_cost_mode="tunnel_only",
            winding_mode="max",
            k_uplift_beta=float(beta),
            action_base=action_base,
        )
        exact = [
            c
            for c in candidates
            if c.exact_mu_root and phase_min <= c.phase_per_pass <= phase_max
        ]
        if not exact:
            continue
        best = min(exact, key=lambda c: c.tau_rel_error)

        lo = evaluate(
            phase=best.phase_per_pass,
            transport=best.transport_strength,
            pinhole=best.hard_pinhole_gamma,
            resistance=0.99 * best.resistance_scale,
            beta=float(beta),
            action_base_locked=action_base,
            n_points=n_points,
        )
        hi = evaluate(
            phase=best.phase_per_pass,
            transport=best.transport_strength,
            pinhole=best.hard_pinhole_gamma,
            resistance=1.01 * best.resistance_scale,
            beta=float(beta),
            action_base_locked=action_base,
            n_points=n_points,
        )
        out.append(
            LocalBest(
                beta=float(beta),
                phase_per_pass=best.phase_per_pass,
                transport_strength=best.transport_strength,
                hard_pinhole_gamma=best.hard_pinhole_gamma,
                resistance_scale=best.resistance_scale,
                mu_rel_error=best.mu_e_rel_error,
                tau_pred_mev=best.tau_mev_pred,
                tau_rel_error=best.tau_rel_error,
                basin_mu_err_lo=lo.mu_rel_error,
                basin_mu_err_hi=hi.mu_rel_error,
            )
        )
    out.sort(key=lambda r: r.tau_rel_error)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--beta-min", type=float, default=130.0)
    p.add_argument("--beta-max", type=float, default=150.0)
    p.add_argument("--beta-step", type=float, default=0.5)
    p.add_argument("--phase-min", type=float, default=0.001)
    p.add_argument("--phase-max", type=float, default=0.01)
    p.add_argument("--phase-steps", type=int, default=13)
    p.add_argument("--transport-min", type=float, default=24.8)
    p.add_argument("--transport-max", type=float, default=25.2)
    p.add_argument("--transport-steps", type=int, default=9)
    p.add_argument("--pinhole-min", type=float, default=20.0)
    p.add_argument("--pinhole-max", type=float, default=25.0)
    p.add_argument("--pinhole-steps", type=int, default=11)
    p.add_argument("--resistance-min", type=float, default=0.19)
    p.add_argument("--resistance-max", type=float, default=0.24)
    p.add_argument("--n-points", type=int, default=24)
    p.add_argument("--top-k", type=int, default=5)
    p.add_argument("--action-base", type=float, default=S3_ACTION_BASE)
    p.add_argument("--beta-winding-integer", type=int, default=5)
    p.add_argument("--beta-family-integers", type=str, default="5,6")
    args = p.parse_args()

    beta_geom = derive_geometric_beta(winding_integer=args.beta_winding_integer)
    beta_family = _parse_int_list(args.beta_family_integers)
    print(f"Tau target = {TAU_MEV:.6f} MeV")
    print(f"action_base (locked) = {args.action_base:.12f}")
    print(
        "geometric beta estimate = "
        f"{beta_geom:.6f} (winding_integer={args.beta_winding_integer})"
    )
    print()

    results = refine_locked_scan(
        beta_min=args.beta_min,
        beta_max=args.beta_max,
        beta_step=args.beta_step,
        phase_min=args.phase_min,
        phase_max=args.phase_max,
        phase_steps=args.phase_steps,
        transport_min=args.transport_min,
        transport_max=args.transport_max,
        transport_steps=args.transport_steps,
        pinhole_min=args.pinhole_min,
        pinhole_max=args.pinhole_max,
        pinhole_steps=args.pinhole_steps,
        resistance_min=args.resistance_min,
        resistance_max=args.resistance_max,
        n_points=args.n_points,
        action_base=args.action_base,
    )

    print(f"Betas with exact mu/e roots in local window = {len(results)}")
    print()
    for i, r in enumerate(results[: args.top_k], start=1):
        print(f"#{i} beta={r.beta:.6f}")
        print(f"  phase_per_pass     = {r.phase_per_pass:.12f}")
        print(f"  transport_strength = {r.transport_strength:.12f}")
        print(f"  hard_pinhole_gamma = {r.hard_pinhole_gamma:.12f}")
        print(f"  resistance_scale   = {r.resistance_scale:.12f}")
        print(f"  mu/e rel. error    = {r.mu_rel_error:.6%}")
        print(f"  tau_pred (MeV)     = {r.tau_pred_mev:.6f}")
        print(f"  tau rel. error     = {r.tau_rel_error:.6%}")
        print(
            "  basin mu/e rel. err (±1% resistance) = "
            f"{r.basin_mu_err_lo:.6%}, {r.basin_mu_err_hi:.6%}"
        )
        print()

    if results:
        best = results[0]
        print(
            f"Best local locked fit beta={best.beta:.6f}, tau={best.tau_pred_mev:.3f} MeV "
            f"({best.tau_rel_error:.3%} error)."
        )
        print()
        print("Integer-anchored geometric beta family check:")
        for n in beta_family:
            b = derive_geometric_beta(winding_integer=n)
            near = min(results, key=lambda r: abs(r.beta - b))
            print(
                f"  n={n}: beta_geom={b:.6f}, nearest_scanned={near.beta:.6f}, "
                f"tau_err={near.tau_rel_error:.6%}"
            )


if __name__ == "__main__":
    main()
