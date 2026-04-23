#!/usr/bin/env python3
"""Targeted sweep over k_uplift_beta under exact mu/e constraint."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.calibrate_muon_ratio import TAU_MEV, calibrate_grid  # noqa: E402
from geometrodynamics.tangherlini.lepton_spectrum import S3_ACTION_BASE  # noqa: E402


@dataclass(frozen=True)
class BetaResult:
    beta: float
    phase_per_pass: float
    transport_strength: float
    hard_pinhole_gamma: float
    resistance_scale: float
    tau_pred_mev: float
    tau_rel_error: float


def sweep_beta(
    beta_min: float,
    beta_max: float,
    beta_steps: int,
    phase_steps: int,
    transport_steps: int,
    pinhole_steps: int,
    n_points: int,
    action_base: float,
) -> list[BetaResult]:
    out: list[BetaResult] = []
    for beta in np.linspace(beta_min, beta_max, beta_steps):
        candidates = calibrate_grid(
            phase_steps=phase_steps,
            transport_min=20.0,
            transport_max=30.0,
            transport_steps=transport_steps,
            resistance_min=0.2,
            resistance_max=0.5,
            pinhole_min=0.0,
            pinhole_max=30.0,
            pinhole_steps=pinhole_steps,
            n_points=n_points,
            depth_cost_mode="tunnel_only",
            winding_mode="max",
            k_uplift_beta=float(beta),
            action_base=action_base,
        )
        exact = [c for c in candidates if c.exact_mu_root]
        if not exact:
            continue
        best = min(exact, key=lambda c: c.tau_rel_error)
        out.append(
            BetaResult(
                beta=float(beta),
                phase_per_pass=best.phase_per_pass,
                transport_strength=best.transport_strength,
                hard_pinhole_gamma=best.hard_pinhole_gamma,
                resistance_scale=best.resistance_scale,
                tau_pred_mev=best.tau_mev_pred,
                tau_rel_error=best.tau_rel_error,
            )
        )
    out.sort(key=lambda r: r.tau_rel_error)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--beta-min", type=float, default=40.0)
    p.add_argument("--beta-max", type=float, default=160.0)
    p.add_argument("--beta-center", type=float, default=None)
    p.add_argument("--beta-span", type=float, default=20.0)
    p.add_argument("--beta-steps", type=int, default=31)
    p.add_argument("--phase-steps", type=int, default=5)
    p.add_argument("--transport-steps", type=int, default=5)
    p.add_argument("--pinhole-steps", type=int, default=5)
    p.add_argument("--n-points", type=int, default=24)
    p.add_argument("--top-k", type=int, default=10)
    p.add_argument("--action-base", type=float, default=S3_ACTION_BASE)
    args = p.parse_args()

    beta_min = args.beta_min
    beta_max = args.beta_max
    if args.beta_center is not None:
        half = max(args.beta_span, 0.0) / 2.0
        beta_min = args.beta_center - half
        beta_max = args.beta_center + half

    results = sweep_beta(
        beta_min=beta_min,
        beta_max=beta_max,
        beta_steps=args.beta_steps,
        phase_steps=args.phase_steps,
        transport_steps=args.transport_steps,
        pinhole_steps=args.pinhole_steps,
        n_points=args.n_points,
        action_base=args.action_base,
    )

    print(f"Tau target = {TAU_MEV:.6f} MeV")
    print(f"action_base (locked) = {args.action_base:.12f}")
    print(f"Betas with exact mu/e roots = {len(results)}")
    print()
    for i, r in enumerate(results[: args.top_k], start=1):
        print(f"#{i} beta={r.beta:.6f}")
        print(f"  tau_pred (MeV)     = {r.tau_pred_mev:.6f}")
        print(f"  tau rel. error     = {r.tau_rel_error:.6%}")
        print(f"  phase_per_pass     = {r.phase_per_pass:.12f}")
        print(f"  transport_strength = {r.transport_strength:.12f}")
        print(f"  hard_pinhole_gamma = {r.hard_pinhole_gamma:.12f}")
        print(f"  resistance_scale   = {r.resistance_scale:.12f}")
        print()

    if results:
        best = results[0]
        print(
            f"Best beta ~{best.beta:.6f} gives tau {best.tau_pred_mev:.3f} MeV "
            f"({best.tau_rel_error:.3%} error)."
        )


if __name__ == "__main__":
    main()
