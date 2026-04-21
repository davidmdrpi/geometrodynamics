#!/usr/bin/env python3
"""Calibrate knotted ladder parameters against the muon/electron ratio.

For each (phase_per_pass, transport_strength) pair on a grid, this script
solves for the resistance_scale that exactly reproduces the experimental
muon/electron mass ratio (depths 3/1), then reports the top candidates by
smallest tau prediction error (depth 5 after electron calibration).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.optimize import minimize_scalar

from geometrodynamics.tangherlini.lepton_spectrum import (
    S3_ACTION_BASE,
    calibrate_electron_predict_heavier,
    compute_knotted_lepton_spectrum,
    compute_tunneling_envelope,
)

ELECTRON_MEV = 0.51099895
MUON_MEV = 105.6583755
TAU_MEV = 1776.86
TARGET_MU_E_RATIO = MUON_MEV / ELECTRON_MEV
RESISTANCE_MODEL = "exponential"


@dataclass(frozen=True)
class CalibrationCandidate:
    phase_per_pass: float
    transport_strength: float
    hard_pinhole_gamma: float
    resistance_scale: float
    mu_e_ratio: float
    mu_e_rel_error: float
    tau_mev_pred: float
    tau_rel_error: float
    exact_mu_root: bool


def _mu_e_ratio(
    phase_per_pass: float,
    transport_strength: float,
    hard_pinhole_gamma: float,
    resistance_scale: float,
    n_points: int,
    depth_cost_mode: str,
    winding_mode: str,
    k_uplift_beta: float,
    action_base: float,
) -> float:
    spec = compute_knotted_lepton_spectrum(
        depths=(1, 3, 5),
        phase_per_pass=phase_per_pass,
        transport_strength=transport_strength,
        resistance_scale=resistance_scale,
        resistance_model=RESISTANCE_MODEL,
        hard_pinhole_gamma=hard_pinhole_gamma,
        n_points=n_points,
        depth_cost_mode=depth_cost_mode,
        winding_mode=winding_mode,
        k_uplift_beta=k_uplift_beta,
        action_base=action_base,
    )
    return spec[3] / spec[1]


def _solve_resistance_for_ratio(
    phase_per_pass: float,
    transport_strength: float,
    hard_pinhole_gamma: float,
    resistance_min: float,
    resistance_max: float,
    n_points: int,
    depth_cost_mode: str,
    winding_mode: str,
    k_uplift_beta: float,
    action_base: float,
) -> float | None:
    def objective(resistance_scale: float) -> float:
        return (
            _mu_e_ratio(
                phase_per_pass,
                transport_strength,
                hard_pinhole_gamma,
                resistance_scale,
                n_points=n_points,
                depth_cost_mode=depth_cost_mode,
                winding_mode=winding_mode,
                k_uplift_beta=k_uplift_beta,
                action_base=action_base,
            )
            - TARGET_MU_E_RATIO
        )

    lo = resistance_min
    hi = resistance_max
    f_lo = objective(lo)
    f_hi = objective(hi)

    # Expand bracket adaptively if needed.
    for _ in range(12):
        if f_lo == 0:
            return lo
        if f_hi == 0:
            return hi
        if np.sign(f_lo) != np.sign(f_hi):
            return float(brentq(objective, lo, hi, xtol=1e-10, rtol=1e-10, maxiter=200))
        hi *= 2.0
        f_hi = objective(hi)

    return None


def calibrate_grid(
    phase_steps: int,
    transport_min: float,
    transport_max: float,
    transport_steps: int,
    resistance_min: float,
    resistance_max: float,
    pinhole_min: float,
    pinhole_max: float,
    pinhole_steps: int,
    n_points: int,
    depth_cost_mode: str,
    winding_mode: str,
    k_uplift_beta: float,
    action_base: float = S3_ACTION_BASE,
) -> list[CalibrationCandidate]:
    phases = np.linspace(1e-3, np.pi / 8.0, phase_steps)
    transports = np.linspace(transport_min, transport_max, transport_steps)
    pinholes = np.linspace(pinhole_min, pinhole_max, pinhole_steps)

    candidates: list[CalibrationCandidate] = []
    for phase in phases:
        for transport in transports:
            for pinhole in pinholes:
                res = _solve_resistance_for_ratio(
                    phase_per_pass=float(phase),
                    transport_strength=float(transport),
                    hard_pinhole_gamma=float(pinhole),
                    resistance_min=resistance_min,
                    resistance_max=resistance_max,
                    n_points=n_points,
                    depth_cost_mode=depth_cost_mode,
                    winding_mode=winding_mode,
                    k_uplift_beta=k_uplift_beta,
                    action_base=action_base,
                )
                exact = res is not None
                if res is None:
                    obj = lambda x: abs(
                        _mu_e_ratio(
                            float(phase),
                            float(transport),
                            float(pinhole),
                            float(x),
                            n_points=n_points,
                            depth_cost_mode=depth_cost_mode,
                            winding_mode=winding_mode,
                            k_uplift_beta=k_uplift_beta,
                            action_base=action_base,
                        )
                        - TARGET_MU_E_RATIO
                    )
                    approx = minimize_scalar(
                        obj,
                        bounds=(resistance_min, resistance_max),
                        method="bounded",
                        options={"xatol": 1e-9, "maxiter": 400},
                    )
                    if not approx.success:
                        continue
                    res = float(approx.x)

                fit = calibrate_electron_predict_heavier(
                    depths=(1, 3, 5),
                    phase_per_pass=float(phase),
                    transport_strength=float(transport),
                    hard_pinhole_gamma=float(pinhole),
                    resistance_scale=float(res),
                    resistance_model=RESISTANCE_MODEL,
                    n_points=n_points,
                    depth_cost_mode=depth_cost_mode,
                    winding_mode=winding_mode,
                    k_uplift_beta=k_uplift_beta,
                    action_base=action_base,
                )
                mu_e = fit.predicted_mev[3] / fit.predicted_mev[1]
                tau_pred = fit.predicted_mev[5]
                tau_rel = abs(tau_pred - TAU_MEV) / TAU_MEV
                mu_rel = abs(mu_e - TARGET_MU_E_RATIO) / TARGET_MU_E_RATIO

                candidates.append(
                    CalibrationCandidate(
                        phase_per_pass=float(phase),
                        transport_strength=float(transport),
                        hard_pinhole_gamma=float(pinhole),
                        resistance_scale=float(res),
                        mu_e_ratio=float(mu_e),
                        mu_e_rel_error=float(mu_rel),
                        tau_mev_pred=float(tau_pred),
                        tau_rel_error=float(tau_rel),
                        exact_mu_root=exact,
                    )
                )

    candidates.sort(key=lambda c: (c.mu_e_rel_error, c.tau_rel_error))
    return candidates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase-steps", type=int, default=21)
    parser.add_argument("--transport-min", type=float, default=20.0)
    parser.add_argument("--transport-max", type=float, default=30.0)
    parser.add_argument("--transport-steps", type=int, default=23)
    parser.add_argument("--resistance-min", type=float, default=0.2)
    parser.add_argument("--resistance-max", type=float, default=0.5)
    parser.add_argument("--pinhole-min", type=float, default=0.0)
    parser.add_argument("--pinhole-max", type=float, default=12.0)
    parser.add_argument("--pinhole-steps", type=int, default=13)
    parser.add_argument("--n-points", type=int, default=32)
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--depth-cost-mode", type=str, default="tunnel_only")
    parser.add_argument("--winding-mode", type=str, default="max")
    parser.add_argument("--k-uplift-beta", type=float, default=0.0)
    parser.add_argument("--action-base", type=float, default=S3_ACTION_BASE)
    args = parser.parse_args()

    candidates = calibrate_grid(
        phase_steps=args.phase_steps,
        transport_min=args.transport_min,
        transport_max=args.transport_max,
        transport_steps=args.transport_steps,
        resistance_min=args.resistance_min,
        resistance_max=args.resistance_max,
        pinhole_min=args.pinhole_min,
        pinhole_max=args.pinhole_max,
        pinhole_steps=args.pinhole_steps,
        n_points=args.n_points,
        depth_cost_mode=args.depth_cost_mode,
        winding_mode=args.winding_mode,
        k_uplift_beta=args.k_uplift_beta,
        action_base=args.action_base,
    )

    print(f"Target mu/e ratio: {TARGET_MU_E_RATIO:.12f}")
    n_exact = sum(1 for c in candidates if c.exact_mu_root)
    print(f"Found {len(candidates)} parameter sets ({n_exact} exact mu/e roots).")
    print()
    for i, cand in enumerate(candidates[: args.top_k], start=1):
        print(f"#{i}")
        print(f"  phase_per_pass     = {cand.phase_per_pass:.12f}")
        print(f"  transport_strength = {cand.transport_strength:.12f}")
        print(f"  hard_pinhole_gamma = {cand.hard_pinhole_gamma:.12f}")
        print(f"  resistance_scale   = {cand.resistance_scale:.12f}")
        print(f"  mu/e ratio         = {cand.mu_e_ratio:.12f}")
        print(f"  mu/e rel. error    = {cand.mu_e_rel_error:.6%}")
        print(f"  exact root         = {cand.exact_mu_root}")
        print(f"  tau_pred (MeV)     = {cand.tau_mev_pred:.6f}")
        print(f"  tau rel. error     = {cand.tau_rel_error:.6%}")
        print()

    if candidates:
        best = candidates[0]
        print(f"action_base (locked) = {args.action_base:.12f}")
        env3 = compute_tunneling_envelope(
            depth=3,
            phase_per_pass=best.phase_per_pass,
            resistance_scale=best.resistance_scale,
            action_base=args.action_base,
            winding_mode=args.winding_mode,
        )
        upper = env3[np.triu_indices(3, k=1)]
        print("k=3 tunneling envelope (exp(-S_ij)):")
        print(f"  min={upper.min():.4e}  mean={upper.mean():.4e}  max={upper.max():.4e}")


if __name__ == "__main__":
    main()
