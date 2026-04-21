#!/usr/bin/env python3
"""Local basin mapping around an exact mu/e root for the lepton surrogate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from scripts.calibrate_muon_ratio import ELECTRON_MEV, MUON_MEV, TAU_MEV
from geometrodynamics.tangherlini.lepton_spectrum import (
    LEPTON_BASELINE_PHASE,
    LEPTON_BASELINE_PINHOLE,
    LEPTON_BASELINE_RESISTANCE,
    LEPTON_BASELINE_TRANSPORT,
    S3_ACTION_BASE,
    TAU_BETA_50PI,
    calibrate_electron_predict_heavier,
)

TARGET_RATIO = MUON_MEV / ELECTRON_MEV


@dataclass
class EvalResult:
    mu_e_ratio: float
    mu_rel_error: float
    tau_pred_mev: float
    tau_rel_error: float


def evaluate(
    *,
    phase: float,
    transport: float,
    pinhole: float,
    resistance: float,
    beta: float,
    action_base_locked: float,
    n_points: int,
) -> EvalResult:
    fit = calibrate_electron_predict_heavier(
        depths=(1, 3, 5),
        phase_per_pass=phase,
        transport_strength=transport,
        hard_pinhole_gamma=pinhole,
        resistance_scale=resistance,
        resistance_model="exponential",
        depth_cost_mode="tunnel_only",
        winding_mode="max",
        k_uplift_beta=beta,
        action_base=action_base_locked,
        n_points=n_points,
    )
    ratio = fit.predicted_mev[3] / fit.predicted_mev[1]
    mu_rel = abs(ratio - TARGET_RATIO) / TARGET_RATIO
    tau = fit.predicted_mev[5]
    tau_rel = abs(tau - TAU_MEV) / TAU_MEV
    return EvalResult(ratio, mu_rel, tau, tau_rel)


def objective(x: np.ndarray, n_points: int, action_base_locked: float) -> float:
    phase, transport, pinhole, resistance, beta = x
    res = evaluate(
        phase=phase,
        transport=transport,
        pinhole=pinhole,
        resistance=resistance,
        beta=beta,
        action_base_locked=action_base_locked,
        n_points=n_points,
    )
    # Strongly enforce mu/e root, then optimize tau.
    return 1e7 * (res.mu_rel_error**2) + res.tau_rel_error**2


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-points", type=int, default=24)
    p.add_argument("--phase", type=float, default=LEPTON_BASELINE_PHASE)
    p.add_argument("--transport", type=float, default=LEPTON_BASELINE_TRANSPORT)
    p.add_argument("--pinhole", type=float, default=LEPTON_BASELINE_PINHOLE)
    p.add_argument("--resistance", type=float, default=LEPTON_BASELINE_RESISTANCE)
    p.add_argument("--beta", type=float, default=TAU_BETA_50PI)
    p.add_argument("--beta-min", type=float, default=150.0)
    p.add_argument("--beta-max", type=float, default=165.0)
    p.add_argument("--action-base", type=float, default=S3_ACTION_BASE)
    args = p.parse_args()
    action_base_locked = float(args.action_base)

    x0 = np.array([
        args.phase,
        args.transport,
        args.pinhole,
        args.resistance,
        args.beta,
    ])

    bounds = [
        (1e-3, np.pi / 8),
        (20.0, 30.0),
        (0.0, 30.0),
        (0.2, 0.5),
        (args.beta_min, args.beta_max),
    ]

    sol = minimize(
        objective,
        x0,
        args=(args.n_points, action_base_locked),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 400},
    )

    best = evaluate(
        phase=float(sol.x[0]),
        transport=float(sol.x[1]),
        pinhole=float(sol.x[2]),
        resistance=float(sol.x[3]),
        beta=float(sol.x[4]),
        action_base_locked=action_base_locked,
        n_points=args.n_points,
    )

    print("Optimized local solution:")
    print(f"  phase_per_pass     = {sol.x[0]:.12f}")
    print(f"  transport_strength = {sol.x[1]:.12f}")
    print(f"  hard_pinhole_gamma = {sol.x[2]:.12f}")
    print(f"  resistance_scale   = {sol.x[3]:.12f}")
    print(f"  k_uplift_beta      = {sol.x[4]:.12f}")
    print(f"  action_base (lock) = {action_base_locked:.12f}")
    print(f"  mu/e ratio         = {best.mu_e_ratio:.12f}")
    print(f"  mu/e rel. error    = {best.mu_rel_error:.6%}")
    print(f"  tau_pred (MeV)     = {best.tau_pred_mev:.6f}")
    print(f"  tau rel. error     = {best.tau_rel_error:.6%}")

    # 1% perturbation basin map for resistance only (action base is locked).
    print("\n1% perturbation map (holding optimized point):")
    base_resistance = float(sol.x[3])
    for frac in (0.99, 1.00, 1.01):
        x = sol.x.copy()
        x[3] = base_resistance * frac
        r = evaluate(
            phase=float(x[0]),
            transport=float(x[1]),
            pinhole=float(x[2]),
            resistance=float(x[3]),
            beta=float(x[4]),
            action_base_locked=action_base_locked,
            n_points=args.n_points,
        )
        print(
            f"  resistance x{frac:.2f}: mu_err={r.mu_rel_error:.6%}, "
            f"tau={r.tau_pred_mev:.3f} MeV, tau_err={r.tau_rel_error:.6%}"
        )


if __name__ == "__main__":
    main()
