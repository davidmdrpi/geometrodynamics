#!/usr/bin/env python3
"""Probe locked manifold with beta hard-fixed to 50π and optimize remaining knobs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geometrodynamics.tangherlini.lepton_spectrum import (  # noqa: E402
    S3_ACTION_BASE,
    TAU_BETA_50PI,
    tau_uplift_2pi_quanta,
)
from scripts.map_basin_k_uplift import evaluate  # noqa: E402


def _objective(x: np.ndarray, n_points: int, beta: float, action_base_locked: float) -> float:
    phase, transport, pinhole, resistance = x
    res = evaluate(
        phase=float(phase),
        transport=float(transport),
        pinhole=float(pinhole),
        resistance=float(resistance),
        beta=float(beta),
        action_base_locked=float(action_base_locked),
        n_points=n_points,
    )
    return 1e7 * (res.mu_rel_error**2) + res.tau_rel_error**2


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-points", type=int, default=24)
    p.add_argument("--phase", type=float, default=0.001)
    p.add_argument("--transport", type=float, default=25.1)
    p.add_argument("--pinhole", type=float, default=22.5)
    p.add_argument("--resistance", type=float, default=0.217875196931)
    p.add_argument("--action-base", type=float, default=S3_ACTION_BASE)
    p.add_argument("--beta", type=float, default=TAU_BETA_50PI)
    args = p.parse_args()

    x0 = np.array([args.phase, args.transport, args.pinhole, args.resistance], dtype=float)
    bounds = [
        (1e-3, 0.01),
        (24.5, 25.5),
        (20.0, 25.0),
        (0.18, 0.25),
    ]

    sol = minimize(
        _objective,
        x0,
        args=(args.n_points, args.beta, args.action_base),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 400},
    )

    best = evaluate(
        phase=float(sol.x[0]),
        transport=float(sol.x[1]),
        pinhole=float(sol.x[2]),
        resistance=float(sol.x[3]),
        beta=float(args.beta),
        action_base_locked=float(args.action_base),
        n_points=args.n_points,
    )

    print("Locked beta probe:")
    print(f"  beta                 = {args.beta:.12f}")
    print(f"  uplift_2pi_quanta    = {tau_uplift_2pi_quanta(args.beta):.6f}")
    print(f"  phase_per_pass       = {sol.x[0]:.12f}")
    print(f"  transport_strength   = {sol.x[1]:.12f}")
    print(f"  hard_pinhole_gamma   = {sol.x[2]:.12f}")
    print(f"  resistance_scale     = {sol.x[3]:.12f}")
    print(f"  mu/e rel. error      = {best.mu_rel_error:.6%}")
    print(f"  tau_pred (MeV)       = {best.tau_pred_mev:.6f}")
    print(f"  tau rel. error       = {best.tau_rel_error:.6%}")


if __name__ == "__main__":
    main()
