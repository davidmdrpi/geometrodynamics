"""Instanton-transition lepton spectrum surrogate.

This module models depth-k lepton states as a k-state instanton transition
matrix, rather than a coupled differential operator.  Off-diagonal couplings
are exponentially suppressed by a topological action S_topological.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Iterable, Sequence

import numpy as np
from scipy.linalg import eigh

from geometrodynamics.constants import CAVITY_GAMMA, R_MID
from geometrodynamics.hopf.connection import hopf_holonomy

@dataclass(frozen=True)
class Crossing:
    """Crossing descriptor.

    ``identified=True`` means a true bulk-identified node; these are allowed
    to modify transition rates.  Projected crossings are ignored.
    """

    pass_a: int
    pass_b: int
    identified: bool
    coupling: float = 0.0


@dataclass(frozen=True)
class LadderFit:
    """Electron-calibrated ladder fit summary for depths (1, 3, 5)."""

    scale_mev: float
    predicted_mev: dict[int, float]
    relative_error: dict[int, float]
    phase_per_pass: float
    resistance_scale: float
    transport_strength: float


def _resistance_sequence(k: int, model: str, scale: float) -> np.ndarray:
    j = np.arange(1, k + 1, dtype=float)
    if model == "none":
        return np.zeros(k, dtype=float)
    if model == "writhe":
        return scale * np.sqrt(j)
    if model == "curvature":
        return scale * j
    if model == "exponential":
        return scale * (np.exp(j) - 1.0)
    raise ValueError(f"Unknown resistance model: {model}")


def _as_crossings(crossings: Iterable[Crossing | dict] | None) -> list[Crossing]:
    if crossings is None:
        return []
    out: list[Crossing] = []
    for c in crossings:
        out.append(c if isinstance(c, Crossing) else Crossing(**c))
    return out


def _instanton_action(
    i: int,
    j: int,
    k: int,
    phase_per_pass: float,
    action_base: float,
    action_slope: float,
    resistance: np.ndarray,
    winding_mode: str,
) -> float:
    del k
    sep = abs(i - j)
    # === BAM GEOMETRIC ACTION ===
    # Use difference in winding multiplicity, not global offset.
    if winding_mode == "delta":
        delta_k = sep
    elif winding_mode == "max":
        delta_k = max(i + 1, j + 1)
    else:
        raise ValueError(f"Unknown winding_mode: {winding_mode}")
    alpha = max(action_slope, 1e-6)
    hopf_term = 0.25 * abs(hopf_holonomy(phase_per_pass)) * delta_k
    closure_term = 0.5 * CAVITY_GAMMA * (1.0 - np.cos(phase_per_pass * sep))
    resistance_term = 0.05 * (resistance[i] + resistance[j])
    s_topological = alpha * delta_k + hopf_term + closure_term + resistance_term
    return s_topological


@lru_cache(maxsize=1)
def _default_core_scale() -> float:
    """Use existing condensate geometry to set a default core scale l."""
    from geometrodynamics.blackhole.condensate import build_schwarzschild_condensate

    return float(build_schwarzschild_condensate(mass=1.0).core_scale)


def _build_instanton_matrix(
    k: int,
    phase_per_pass: float,
    transport_strength: float,
    resistance_model: str,
    resistance_scale: float,
    action_base: float,
    action_slope: float,
    hard_pinhole_gamma: float,
    hard_pinhole_depths: Sequence[int],
    crossings: list[Crossing],
    depth_cost_mode: str,
    winding_mode: str,
) -> np.ndarray:
    resistance = _resistance_sequence(k, resistance_model, resistance_scale)
    l_core = _default_core_scale()
    if action_base <= 0.0:
        # Requested geometric baseline: S0 = π * (R_MID / l)
        action_base = float(np.pi * (R_MID / l_core))

    if depth_cost_mode not in {"both", "diag_only", "tunnel_only"}:
        raise ValueError(f"Unknown depth_cost_mode: {depth_cost_mode}")

    diag_action = action_base + resistance_scale * k
    tunnel_action_base = action_base if depth_cost_mode in {"both", "tunnel_only"} else 0.0

    h = np.zeros((k, k), dtype=float)

    # Flattened diagonal: compact S^3 bulk with linear winding cost.
    for i in range(k):
        h[i, i] = diag_action

    # Instanton transitions: exp(-S_topological) suppression.
    for i in range(k):
        for j in range(i + 1, k):
            s_ij = _instanton_action(
                i=i,
                j=j,
                k=k,
                phase_per_pass=phase_per_pass,
                action_base=tunnel_action_base,
                action_slope=action_slope,
                resistance=resistance,
                winding_mode=winding_mode,
            )
            amp = transport_strength * np.exp(-s_ij)
            # Hopf/Möbius phase proxy projected to scalar real sector.
            amp *= np.cos(phase_per_pass * abs(i - j))
            h[i, j] -= amp
            h[j, i] -= amp

    # True identified crossings lower (or raise) local action barrier.
    for c in crossings:
        if not c.identified:
            continue
        ia = (c.pass_a - 1) % k
        ib = (c.pass_b - 1) % k
        g = float(c.coupling)
        h[ia, ib] -= g
        h[ib, ia] -= g
        h[ia, ia] += abs(g)
        h[ib, ib] += abs(g)

    # Optional Path-A pinhole as a localized delta-like barrier at mid pass.
    if hard_pinhole_gamma != 0.0 and k in hard_pinhole_depths:
        m = (k - 1) // 2
        h[m, m] += hard_pinhole_gamma

    return h


def _build_generation_block(
    depths: Sequence[int],
    phase_per_pass: float,
    transport_strength: float,
    resistance_model: str,
    resistance_scale: float,
    hard_pinhole_gamma: float,
    action_base: float,
    action_slope: float,
    depth_cost_mode: str,
    winding_mode: str,
    k_uplift_beta: float,
) -> np.ndarray:
    """Strongly mixed generation-level block for selected depths."""
    n = len(depths)
    h = np.zeros((n, n), dtype=float)
    if action_base <= 0.0:
        l_core = _default_core_scale()
        action_base = float(np.pi * (R_MID / l_core))
    if resistance_model == "none":
        res_diag = np.zeros(n, dtype=float)
    elif resistance_model == "writhe":
        res_diag = resistance_scale * np.sqrt(np.asarray(depths, dtype=float))
    elif resistance_model == "curvature":
        res_diag = resistance_scale * np.asarray(depths, dtype=float)
    elif resistance_model == "exponential":
        res_diag = resistance_scale * (np.exp(np.asarray(depths, dtype=float)) - 1.0)
    else:
        raise ValueError(f"Unknown resistance model: {resistance_model}")

    for a, k in enumerate(depths):
        base = action_base if depth_cost_mode in {"both", "diag_only", "tunnel_only"} else 0.0
        pin = hard_pinhole_gamma if k in (3, 5) else 0.0
        uplift = k_uplift_beta * max(0.0, k - 3.0) ** 2
        h[a, a] = base + resistance_scale * (k**2) + res_diag[a] + pin + uplift
    alpha = max(action_slope, 1e-6)
    for a in range(n):
        for b in range(a + 1, n):
            if winding_mode == "delta":
                dk = abs(depths[a] - depths[b])
                alpha_eff = alpha
            elif winding_mode == "max":
                dk = max(depths[a], depths[b])
                alpha_eff = 0.35 * alpha
            else:
                raise ValueError(f"Unknown winding_mode: {winding_mode}")
            amp = transport_strength * np.exp(-alpha_eff * dk)
            amp *= np.cos(phase_per_pass * dk)
            if depth_cost_mode == "diag_only":
                amp = 0.0
            h[a, b] = -amp
            h[b, a] = -amp
    return h


def compute_knotted_lepton_spectrum(
    depths: Sequence[int] = (1, 3, 5),
    l: int = 1,
    n_points: int = 64,
    rs: float = 1.0,
    r_outer: float = 50.0,
    phase_per_pass: float = 0.0,
    transport_strength: float = 4.0,
    resistance_model: str = "writhe",
    resistance_scale: float = 0.5,
    crossings: Iterable[Crossing | dict] | None = None,
    mode_selection: str = "depth_index",
    hard_pinhole_gamma: float = 0.0,
    hard_pinhole_depths: Sequence[int] = (3, 5),
    action_base: float = 0.0,
    action_slope: float = 0.5,
    depth_cost_mode: str = "both",
    winding_mode: str = "delta",
    hierarchy_block: bool = True,
    k_uplift_beta: float = 0.0,
) -> dict[int, float]:
    """Compute depth spectrum from instanton transition matrices.

    Parameters kept for API compatibility (`l`, `n_points`, `rs`, `r_outer`) are
    unused in this instanton formulation.
    """
    del l, n_points, rs, r_outer
    xs = _as_crossings(crossings)

    if hierarchy_block and len(depths) >= 2:
        dlist = tuple(depths)
        h = _build_generation_block(
            depths=dlist,
            phase_per_pass=phase_per_pass,
            transport_strength=transport_strength,
            resistance_model=resistance_model,
            resistance_scale=resistance_scale,
            hard_pinhole_gamma=hard_pinhole_gamma,
            action_base=action_base,
            action_slope=action_slope,
            depth_cost_mode=depth_cost_mode,
            winding_mode=winding_mode,
            k_uplift_beta=k_uplift_beta,
        )
        w = np.sort(eigh(h, eigvals_only=True))
        positive = w[w > 0]
        if positive.size < len(dlist):
            positive = np.pad(positive, (0, len(dlist) - positive.size), mode="edge")
        return {k: float(positive[idx]) for idx, k in enumerate(sorted(dlist))}

    out: dict[int, float] = {}
    for k in depths:
        if k < 1:
            raise ValueError("depth values must be >= 1")

        h = _build_instanton_matrix(
            k=k,
            phase_per_pass=phase_per_pass,
            transport_strength=transport_strength,
            resistance_model=resistance_model,
            resistance_scale=resistance_scale,
            action_base=action_base,
            action_slope=action_slope,
            hard_pinhole_gamma=hard_pinhole_gamma,
            hard_pinhole_depths=hard_pinhole_depths,
            crossings=xs,
            depth_cost_mode=depth_cost_mode,
            winding_mode=winding_mode,
        )
        w = np.sort(eigh(h, eigvals_only=True))
        positive = w[w > 0]
        if positive.size == 0:
            raise RuntimeError("No positive eigenvalues found")

        if mode_selection == "ground":
            idx = 0
        elif mode_selection == "depth_index":
            idx = min(k - 1, positive.size - 1)
        else:
            raise ValueError(f"Unknown mode selection: {mode_selection}")
        out[k] = float(positive[idx])

    return out


def compute_tunneling_envelope(
    depth: int,
    *,
    phase_per_pass: float = 0.0,
    resistance_model: str = "writhe",
    resistance_scale: float = 0.5,
    action_base: float = 0.0,
    winding_mode: str = "delta",
) -> np.ndarray:
    """Return matrix of effective instanton amplitudes exp(-S_ij)."""
    if depth < 2:
        return np.zeros((depth, depth), dtype=float)
    l_core = _default_core_scale()
    if action_base <= 0.0:
        action_base = float(np.pi * (R_MID / l_core))
    resistance = _resistance_sequence(depth, resistance_model, resistance_scale)
    amps = np.zeros((depth, depth), dtype=float)
    for i in range(depth):
        for j in range(i + 1, depth):
            s_ij = _instanton_action(
                i=i,
                j=j,
                k=depth,
                phase_per_pass=phase_per_pass,
                action_base=action_base,
                action_slope=0.0,
                resistance=resistance,
                winding_mode=winding_mode,
            )
            a = float(np.exp(-s_ij))
            amps[i, j] = a
            amps[j, i] = a
    return amps


def calibrate_electron_predict_heavier(
    *,
    depths: Sequence[int] = (1, 3, 5),
    electron_mass_mev: float = 0.51099895,
    observed_mev: dict[int, float] | None = None,
    **kwargs,
) -> LadderFit:
    """Calibrate only to electron depth (1), predict higher depths."""
    if 1 not in depths:
        raise ValueError("depths must include 1 for electron calibration")

    observed = observed_mev or {1: 0.51099895, 3: 105.6583755, 5: 1776.86}
    eig = compute_knotted_lepton_spectrum(depths=depths, **kwargs)
    scale = electron_mass_mev / eig[1]
    predicted = {k: scale * v for k, v in eig.items()}

    rel = {
        k: abs(predicted[k] - target) / target
        for k, target in observed.items()
        if k in predicted and target > 0
    }

    return LadderFit(
        scale_mev=scale,
        predicted_mev=predicted,
        relative_error=rel,
        phase_per_pass=float(kwargs.get("phase_per_pass", 0.0)),
        resistance_scale=float(kwargs.get("resistance_scale", 0.5)),
        transport_strength=float(kwargs.get("transport_strength", 4.0)),
    )


def tune_transport_and_resistance(
    *,
    phase_grid: Sequence[float] = tuple(np.linspace(0.0, np.pi, 13)),
    resistance_grid: Sequence[float] = tuple(np.linspace(0.1, 2.0, 20)),
    transport_strength: float = 4.0,
    depths: Sequence[int] = (1, 3, 5),
    **kwargs,
) -> LadderFit:
    """Grid-search phase/resistance after electron-only calibration."""
    best: LadderFit | None = None
    best_obj = float("inf")

    for phi in phase_grid:
        for kappa in resistance_grid:
            fit = calibrate_electron_predict_heavier(
                depths=depths,
                phase_per_pass=float(phi),
                resistance_scale=float(kappa),
                transport_strength=transport_strength,
                **kwargs,
            )
            obj = fit.relative_error.get(3, 1.0) + fit.relative_error.get(5, 1.0)
            if obj < best_obj:
                best_obj = obj
                best = fit

    if best is None:
        raise RuntimeError("Tuning grid produced no candidate fit")
    return best
