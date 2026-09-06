"""Finite-width pointer coupled to complete classical rotor histories.

Public freeze: docs/pointer_spread_prereg.md, baf856d. The rotor coupling,
positive closure weight and preparation semantics are explicit additional
assumptions. This module supplies no BAM field or intervention derivation.
"""

from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.special import ndtr, ndtri, roots_hermitenorm
from scipy.stats import qmc

from geometrodynamics.bulk.source_readout import (
    READOUT_AXIS, FIXED_ANALYZER, FUTURE_ANALYZERS, gaussian_tail,
)


PUBLIC_PREREG = "baf856df06b26de6c9750697ba5df5ad521000c5"
SPREADS = (0., .01, .025, .05, .1, .25, .5)
SEEDS = (20260906, 20260907, 20260908, 20260909)


@dataclass(frozen=True)
class PointerModel:
    source_sigma: float = .1
    beta: float = 64.
    duration: float = 1.
    future_time: float = 2.
    strength: float = 1.
    mass: float = 1.
    noise: float = .15
    threshold: float = .6
    axis: tuple = tuple(READOUT_AXIS)

    def __post_init__(self):
        values = (self.source_sigma, self.beta, self.duration, self.future_time,
                  self.strength, self.mass, self.noise, self.threshold)
        if not all(np.isfinite(values)):
            raise ValueError("model parameters must be finite")
        if (self.source_sigma < 0 or self.beta < 0 or self.duration <= 0
                or self.future_time < self.duration or self.mass <= 0
                or self.noise <= 0 or self.threshold < 0):
            raise ValueError("invalid model scales or boundary times")
        axis = np.asarray(self.axis)
        if axis.shape != (3,) or not np.all(np.isfinite(axis)) or np.linalg.norm(axis) == 0:
            raise ValueError("readout axis must be a finite nonzero three-vector")


def qmultiply(left, right):
    """Quaternion product on the last axis, with NumPy broadcasting."""
    left, right = np.asarray(left), np.asarray(right)
    scalar = left[..., :1] * right[..., :1] - np.sum(
        left[..., 1:] * right[..., 1:], axis=-1, keepdims=True)
    vector = (left[..., :1] * right[..., 1:] + right[..., :1] * left[..., 1:]
              + np.cross(left[..., 1:], right[..., 1:]))
    return np.concatenate((scalar, vector), axis=-1)


def minimal_lift(start, end):
    """Short-geodesic lift. Antipodal endpoints have no assigned holonomy."""
    start, end = np.asarray(start), np.asarray(end)
    scalar = 1 + np.sum(start * end, axis=-1, keepdims=True)
    quaternion = np.concatenate((scalar, np.cross(start, end)), axis=-1)
    norm = np.linalg.norm(quaternion, axis=-1, keepdims=True)
    if np.any(norm < 1e-12):
        raise ValueError("antipodal geodesic puncture: lift undefined")
    return quaternion / norm


def geodesic_drift(x, p, transport, time):
    """Exact free Hamiltonian flow and Levi-Civita transport, including winding.

    transport acts by Ad_C. No renormalization/projection is used in the
    flow; sphere and cotangent constraints are measured independently.
    """
    speed2 = np.sum(p * p, axis=-1)
    speed = np.sqrt(speed2)
    angle = speed * time
    cosine = np.cos(angle)[:, None]
    sine_over_speed = (time * np.sinc(angle / math.pi))[:, None]
    rotation = np.column_stack((np.cos(angle / 2),
        (.5 * time * np.sinc(angle / (2 * math.pi)))[:, None] * np.cross(x, p)))
    return (cosine * x + sine_over_speed * p,
            cosine * p - (speed2[:, None] * sine_over_speed) * x,
            qmultiply(rotation, transport))


def evolve(x0, p0, pointer_momentum, model=PointerModel(), steps=64):
    """Finite pulse with reciprocal force; early record precedes closure.

    Strang drift/kick/drift uses midpoint h(t). Each constituent is an exact
    canonical flow. Q0 is integrated later as one common Gaussian kernel.
    """
    if not isinstance(steps, int) or steps < 2:
        raise ValueError("at least two pulse steps required")
    x, p = np.array(x0, dtype=float, copy=True), np.array(p0, dtype=float, copy=True)
    if x.ndim != 2 or x.shape[1] != 3 or p.shape != x.shape:
        raise ValueError("x0 and p0 must have shape (n,3)")
    if not np.allclose(np.sum(x * x, axis=1), 1., atol=1e-10, rtol=0):
        raise ValueError("source position must lie on the unit sphere")
    if np.max(np.abs(np.sum(x * p, axis=1))) > 1e-10:
        raise ValueError("source momentum must be tangent")
    P = np.broadcast_to(np.asarray(pointer_momentum, dtype=float), (len(x),))
    if not np.all(np.isfinite(p)) or not np.all(np.isfinite(P)):
        raise ValueError("momenta must be finite")
    m = np.asarray(model.axis) / np.linalg.norm(model.axis)
    C = np.tile([1., 0., 0., 0.], (len(x), 1))
    record = model.duration * P / model.mass
    dt = model.duration / steps
    for k in range(steps):
        x, p, C = geodesic_drift(x, p, C, .5 * dt)
        h = (2 * model.strength / model.duration
             * math.sin(math.pi * (k + .5) / steps) ** 2)
        f = x @ m
        p -= (dt * h * P)[:, None] * (m - f[:, None] * x)
        record = record + dt * h * f
        x, p, C = geodesic_drift(x, p, C, .5 * dt)
    x_read, p_read = x.copy(), p.copy()
    x, p, C = geodesic_drift(x, p, C, model.future_time - model.duration)
    diagnostics = {
        "sphere_residual": float(max(np.max(np.abs(np.sum(z * z, axis=1) - 1))
                                      for z in (x_read, x))),
        "tangency_residual": float(max(np.max(np.abs(np.sum(z * v, axis=1)))
                                        for z, v in ((x_read, p_read), (x, p)))),
        "transport_norm_residual": float(np.max(np.abs(np.sum(C * C, axis=1) - 1))),
        "axial_momentum_residual": float(np.max(np.abs(
            np.cross(x, p) @ m - np.cross(x0, p0) @ m))),
    }
    return {"x0": np.asarray(x0), "p0": np.asarray(p0), "x_read": x_read,
            "p_read": p_read, "x_final": x, "p_final": p, "transport": C,
            "record_shift": record, "P": P, "diagnostics": diagnostics}


def free_path(x0, p0, future_time=2.):
    C = np.tile([1., 0., 0., 0.], (len(x0), 1))
    x, p, C = geodesic_drift(x0, p0, C, future_time)
    return {"x0": x0, "x_final": x, "p_final": p, "transport": C}


def history_weights(paths, a, b, beta=64.):
    """Score the complete x0--source path--x_B--u--w--x0 itinerary.

    The actual source transport is retained. The positive frame score is a
    chosen history law, not a real-time probability derived from H.
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    x0, end, C = paths["x0"], paths["x_final"], paths["transport"]
    energies, axis_error = [], 0.
    for sa in (-1, 1):
        for sb in (-1, 1):
            u, w = sa * a, -sb * b
            G = qmultiply(minimal_lift(w, x0), qmultiply(
                minimal_lift(u, w), qmultiply(minimal_lift(end, u), C)))
            norm2 = np.sum(G * G, axis=1)
            energy = .5 * np.sum(G[:, 1:] ** 2, axis=1) / norm2
            energies.append(energy)
            axis_error = max(axis_error, float(np.max(np.linalg.norm(
                np.cross(G[:, 1:], x0), axis=1) / np.sqrt(norm2))))
    energies = np.column_stack(energies)
    return {"weight": np.mean(np.exp(-beta * energies), axis=1),
            "energies": energies, "loop_axis_residual": axis_error}


@lru_cache(maxsize=24)
def source_nodes(power=11, seed=20260906, source_sigma=.1):
    """Four source coordinates, with explicit antipodal phase-space partners."""
    u = qmc.Sobol(4, scramble=True, seed=seed).random_base2(power)
    z, phi = 2 * u[:, 0] - 1, 2 * math.pi * u[:, 1]
    radius = np.sqrt(1 - z * z)
    x = np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))
    etheta = np.column_stack((z * np.cos(phi), z * np.sin(phi), -radius))
    ephi = np.column_stack((-np.sin(phi), np.cos(phi), np.zeros(len(phi))))
    # Clipping protects the inverse-normal transform at representable endpoints.
    normal = ndtri(np.clip(u[:, 2:], np.finfo(float).eps, 1 - np.finfo(float).eps))
    p = source_sigma * (normal[:, :1] * etheta + normal[:, 1:] * ephi)
    x, p = np.concatenate((x, -x)), np.concatenate((p, -p))
    x.setflags(write=False)
    p.setflags(write=False)
    return x, p


def pointer_nodes(sigma, count=16):
    if not np.isfinite(sigma) or sigma < 0:
        raise ValueError("pointer width must be finite and nonnegative")
    if sigma == 0:
        return np.array([0.]), np.array([1.])
    nodes, weights = roots_hermitenorm(count)
    return sigma * nodes, weights / math.sqrt(2 * math.pi)


def normalized_measures(weight, pointer_prior, free_weight):
    """Two explicit preparation semantics plus the frozen-posterior control.

    weight has shape (nP,nSource). Each row is a conditional source measure
    at a fixed externally prepared P. No empirical tuning enters its weights.
    """
    weight, prior = np.asarray(weight), np.asarray(pointer_prior)
    row_sum = weight.sum(axis=1)
    if np.any(row_sum <= 0) or np.any(weight < 0):
        raise ValueError("positive resolved history weights required")
    fixed = prior[:, None] * weight / row_sum[:, None]
    joint = prior[:, None] * weight
    joint /= joint.sum()
    frozen = prior[:, None] * np.asarray(free_weight)[None, :] / np.sum(free_weight)
    ess = row_sum ** 2 / np.sum(weight ** 2, axis=1)
    return {"fixed": fixed, "joint": joint, "frozen": frozen,
            "Z_by_P": row_sum / weight.shape[1], "ess_by_P": ess}


def record_law(shift, P, weights, model=PointerModel()):
    w = np.asarray(weights).ravel()
    shift, P = np.asarray(shift).ravel(), np.asarray(P).ravel()
    mean, pmean = float(w @ shift), float(w @ P)
    return {"tail": float(w @ gaussian_tail(shift, model.threshold, model.noise)),
            "mean": mean, "variance": float(w @ ((shift - mean) ** 2)) + model.noise ** 2,
            "positive_probability": float(w @ ndtr(shift / model.noise)),
            "P_mean": pmean, "P_variance": float(w @ ((P - pmean) ** 2))}


def simulate_spread(sigma, *, model=PointerModel(), power=11, seed=20260906,
                    hermite=16, steps=64):
    """One deterministic numerical replicate, with all three laws reported."""
    source_x, source_p = source_nodes(power, seed, model.source_sigma)
    pn, pw = pointer_nodes(sigma, hermite)
    ns = len(source_x)
    x0, p0 = np.tile(source_x, (len(pn), 1)), np.tile(source_p, (len(pn), 1))
    P = np.repeat(pn, ns)
    paths = evolve(x0, p0, P, model, steps)
    free = free_path(source_x, source_p, model.future_time)
    identity = np.tile([1., 0., 0., 0.], (ns, 1))
    xr, pr, _ = geodesic_drift(source_x, source_p, identity, model.duration)
    dx = paths["x_read"].reshape(len(pn), ns, 3) - xr
    dp = paths["p_read"].reshape(len(pn), ns, 3) - pr
    diagnostics = dict(paths["diagnostics"])
    diagnostics.update({
        "reference_position_recoil_rms": float(np.sqrt(pw @ np.mean(np.sum(dx * dx, axis=2), axis=1))),
        "reference_momentum_recoil_rms": float(np.sqrt(pw @ np.mean(np.sum(dp * dp, axis=2), axis=1))),
        "record_antipodal_residual": float(np.max(np.abs(
            paths["record_shift"].reshape(len(pn), ns)
            + np.roll(paths["record_shift"].reshape(len(pn), ns)[::-1], ns // 2, axis=1))))})
    choices = []
    for choice, b in enumerate(FUTURE_ANALYZERS):
        scored = history_weights(paths, FIXED_ANALYZER, b, model.beta)
        weight = scored["weight"].reshape(len(pn), ns)
        free_weight = history_weights(free, FIXED_ANALYZER, b, model.beta)["weight"]
        measures = normalized_measures(weight, pw, free_weight)
        row = {"choice": choice, "Z_by_P": measures["Z_by_P"].tolist(),
               "min_effective_source_samples": float(np.min(measures["ess_by_P"]))}
        for kind in ("fixed", "joint", "frozen"):
            row[kind] = record_law(paths["record_shift"], P, measures[kind], model)
        choices.append(row)
        diagnostics["loop_axis_residual"] = max(diagnostics.get("loop_axis_residual", 0.), scored["loop_axis_residual"])
        diagnostics["weight_parity_residual"] = max(diagnostics.get("weight_parity_residual", 0.),
            float(np.max(np.abs(weight - np.roll(weight[::-1], ns // 2, axis=1)))))
    return {"sigma_pointer": sigma, "source_sigma": model.source_sigma, "seed": seed,
            "source_power": power, "source_samples": ns, "hermite_nodes": len(pn),
            "steps": steps, "pointer_nodes": pn.tolist(), "pointer_prior": pw.tolist(),
            "choices": choices, "diagnostics": diagnostics,
            "contrasts": {kind: choices[1][kind]["tail"] - choices[0][kind]["tail"]
                          for kind in ("fixed", "joint", "frozen")}}
