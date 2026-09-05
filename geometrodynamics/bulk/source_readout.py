"""Conditional source-record laws and a classical canonical pointer.

Preregistered in docs/source_readout_prereg.md. These synthetic functions
of the triangle coordinate are not a BAM field observable or intervention.
"""

import math

import numpy as np
from scipy.special import ndtr, roots_legendre

from geometrodynamics.bulk.closure_equilibrium import rotor_potential


READOUT_AXIS = np.array([1., 2., 0.]) / math.sqrt(5.)
FIXED_ANALYZER = np.array([0., 0., 1.])
FUTURE_ANALYZERS = (np.array([1., 0., 0.]), np.array([0., 1., 0.]))


def _unit(v):
    v = np.asarray(v, dtype=float)
    if v.shape != (3,) or not np.all(np.isfinite(v)) or np.linalg.norm(v) == 0:
        raise ValueError("expected a finite nonzero three-vector")
    return v / np.linalg.norm(v)


def _basis(a, b):
    a, b = _unit(a), _unit(b)
    normal = np.cross(a, b)
    if np.linalg.norm(normal) < 1e-10:
        raise ValueError("analyzers must be non-collinear")
    normal = _unit(normal)
    return a, b, np.cross(normal, a), normal


def source_circle(a, b, n=32768):
    """Positive, equal-sector, sharp phase-coarea law on the closure circle.

    Outcome summation makes relabelling the partner sign immaterial. The
    common cross-product denominator cancels in the normalized weights.
    """
    a, b, tangent, _ = _basis(a, b)
    phi = 2 * math.pi * (np.arange(n) + .5) / n
    x = np.cos(phi)[:, None] * a + np.sin(phi)[:, None] * tangent
    density = np.zeros(n)
    for sa in (-1, 1):
        for sb in (-1, 1):
            u, w = sa * a, -sb * b
            density += np.abs(1 + u @ w + x @ (u + w))
    return x, density / density.sum()


def source_sphere(a, b, beta=64., n_normal=128, n_azimuth=512):
    """Independent whole-sphere quadrature of the outcome-summed Gibbs law.

    Resolves both sides of the closure plane with Gauss-Legendre nodes;
    no sharp-closure formula or sector-normalized mixture is substituted.
    """
    if not np.isfinite(beta) or beta < 0:
        raise ValueError("beta must be finite and nonnegative")
    a, b, tangent, normal = _basis(a, b)
    z, wz = roots_legendre(n_normal)
    phi = 2 * math.pi * (np.arange(n_azimuth) + .5) / n_azimuth
    circle = np.cos(phi)[:, None] * a + np.sin(phi)[:, None] * tangent
    x = (np.sqrt(1 - z[:, None, None] ** 2) * circle[None, :, :]
         + z[:, None, None] * normal).reshape(-1, 3)
    density = np.zeros(len(x))
    for sa in (-1, 1):
        for sb in (-1, 1):
            density += np.exp(-beta * rotor_potential(x, sa * a, -sb * b))
    weights = density * np.repeat(wz, n_azimuth)
    return x, weights / weights.sum()


def gaussian_tail(values, threshold=.6, noise=.15):
    """P(|F+Z|>threshold), with the SAME independent noise at either setting."""
    if noise < 0 or threshold < 0:
        raise ValueError("noise and threshold must be nonnegative")
    values = np.asarray(values)
    if noise == 0:
        return (np.abs(values) > threshold).astype(float)
    return ndtr((-threshold - values) / noise) + ndtr((values - threshold) / noise)


def record_statistics(x, weights, axis=READOUT_AXIS, noise=0., threshold=.6):
    """Full-law diagnostics, without asserting an operational BAM record."""
    f = np.asarray(x) @ _unit(axis)
    w = np.asarray(weights)
    mean = float(w @ f)
    variance = float(w @ (f - mean) ** 2)
    positive = ndtr(f / noise) if noise > 0 else (f > 0)
    return {"mean": mean, "variance": variance,
            "record_variance": variance + noise ** 2,
            "positive_probability": float(w @ positive),
            "tail_probability": float(w @ gaussian_tail(f, threshold, noise))}


def orthogonal_noiseless_tail(threshold=.6):
    """Independent antiderivative of (1+|sin phi|+|cos phi|)/(2pi+8)."""
    tails = []
    for amplitude in (1 / math.sqrt(5), 2 / math.sqrt(5)):
        v = threshold / amplitude
        if v >= 1:
            tails.append(0.)
        elif v <= 0:
            tails.append(1.)
        else:
            tails.append(2 * (math.pi / 2 - math.asin(v)
                             + math.sqrt(1 - v * v) + 1 - v) / (math.pi + 4))
    return tails


def chart_observable(chi, phi, axis=READOUT_AXIS):
    """F=m.x and its coordinate derivatives on a regular polar chart."""
    m = _unit(axis)
    x = np.array([np.sin(chi) * np.cos(phi), np.sin(chi) * np.sin(phi), np.cos(chi)])
    dchi = np.array([np.cos(chi) * np.cos(phi), np.cos(chi) * np.sin(phi), -np.sin(chi)])
    dphi = np.array([-np.sin(chi) * np.sin(phi), np.sin(chi) * np.cos(phi), 0.])
    return m @ x, m @ dchi, m @ dphi


def pointer_kick(state, strength=1., axis=READOUT_AXIS):
    """Exact flow of g P F in (chi,p_chi,phi,p_phi,Q,P), including recoil.

    P=0 preserves the source phase-space state. This ideal preparation and
    interaction are added classical assumptions, not derived BAM equipment.
    """
    out = np.array(state, dtype=float, copy=True)
    f, dchi, dphi = chart_observable(out[0], out[2], axis)
    out[1] -= strength * out[5] * dchi
    out[3] -= strength * out[5] * dphi
    out[4] += strength * f
    return out
