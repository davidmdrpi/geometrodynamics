"""
SU(3) color algebra for the flux-tube network.

Provides color/anti-color bookkeeping, singlet checks, generator
matrices, and gluon wavefront projection amplitudes.
"""

from __future__ import annotations

from math import sqrt
from typing import TYPE_CHECKING, List

import numpy as np

if TYPE_CHECKING:
    from geometrodynamics.qcd.network import Branch

COLORS = {"r", "g", "b"}
ANTI_COLORS = {"r̄", "ḡ", "b̄"}
ANTI_OF = {"r": "r̄", "g": "ḡ", "b": "b̄", "r̄": "r", "ḡ": "g", "b̄": "b"}
BASE_OF = {"r": "r", "g": "g", "b": "b", "r̄": "r", "ḡ": "g", "b̄": "b"}
IDX_OF = {"r": 0, "g": 1, "b": 2, "r̄": 0, "ḡ": 1, "b̄": 2}


def is_singlet(colors: List[str]) -> bool:
    """Check whether a list of color labels forms a color singlet."""
    if not colors:
        return True
    charge = {"r": 0, "g": 0, "b": 0}
    for c in colors:
        if c in COLORS:
            charge[c] += 1
        elif c in ANTI_COLORS:
            charge[BASE_OF[c]] -= 1
        else:
            raise ValueError(f"Unknown color: {c!r}")
    if all(v == 0 for v in charge.values()):
        return True
    if set(colors) in ({"r", "g", "b"}, {"r̄", "ḡ", "b̄"}):
        return True
    return False


def color_charge_vector(colors: List[str]) -> np.ndarray:
    """Net color charge as a 3-vector (r, g, b components)."""
    q = np.zeros(3)
    for c in colors:
        if c in COLORS:
            q[IDX_OF[c]] += 1.0
        elif c in ANTI_COLORS:
            q[IDX_OF[c]] -= 1.0
    return q


def gluon_generator_matrix(k: int) -> np.ndarray:
    """Return the k-th SU(3) generator (Gell-Mann basis, k=0..7)."""
    G = np.zeros((3, 3), dtype=complex)
    pairs = [(0, 1), (1, 0), (0, 2), (2, 0), (1, 2), (2, 1)]
    if k < 6:
        i, j = pairs[k]
        G[i, j] = 1.0
    elif k == 6:
        G[0, 0] = +1 / sqrt(2)
        G[1, 1] = -1 / sqrt(2)
    elif k == 7:
        G[0, 0] = G[1, 1] = +1 / sqrt(6)
        G[2, 2] = -2 / sqrt(6)
    return G


def generator_in_rep(gen_index: int, rep: str) -> np.ndarray:
    """Generator in fundamental or anti-fundamental representation."""
    G = gluon_generator_matrix(gen_index)
    return G if rep == "fund" else -np.conj(G)


def gluon_wavefront_amplitude(branch: "Branch", gen_index: int) -> float:
    """Projection amplitude for generator k on a given branch."""
    Ts = generator_in_rep(gen_index, branch.rep_src)
    Td = generator_in_rep(gen_index, branch.rep_dst)
    i = IDX_OF.get(branch.color_pair[0], 0)
    j = IDX_OF.get(branch.color_pair[1], 0)
    if gen_index >= 6:
        amp = abs(float(Ts[i, i].real) - float(Td[j, j].real))
    else:
        amp = abs(Ts[i, j]) + abs(Td[i, j])
    return float(amp)


def seed_gluon_wavefront(
    branch: "Branch", gen_index: int, amplitude: float = 1.0
) -> np.ndarray:
    """Gaussian gluon wavefront seeded at the branch midpoint."""
    proj = gluon_wavefront_amplitude(branch, gen_index)
    N = len(branch.psi) if branch.psi is not None else 100
    s = np.linspace(0, branch.length, N)
    mid, sig = branch.length / 2.0, branch.length / 8.0
    return amplitude * proj * np.exp(-0.5 * ((s - mid) / sig) ** 2)
