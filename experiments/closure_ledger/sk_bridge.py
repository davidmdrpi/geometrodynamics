"""
S(k) bridge for the closure-phase ledger (Layer 2 wiring).

Implements the candidate maps listed in `blockers.py` and computes the
per-generation radial bulk phase

    Φ_radial(k) = Σ_(l, n) ∈ S(k)  ∫ √max(ω(l, n)² − V_eff(r, l), 0) dr*

where V_eff(r, l) is the 5D Tangherlini effective potential supplied by
`geometrodynamics.tangherlini.radial.V_tangherlini`, ω(l, n) is the
n-th eigenfrequency of the Chebyshev radial solver at angular index l,
and the integration runs over the classically-allowed region of the
tortoise grid.

Wired candidates:
    A_lowest_radial_per_l   — S(k) = {(l, 0) : l = 1, 3, ..., k}
    B_fixed_total_quantum_number
                            — S(k) = {(l, (k−l)/2) : l = 1, 3, ..., k}
                              (uses excited modes; relies on the solver
                              returning at least n+1 modes for each l)

`C_closure_coherent_superposition` is intentionally not wired: it
requires the S³ closure-phase condition on ω(l, n) to be defined
explicitly, which is still a thesis-level open question.

`none` returns an empty membership; the radial channel falls back to
the Layer-1 "missing" status, reproducing the pre-bridge ledger.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


SK_CANDIDATES = (
    "A_lowest_radial_per_l",
    "B_fixed_total_quantum_number",
    "none",
)
DEFAULT_SK_CANDIDATE = "A_lowest_radial_per_l"


@dataclass(frozen=True)
class Mode:
    l: int
    n: int


def s_k_membership(k: int, candidate: str) -> list[Mode]:
    """Return the (l, n) modes that constitute S(k) under the chosen candidate."""
    if candidate == "none":
        return []
    if candidate == "A_lowest_radial_per_l":
        # Sum of odd-l radial ground states up to angular harmonic l = k.
        return [Mode(l=ll, n=0) for ll in range(1, k + 1, 2)]
    if candidate == "B_fixed_total_quantum_number":
        # 2n + l = k, l odd ≥ 1, n ≥ 0  ⟹  l ∈ {1, 3, ..., k}, n = (k−l)/2.
        return [Mode(l=ll, n=(k - ll) // 2) for ll in range(1, k + 1, 2)]
    if candidate == "C_closure_coherent_superposition":
        raise NotImplementedError(
            "C_closure_coherent_superposition requires the S³ closure-phase "
            "condition on ω(l, n) to be defined; still a thesis-level open "
            "question. Pass candidate='none' or a wired candidate."
        )
    raise ValueError(f"Unknown S(k) candidate: {candidate}")


@dataclass
class ModePhase:
    l: int
    n: int
    omega: Optional[float]
    phi: Optional[float]
    status: str           # "computed" | "mode_unavailable" | "import_failed"
    error: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "l": self.l,
            "n": self.n,
            "omega": self.omega,
            "phi": self.phi,
            "status": self.status,
            "error": self.error,
        }


def phi_radial_for_mode(l: int, n: int, *, N: int = 80) -> ModePhase:
    """
    Compute Φ(l, n) by integrating √max(ω² − V_eff, 0) over the tortoise grid.

    Returns a ModePhase with status:
        "computed"          — phi is the integrated radial action.
        "mode_unavailable"  — solver returned fewer than n+1 modes.
        "import_failed"     — Tangherlini infrastructure not importable.
    """
    try:
        import numpy as np
        from geometrodynamics.tangherlini.radial import (
            solve_radial_modes,
            r_to_rstar,
            V_tangherlini,
        )
        from geometrodynamics.constants import R_MID
    except ImportError as exc:
        return ModePhase(
            l=l, n=n, omega=None, phi=None,
            status="import_failed", error=str(exc),
        )

    omegas, _funcs, r_grid = solve_radial_modes(l, N=N, n_modes=n + 1)
    if len(omegas) <= n:
        return ModePhase(
            l=l, n=n, omega=None, phi=None,
            status="mode_unavailable",
            error=(
                f"solver returned {len(omegas)} modes for l={l}; "
                f"need at least n+1={n + 1}"
            ),
        )

    omega_n = float(omegas[n])
    rs = float(R_MID)
    Vg = V_tangherlini(r_grid, l, rs)
    integrand = np.sqrt(np.maximum(omega_n ** 2 - Vg, 0.0))
    rstar = np.array([r_to_rstar(float(r), rs) for r in r_grid])
    order = np.argsort(rstar)
    phi = float(np.trapezoid(integrand[order], rstar[order]))
    return ModePhase(l=l, n=n, omega=omega_n, phi=phi, status="computed")


@dataclass
class RadialPhaseResult:
    candidate: str
    k: int
    modes: list[ModePhase]
    total_phi: Optional[float]
    status: str            # "computed" | "partial" | "missing"

    def to_dict(self) -> dict:
        return {
            "candidate": self.candidate,
            "k": self.k,
            "modes": [m.to_dict() for m in self.modes],
            "total_phi": self.total_phi,
            "status": self.status,
        }


def phi_radial_from_sk(
    k: int, candidate: str, *, N: int = 80,
) -> RadialPhaseResult:
    """
    Compute Φ_radial(k) = Σ_(l,n)∈S(k) Φ(l, n).

    Returns a RadialPhaseResult with overall status:
        "computed" — every mode integral resolved.
        "partial"  — at least one mode resolved, at least one failed.
        "missing"  — no modes resolved (or candidate == "none").
    """
    if candidate == "none":
        return RadialPhaseResult(
            candidate=candidate, k=k, modes=[],
            total_phi=None, status="missing",
        )
    members = s_k_membership(k, candidate)
    mode_results = [phi_radial_for_mode(m.l, m.n, N=N) for m in members]
    if mode_results and all(mr.status == "computed" for mr in mode_results):
        total = sum(mr.phi for mr in mode_results)
        return RadialPhaseResult(
            candidate=candidate, k=k, modes=mode_results,
            total_phi=float(total), status="computed",
        )
    if any(mr.status == "computed" for mr in mode_results):
        return RadialPhaseResult(
            candidate=candidate, k=k, modes=mode_results,
            total_phi=None, status="partial",
        )
    return RadialPhaseResult(
        candidate=candidate, k=k, modes=mode_results,
        total_phi=None, status="missing",
    )
