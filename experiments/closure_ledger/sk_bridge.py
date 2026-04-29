"""
S(k) bridge for the closure-phase ledger (Layer 2 wiring).

Implements the candidate maps listed in `blockers.py` and computes the
per-generation radial bulk phase

    Φ_radial(k) = Σ_(l, n) ∈ S(k)  ω(l, n)-weighted WKB radial action
                ≡ Σ_(l, n) ∈ S(k)  ∫ √max(ω(l, n)² − V_eff(r, l), 0) dr*

where V_eff(r, l) is the 5D Tangherlini effective potential supplied by
`geometrodynamics.tangherlini.radial.V_tangherlini`, ω(l, n) is the
n-th eigenfrequency of the Chebyshev radial solver at angular index l,
and the integration runs over the classically-allowed region of the
tortoise grid.

This is a **WKB radial-action convention** — not "the" radial phase.
The convention is one-sided (the integration covers the entire
classically-allowed region of the tortoise grid; for symmetric soft
turning points it equals 2 ∫_{r₁}^{r₂} p dr*), and applies no Maslov
correction by default. The Maslov shift can be added via the
`maslov_correction` parameter to `phi_radial_for_mode`. Candidate
results are sensitive to both the convention and the correction; see
`phi_convergence_table` for grid-stability evidence and the test suite
for the Maslov-additivity check.

Wired candidates:
    A_lowest_radial_per_l       — S(k) = {(l, 0) : l = 1, 3, ..., k}
    B1_single_angular_mode      — S(k) = {(l = k, n = 0)}
    B2_single_radial_excitation — S(k) = {(l = 1, n = (k − 1) / 2)}

`C_eigenvector_weighted` is intentionally not wired: it requires
deriving S(k) from the lepton instanton-surrogate eigenvectors / the
generation block of the Hamiltonian rather than imposing it by hand.
That is the principled bridge from the surrogate to the Tangherlini
radial operator and is still a thesis-level open question.

`none` returns an empty membership; the radial channel falls back to
the Layer-1 "missing" status, reproducing the pre-bridge ledger.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


WKB_CONVENTION = "WKB_one_sided_radial_action_no_maslov"

SK_CANDIDATES = (
    "A_lowest_radial_per_l",
    "B1_single_angular_mode",
    "B2_single_radial_excitation",
    "none",
)
DEFAULT_SK_CANDIDATE = "A_lowest_radial_per_l"

# All wired candidates in the order they were proposed (used by the
# comparison runner).
WIRED_CANDIDATES = (
    "A_lowest_radial_per_l",
    "B1_single_angular_mode",
    "B2_single_radial_excitation",
)


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
    if candidate == "B1_single_angular_mode":
        # Single angular mode per generation: l = k, n = 0.
        return [Mode(l=k, n=0)]
    if candidate == "B2_single_radial_excitation":
        # Single radial excitation in the l = 1 ladder: n = (k − 1) / 2.
        if (k - 1) % 2 != 0:
            raise ValueError(
                f"B2_single_radial_excitation requires odd k; got k={k}"
            )
        return [Mode(l=1, n=(k - 1) // 2)]
    if candidate == "C_eigenvector_weighted":
        raise NotImplementedError(
            "C_eigenvector_weighted requires the lepton instanton-surrogate "
            "eigenvectors (or generation-block Hamiltonian) to define S(k); "
            "still a thesis-level open question. Pass candidate='none' or "
            "a wired candidate."
        )
    raise ValueError(f"Unknown S(k) candidate: {candidate}")


@dataclass
class ModePhase:
    l: int
    n: int
    omega: Optional[float]
    phi: Optional[float]
    status: str           # "computed" | "mode_unavailable" | "import_failed"
    convention: str = WKB_CONVENTION
    maslov_correction: float = 0.0
    N_grid: int = 0
    error: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "l": self.l,
            "n": self.n,
            "omega": self.omega,
            "phi": self.phi,
            "status": self.status,
            "convention": self.convention,
            "maslov_correction": self.maslov_correction,
            "N_grid": self.N_grid,
            "error": self.error,
        }


def phi_radial_for_mode(
    l: int,
    n: int,
    *,
    N: int = 80,
    maslov_correction: float = 0.0,
) -> ModePhase:
    """
    Compute Φ(l, n) by integrating √max(ω² − V_eff, 0) over the tortoise grid.

    Convention: one-sided WKB radial action, no Maslov correction by
    default. Pass `maslov_correction` (e.g. π/4 for one soft turning
    point, π/2 for two soft turning points) to add a constant offset to
    the returned phi.

    Returns a ModePhase with status:
        "computed"          — phi is the integrated radial action + Maslov.
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
            maslov_correction=maslov_correction, N_grid=N,
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
            maslov_correction=maslov_correction, N_grid=N,
        )

    omega_n = float(omegas[n])
    rs = float(R_MID)
    Vg = V_tangherlini(r_grid, l, rs)
    integrand = np.sqrt(np.maximum(omega_n ** 2 - Vg, 0.0))
    rstar = np.array([r_to_rstar(float(r), rs) for r in r_grid])
    order = np.argsort(rstar)
    phi_wkb = float(np.trapezoid(integrand[order], rstar[order]))
    return ModePhase(
        l=l, n=n, omega=omega_n,
        phi=phi_wkb + maslov_correction,
        status="computed",
        maslov_correction=maslov_correction, N_grid=N,
    )


def phi_convergence_table(
    l: int, n: int, *, Ns: tuple[int, ...] = (40, 60, 80, 100, 120),
) -> list[dict]:
    """
    Tabulate Φ(l, n) at a sequence of Chebyshev grid sizes.

    Used to confirm that the WKB radial-action integral is grid-stable.
    Returns a list of {"N": int, "phi": float, "status": str} dicts.
    """
    rows: list[dict] = []
    for N in Ns:
        mp = phi_radial_for_mode(l, n, N=N)
        rows.append({"N": N, "phi": mp.phi, "status": mp.status})
    return rows


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
    k: int, candidate: str, *,
    N: int = 80, maslov_correction: float = 0.0,
) -> RadialPhaseResult:
    """
    Compute Φ_radial(k) = Σ_(l,n)∈S(k) Φ(l, n).

    `maslov_correction` is added per mode (so a candidate with |S(k)|
    modes accumulates |S(k)| · maslov_correction). Default 0.

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
    mode_results = [
        phi_radial_for_mode(m.l, m.n, N=N, maslov_correction=maslov_correction)
        for m in members
    ]
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
