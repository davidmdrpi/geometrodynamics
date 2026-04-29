"""
S(k) bridge for the closure-phase ledger (Layer 2 wiring).

Implements the candidate maps listed in `blockers.py` and computes the
per-generation radial bulk phase

    Φ_radial(k) = Σ_(l, n) ∈ S(k)  w(l, n; k) · ∫ √max(ω² − V_eff, 0) dr*

where V_eff(r, l) is the 5D Tangherlini effective potential supplied by
`geometrodynamics.tangherlini.radial.V_tangherlini`, ω(l, n) is the
n-th eigenfrequency of the Chebyshev radial solver at angular index l,
and w(l, n; k) is a candidate-specific weight (1.0 for hand-imposed
candidates A/B1/B2; |v_species,i|² for the eigenvector-weighted
candidates C1/C2). Integration runs over the classically-allowed
region of the tortoise grid.

This is a **WKB radial-action convention** — not "the" radial phase.
The convention is one-sided (the integration covers the entire
classically-allowed region of the tortoise grid; for symmetric soft
turning points it equals 2 ∫_{r₁}^{r₂} p dr*), and applies no Maslov
correction by default. The Maslov shift can be added per mode via the
`maslov_correction` parameter to `phi_radial_for_mode`. Candidate
results are sensitive to both the convention and the correction.

Wired candidates:
    A_lowest_radial_per_l       — S(k) = {(l, 0) : l = 1, 3, ..., k}
    B1_single_angular_mode      — S(k) = {(l = k, n = 0)}
    B2_single_radial_excitation — S(k) = {(l = 1, n = (k − 1) / 2)}
    C1_eigenvector_weighted_B1  — Σ_i |v_species,i|² Φ(l = k_i, n = 0)
    C2_eigenvector_weighted_B2  — Σ_i |v_species,i|² Φ(l = 1, n = (k_i−1)/2)

C1 and C2 derive their weights from the **locked lepton generation
block** of `geometrodynamics.tangherlini.lepton_spectrum` (no fitted
parameters introduced here). The species-to-eigenvector assignment is
the standard ascending-eigenvalue ordering: smallest → e, next → μ,
largest → τ.

`none` returns an empty membership; the radial channel falls back to
the Layer-1 "missing" status, reproducing the pre-bridge ledger.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Optional


WKB_CONVENTION = "WKB_one_sided_radial_action_no_maslov"

LEPTON_DEPTHS = (1, 3, 5)   # fixed depth basis for C1/C2 eigenvector weights

SK_CANDIDATES = (
    "A_lowest_radial_per_l",
    "B1_single_angular_mode",
    "B2_single_radial_excitation",
    "C1_eigenvector_weighted_B1",
    "C2_eigenvector_weighted_B2",
    "none",
)
DEFAULT_SK_CANDIDATE = "A_lowest_radial_per_l"

WIRED_CANDIDATES = (
    "A_lowest_radial_per_l",
    "B1_single_angular_mode",
    "B2_single_radial_excitation",
    "C1_eigenvector_weighted_B1",
    "C2_eigenvector_weighted_B2",
)


@dataclass(frozen=True)
class Mode:
    l: int
    n: int


def s_k_membership(k: int, candidate: str) -> list[Mode]:
    """
    Return the (l, n) modes that constitute S(k) under the chosen candidate.

    For the eigenvector-weighted candidates C1/C2, the membership returned
    here is the full depth-basis mode set (one per LEPTON_DEPTHS entry);
    use `s_k_weighted_modes` to also receive the |v_species,i|² weights.
    """
    if candidate == "none":
        return []
    if candidate == "A_lowest_radial_per_l":
        return [Mode(l=ll, n=0) for ll in range(1, k + 1, 2)]
    if candidate == "B1_single_angular_mode":
        return [Mode(l=k, n=0)]
    if candidate == "B2_single_radial_excitation":
        if (k - 1) % 2 != 0:
            raise ValueError(
                f"B2_single_radial_excitation requires odd k; got k={k}"
            )
        return [Mode(l=1, n=(k - 1) // 2)]
    if candidate == "C1_eigenvector_weighted_B1":
        # All B1 modes across the depth basis; weights selected by k.
        return [Mode(l=ll, n=0) for ll in LEPTON_DEPTHS]
    if candidate == "C2_eigenvector_weighted_B2":
        # All B2 modes across the depth basis; weights selected by k.
        return [Mode(l=1, n=(ll - 1) // 2) for ll in LEPTON_DEPTHS]
    raise ValueError(f"Unknown S(k) candidate: {candidate}")


@lru_cache(maxsize=1)
def _locked_lepton_eigenvectors() -> tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]:
    """
    Diagonalize the locked lepton generation block and return its
    (eigenvalues_ascending, eigenvectors_columns) as plain Python tuples.

    The block is built from the canonical baseline parameters in
    `geometrodynamics.tangherlini.lepton_spectrum`:

        depths            = LEPTON_BASELINE_DEPTHS                 = (1, 3, 5)
        phase_per_pass    = LEPTON_BASELINE_PHASE                  = 0.001
        transport_strength= LEPTON_BASELINE_TRANSPORT              = 25.1
        hard_pinhole_gamma= LEPTON_BASELINE_PINHOLE                = 22.5
        resistance_scale  = LEPTON_BASELINE_RESISTANCE             = 0.217869...
        resistance_model  = "exponential"
        depth_cost_mode   = "tunnel_only"
        winding_mode      = "max"
        action_base       = S3_ACTION_BASE   (= 2π)
        k_uplift_beta     = TAU_BETA_50PI    (= 50π)
        action_slope      = 0.5              (default of compute_knotted_lepton_spectrum)

    Eigenvalues are sorted ascending; columns of the returned matrix
    are the corresponding eigenvectors. With these parameters every
    eigenvalue is positive and species_index → (e, μ, τ) follows the
    same convention as `compute_knotted_lepton_spectrum`.

    Cached so repeated calls are deterministic and free.
    """
    import numpy as np
    from scipy.linalg import eigh
    from geometrodynamics.tangherlini.lepton_spectrum import (
        _build_generation_block,
        LEPTON_BASELINE_DEPTHS,
        LEPTON_BASELINE_PHASE,
        LEPTON_BASELINE_TRANSPORT,
        LEPTON_BASELINE_PINHOLE,
        LEPTON_BASELINE_RESISTANCE,
        S3_ACTION_BASE,
        TAU_BETA_50PI,
    )

    h = _build_generation_block(
        depths=LEPTON_BASELINE_DEPTHS,
        phase_per_pass=LEPTON_BASELINE_PHASE,
        transport_strength=LEPTON_BASELINE_TRANSPORT,
        resistance_model="exponential",
        resistance_scale=LEPTON_BASELINE_RESISTANCE,
        hard_pinhole_gamma=LEPTON_BASELINE_PINHOLE,
        action_base=S3_ACTION_BASE,
        action_slope=0.5,
        depth_cost_mode="tunnel_only",
        winding_mode="max",
        k_uplift_beta=TAU_BETA_50PI,
    )
    w, V = eigh(h)
    return (
        tuple(float(x) for x in w),
        tuple(tuple(float(x) for x in col) for col in V.T),
    )


def lepton_eigenvector_weights() -> dict[int, list[float]]:
    """
    Return |v_species,i|² weights for each species, keyed by depth k.

    Result: {k: [w_for_depth_basis_index_0, w_for_idx_1, w_for_idx_2]}
    where the species (and hence the eigenvector) is selected by the
    species_index_for_k convention: k=1→e, k=3→μ, k=5→τ.

    Each row sums to 1 by orthonormality of the eigenvectors.
    """
    _, eigenvectors = _locked_lepton_eigenvectors()
    weights_by_k: dict[int, list[float]] = {}
    for sp_idx, k in enumerate(LEPTON_DEPTHS):
        v = eigenvectors[sp_idx]
        weights_by_k[k] = [float(c) ** 2 for c in v]
    return weights_by_k


def s_k_weighted_modes(k: int, candidate: str) -> list[tuple[Mode, float]]:
    """
    Return the (mode, weight) list for the candidate.

    For unweighted candidates (A/B1/B2), every weight is 1.0.
    For C1/C2, the weights are the squared eigenvector amplitudes of
    the species selected by k (k=1 → electron, k=3 → muon, k=5 → tau).
    """
    if candidate in ("C1_eigenvector_weighted_B1", "C2_eigenvector_weighted_B2"):
        if k not in LEPTON_DEPTHS:
            raise ValueError(
                f"C1/C2 require k ∈ {LEPTON_DEPTHS}; got k={k}"
            )
        weights = lepton_eigenvector_weights()[k]
        modes = s_k_membership(k, candidate)
        return list(zip(modes, weights))
    return [(m, 1.0) for m in s_k_membership(k, candidate)]


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
    weight: float = 1.0
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
            "weight": self.weight,
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
    weighted = s_k_weighted_modes(k, candidate)
    mode_results: list[ModePhase] = []
    for mode, weight in weighted:
        mp = phi_radial_for_mode(
            mode.l, mode.n, N=N, maslov_correction=maslov_correction,
        )
        mp.weight = float(weight)
        mode_results.append(mp)
    if mode_results and all(mr.status == "computed" for mr in mode_results):
        total = sum(mr.weight * mr.phi for mr in mode_results)
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
