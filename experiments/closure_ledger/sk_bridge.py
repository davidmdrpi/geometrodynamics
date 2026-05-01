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
    C1_maslov_standard          — C1 mode set with −π/2 per turning point
    B2_maslov_standard          — B2 mode set with −π/2 per turning point
    C2_maslov_standard          — C2 mode set with −π/2 per turning point
    D0_overlap_phase            — v^T Φ v with Φ_ij = π·⟨u_i|u_j⟩
    D1_potential_difference_phase — v^T Φ v with Φ_ij = ⟨u_i|V_j−V_i|u_j⟩
                                  (off-diagonal mirrored to make Hermitian)
    D2_symmetrized_momentum_phase — v^T Φ v with
                                  Φ_ij = ⟨u_i|√(ω̄²−V̄)_+|u_j⟩, ω̄/V̄ symmetric

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
    "C1_maslov_standard",
    "B2_maslov_standard",
    "C2_maslov_standard",
    "D0_overlap_phase",
    "D1_potential_difference_phase",
    "D2_symmetrized_momentum_phase",
    "none",
)
DEFAULT_SK_CANDIDATE = "A_lowest_radial_per_l"

WIRED_CANDIDATES = (
    "A_lowest_radial_per_l",
    "B1_single_angular_mode",
    "B2_single_radial_excitation",
    "C1_eigenvector_weighted_B1",
    "C2_eigenvector_weighted_B2",
    "C1_maslov_standard",
    "B2_maslov_standard",
    "C2_maslov_standard",
    "D0_overlap_phase",
    "D1_potential_difference_phase",
    "D2_symmetrized_momentum_phase",
)

# Candidates that derive the per-mode Maslov correction from a turning-point
# count rather than from a caller-supplied constant offset. The "standard"
# policy applies −π/2 per detected turning point (sign change of ω² − V_eff
# inside the integration grid).
MASLOV_STANDARD_CANDIDATES = (
    "C1_maslov_standard",
    "B2_maslov_standard",
    "C2_maslov_standard",
)

# Operator-valued candidates: the radial phase is v^T Φ v where Φ is a
# 3×3 Hermitian matrix indexed by the depth-basis B1 modes (l ∈ {1,3,5},
# n=0). Off-diagonal entries — the new degree of freedom relative to the
# scalar-per-mode C-family — are determined by an explicit operator
# kernel; v are the locked lepton eigenvectors used by C1/C2.
D_OPERATOR_VARIANTS = (
    "D0_overlap_phase",
    "D1_potential_difference_phase",
    "D2_symmetrized_momentum_phase",
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
    if candidate == "C1_maslov_standard":
        # Same membership as C1; differs only in the Maslov policy applied
        # by phi_radial_from_sk.
        return [Mode(l=ll, n=0) for ll in LEPTON_DEPTHS]
    if candidate == "B2_maslov_standard":
        # Same membership as B2; differs only in the Maslov policy.
        if (k - 1) % 2 != 0:
            raise ValueError(
                f"B2_maslov_standard requires odd k; got k={k}"
            )
        return [Mode(l=1, n=(k - 1) // 2)]
    if candidate == "C2_maslov_standard":
        # Same membership as C2; differs only in the Maslov policy.
        return [Mode(l=1, n=(ll - 1) // 2) for ll in LEPTON_DEPTHS]
    if candidate in D_OPERATOR_VARIANTS:
        # Operator candidates use the B1 depth-basis modes as the row/column
        # index of a 3×3 Hermitian Φ matrix; the scalar per-mode list returned
        # here matches C1's diagonal labels for compatibility with the
        # mode-map comparison rendering.
        return [Mode(l=ll, n=0) for ll in LEPTON_DEPTHS]
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
    if candidate in (
        "C1_eigenvector_weighted_B1",
        "C2_eigenvector_weighted_B2",
        "C1_maslov_standard",
        "C2_maslov_standard",
    ):
        if k not in LEPTON_DEPTHS:
            raise ValueError(
                f"C-family candidates require k ∈ {LEPTON_DEPTHS}; got k={k}"
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
    maslov_policy: str = "none"   # "none" | "standard"
    n_turning_points: int = 0
    N_grid: int = 0
    weight: float = 1.0
    error: Optional[str] = None
    # Operator-valued candidates only: the row/column pair of the 3×3 Φ
    # matrix this entry contributes to. For diagonal-only candidates these
    # remain None.
    pair_l: Optional[int] = None
    pair_n: Optional[int] = None
    operator_variant: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "l": self.l,
            "n": self.n,
            "omega": self.omega,
            "phi": self.phi,
            "status": self.status,
            "convention": self.convention,
            "maslov_correction": self.maslov_correction,
            "maslov_policy": self.maslov_policy,
            "n_turning_points": self.n_turning_points,
            "N_grid": self.N_grid,
            "weight": self.weight,
            "error": self.error,
            "pair_l": self.pair_l,
            "pair_n": self.pair_n,
            "operator_variant": self.operator_variant,
        }


def _count_turning_points(omega2_minus_V) -> int:
    """
    Count classical turning points along the integration grid.

    A turning point is detected as a sign change of (ω² − V_eff) between
    adjacent grid samples; zeros count as their adjacent sign for this
    purpose. This counts the genuinely soft turning points where the
    classical momentum vanishes inside the grid; a grid endpoint where
    (ω² − V_eff) > 0 is a hard wall and is NOT counted.
    """
    import numpy as np
    arr = np.asarray(omega2_minus_V)
    sign = np.sign(arr)
    # Treat any zero as belonging to its previous sign so isolated zeros
    # don't double-count.
    for i in range(1, len(sign)):
        if sign[i] == 0:
            sign[i] = sign[i - 1]
    if sign[0] == 0:
        # Pick whichever non-zero sign appears first.
        for s in sign:
            if s != 0:
                sign[0] = s
                break
    changes = int(np.sum(np.diff(sign) != 0))
    return changes


def phi_radial_for_mode(
    l: int,
    n: int,
    *,
    N: int = 80,
    maslov_correction: float = 0.0,
    maslov_policy: str = "none",
) -> ModePhase:
    """
    Compute Φ(l, n) by integrating √max(ω² − V_eff, 0) over the tortoise grid.

    Convention: one-sided WKB radial action, no Maslov correction by
    default.

    `maslov_policy` controls how the Maslov phase is determined:
        "none"     — use the caller-supplied `maslov_correction` verbatim
                     (default; preserves the prior behaviour).
        "standard" — derive the correction as `−π/2 · n_turning_points`,
                     where `n_turning_points` is the number of sign changes
                     of (ω² − V_eff) inside the tortoise grid. The
                     `maslov_correction` argument is ignored under this
                     policy and the computed value is reported on the
                     returned ModePhase.

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
            maslov_correction=maslov_correction,
            maslov_policy=maslov_policy, N_grid=N,
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
            maslov_correction=maslov_correction,
            maslov_policy=maslov_policy, N_grid=N,
        )

    import math as _math
    omega_n = float(omegas[n])
    rs = float(R_MID)
    Vg = V_tangherlini(r_grid, l, rs)
    diff = omega_n ** 2 - Vg
    integrand = np.sqrt(np.maximum(diff, 0.0))
    rstar = np.array([r_to_rstar(float(r), rs) for r in r_grid])
    order = np.argsort(rstar)
    phi_wkb = float(np.trapezoid(integrand[order], rstar[order]))

    n_tp = _count_turning_points(diff[order])
    if maslov_policy == "standard":
        applied_correction = -(_math.pi / 2.0) * n_tp
    elif maslov_policy == "none":
        applied_correction = float(maslov_correction)
    else:
        raise ValueError(
            f"Unknown maslov_policy: {maslov_policy!r}; "
            "expected 'none' or 'standard'."
        )

    return ModePhase(
        l=l, n=n, omega=omega_n,
        phi=phi_wkb + applied_correction,
        status="computed",
        maslov_correction=applied_correction,
        maslov_policy=maslov_policy,
        n_turning_points=n_tp,
        N_grid=N,
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


@lru_cache(maxsize=64)
def _operator_phase_matrix(
    variant: str,
    *,
    N: int = 80,
    phase_scale: float = 3.141592653589793,
) -> tuple[
    tuple[tuple[float, ...], ...],   # Phi (3×3)
    tuple[float, ...],                # omegas[i] for each depth-basis mode
    tuple[tuple[int, int], ...],      # mode labels (l_i, n_i)
]:
    """
    Build the 3×3 Hermitian radial phase matrix Φ for an operator-valued
    candidate. Indexed by the B1 depth-basis modes (l ∈ {1,3,5}, n=0).

    Eigenfunctions u_i are L²-normalized on the tortoise grid before
    matrix-element evaluation. The grid is shared across all l in
    `solve_radial_modes` (it depends only on N, rs, r_outer), so all
    three eigenfunctions share a common discretization.

    Variants:
        D0_overlap_phase            — Φ_ij = phase_scale · ⟨u_i | u_j⟩
                                      (positive on diagonal; off-diagonal
                                      measures orthogonality breaking
                                      between different-l eigenfunctions)
        D1_potential_difference_phase — Φ_ij = ⟨u_i | V_j − V_i | u_j⟩,
                                        Hermitized by mirroring i<j to j<i
                                        (the canonical quark transport
                                        matrix-element form, reused as the
                                        upper triangle and reflected to make
                                        the matrix real symmetric).
                                        Diagonal Φ_ii = 0 by construction.
        D2_symmetrized_momentum_phase — Φ_ij = ⟨u_i | √max(ω̄² − V̄, 0) | u_j⟩
                                        with ω̄² = (ω_i² + ω_j²)/2 and
                                        V̄(r) = (V_i(r) + V_j(r))/2. Symmetric
                                        by construction (kernel symmetric in
                                        i,j; both u real).
    """
    import math as _math
    import numpy as np
    from geometrodynamics.tangherlini.radial import (
        solve_radial_modes,
        r_to_rstar,
        V_tangherlini,
    )
    from geometrodynamics.constants import R_MID

    if variant not in D_OPERATOR_VARIANTS:
        raise ValueError(f"Unknown operator variant: {variant!r}")

    rs = float(R_MID)
    omegas: list[float] = []
    us: list = []
    Vs: list = []
    rstar_sorted = None
    for ll in LEPTON_DEPTHS:
        oms, funcs, rg = solve_radial_modes(ll, N=N, n_modes=1)
        omega = float(oms[0])
        u = np.asarray(funcs[0]["u_half"], dtype=float)
        Vg = np.asarray(V_tangherlini(rg, ll, rs), dtype=float)
        rstar = np.array([r_to_rstar(float(r), rs) for r in rg])
        order = np.argsort(rstar)
        if rstar_sorted is None:
            rstar_sorted = rstar[order]
        u_sorted = u[order]
        norm2 = float(np.trapezoid(u_sorted ** 2, rstar_sorted))
        u_normed = u_sorted / _math.sqrt(norm2)
        omegas.append(omega)
        us.append(u_normed)
        Vs.append(Vg[order])

    Phi = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            if variant == "D0_overlap_phase":
                Phi[i, j] = phase_scale * float(
                    np.trapezoid(us[i] * us[j], rstar_sorted)
                )
            elif variant == "D1_potential_difference_phase":
                if i == j:
                    Phi[i, j] = 0.0
                else:
                    Phi[i, j] = float(np.trapezoid(
                        us[i] * (Vs[j] - Vs[i]) * us[j], rstar_sorted,
                    ))
            elif variant == "D2_symmetrized_momentum_phase":
                omega_bar2 = (omegas[i] ** 2 + omegas[j] ** 2) / 2.0
                V_bar = (Vs[i] + Vs[j]) / 2.0
                kernel = np.sqrt(np.maximum(omega_bar2 - V_bar, 0.0))
                Phi[i, j] = float(np.trapezoid(
                    us[i] * kernel * us[j], rstar_sorted,
                ))
            if i != j:
                Phi[j, i] = Phi[i, j]

    Phi_t = tuple(tuple(float(x) for x in row) for row in Phi)
    omegas_t = tuple(omegas)
    labels_t = tuple((int(ll), 0) for ll in LEPTON_DEPTHS)
    return Phi_t, omegas_t, labels_t


def _operator_mode_phases(
    k: int, variant: str, *, N: int = 80,
) -> list[ModePhase]:
    """
    Materialize the contraction v^T Φ v as a list of 9 ModePhase entries
    (one per (i, j) pair of the 3×3 Φ matrix). Each entry carries
    weight = v_species,i · v_species,j and phi = Φ_ij; the candidate's
    total radial phase is Σ weight·phi over the 9 entries, which equals
    v^T Φ v exactly.
    """
    Phi, omegas, labels = _operator_phase_matrix(variant, N=N)
    weights = lepton_eigenvector_weights()  # placeholder to ensure cache warm
    # The C-family weights helper returns |v_i|² but we need signed v_i
    # here; pull the eigenvectors directly.
    _, eigenvectors = _locked_lepton_eigenvectors()
    if k not in LEPTON_DEPTHS:
        raise ValueError(
            f"Operator candidates require k ∈ {LEPTON_DEPTHS}; got k={k}"
        )
    sp_idx = LEPTON_DEPTHS.index(k)
    v = eigenvectors[sp_idx]   # signed amplitudes over the depth basis
    entries: list[ModePhase] = []
    for i, (l_i, n_i) in enumerate(labels):
        for j, (l_j, n_j) in enumerate(labels):
            entries.append(ModePhase(
                l=l_i, n=n_i,
                omega=float(omegas[i]),
                phi=float(Phi[i][j]),
                status="computed",
                convention="operator_radial_phase",
                maslov_correction=0.0,
                maslov_policy="none",
                n_turning_points=0,
                N_grid=N,
                weight=float(v[i] * v[j]),
                pair_l=l_j, pair_n=n_j,
                operator_variant=variant,
            ))
    return entries


def phi_radial_from_sk(
    k: int, candidate: str, *,
    N: int = 80, maslov_correction: float = 0.0,
) -> RadialPhaseResult:
    """
    Compute Φ_radial(k) = Σ_(l,n)∈S(k) w(l, n; k) · Φ(l, n).

    `maslov_correction` is added per mode for candidates that use the
    "none" Maslov policy (so a candidate with |S(k)| modes accumulates
    |S(k)| · maslov_correction). Candidates listed in
    MASLOV_STANDARD_CANDIDATES instead derive the per-mode correction
    as −π/2 per detected classical turning point.

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
    if candidate in D_OPERATOR_VARIANTS:
        try:
            mode_results = _operator_mode_phases(k, candidate, N=N)
        except ImportError as exc:
            return RadialPhaseResult(
                candidate=candidate, k=k, modes=[],
                total_phi=None, status="missing",
            )
        total = sum(mr.weight * mr.phi for mr in mode_results)
        return RadialPhaseResult(
            candidate=candidate, k=k, modes=mode_results,
            total_phi=float(total), status="computed",
        )
    policy = (
        "standard" if candidate in MASLOV_STANDARD_CANDIDATES else "none"
    )
    weighted = s_k_weighted_modes(k, candidate)
    mode_results: list[ModePhase] = []
    for mode, weight in weighted:
        mp = phi_radial_for_mode(
            mode.l, mode.n, N=N,
            maslov_correction=maslov_correction,
            maslov_policy=policy,
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
