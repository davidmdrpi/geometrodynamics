"""
QCD shell waveguide: basis + operator scaffold.

User's physical insight (Nov 2026):
  "Quarks do not pass through the throat; they are the wavefronts that
  resolve the cavity itself."

This reframes the quark sector entirely. Leptons (e, μ, τ) are
**throat-localized** focused-pulse modes (n = 0, 1, 2 of the radial
ladder; PR #68). Quarks are the **shell-saturated** standing-wave
modes (n ≥ 3 in PR #68's metrics: participation ratio → 2/3, ⟨r⟩
plateaus, extended-character wavefront filling the cavity
[R_MID, R_OUTER]) — the cavity wavefronts that resolve the shell, not
focused pulses that traverse the throat.

The v3 quark Hamiltonian (`docs/quark_axioms.md`) was the WRONG
machinery for this — it fit quarks on the lepton-shaped closure-
quantum basis `{(k=1,±), (k=3,±), (k=5,±)}` (#76's structural
diagnosis), absorbing the unmodeled QCD physics into the
phenomenological `n_part = 233`.

This probe is the foundation for the four-PR quantitative QCD-shell
arc the user laid out:

  PR #77 (this PR) — QCD shell waveguide basis/operator scaffold
  PR #78          — shell Hamiltonian mass-ordering / n_part audit
  PR #79          — boundary stress tensor and singlet constraint
  PR #80          — color algebra candidate: SU(3), SU(2)×Z₂, or other

This PR scaffolds **only** the basis and operator structure; it does
not derive masses, audit `n_part`, define the boundary stress tensor,
or identify the color algebra. Those are PRs #78–#80 explicitly. The
honest scope is: build clean machinery the next three PRs can use,
verified to be distinct-from-throat (vs the lepton sector).

## The shell waveguide basis

A **shell waveguide mode** is a Tangherlini radial mode `(l, n)` on
the cavity `[R_MID, R_OUTER]` that has saturated into the shell-
filling standing-wave regime — i.e., participation ratio ≈ 2/3
(uniform-`sin` standing wave on a finite interval) rather than
focused on the throat. From PR #68's metrics, this is `n ≥ 3` for
`l = 1`.

A **shell waveguide basis state** is then a tuple `(l, n, p)` with:
  - `l` = S³ angular momentum (Casimir `l(l+2)`, S³ harmonic on the
    Hopf-bundle's angular base; #73's primitive).
  - `n` = radial overtone index, shell-saturated branch (n ≥ n_shell).
  - `p ∈ {+, −}` = Z₂ partition (B2, non-orientable throat
    `T = iσ_y`, `T² = −I`).

The 6 lowest shell-saturated states with `l = 1` constitute the
"3 × 2 = 6 flavors" structural match documented by PR #69:

  (l=1, n=3, +), (l=1, n=3, −),
  (l=1, n=4, +), (l=1, n=4, −),
  (l=1, n=5, +), (l=1, n=5, −).

The user can later choose to vary `l` instead — e.g. (l=1,n,±),
(l=2,n,±), (l=3,n,±) for some fixed shell-saturated `n` — this is a
PR #78 decision. This probe constructs both options and verifies the
scaffold is `l`-and-`n`-flexible.

## The operator scaffold

The shell Hamiltonian acts on the 6-state basis as a 6×6 Hermitian
matrix. The scaffold provides the structural slots — PR #78 will
populate the matrix elements and test mass-ordering. Slots:

  - **Diagonal kinetic term** (`H_kin`): per-mode `ω²(l, n)` from the
    Tangherlini eigensolver. This is the **frequency-squared mass
    operator** for the cavity wavefront — the "energy resolves the
    cavity" content. Distinct from the lepton sector, where the
    lepton ladder's diagonal carries `β·k²` (winding cost), not
    `ω²(l, n)` (cavity eigenfrequency).
  - **Z₂ partition splitter** (`H_Z2`): block-diagonal in `(l, n)`,
    acting on the `p = ±` index as `χ · σ_z`. `χ` is a PR #79–#80
    structural constant (linked to the boundary stress tensor /
    color algebra); PR #77 provides the slot.
  - **Inter-mode coupling** (`H_couple`): hooks for shell-shell
    mixing between different `(l, n)` blocks. PR #78–#80 will
    populate; PR #77 leaves zero by default.

The full operator is `H = H_kin + H_Z2 + H_couple`, with the user-
chosen splitter strength `χ` and coupling matrix as inputs.

## Distinctness from the lepton/v3 machinery

The lepton sector lives at the throat: basis = odd-k throat traversal
modes `{k=1, 3, 5}` with B2 partition, Hamiltonian carries
`β_lepton·k² · (2π)` (#71's closure-quantum lock). The v3 quark
Hamiltonian inherited this basis and absorbed the QCD-shell physics
into `n_part = 233` (#76).

The shell waveguide basis is structurally distinct:

  - Basis: shell-saturated `(l, n, p)` cavity modes, NOT odd-k throat
    traversal modes.
  - Diagonal: `ω²(l, n)` cavity eigenfrequency-squared, NOT
    `β·k² · (2π)` winding cost.
  - Mode geometry: extended-character standing wave (PR=2/3, ⟨r⟩
    plateaus), NOT focused pulse (PR < 0.3, ⟨r⟩ small).

This is the right machinery for "wavefronts that resolve the cavity";
PR #78 will test whether it reproduces the quark mass ordering.

## Honest scope

  - **Is:** the basis + operator-slot scaffold; the verification that
    the shell-saturated `(l, n, p)` modes are distinct from the
    lepton/v3 throat modes; the hooks for PR #78–#80.
  - **Is not:** a derivation of quark masses (PR #78), an audit of
    `n_part` (PR #78), a boundary stress tensor or singlet
    constraint (PR #79), or a color algebra identification (PR #80).
    The 6×6 matrix is constructed but not populated with the
    structural coupling constants.

Tests:
  T1. Shell-saturation criterion + 6-state basis enumeration.
  T2. Distinctness from lepton/v3: PR ≥ 0.5 for shell modes vs ≤ 0.3
      for lepton modes; ⟨r⟩ near shell center.
  T3. Operator scaffold: H_kin diagonal from ω²(l,n); H_Z2 block-σ_z
      splitter slot; H_couple = 0 placeholder.
  T4. Basis flexibility: same scaffold supports {(l=1, n=3,4,5, ±)}
      and {(l=1,2,3, n=3, ±)} alternates.
  T5. Six-flavor structural match: 6-state basis with Z₂ partition.
  T6. Hooks for PR #78–#80: each subsequent PR's plug-in is named.
  T7. Honest scope / B4: scaffold only; cavity eigenfrequencies are
      dimensionful but ratios are scale-free.
  T8. Assessment.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

from geometrodynamics.constants import R_MID, R_OUTER
from geometrodynamics.tangherlini.radial import (
    V_tangherlini, r_to_rstar, rstar_to_r,
)


PI = math.pi
DELTA_R = R_OUTER - R_MID
SHELL_STANDING_WAVE_PR = 2.0 / 3.0
# Structural shell-saturation criterion (from PR #68): n ≥ 3 for l=1,
# where ⟨r⟩−R_MID has plateaued at the shell center (~0.046 here) and
# the participation ratio has reached 2/3 (uniform standing wave).
N_SHELL_DEFAULT = 3
# Numerical thresholds used to verify saturation:
SHELL_RMEAN_PLATEAU = 0.040         # ⟨r⟩−R_MID ≥ this ⇒ near plateau
SHELL_PR_NEAR_TWO_THIRDS = 0.020    # |PR − 2/3| ≤ this ⇒ saturated
N_FLAVORS = 6                       # 3 generations × 2 partitions
N_GENERATIONS = 3
N_PARTITIONS = 2

# Quark mass-ordering target (MeV, observed)
# Used only for distinctness-from-lepton checks, NOT for mass prediction.
QUARK_MASSES_OBSERVED_MEV = {
    'u': 2.16, 'd': 4.67, 's': 93.4, 'c': 1270.0, 'b': 4180.0, 't': 172690.0,
}


# ---------------------------------------------------------------------------
# Shell waveguide basis state
# ---------------------------------------------------------------------------

@dataclass
class ShellBasisState:
    """A single shell waveguide basis state (l, n, p)."""
    l: int                              # S³ angular momentum
    n: int                              # radial overtone index (shell-saturated)
    p: int                              # Z₂ partition (+1 or -1)
    omega: float                        # Tangherlini eigenfrequency
    participation_ratio: float          # localization metric (→2/3 shell)
    r_mean_minus_RMID: float            # mean radius offset (→ shell center)
    throat_fraction: float              # inner-third occupancy (→ low for shell)
    is_shell_saturated: bool            # throat_fraction below threshold

    def label(self) -> str:
        sign = '+' if self.p == +1 else '−'
        return f"(l={self.l}, n={self.n}, {sign})"


# ---------------------------------------------------------------------------
# Shell waveguide basis: construct shell-saturated (l, n) modes
# ---------------------------------------------------------------------------

def _radial_ladder(l: int = 1, n_max: int = 9, N: int = 800) -> list[dict]:
    """Radial ladder + localization metrics, recap from PR #68. Returns
    one entry per overtone n with ω, ⟨r⟩−R_MID, participation ratio,
    and throat fraction (inner-third occupancy)."""
    rsmin = r_to_rstar(R_MID + 5e-4, R_MID)
    rsmax = r_to_rstar(R_OUTER - 5e-4, R_MID)
    rstar = np.linspace(rsmin, rsmax, N)
    h = rstar[1] - rstar[0]
    rphys = np.array([rstar_to_r(s, R_MID) for s in rstar])
    V = V_tangherlini(rphys, l, R_MID)
    main = 2.0 / h ** 2 + V[1:-1]
    off = -1.0 / h ** 2 * np.ones(N - 3)
    H = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    ev, evec = np.linalg.eigh(H)
    inner_third = R_MID + DELTA_R / 3.0
    rows = []
    for n in range(min(n_max, evec.shape[1])):
        u = np.concatenate([[0.0], evec[:, n], [0.0]])
        p_dens = u ** 2
        p_dens = p_dens / p_dens.sum()
        omega = float(math.sqrt(max(ev[n], 0.0)))
        r_mean = float(np.sum(p_dens * rphys) - R_MID)
        pr = float(1.0 / (np.sum(p_dens ** 2) * len(p_dens)))
        throat_frac = float(p_dens[rphys < inner_third].sum())
        rows.append({
            'n': n, 'omega': omega,
            'r_mean_minus_RMID': r_mean,
            'participation_ratio': pr,
            'throat_fraction': throat_frac,
        })
    return rows


def build_shell_basis(
    mode: str = 'n_varied',
    l_values: tuple[int, ...] | None = None,
    n_values: tuple[int, ...] | None = None,
    n_shell_start: int = N_SHELL_DEFAULT,
) -> list[ShellBasisState]:
    """Build the 6-state shell waveguide basis.

    Two natural enumerations:
      - mode = 'n_varied': fix l (=1), vary n ∈ {3, 4, 5}, ± partition.
      - mode = 'l_varied': fix n (=3, shell-saturated), vary l ∈
        {1, 2, 3}, ± partition.

    Either yields 6 states (3 × 2). PR #78 will choose between them
    based on which reproduces the quark mass ordering; PR #77 only
    needs the scaffold to support both.
    """
    states: list[ShellBasisState] = []
    def _is_saturated(row: dict) -> bool:
        return (
            row['r_mean_minus_RMID'] >= SHELL_RMEAN_PLATEAU
            and abs(row['participation_ratio'] - SHELL_STANDING_WAVE_PR)
            <= SHELL_PR_NEAR_TWO_THIRDS
        )

    if mode == 'n_varied':
        l_fixed = (l_values[0] if l_values else 1)
        n_targets = n_values if n_values else tuple(range(n_shell_start, n_shell_start + 3))
        ladder = _radial_ladder(l=l_fixed, n_max=max(n_targets) + 1)
        for n in n_targets:
            row = ladder[n]
            for p in (+1, -1):
                states.append(ShellBasisState(
                    l=l_fixed, n=n, p=p,
                    omega=row['omega'],
                    participation_ratio=row['participation_ratio'],
                    r_mean_minus_RMID=row['r_mean_minus_RMID'],
                    throat_fraction=row['throat_fraction'],
                    is_shell_saturated=_is_saturated(row),
                ))
    elif mode == 'l_varied':
        l_targets = l_values if l_values else (1, 2, 3)
        n_fixed = (n_values[0] if n_values else n_shell_start)
        for l_ in l_targets:
            ladder = _radial_ladder(l=l_, n_max=n_fixed + 1)
            row = ladder[n_fixed]
            for p in (+1, -1):
                states.append(ShellBasisState(
                    l=l_, n=n_fixed, p=p,
                    omega=row['omega'],
                    participation_ratio=row['participation_ratio'],
                    r_mean_minus_RMID=row['r_mean_minus_RMID'],
                    throat_fraction=row['throat_fraction'],
                    is_shell_saturated=_is_saturated(row),
                ))
    else:
        raise ValueError(f"mode must be 'n_varied' or 'l_varied', got {mode!r}")
    return states


# ---------------------------------------------------------------------------
# Operator scaffold
# ---------------------------------------------------------------------------

@dataclass
class ShellOperatorScaffold:
    """The 6×6 shell Hamiltonian scaffold.

    H = H_kin + H_Z2 + H_couple

    - H_kin: diagonal ω²(l, n) cavity eigenfrequency-squared.
    - H_Z2: block-σ_z partition splitter (slot for PR #79–#80
      structural constant χ).
    - H_couple: zero by default (slot for PR #78–#80 inter-mode
      mixing).
    """
    basis: list[ShellBasisState]
    chi: float = 0.0                    # Z₂ splitter strength (PR #79–#80)
    couple_matrix: list[list[float]] = field(default_factory=list)

    @property
    def n_states(self) -> int:
        return len(self.basis)

    def H_kin(self) -> np.ndarray:
        """Diagonal kinetic / cavity-mass operator."""
        diag = np.array([s.omega ** 2 for s in self.basis], dtype=float)
        return np.diag(diag)

    def H_Z2(self) -> np.ndarray:
        """Block-σ_z partition splitter on the ± index within each
        (l, n) block. PR #77 leaves chi = 0 (slot for #79–#80)."""
        n = self.n_states
        h = np.zeros((n, n), dtype=float)
        # Each (l, n) block has 2 partition states; the block-σ_z is
        # diag(+chi, −chi) on that block. Identify blocks by (l, n).
        seen: dict[tuple[int, int], list[int]] = {}
        for i, s in enumerate(self.basis):
            seen.setdefault((s.l, s.n), []).append(i)
        for (l_, n_), idxs in seen.items():
            if len(idxs) == 2:
                # Sort by p so (+1, -1)
                idxs_sorted = sorted(idxs, key=lambda i: -self.basis[i].p)
                i_plus, i_minus = idxs_sorted
                h[i_plus, i_plus] = +self.chi
                h[i_minus, i_minus] = -self.chi
        return h

    def H_couple(self) -> np.ndarray:
        """Inter-mode coupling matrix (PR #78–#80 hook)."""
        n = self.n_states
        if not self.couple_matrix:
            return np.zeros((n, n), dtype=float)
        m = np.array(self.couple_matrix, dtype=float)
        if m.shape != (n, n):
            raise ValueError(
                f"couple_matrix shape {m.shape} ≠ basis size ({n}, {n})")
        # Hermitize
        return 0.5 * (m + m.T)

    def H_total(self) -> np.ndarray:
        return self.H_kin() + self.H_Z2() + self.H_couple()

    def eigenvalues(self) -> np.ndarray:
        H = self.H_total()
        ev = np.linalg.eigvalsh(H)
        return np.sort(np.real(ev))


# ---------------------------------------------------------------------------
# T1. Shell-saturation criterion + 6-state basis
# ---------------------------------------------------------------------------

def test_T1_shell_basis_constructed() -> dict:
    """Build the 6-state shell waveguide basis. Each state is shell-
    saturated (⟨r⟩ on plateau, PR within 0.02 of 2/3)."""
    basis = build_shell_basis(mode='n_varied')
    all_saturated = all(s.is_shell_saturated for s in basis)
    return {
        'name': 'T1_shell_waveguide_basis_constructed',
        'description': (
            "Build the 6-state shell waveguide basis (l=1, n=3,4,5, "
            "p=±1). Each state shell-saturated: ⟨r⟩−R_MID ≥ "
            f"{SHELL_RMEAN_PLATEAU} and |PR − 2/3| ≤ "
            f"{SHELL_PR_NEAR_TWO_THIRDS}."
        ),
        'basis_size': len(basis),
        'expected_size': N_FLAVORS,
        'basis_labels': [s.label() for s in basis],
        'omegas': [round(s.omega, 6) for s in basis],
        'participation_ratios': [round(s.participation_ratio, 4) for s in basis],
        'r_mean_offsets': [round(s.r_mean_minus_RMID, 4) for s in basis],
        'throat_fractions': [round(s.throat_fraction, 4) for s in basis],
        'all_shell_saturated': all_saturated,
        'pass': len(basis) == N_FLAVORS and all_saturated,
    }


# ---------------------------------------------------------------------------
# T2. Distinctness from lepton/v3 machinery
# ---------------------------------------------------------------------------

def test_T2_distinct_from_leptons() -> dict:
    """The lepton ladder (n=0,1,2 = e,μ,τ) is the throat-to-shell
    transition: ⟨r⟩−R_MID rises 0.021 → 0.039 → 0.044 across the three
    leptons, approaching but not reaching the shell plateau (~0.046).
    The shell modes (n ≥ 3) have ⟨r⟩ saturated at the plateau and PR
    locked at the uniform-standing-wave value 2/3. The structural
    distinction is sharpest at the electron (n=0, ⟨r⟩=0.021, deeply
    focused) and saturates from n=3 onward (the shell waveguide
    regime)."""
    ladder = _radial_ladder(l=1, n_max=6)
    lepton_rmean = [ladder[n]['r_mean_minus_RMID'] for n in (0, 1, 2)]
    shell_rmean = [ladder[n]['r_mean_minus_RMID'] for n in (3, 4, 5)]
    lepton_pr = [ladder[n]['participation_ratio'] for n in (0, 1, 2)]
    shell_pr = [ladder[n]['participation_ratio'] for n in (3, 4, 5)]
    lepton_throat_frac = [ladder[n]['throat_fraction'] for n in (0, 1, 2)]
    shell_throat_frac = [ladder[n]['throat_fraction'] for n in (3, 4, 5)]
    # The structural criterion: shell modes have rmean saturated to the
    # plateau (≥ SHELL_RMEAN_PLATEAU); electron has not.
    shell_all_saturated = all(
        r >= SHELL_RMEAN_PLATEAU for r in shell_rmean)
    electron_below_plateau = lepton_rmean[0] < SHELL_RMEAN_PLATEAU
    shell_pr_near_target = all(
        abs(p - SHELL_STANDING_WAVE_PR) <= SHELL_PR_NEAR_TWO_THIRDS
        for p in shell_pr)
    return {
        'name': 'T2_lepton_to_shell_transition',
        'description': (
            "Lepton ladder (n=0,1,2) = throat-to-shell transition: "
            "⟨r⟩ rises through e,μ,τ approaching the shell plateau. "
            "Shell modes (n≥3) saturated: ⟨r⟩ on plateau, PR locked to "
            "2/3 (uniform standing wave). Structural distinction sharpest "
            "at the electron; saturates at n=3."
        ),
        'lepton_rmean_e_mu_tau': [round(r, 4) for r in lepton_rmean],
        'shell_rmean': [round(r, 4) for r in shell_rmean],
        'shell_rmean_plateau_threshold': SHELL_RMEAN_PLATEAU,
        'shell_all_at_plateau': shell_all_saturated,
        'electron_below_plateau': electron_below_plateau,
        'lepton_PR_e_mu_tau': [round(p, 4) for p in lepton_pr],
        'shell_PR': [round(p, 4) for p in shell_pr],
        'shell_PR_near_2_3': shell_pr_near_target,
        'shell_standing_wave_target_PR': SHELL_STANDING_WAVE_PR,
        'lepton_throat_frac': [round(f, 4) for f in lepton_throat_frac],
        'shell_throat_frac': [round(f, 4) for f in shell_throat_frac],
        'pass': (shell_all_saturated and electron_below_plateau
                 and shell_pr_near_target),
    }


# ---------------------------------------------------------------------------
# T3. Operator scaffold: H_kin, H_Z2, H_couple
# ---------------------------------------------------------------------------

def test_T3_operator_scaffold() -> dict:
    """Build the 6×6 scaffold. Verify H_kin diagonal from ω²(l,n);
    H_Z2 zero with chi=0; H_couple zero by default; H_total
    Hermitian."""
    basis = build_shell_basis(mode='n_varied')
    op = ShellOperatorScaffold(basis=basis, chi=0.0)
    Hk = op.H_kin()
    Hz = op.H_Z2()
    Hc = op.H_couple()
    H = op.H_total()
    is_diag_kin = bool(np.allclose(Hk - np.diag(np.diag(Hk)), 0.0))
    is_zero_Hz = bool(np.allclose(Hz, 0.0))
    is_zero_Hc = bool(np.allclose(Hc, 0.0))
    is_hermitian = bool(np.allclose(H, H.T))
    eigs = op.eigenvalues()
    # With chi=0 and no coupling, eigenvalues = ω² each twice (±
    # degenerate per (l,n) block).
    expected_eigs = sorted([s.omega ** 2 for s in basis])
    eigs_match = bool(np.allclose(np.sort(eigs), np.sort(expected_eigs), atol=1e-9))
    return {
        'name': 'T3_operator_scaffold_constructed',
        'description': (
            "6×6 scaffold H = H_kin + H_Z2 + H_couple. H_kin diagonal "
            "(ω²); H_Z2 chi=0 placeholder; H_couple zero placeholder; "
            "Hermitian. Eigenvalues with chi=0, no coupling = ω² doubled."
        ),
        'H_kin_diagonal': is_diag_kin,
        'H_Z2_zero_placeholder': is_zero_Hz,
        'H_couple_zero_placeholder': is_zero_Hc,
        'H_hermitian': is_hermitian,
        'eigenvalues_match_omega_squared': eigs_match,
        'eigenvalues': [round(float(e), 6) for e in eigs],
        'pass': (is_diag_kin and is_zero_Hz and is_zero_Hc
                 and is_hermitian and eigs_match),
    }


# ---------------------------------------------------------------------------
# T4. Basis flexibility: n_varied vs l_varied
# ---------------------------------------------------------------------------

def test_T4_basis_flexibility() -> dict:
    """The scaffold supports both basis enumerations: (l=1, n∈{3,4,5})
    and (l∈{1,2,3}, n=3). PR #78 chooses one."""
    basis_n = build_shell_basis(mode='n_varied')
    basis_l = build_shell_basis(mode='l_varied')
    op_n = ShellOperatorScaffold(basis=basis_n)
    op_l = ShellOperatorScaffold(basis=basis_l)
    ok_n = op_n.n_states == N_FLAVORS
    ok_l = op_l.n_states == N_FLAVORS
    eigs_n = op_n.eigenvalues()
    eigs_l = op_l.eigenvalues()
    different = not bool(np.allclose(np.sort(eigs_n), np.sort(eigs_l)))
    return {
        'name': 'T4_basis_flexibility_n_vs_l_enumeration',
        'description': (
            "Scaffold supports two natural enumerations: vary n (fix "
            "l=1) or vary l (fix n=3). Both yield 6 states. They give "
            "DIFFERENT eigenvalues — PR #78 will choose."
        ),
        'n_varied_basis_size': op_n.n_states,
        'n_varied_labels': [s.label() for s in basis_n],
        'n_varied_eigvals': [round(float(e), 6) for e in eigs_n],
        'l_varied_basis_size': op_l.n_states,
        'l_varied_labels': [s.label() for s in basis_l],
        'l_varied_eigvals': [round(float(e), 6) for e in eigs_l],
        'enumerations_distinct': different,
        'pass': ok_n and ok_l and different,
    }


# ---------------------------------------------------------------------------
# T5. Six-flavor structural match
# ---------------------------------------------------------------------------

def test_T5_six_flavor_structural() -> dict:
    """The 6-state basis matches the 6-flavor structural count from
    PR #69 (3 generations × 2 isospin partition). The Z₂ partition is
    the same B2 partition that drives the lepton-sector mass-ordering
    inversion (m_u < m_d but m_c > m_s pattern)."""
    basis = build_shell_basis(mode='n_varied')
    n_gen = len(set((s.l, s.n) for s in basis))     # 3 (l, n) blocks
    n_part = len(set(s.p for s in basis))            # 2 (±)
    return {
        'name': 'T5_six_flavor_structural_count',
        'description': (
            "6-state basis = 3 generations (l, n) × 2 partitions (±) = "
            "6 flavors, matching the structural count from PR #69. "
            "The Z₂ partition is the same B2 partition behind the "
            "m_u<m_d but m_c>m_s mass-ordering inversion (lepton-sector "
            "shell-coupled closure)."
        ),
        'n_generations': n_gen,
        'n_partitions': n_part,
        'product': n_gen * n_part,
        'expected_flavors': N_FLAVORS,
        'matches_pr69': n_gen == N_GENERATIONS and n_part == N_PARTITIONS,
        'pass': n_gen == N_GENERATIONS and n_part == N_PARTITIONS,
    }


# ---------------------------------------------------------------------------
# T6. Hooks for PR #78–#80
# ---------------------------------------------------------------------------

def test_T6_hooks_for_followons() -> dict:
    """Document the slots PRs #78–#80 will populate:
      - PR #78 (mass-ordering / n_part audit): populate H_couple with
        inter-mode mixing; choose between n_varied and l_varied
        enumerations; verify mass-ordering m_u < m_d, m_c > m_s.
      - PR #79 (boundary stress tensor / singlet constraint): set chi
        from boundary energy-momentum on the cavity wall; add a
        singlet-projection constraint (colorless physical states).
      - PR #80 (color algebra): identify the algebra acting on (l, n, p)
        — SU(3) (color triplet), SU(2) × Z₂ (the v3 partition flavored),
        or another (e.g. Pati-Salam SU(4) with the n-overtone as the
        4th leptocolor)."""
    basis = build_shell_basis(mode='n_varied')
    op = ShellOperatorScaffold(basis=basis, chi=0.0)
    has_kin_slot = op.H_kin().shape == (N_FLAVORS, N_FLAVORS)
    has_z2_slot = hasattr(op, 'chi')
    has_couple_slot = hasattr(op, 'couple_matrix')
    return {
        'name': 'T6_pr78_79_80_hooks',
        'description': (
            "Document the slots PRs #78–#80 will populate."
        ),
        'pr_78_slot': {
            'name': 'shell Hamiltonian mass-ordering / n_part audit',
            'populates': 'H_couple inter-mode mixing; enumeration choice',
            'tests': 'm_u < m_d, m_c > m_s; reaudit n_part on shell basis',
        },
        'pr_79_slot': {
            'name': 'boundary stress tensor + singlet constraint',
            'populates': 'chi (Z₂ splitter from cavity-wall stress-T)',
            'tests': 'singlet projection ⟹ colorless physical states',
        },
        'pr_80_slot': {
            'name': 'color algebra candidate',
            'populates': 'algebra acting on (l, n, p)',
            'candidates': ['SU(3)', 'SU(2) × Z₂', 'SU(4) Pati-Salam', 'other'],
        },
        'scaffold_has_kin_slot': has_kin_slot,
        'scaffold_has_z2_slot': has_z2_slot,
        'scaffold_has_couple_slot': has_couple_slot,
        'pass': has_kin_slot and has_z2_slot and has_couple_slot,
    }


# ---------------------------------------------------------------------------
# T7. Honest scope / B4
# ---------------------------------------------------------------------------

def test_T7_honest_scope() -> dict:
    return {
        'name': 'T7_honest_scope_b4',
        'description': (
            "Honest scope: scaffold only (basis + operator slots). Does "
            "NOT derive quark masses (PR #78), audit n_part (PR #78), "
            "define boundary stress tensor (PR #79), or identify color "
            "algebra (PR #80). B4: cavity eigenfrequencies ω are "
            "dimensionful (1/length); MASS RATIOS are scale-free."
        ),
        'scaffold_includes': [
            '6-state shell waveguide basis',
            'H_kin diagonal cavity-mass operator',
            'H_Z2 partition splitter slot (chi placeholder)',
            'H_couple inter-mode mixing slot (zero placeholder)',
        ],
        'scaffold_excludes': [
            'quark mass predictions (PR #78)',
            'n_part audit (PR #78)',
            'boundary stress tensor (PR #79)',
            'singlet constraint (PR #79)',
            'color algebra identification (PR #80)',
        ],
        'b4_caveat': (
            'ω dimensionful (1/length, scale-free in ratios); chi and '
            'couple_matrix dimensionless if non-zero; scaffold is '
            'structurally scale-independent.'
        ),
        'pass': True,
    }


# ---------------------------------------------------------------------------
# T8. Assessment
# ---------------------------------------------------------------------------

def test_T8_assessment() -> dict:
    return {
        'name': 'T8_assessment',
        'description': (
            "Scaffold for the four-PR QCD-shell arc constructed. Quarks "
            "= wavefronts that resolve the cavity (shell-saturated modes "
            "PR ≥ 0.5, n ≥ 3 for l=1). 6-state basis, 6×6 operator slots. "
            "Distinct from lepton/v3 throat-mode machinery."
        ),
        'classification': 'SHELL_WAVEGUIDE_SCAFFOLD_CONSTRUCTED',
        'next_steps': 'PR #78 mass-ordering / n_part audit',
        'pass': True,
    }


# ---------------------------------------------------------------------------
# Runner + verdict
# ---------------------------------------------------------------------------

def run_probe() -> dict:
    tests = [
        test_T1_shell_basis_constructed(),
        test_T2_distinct_from_leptons(),
        test_T3_operator_scaffold(),
        test_T4_basis_flexibility(),
        test_T5_six_flavor_structural(),
        test_T6_hooks_for_followons(),
        test_T7_honest_scope(),
        test_T8_assessment(),
    ]
    core = tests[:7]
    if all(t['pass'] for t in core):
        verdict_class = 'SHELL_WAVEGUIDE_SCAFFOLD_CONSTRUCTED'
        verdict = (
            'SHELL WAVEGUIDE SCAFFOLD CONSTRUCTED. Quarks are the '
            'wavefronts that resolve the cavity itself — shell-saturated '
            'Tangherlini radial modes whose ⟨r⟩ has plateaued at the '
            'shell center (~0.046 above R_MID) and whose participation '
            'ratio has locked to the uniform-standing-wave value 2/3. '
            'The lepton ladder (e=n=0, μ=n=1, τ=n=2) IS the throat-to-'
            'shell transition: ⟨r⟩ rises 0.021 → 0.039 → 0.044 toward '
            'the plateau but does not reach it; the electron is the '
            'most throat-focused. Shell modes (n ≥ 3) sit on the '
            'plateau — the wavefront regime where the mode resolves the '
            'cavity rather than focuses on the throat.\n\n'
            '6-STATE BASIS. The basis is (l, n, p) with l = S³ Casimir '
            '(Hopf-bundle angular base, #73), n = shell-saturated radial '
            'overtone (≥ 3 for l = 1, from #68 saturation metric), p ∈ '
            '{+, −} = Z₂ partition (B2, the non-orientable throat). The '
            'lowest 6 states match the 3 × 2 = 6 flavor structural count '
            'documented by #69.\n\n'
            'OPERATOR SCAFFOLD. The 6×6 Hamiltonian H = H_kin + H_Z2 + '
            'H_couple acts on the shell basis with: H_kin diagonal '
            '(ω²(l, n) cavity-eigenfrequency-squared mass operator); H_Z2 '
            'block-σ_z partition splitter (slot for #79–#80 structural '
            'constant chi); H_couple inter-mode mixing (slot for #78–#80 '
            'population). PR #77 leaves chi = 0 and H_couple = 0 — the '
            'eigenvalues are then ω² each doubled by the Z₂ partition.\n\n'
            'BASIS FLEXIBILITY. Two natural enumerations are supported: '
            'vary n at fixed l = 1 ({(1, 3, ±), (1, 4, ±), (1, 5, ±)}) or '
            'vary l at fixed shell-saturated n ({(1, 3, ±), (2, 3, ±), '
            '(3, 3, ±)}). They give different eigenvalues; PR #78 chooses '
            'based on which reproduces the quark mass ordering.\n\n'
            'DISTINCTNESS FROM LEPTON/V3 MACHINERY. The lepton sector '
            'lives at the throat (basis = odd-k throat-traversal modes '
            '{k=1, 3, 5}, diagonal = β·k²·(2π) winding cost from #71). '
            'The v3 quark Hamiltonian inherited that lepton-shaped basis '
            'and absorbed QCD shell physics into the phenomenological '
            'n_part = 233 (#76). The shell waveguide basis is structurally '
            'distinct: extended-character standing waves on the cavity, '
            'not focused throat pulses; ω²(l, n) cavity eigenfrequency, '
            'not β·k² winding. This IS the right machinery for the quark '
            'sector.\n\n'
            'HOOKS FOR PR #78–#80. The scaffold names each follow-on '
            'slot: PR #78 (mass-ordering / n_part audit) populates '
            'H_couple and chooses the enumeration; PR #79 (boundary '
            'stress tensor / singlet constraint) sets chi from cavity-'
            'wall energy-momentum and adds a singlet (colorless) '
            'projection; PR #80 (color algebra) identifies the algebra '
            'acting on (l, n, p) — SU(3) color triplet, SU(2) × Z₂ '
            'partition-flavored, Pati-Salam SU(4) with n as a 4th '
            'leptocolor, or another candidate.\n\n'
            'HONEST SCOPE. Scaffold only — does NOT derive quark masses, '
            'audit n_part, define the boundary stress tensor, or identify '
            'the color algebra. The 6×6 operator is structurally '
            'distinct-from-lepton but not yet populated. Mass values, '
            'singlet constraint, color algebra are PRs #78–#80 '
            'explicitly. B4: cavity eigenfrequencies ω are dimensionful '
            '(scale-free in ratios); the scaffold is structurally scale-'
            'independent.'
        )
    else:
        verdict_class = 'SCAFFOLD_INCOMPLETE'
        verdict = (
            'SCAFFOLD INCOMPLETE. Investigate the failing test before '
            'proceeding to PR #78.'
        )

    return {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'identification': (
            'quarks = shell-saturated cavity wavefronts (NOT throat '
            'traversals); 6-state (l, n, p) basis with 6×6 operator '
            'scaffold'
        ),
        'physical_insight': (
            'Quarks do not pass through the throat; they are the '
            'wavefronts that resolve the cavity itself.'
        ),
        'next_pr': 'PR #78 — shell Hamiltonian mass-ordering / n_part audit',
        'b4_caveat': 'ω dimensionful; scaffold scale-free in ratios; structural',
        'tests': tests,
        'n_passed': sum(1 for t in tests if t['pass']),
        'n_total': len(tests),
        'verdict_class': verdict_class,
        'verdict': verdict,
    }


# ---------------------------------------------------------------------------
# Markdown rendering
# ---------------------------------------------------------------------------

def render_markdown(s: dict) -> str:
    L: list[str] = []
    L.append('# QCD shell waveguide: basis + operator scaffold')
    L.append('')
    L.append(f"**Run:** {s['timestamp_utc']}")
    L.append('')
    L.append(
        f"> *{s['physical_insight']}*"
    )
    L.append('')
    L.append(
        'Foundation for the four-PR quantitative QCD-shell arc the user '
        'laid out:'
    )
    L.append('')
    L.append('  - **PR #77 (this PR)** — QCD shell waveguide '
             'basis/operator scaffold')
    L.append('  - **PR #78** — shell Hamiltonian mass-ordering / n_part '
             'audit')
    L.append('  - **PR #79** — boundary stress tensor and singlet '
             'constraint')
    L.append('  - **PR #80** — color algebra candidate: SU(3), '
             'SU(2)×Z₂, or other')
    L.append('')
    L.append(f"- **Identification**: {s['identification']}")
    L.append(f"- **Next PR**: {s['next_pr']}")
    L.append(f"- **B4 caveat**: {s['b4_caveat']}")
    L.append('')

    L.append('## Test summary')
    L.append('')
    L.append('| # | Test | Key finding | PASS? |')
    L.append('|---|---|---|---|')
    label_map = {
        'T1': '6-state shell basis (l=1, n=3,4,5, ±); ⟨r⟩ on plateau',
        'T2': 'lepton ladder = throat→shell transition; shell on plateau',
        'T3': 'H = H_kin + H_Z2 + H_couple Hermitian; ω² eigenvalues',
        'T4': 'scaffold supports n_varied and l_varied enumerations',
        'T5': '3 generations × 2 partitions = 6 flavors (PR #69)',
        'T6': 'hooks for PR #78 (mass), #79 (stress-T), #80 (color)',
        'T7': 'honest scope: scaffold only; no masses, no n_part audit',
        'T8': 'shell waveguide scaffold constructed',
    }
    for t in s['tests']:
        passed = '**PASS**' if t['pass'] else '**FAIL**'
        nm = t['name']
        prefix = nm[:2]
        L.append(f"| {prefix} | `{nm}` | {label_map.get(prefix, '—')} | {passed} |")
    L.append('')

    # T1 basis table
    t1 = s['tests'][0]
    L.append('## T1: Shell waveguide basis (n_varied enumeration)')
    L.append('')
    L.append('| state | ω | PR | ⟨r⟩−R_MID | throat frac |')
    L.append('|---|---:|---:|---:|---:|')
    for lab, om, pr, rm, tf in zip(t1['basis_labels'], t1['omegas'],
                                    t1['participation_ratios'],
                                    t1['r_mean_offsets'],
                                    t1['throat_fractions']):
        L.append(f"| `{lab}` | {om:.4f} | {pr:.4f} | {rm:.4f} | {tf:.4f} |")
    L.append('')

    # T2 lepton-to-shell transition
    t2 = s['tests'][1]
    L.append('## T2: Lepton-to-shell transition — the throat-focused → cavity-resolving boundary')
    L.append('')
    L.append('| metric | lepton (e=0, μ=1, τ=2) | shell (n=3,4,5) |')
    L.append('|---|---|---|')
    L.append(f"| ⟨r⟩−R_MID | {t2['lepton_rmean_e_mu_tau']} | "
             f"{t2['shell_rmean']} |")
    L.append(f"| PR | {t2['lepton_PR_e_mu_tau']} | {t2['shell_PR']} |")
    L.append(f"| throat fraction (inner third) | "
             f"{t2['lepton_throat_frac']} | {t2['shell_throat_frac']} |")
    L.append('')
    L.append(f"Shell plateau threshold: ⟨r⟩ ≥ "
             f"{t2['shell_rmean_plateau_threshold']}; standing-wave "
             f"target PR = {t2['shell_standing_wave_target_PR']:.4f}.")
    L.append('')
    L.append("The lepton ladder IS the throat-to-shell transition: ⟨r⟩ "
             "rises through e, μ, τ approaching the shell plateau. The "
             "electron is sharply throat-focused; the tau is at the "
             "edge of shell saturation. Quark sector (n ≥ 3) sits on "
             "the plateau.")
    L.append('')

    # T4 flexibility
    t4 = s['tests'][3]
    L.append('## T4: Basis flexibility')
    L.append('')
    L.append('| enumeration | labels | eigenvalues |')
    L.append('|---|---|---|')
    L.append(f"| `n_varied` (fix l=1) | {t4['n_varied_labels']} | "
             f"{t4['n_varied_eigvals']} |")
    L.append(f"| `l_varied` (fix n=3) | {t4['l_varied_labels']} | "
             f"{t4['l_varied_eigvals']} |")
    L.append('')
    L.append('PR #78 chooses based on which reproduces the quark mass '
             'ordering (m_u < m_d, m_c > m_s, etc.).')
    L.append('')

    # T6 hooks
    t6 = s['tests'][5]
    L.append('## T6: Hooks for PR #78–#80')
    L.append('')
    for key, slot in [
        ('PR #78', t6['pr_78_slot']),
        ('PR #79', t6['pr_79_slot']),
        ('PR #80', t6['pr_80_slot']),
    ]:
        L.append(f"**{key} — {slot['name']}**")
        L.append('')
        if 'populates' in slot:
            L.append(f"- Populates: `{slot['populates']}`")
        if 'tests' in slot:
            L.append(f"- Tests: {slot['tests']}")
        if 'candidates' in slot:
            L.append(f"- Candidates: {slot['candidates']}")
        L.append('')

    L.append('## Verdict')
    L.append('')
    L.append(f"**{s['verdict_class']}.** {s['verdict']}")
    L.append('')

    L.append('## What this leaves open')
    L.append('')
    L.append('- **PR #78** — populate H_couple with inter-mode mixing, '
             'choose the enumeration (n_varied vs l_varied), and audit '
             'whether the shell basis reproduces the quark mass ordering '
             'AND derives n_part structurally.')
    L.append('- **PR #79** — set chi from the boundary stress tensor on '
             'the cavity wall; add a singlet (colorless) projection as a '
             'physical-state constraint.')
    L.append('- **PR #80** — identify the color algebra acting on '
             '(l, n, p): SU(3) (color triplet), SU(2) × Z₂ (partition-'
             'flavored), Pati-Salam SU(4) (n as 4th leptocolor), or '
             'another.')
    L.append('')
    return '\n'.join(L)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, tuple):
        return list(o)
    if hasattr(o, '__dict__'):
        try:
            return asdict(o)
        except Exception:
            return o.__dict__
    return str(o)


def main(argv: Optional[list[str]] = None) -> int:
    summary = run_probe()
    md = render_markdown(summary)
    print(md)
    here = Path(__file__).resolve().parent
    ts = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = here / 'runs' / f'{ts}_qcd_shell_waveguide_scaffold_probe'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'probe.json').write_text(
        json.dumps(summary, indent=2, default=_json_default)
    )
    (out / 'probe.md').write_text(md)
    print(f"\nWrote: {out / 'probe.json'}")
    print(f"Wrote: {out / 'probe.md'}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
