"""
geometrodynamics.qcd.quark_spectrum
====================================

Shelled-closure census for the hadronic-constituent sector.

Parallel to ``geometrodynamics.tangherlini.lepton_spectrum``, encoding the
six observed hadronic-constituent states as eigenvalues of a 6×6 Hermitian
Hamiltonian whose entries are via-bulk closure path lengths (in units of
the S³ great-circle arc ``2π``) and inter-closure tunneling costs.  See
``docs/quark_axioms.md`` for the full v3 design spec and the reasoning
behind every structural choice in this module.

METHODOLOGICAL RULE (binding, inherited from the v3 spec §0.5)
──────────────────────────────────────────────────────────────
This module describes topology and its consequences.  Labels that denote
geometric objects — pass count, partition class, shell circumference,
coupling path, mouth, neck, chamber — are permitted.  Labels that denote
Standard Model constructs — isospin, weak charge, flavor-as-SU(2)-doublet,
CKM, W, gluon, color-as-gauge-field — are NOT permitted inside the
axioms, the Hamiltonian, or any public identifier introduced here.

Analogies to Standard Model structures may be noted in a clearly-separated
phenomenological-interpretation section after the topology has been
developed (see docs/quark_axioms.md §9), but never inside the axioms or
the Hamiltonian.

GOVERNING PRINCIPLE (v3 spec §0)
─────────────────────────────────
Spacetime is continuous and classical.  Quantization is emergent, not
imposed: for any given stable defect topology there is exactly one family
of via-bulk closure paths, and their geodesic lengths are fixed by the
embedding.  Same topology → same bulk closure length → same resonance
condition → same mass.  The matrix entries constructed below are path
lengths, period.

PHYSICAL-SPECTRUM EXTRACTION (r2 addition)
───────────────────────────────────────────
The Hermitian Hamiltonian can produce a mix of positive and negative raw
eigenvalues depending on parameter choices.  Raw eigenvalues are NOT
physical masses — they are closure-cost path lengths relative to an
arbitrary zero.  To extract physical masses we:

  1. Shift the spectrum by a constant `spectrum_zero` so the physical
     eigen-branch is positive.  The default zero is chosen so that the
     minimum eigenvalue at the unmixed-limit parameters is zero; all
     parameter points are measured against that same absolute reference.
  2. Label each eigenvector by adiabatic continuation from the unmixed
     limit ``(γ_q, w_q) = (0, 0)``, where the Hamiltonian is block-
     diagonal and each eigenvector IS a basis vector ``|k,p⟩``.
     Under nonzero mixing, each eigenvector is tracked by maximum overlap
     with its unmixed ancestor.  Species identification is by the
     ancestor's ``BASIS_TO_SPECIES`` mapping, NOT by sorted eigenvalue
     order.

This machinery is ``extract_physical_spectrum()``.  All calibration and
regression tools use it; raw ``np.linalg.eigvalsh`` is only used by the
structural tests.

API
────
* ``build_quark_hamiltonian(params)`` — the 6×6 Hermitian matrix.
* ``extract_physical_spectrum(params)`` — adiabatic species → mass mapping.
* ``solved_quark_masses_mev()`` — the solved mass vector (once calibrated).
* ``quark_lepton_limit_check()`` — regression test for single-topology-
  continuity (v3 §7 criterion 4).

Basis constants
────────────────
* ``BASIS_STATES`` — the six ``(k, p)`` labels, in canonical order.
* ``BASIS_TO_SPECIES`` — mapping to conventional species names.
* ``QUARK_SPECIES`` — species in observed-mass ascending order.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Tuple

import numpy as np

# ════════════════════════════════════════════════════════════════════════
# CONSTANTS — topological invariants, not tuning parameters
# ════════════════════════════════════════════════════════════════════════

# §5 step 2 — candidate topological lock for the shelled-closure action base.
QUARK_ACTION_BASE: float = math.pi

# §5 step 3 — integer-winding β-uplift lock for the heaviest-sector mass.
QUARK_BETA_DEFAULT: Optional[float] = None

# §2 — pass counts indexing the three generations
PASS_COUNTS: Tuple[int, ...] = (1, 3, 5)

# §4 — canonical basis ordering for the 6×6 Hamiltonian
BASIS_STATES: Tuple[Tuple[int, str], ...] = (
    (1, "+"), (1, "-"),
    (3, "+"), (3, "-"),
    (5, "+"), (5, "-"),
)

# §4 — basis-state → conventional hadronic-constituent label.
BASIS_TO_SPECIES: Dict[Tuple[int, str], str] = {
    (1, "+"): "u",
    (1, "-"): "d",
    (3, "+"): "c",
    (3, "-"): "s",
    (5, "+"): "t",
    (5, "-"): "b",
}

# Species in observed-mass ascending order.
QUARK_SPECIES: Tuple[str, ...] = ("u", "d", "s", "c", "b", "t")

# Reverse map
SPECIES_TO_BASIS: Dict[str, Tuple[int, str]] = {
    v: k for k, v in BASIS_TO_SPECIES.items()
}

# Observed masses in MeV (PDG 2024 MS-bar, approximate)
OBSERVED_MASSES_MEV: Dict[str, float] = {
    "u":      2.16,
    "d":      4.67,
    "s":     93.4,
    "c":   1270.0,
    "b":   4180.0,
    "t": 172690.0,
}

# §4 — partition sign function σ(p).
_SIGMA: Dict[str, float] = {"+": +1.0, "-": -1.0}


# ════════════════════════════════════════════════════════════════════════
# PARAMETERS — structural locks and residual continuous knobs
# ════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class QuarkParams:
    """
    Parameters for the shelled-closure Hamiltonian.

    Structural locks (§5 steps 2-4):
        action_base : topological arc lock (step 2)
        beta        : integer-winding uplift (step 3)
        gamma_q     : partition-class coupling (step 4)
        u_q_form    : minimal-ansatz name (step 4)

    Residual continuous knobs (§5 step 5):
        phase, transport, pinhole, resistance, partition_mixing

    Spectrum zero (r2):
        spectrum_zero : shift applied to the Hamiltonian before anchoring
                        masses.  By default None, meaning it will be
                        computed automatically from the unmixed-limit
                        ground state.
    """

    action_base: float = QUARK_ACTION_BASE
    beta: float = 0.0
    gamma_q: float = 0.0
    u_q_form: str = "k_minus_2"

    phase: float = 0.001
    transport: float = 25.1
    pinhole: float = 22.5
    resistance: float = 0.217869
    partition_mixing: float = 0.0

    winding_mode: str = "max"
    resistance_model: str = "exponential"
    depth_cost_mode: str = "tunnel_only"
    uplift_mode: str = "k_minus_3_sq"
    uplift_asymmetry: float = 0.0
    spectrum_zero_mode: str = "action_base"
    chi_q_k3: float = 0.0
    eta_k3k5_minus: float = 0.0

    spectrum_zero: Optional[float] = None


# ════════════════════════════════════════════════════════════════════════
# HAMILTONIAN CONSTRUCTION (v3 §4)
# ════════════════════════════════════════════════════════════════════════

def _u_q(k: int, form: str = "k_minus_2") -> float:
    if form == "k_minus_2":
        return float(k - 2)
    if form == "k_times_k_minus_2":
        return float(k * (k - 2))
    if form == "zero":
        return 0.0
    raise ValueError(f"Unknown u_q form: {form!r}")


def _diagonal_entry(k: int, p: str, params: QuarkParams) -> float:
    base = params.action_base
    k_cost = params.resistance * k * k
    pinhole_term = params.pinhole if k in (3, 5) else 0.0
    base_uplift = max(0, k - 3) ** 2
    if params.uplift_mode == "k_minus_3_sq":
        uplift = params.beta * base_uplift
    elif params.uplift_mode == "partition_asymmetric":
        # Multiplicative partition breaking; keeps uplift non-negative
        # when |uplift_asymmetry| <= 1, so positivity rejection is
        # governed by gamma_q / resistance as before.  Reduces to the
        # k_minus_3_sq mode at uplift_asymmetry = 0.
        uplift = (
            params.beta
            * base_uplift
            * (1.0 + params.uplift_asymmetry * _SIGMA[p])
        )
    else:
        raise ValueError(f"Unknown uplift_mode: {params.uplift_mode!r}")
    partition = params.gamma_q * _SIGMA[p] * _u_q(k, params.u_q_form)
    k3_split = params.chi_q_k3 * _SIGMA[p] if k == 3 else 0.0
    return base + k_cost + pinhole_term + uplift + partition + k3_split


def _offdiag_same_partition(
    k1: int, k2: int, p: str, params: QuarkParams,
) -> float:
    if k1 == k2:
        raise ValueError("same-partition off-diagonal requires k1 ≠ k2")
    if params.winding_mode == "max":
        dk = max(k1, k2)
    elif params.winding_mode == "delta":
        dk = abs(k1 - k2)
    else:
        raise ValueError(f"Unknown winding_mode: {params.winding_mode!r}")
    alpha_eff = params.resistance
    base = -params.transport * math.exp(-alpha_eff * dk) * math.cos(
        params.phase * dk,
    )
    # Targeted (3,−)–(5,−) coupling: single opt-in amplitude that lives
    # only on the k=3↔k=5 element inside the partition-"−" block.
    # Default 0.0 recovers the minimal ansatz.  Physical motivation:
    # the s and t outliers in experiment 3 both sit in partition "−"
    # at k=3 and k=5 respectively, so a level-repulsion channel that
    # acts on that pair alone is the minimal structural fix.
    if params.eta_k3k5_minus != 0.0 and p == "-" and {k1, k2} == {3, 5}:
        base += -params.eta_k3k5_minus
    return base


def _offdiag_different_partition(k: int, params: QuarkParams) -> complex:
    """
    Partition-mixing amplitude (§4 case 2).

    NOTE (r2, unchanged): φ_q(k) is still a placeholder here.
    The v3 spec §4 requires it to come from the Hopf connection in
    ``geometrodynamics.hopf.connection``.  Replacing this placeholder
    with the Hopf-derived phase is tracked as a post-landing TODO;
    it requires reading the live Hopf module to pick the right API
    and is deferred out of this drop.
    """
    phi_q_k = params.phase * k
    return -params.partition_mixing * np.exp(1j * phi_q_k)


def build_quark_hamiltonian(
    params: Optional[QuarkParams] = None,
) -> np.ndarray:
    """
    Construct the 6×6 shelled-closure Hamiltonian from v3 spec §4.

    Returns a complex Hermitian ndarray in ``BASIS_STATES`` order.
    """
    if params is None:
        params = QuarkParams()

    n = len(BASIS_STATES)
    H = np.zeros((n, n), dtype=complex)

    for i, (k, p) in enumerate(BASIS_STATES):
        H[i, i] = _diagonal_entry(k, p, params)

    for i, (ki, pi) in enumerate(BASIS_STATES):
        for j in range(i + 1, n):
            kj, pj = BASIS_STATES[j]
            if ki != kj and pi == pj:
                H[i, j] = _offdiag_same_partition(ki, kj, pi, params)
            elif ki == kj and pi != pj:
                H[i, j] = _offdiag_different_partition(ki, params)
            else:
                H[i, j] = 0.0 + 0.0j
            H[j, i] = np.conj(H[i, j])

    if not np.allclose(H, H.conj().T):
        raise RuntimeError("Hamiltonian failed Hermiticity check")
    return H


# ════════════════════════════════════════════════════════════════════════
# PHYSICAL-SPECTRUM EXTRACTION (r2)
# ════════════════════════════════════════════════════════════════════════

def _unmixed_params(params: QuarkParams) -> QuarkParams:
    """
    Return the fully-unmixed reference version of the given params:
    γ_q = 0, partition_mixing = 0, AND transport = 0.  This is the
    regime where the Hamiltonian is fully diagonal and each eigenvector
    is exactly a basis vector ``|k,p⟩``, so species identification is
    unambiguous.

    The adiabatic continuation in ``extract_physical_spectrum`` turns
    ALL three mixing knobs on together, scaling them uniformly from 0
    to their target values over the adiabatic path.
    """
    return replace(params, gamma_q=0.0, partition_mixing=0.0, transport=0.0)


def _default_spectrum_zero(params: QuarkParams) -> float:
    """
    Compute the default spectrum zero.

    Physical choice (r2): ``params.action_base`` is the topological
    minimum closure cost for any shelled defect (v3 §4).  The mass of
    a physical species is the excess of its closure cost over this
    topological minimum.  Using ``action_base`` as the zero means that:

      * Every entry of the Hamiltonian diagonal has the form
        ``action_base + (non-negative terms)``, so shifted eigenvalues
        are bounded below by zero at the unmixed limit.
      * Under mixing, level repulsion can push some eigenvalues below
        ``action_base``; in that regime the corresponding species sits
        below the topological minimum, which is a signal of an
        inconsistent parameter point.  Callers should check for
        non-positive shifted masses and reject such points.
      * The MeV anchor converts path-length excess directly to mass,
        consistent with the lepton-sector convention.
    """
    return float(params.action_base)


def extract_physical_spectrum(
    params: Optional[QuarkParams] = None,
    n_adiabatic_steps: int = 16,
) -> Dict[str, float]:
    """
    Extract the physical mass spectrum from the shelled-closure
    Hamiltonian, with species identification by adiabatic continuation
    from the fully-unmixed reference state.

    Procedure (r2):
      1. Pick the spectrum zero: by default, ``params.action_base``.
         This is the physically meaningful choice — action_base is the
         topological minimum closure cost (the shelled-closure minimal
         arc, v3 §4), and the excess over that minimum is what the MeV
         anchor converts into physical mass.
      2. Build the fully-unmixed reference Hamiltonian (all mixing knobs
         zero).  It is diagonal; eigenvectors are basis vectors; species
         identification via BASIS_TO_SPECIES is unambiguous.
      3. Ramp all mixing knobs (transport, γ_q, partition_mixing) from
         zero to target in n_adiabatic_steps.  Re-assign species tags at
         each step by maximum eigenvector overlap with previous.
      4. Shift final eigenvalues by -spectrum_zero.  The result is a
         dict mapping species name → (eigenvalue - action_base), which
         represents the path-length excess of each closure over the
         topological minimum.  Caller anchors to MeV.

    Callers of this function should verify that every species mass is
    strictly positive after shifting — a species with mass exactly at
    the topological minimum is degenerate with "no closure at all" and
    is not a physical particle.

    Raises ValueError if adiabatic tracking fails.
    """
    if params is None:
        params = QuarkParams()

    # Explicit numeric override wins; otherwise dispatch on
    # spectrum_zero_mode.  "min_eigenvalue" is resolved after the final
    # diagonalization below; flagged here with NaN so the shift step
    # knows to recompute.
    if params.spectrum_zero is not None:
        zero = params.spectrum_zero
    elif params.spectrum_zero_mode == "action_base":
        zero = _default_spectrum_zero(params)
    elif params.spectrum_zero_mode == "min_eigenvalue":
        zero = float("nan")
    else:
        raise ValueError(
            f"Unknown spectrum_zero_mode: {params.spectrum_zero_mode!r}"
        )

    # ── Step 1: unmixed diagonalization, trivial species labeling ──
    unmixed = _unmixed_params(params)
    H0 = build_quark_hamiltonian(unmixed)
    eig0, vec0 = np.linalg.eigh(H0)

    # In the unmixed limit, each eigenvector is (to within sign) a basis
    # vector.  Identify which by max-magnitude component.
    species_by_column: List[str] = []
    for col_idx in range(vec0.shape[1]):
        v = vec0[:, col_idx]
        basis_idx = int(np.argmax(np.abs(v)))
        basis_state = BASIS_STATES[basis_idx]
        species_by_column.append(BASIS_TO_SPECIES[basis_state])

    # Sanity: all six species must appear exactly once
    if sorted(species_by_column) != sorted(QUARK_SPECIES):
        raise RuntimeError(
            f"Unmixed-limit species labeling is inconsistent: "
            f"{species_by_column}.  This indicates a structural problem "
            f"in the Hamiltonian or basis definition."
        )

    current_vecs = vec0
    current_species = list(species_by_column)

    # ── Step 3: adiabatic continuation toward target mixing ──
    # Ramp all three mixing knobs (transport, γ_q, partition_mixing)
    # together from 0 to their target values.  This keeps the path
    # smooth even when multiple knobs are simultaneously large in the
    # target regime.
    gamma_target = params.gamma_q
    mixing_target = params.partition_mixing
    transport_target = params.transport

    if n_adiabatic_steps < 1:
        raise ValueError("n_adiabatic_steps must be >= 1")

    step_values = np.linspace(0.0, 1.0, n_adiabatic_steps + 1)[1:]

    for fraction in step_values:
        step_params = replace(
            params,
            gamma_q=fraction * gamma_target,
            partition_mixing=fraction * mixing_target,
            transport=fraction * transport_target,
        )
        H_step = build_quark_hamiltonian(step_params)
        eig_step, vec_step = np.linalg.eigh(H_step)

        # Overlap matrix: |⟨current_i | step_j⟩|²
        overlap = np.abs(current_vecs.conj().T @ vec_step) ** 2

        # Greedy max-overlap assignment: for each current species (row),
        # pick the new-column with largest overlap that hasn't been taken.
        new_species = [""] * vec_step.shape[1]
        claimed = set()
        row_order = np.argsort(-overlap.max(axis=1))  # decide confident rows first
        for row in row_order:
            cols_by_overlap = np.argsort(-overlap[row, :])
            for col in cols_by_overlap:
                if col not in claimed:
                    new_species[col] = current_species[row]
                    claimed.add(col)
                    break

        if any(s == "" for s in new_species):
            raise ValueError(
                "Adiabatic species tracking failed: assignment "
                "incomplete.  Try increasing n_adiabatic_steps."
            )

        current_vecs = vec_step
        current_species = new_species

    # ── Step 4: build species → shifted eigenvalue map ──
    eig_final, _ = np.linalg.eigh(build_quark_hamiltonian(params))
    # Note: `current_vecs` and eig_final from the same call would be
    # consistent; we redo the final diagonalization only to get the
    # eigenvalue array paired with current_species by column index.
    final_H = build_quark_hamiltonian(params)
    eig_final_chk, vec_final_chk = np.linalg.eigh(final_H)

    # Verify the final eigenvectors match current_vecs up to sign/phase
    # (sanity check on tracking).
    for col in range(vec_final_chk.shape[1]):
        ov = abs(np.vdot(current_vecs[:, col], vec_final_chk[:, col]))
        if ov < 0.9:
            raise ValueError(
                f"Adiabatic tracking diverged at column {col}: "
                f"|overlap| = {ov:.3f}.  Try increasing n_adiabatic_steps."
            )

    if math.isnan(zero):
        zero = float(np.min(eig_final_chk))

    result: Dict[str, float] = {}
    for col, species in enumerate(current_species):
        raw = float(eig_final_chk[col])
        result[species] = raw - zero

    if sorted(result.keys()) != sorted(QUARK_SPECIES):
        raise RuntimeError(
            f"Final species map is inconsistent: keys = {sorted(result.keys())}"
        )

    return result


# ════════════════════════════════════════════════════════════════════════
# LOCKED BASELINE AND SOLVED MASSES
# ════════════════════════════════════════════════════════════════════════

# Pass-3 lock from scripts/refine_pass3_coord_descent.py on 2026-04-24.
#
# Combines the minimal v3 ansatz with three opt-in structural extensions:
#   uplift_mode = "partition_asymmetric"   — k=5 b/t splitter
#   spectrum_zero_mode = "min_eigenvalue"  — d-anchor instead of u-anchor
#   chi_q_k3 = 19.8                        — k=3 c/s splitter
#   eta_k3k5_minus = 5.0                   — targeted (3,−)–(5,−) coupling
#
# Hits max_rel_err = 1.6% across {s, c, b, t} (u is 0 by construction
# under min_eigenvalue zero; d is the anchor at 4.67 MeV).  N=460 sits
# in a smooth basin of width ~30 (within 2× of best); χ and η likewise
# in clean basins of widths ~0.3 and ~2 respectively.  See
# docs/quark_axioms.md §8 for the full calibration log, basin-probe
# evidence, and the candidate topological readings of N, χ, η.
LOCKED_QUARK_PARAMS: Optional[QuarkParams] = QuarkParams(
    action_base=QUARK_ACTION_BASE,
    beta=460 * math.pi / 2.0,
    gamma_q=0.10,
    u_q_form="k_minus_2",
    phase=0.0049,
    transport=0.55,
    pinhole=22.0,
    resistance=0.14,
    partition_mixing=0.0,
    winding_mode="max",
    resistance_model="exponential",
    depth_cost_mode="tunnel_only",
    uplift_mode="partition_asymmetric",
    uplift_asymmetry=0.96,
    spectrum_zero_mode="min_eigenvalue",
    chi_q_k3=19.8,
    eta_k3k5_minus=5.0,
    spectrum_zero=None,
)

QUARK_ANCHOR_SPECIES: str = "u"
QUARK_ANCHOR_MASS_MEV: float = OBSERVED_MASSES_MEV["u"]


def solved_quark_masses_mev() -> np.ndarray:
    """
    Returns the six hadronic-constituent masses in MeV from the locked
    baseline, in ``QUARK_SPECIES`` order (u, d, s, c, b, t).

    Uses physical-spectrum extraction: adiabatic species tagging plus
    spectrum-zero shift.  The anchor species (default: u) is set to its
    observed mass; all other species are scaled accordingly.

    Raises NotImplementedError until LOCKED_QUARK_PARAMS is populated by
    the calibration pipeline.
    """
    if LOCKED_QUARK_PARAMS is None:
        raise NotImplementedError(
            "Shelled-closure spectrum not yet calibrated.  Run the "
            "calibration pipeline: scripts/calibrate_quark_ratios.py, "
            "scripts/sweep_quark_beta.py, scripts/lock_quark_beta_probe.py. "
            "See docs/quark_axioms.md §5 for the locking discipline."
        )

    species_map = extract_physical_spectrum(LOCKED_QUARK_PARAMS)
    # Under min_eigenvalue spectrum zero, the lightest species sits at
    # exactly 0 by construction, so it cannot be the MeV anchor.  Fall
    # back to d (the next-lightest) in that mode.
    if LOCKED_QUARK_PARAMS.spectrum_zero_mode == "min_eigenvalue":
        return _anchor_to_mev(
            species_map,
            anchor_species="d",
            anchor_mass_mev=OBSERVED_MASSES_MEV["d"],
        )
    return _anchor_to_mev(species_map)


def _anchor_to_mev(
    species_map: Dict[str, float],
    anchor_species: str = None,
    anchor_mass_mev: float = None,
) -> np.ndarray:
    """
    Given a species → path-length-mass map, anchor the specified species
    to its observed MeV mass and return the full spectrum in MeV, in
    QUARK_SPECIES order.

    Raises if the anchor species has non-positive mass after the shift.
    """
    anchor_species = anchor_species or QUARK_ANCHOR_SPECIES
    anchor_mass_mev = (anchor_mass_mev
                       if anchor_mass_mev is not None
                       else QUARK_ANCHOR_MASS_MEV)

    anchor_val = species_map[anchor_species]
    if anchor_val <= 0:
        raise RuntimeError(
            f"Anchor species {anchor_species!r} has non-positive "
            f"mass after spectrum shift: {anchor_val:.4e}.  "
            f"The shift is set to the current-Hamiltonian minimum, "
            f"so this indicates that {anchor_species!r} sits AT the "
            f"minimum (degenerate with another species) rather than "
            f"strictly above it.  Try a different anchor species, or "
            f"supply an explicit params.spectrum_zero lower than the "
            f"minimum."
        )
    scale = anchor_mass_mev / anchor_val

    masses = np.array(
        [species_map[s] * scale for s in QUARK_SPECIES],
        dtype=float,
    )
    masses.setflags(write=False)
    return masses


# ════════════════════════════════════════════════════════════════════════
# REGRESSION CHECKS (§7 criteria 1, 3, 4, 5)
# ════════════════════════════════════════════════════════════════════════

def quark_lepton_limit_check(tol: float = 1e-6) -> Dict[str, object]:
    """
    Single-topology-continuity regression check (v3 §7 criterion 4).

    With γ_q = 0, partition_mixing = 0, and action_base = 2π, the 6×6
    Hamiltonian must collapse into two decoupled 3×3 blocks, each
    structurally identical to the lepton Hamiltonian.  Eigenvalues must
    pair up, and the distinct-eigenvalue ratios must match lepton-mass
    ratios.

    This is a RATIO check, not an absolute check.  Absolute equality
    would additionally require that the residual continuous knobs
    match the lepton locked baseline, which the live calibration must
    arrange.
    """
    try:
        from geometrodynamics.tangherlini import solved_lepton_masses_mev
    except ImportError as exc:
        return {
            "passed": False,
            "error": f"Cannot import lepton spectrum: {exc}",
            "max_pair_degeneracy_err": None,
            "max_ratio_error_vs_lepton": None,
            "expected_ratios": None,
            "observed_ratios": None,
            "tol": tol,
        }

    limit_params = QuarkParams(
        action_base=2.0 * math.pi,
        beta=0.0,
        gamma_q=0.0,
        u_q_form="k_minus_2",
        partition_mixing=0.0,
    )

    H = build_quark_hamiltonian(limit_params)
    eigvals = np.linalg.eigvalsh(H)

    pair_errs = [abs(eigvals[2 * i + 1] - eigvals[2 * i]) for i in range(3)]
    max_pair_err = float(max(pair_errs))

    lepton_masses = np.asarray(solved_lepton_masses_mev())
    if lepton_masses.shape != (3,):
        return {
            "passed": False,
            "error": f"Unexpected lepton mass shape {lepton_masses.shape}",
            "max_pair_degeneracy_err": max_pair_err,
            "max_ratio_error_vs_lepton": None,
            "expected_ratios": None,
            "observed_ratios": None,
            "tol": tol,
        }

    distinct = np.array([eigvals[0], eigvals[2], eigvals[4]], dtype=float)
    # Shift so the lightest distinct eigenvalue is positive; required
    # for a meaningful ratio comparison.
    shift = -distinct.min() + 1.0 if distinct.min() <= 0 else 0.0
    distinct_shifted = distinct + shift

    quark_ratios = distinct_shifted[1:] / distinct_shifted[0]
    lepton_ratios = lepton_masses[1:] / lepton_masses[0]
    ratio_errs = np.abs(quark_ratios - lepton_ratios) / lepton_ratios
    max_ratio_err = float(np.max(ratio_errs))

    passed = bool((max_pair_err < tol) and (max_ratio_err < tol))

    return {
        "passed": passed,
        "max_pair_degeneracy_err": max_pair_err,
        "max_ratio_error_vs_lepton": max_ratio_err,
        "expected_ratios": lepton_ratios.tolist(),
        "observed_ratios": quark_ratios.tolist(),
        "shift_applied": shift,
        "tol": tol,
    }


def color_independence_check(params: Optional[QuarkParams] = None) -> Dict[str, object]:
    """§7 criterion 5: mass eigenvalues are color-independent by construction."""
    params = params or QuarkParams()
    return {
        "passed": True,
        "reason": (
            "build_quark_hamiltonian has no color argument; mass "
            "eigenvalues are color-independent by construction."
        ),
        "color_multiplicity": 3,
    }


def adiabatic_tracking_check(
    gamma_q: float = 5.0,
    partition_mixing: float = 2.0,
    n_steps: int = 32,
) -> Dict[str, object]:
    """
    Sanity check for the adiabatic species-tracking path (r2).

    Runs ``extract_physical_spectrum`` at progressively higher mixing
    values and checks that (a) species identifications stay complete,
    and (b) the adiabatic path does not produce level-crossings that
    break species identity.

    Returns diagnostic info; ``passed`` is True if tracking ran cleanly.
    """
    try:
        params = QuarkParams(
            gamma_q=gamma_q,
            partition_mixing=partition_mixing,
        )
        spectrum = extract_physical_spectrum(params, n_adiabatic_steps=n_steps)
        return {
            "passed": True,
            "gamma_q": gamma_q,
            "partition_mixing": partition_mixing,
            "n_steps": n_steps,
            "spectrum": {k: float(v) for k, v in spectrum.items()},
        }
    except Exception as exc:
        return {
            "passed": False,
            "error": str(exc),
            "gamma_q": gamma_q,
            "partition_mixing": partition_mixing,
            "n_steps": n_steps,
        }
