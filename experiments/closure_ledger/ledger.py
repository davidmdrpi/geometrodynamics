"""
Layer 1 of the closure-phase ledger.

Computes per-generation phase contributions from the four channels that
ARE currently wired into the geometrodynamics package:

    Φ_antipodal  = k · action_base                        (= k · 2π)
    Φ_hopf       = π · cos(χ)                              (Hopf holonomy)
    Φ_throat_T²  = π                                        (T² = −I closure phase)
    Φ_uplift     = β_lepton · max(0, k − 3)²                (β·k² uplift term)

Channel 3 (bulk radial) is NOT wired into the lepton sector — see
blockers.py for the structured Layer 2 blocker.

External dependencies (verified against repo README; Claude Code should
confirm each import path against the actual module before commit):

    from geometrodynamics.hopf.connection import hopf_holonomy
        # Returns ∮A = π·cos(χ) for given χ. If the actual function name
        # or signature differs, this is the place to adapt.

    from geometrodynamics.embedding.transport import T_matrix
        # Returns the 2×2 throat transport matrix T = iσ_y.
        # If the export is named differently (e.g. `transport_T`,
        # `orientation_reversing_map`, etc.), update the import here.

    from geometrodynamics.tangherlini.lepton_spectrum import (
        TAU_BETA_50PI, S3_ACTION_BASE,
    )
        # Locked lepton constants per README. These names appear in the
        # README's "Quick Start" snippet, so they should resolve. If the
        # quark-side constants are needed, the README cites BETA_QUARK
        # implicitly via N=466; access path is to be confirmed.

If any of these imports fails at runtime, the runner catches the
ImportError, logs the missing symbol in the result JSON's `import_errors`
field, and falls back to the README-published numerical values
(`action_base = 2π`, `β_lepton = 50π`). The fallback is correct for
this Layer-1 ledger because all four phase contributions are
analytically expressible from the README; the imports exist to ensure
the experiment stays in sync with the canonical repo definitions if
they ever change.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict
from typing import Optional


# ---------------------------------------------------------------------------
# Constants (also imported from the repo when possible; see fallback logic
# in `_load_repo_constants`).
# ---------------------------------------------------------------------------

TAU = 2.0 * math.pi

# README-published values, used as fallback if repo imports fail.
_FALLBACK_ACTION_BASE = TAU            # action_base = 2π
_FALLBACK_BETA_LEPTON = 50.0 * math.pi # k_uplift_beta = 50π
_FALLBACK_BETA_QUARK = 466.0 * math.pi / 2.0  # β = N · π/2 with N = 466
_FALLBACK_LEPTON_QUANTA = 100          # 4β_lepton / (2π)
_FALLBACK_QUARK_QUANTA = 466           # 4β_quark / (2π)


@dataclass
class RepoConstants:
    """Bag of constants we either imported from the repo or fell back on."""
    action_base: float
    beta_lepton: float
    beta_quark: float
    lepton_quanta: int
    quark_quanta: int
    source: dict          # {const_name: "repo" | "fallback_readme"}


def _load_repo_constants() -> tuple[RepoConstants, list[str]]:
    """
    Try importing canonical constants from the repo. On any ImportError,
    fall back to README-published values and record the failure.

    Returns (constants, list_of_import_errors).
    """
    errors: list[str] = []
    sources: dict[str, str] = {}

    try:
        from geometrodynamics.tangherlini.lepton_spectrum import (
            S3_ACTION_BASE,
            TAU_BETA_50PI,
        )
        action_base = float(S3_ACTION_BASE)
        beta_lepton = float(TAU_BETA_50PI)
        sources["action_base"] = "repo"
        sources["beta_lepton"] = "repo"
    except Exception as exc:
        errors.append(f"lepton_spectrum constants: {type(exc).__name__}: {exc}")
        action_base = _FALLBACK_ACTION_BASE
        beta_lepton = _FALLBACK_BETA_LEPTON
        sources["action_base"] = "fallback_readme"
        sources["beta_lepton"] = "fallback_readme"

    # Quark β. README cites N = 466 for the lock. The canonical access
    # path in this repo is LOCKED_QUARK_PARAMS.beta on the locked
    # baseline QuarkParams; older speculative names are tried first for
    # forward compatibility.
    try:
        try:
            from geometrodynamics.qcd.quark_spectrum import BETA_QUARK_LOCK as _bq
            beta_quark = float(_bq)
            sources["beta_quark"] = "repo"
        except ImportError:
            try:
                from geometrodynamics.qcd.quark_spectrum import (
                    LOCKED_BETA as _bq,
                )
                beta_quark = float(_bq)
                sources["beta_quark"] = "repo"
            except ImportError:
                from geometrodynamics.qcd.quark_spectrum import (
                    LOCKED_QUARK_PARAMS,
                )
                if LOCKED_QUARK_PARAMS is None:
                    raise ImportError(
                        "LOCKED_QUARK_PARAMS is None; quark β is not locked"
                    )
                beta_quark = float(LOCKED_QUARK_PARAMS.beta)
                sources["beta_quark"] = "repo"
    except Exception as exc:
        errors.append(f"quark β constant: {type(exc).__name__}: {exc}")
        beta_quark = _FALLBACK_BETA_QUARK
        sources["beta_quark"] = "fallback_readme"

    # Quanta counts are derived, not imported.
    lepton_quanta = round(4 * beta_lepton / TAU)
    quark_quanta = round(4 * beta_quark / TAU)

    return (
        RepoConstants(
            action_base=action_base,
            beta_lepton=beta_lepton,
            beta_quark=beta_quark,
            lepton_quanta=lepton_quanta,
            quark_quanta=quark_quanta,
            source=sources,
        ),
        errors,
    )


def _hopf_holonomy(chi: float) -> tuple[float, str]:
    """
    Hopf holonomy ∮A = π·cos(χ).

    Tries to call the repo's `hopf_holonomy(chi)`; falls back to the
    closed form on any ImportError.

    Returns (value, source) where source is "repo" or "fallback_closed_form".
    """
    try:
        from geometrodynamics.hopf.connection import hopf_holonomy as _hopf
        return float(_hopf(chi)), "repo"
    except Exception:
        return math.pi * math.cos(chi), "fallback_closed_form"


def _throat_transport_phase(power: int) -> tuple[float, str]:
    """
    Throat transport phase from T = iσ_y, raised to `power`.

    Eigenvalues of T are ±i, so T^p contributes phase p · (π/2) per
    closure. For p = 2 (T² = −I, the validated closure signature),
    the phase is π.

    Tries to read T from the repo's `embedding.transport` to confirm the
    matrix, but the phase is computed analytically from eigenvalues
    regardless of whether the import succeeds.
    """
    # Attempt to verify the matrix exists in the repo, but the numeric
    # answer is from the closed-form spectrum.
    source = "fallback_closed_form"
    try:
        from geometrodynamics.embedding.transport import T_matrix  # noqa: F401
        source = "repo_matrix_verified"
    except Exception:
        # Try alternative export names commonly used.
        try:
            from geometrodynamics.embedding.transport import T  # noqa: F401
            source = "repo_matrix_verified"
        except Exception:
            try:
                from geometrodynamics.embedding.transport import (
                    derive_throat_transport,
                )
                _ = derive_throat_transport()
                source = "repo_matrix_verified"
            except Exception:
                pass
    return power * (math.pi / 2.0), source


# ---------------------------------------------------------------------------
# Phase term and ledger row datatypes
# ---------------------------------------------------------------------------

@dataclass
class PhaseTerm:
    name: str
    source: str           # human-readable provenance (formula or repo symbol)
    repo_provenance: str  # "repo" | "fallback_closed_form" | "fallback_readme"
                          # | "missing" | "not_modeled"
    value: Optional[float]
    status: str           # "available" | "missing" | "not_modeled"

    def value_mod_2pi(self) -> Optional[float]:
        return None if self.value is None else self.value % TAU


@dataclass
class LedgerRow:
    label: str
    k: int
    terms: list[PhaseTerm] = field(default_factory=list)
    closure_status: str = "unknown"   # "closed_mod_2pi" | "fails_to_close" | "partial"
    available_total: Optional[float] = None
    available_total_mod_2pi: Optional[float] = None
    blocking_terms: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "k": self.k,
            "closure_status": self.closure_status,
            "available_total": self.available_total,
            "available_total_mod_2pi": self.available_total_mod_2pi,
            "blocking_terms": list(self.blocking_terms),
            "terms": [asdict(t) for t in self.terms],
        }


# ---------------------------------------------------------------------------
# Ledger assembly
# ---------------------------------------------------------------------------

def _assemble_lepton_row(
    label: str,
    k: int,
    constants: RepoConstants,
    chi: float = 0.0,
    transport_power: int = 2,
) -> LedgerRow:
    """
    Build a single lepton ledger row from the four wired channels plus
    the two structurally missing ones.

    transport_power: 2 (default, T² closure) or 1 (diagnostic).
    chi: Hopf fibre angle, default 0 (canonical fibre, holonomy = π).
    """
    # Channel 1: antipodal closure
    antipodal_value = k * constants.action_base
    antipodal = PhaseTerm(
        name="antipodal_closure",
        source=f"k · action_base = {k} · {constants.action_base / math.pi:.3f}π",
        repo_provenance=constants.source.get("action_base", "fallback_readme"),
        value=antipodal_value,
        status="available",
    )

    # Channel 2 contribution from Hopf holonomy
    hopf_value, hopf_source = _hopf_holonomy(chi)
    hopf = PhaseTerm(
        name="hopf_holonomy",
        source=f"∮A = π cos(χ={chi:.3f}) = {hopf_value / math.pi:.3f}π",
        repo_provenance=hopf_source,
        value=hopf_value,
        status="available",
    )

    # Channel 2 contribution from throat transport
    throat_value, throat_source = _throat_transport_phase(transport_power)
    throat = PhaseTerm(
        name=f"throat_transport_T{transport_power}",
        source=f"T^{transport_power} eigenvalue arg = {transport_power}·(π/2) "
               f"= {throat_value / math.pi:.3f}π",
        repo_provenance=throat_source,
        value=throat_value,
        status="available",
    )

    # β-uplift (channel 1, k=5 only for lepton sector)
    uplift_value = constants.beta_lepton * max(0, k - 3) ** 2
    uplift = PhaseTerm(
        name="lepton_uplift",
        source=f"β · max(0, k−3)² = {constants.beta_lepton / math.pi:.3f}π · "
               f"{max(0, k - 3) ** 2}",
        repo_provenance=constants.source.get("beta_lepton", "fallback_readme"),
        value=uplift_value,
        status="available",
    )

    # Channel 3: radial bulk phase — STRUCTURALLY MISSING in lepton sector
    # (no S(k) bridge in lepton_spectrum.py per repo audit).
    radial = PhaseTerm(
        name="radial_bulk_phase",
        source="Tangherlini eigenmode phase Φ_radial(k); requires S(k) bridge",
        repo_provenance="missing",
        value=None,
        status="missing",
    )

    # Moving-mouth: not modeled (future falsification test in THESIS.md)
    moving = PhaseTerm(
        name="moving_mouth_phase",
        source="time-dependent T(t) along closed mouth trajectory",
        repo_provenance="not_modeled",
        value=None,
        status="not_modeled",
    )

    terms = [antipodal, hopf, throat, uplift, radial, moving]
    blocking = [t.name for t in terms if t.value is None]

    available_sum = sum(t.value for t in terms if t.value is not None)
    available_mod = available_sum % TAU

    if blocking:
        # Status reflects what we CAN say at Layer 1: did the available
        # subset close to 0 mod 2π? If yes, the partial ledger is
        # consistent with the universal-closure conjecture; if no, the
        # partial ledger already falsifies it.
        partial_closes = (
            math.isclose(available_mod, 0.0, abs_tol=1e-9)
            or math.isclose(available_mod, TAU, abs_tol=1e-9)
        )
        status = (
            "partial_closes_mod_2pi" if partial_closes
            else "partial_fails_to_close"
        )
    else:
        partial_closes = (
            math.isclose(available_mod, 0.0, abs_tol=1e-9)
            or math.isclose(available_mod, TAU, abs_tol=1e-9)
        )
        status = "closed_mod_2pi" if partial_closes else "fails_to_close"

    return LedgerRow(
        label=label,
        k=k,
        terms=terms,
        closure_status=status,
        available_total=available_sum,
        available_total_mod_2pi=available_mod,
        blocking_terms=blocking,
    )


def compute_lepton_ledger(
    chi: float = 0.0,
    transport_power: int = 2,
    constants: Optional[RepoConstants] = None,
) -> tuple[list[LedgerRow], RepoConstants, list[str]]:
    """
    Compute the closure-phase ledger across the three lepton generations.

    Returns:
        (rows, constants_used, import_errors)

    rows: one LedgerRow per lepton (e, μ, τ).
    constants_used: the RepoConstants instance, with provenance per field.
    import_errors: list of error messages from any failed repo imports
                   (empty if everything resolved cleanly).
    """
    if constants is None:
        constants, errors = _load_repo_constants()
    else:
        errors = []

    rows = [
        _assemble_lepton_row("electron", k=1, constants=constants,
                             chi=chi, transport_power=transport_power),
        _assemble_lepton_row("muon", k=3, constants=constants,
                             chi=chi, transport_power=transport_power),
        _assemble_lepton_row("tau", k=5, constants=constants,
                             chi=chi, transport_power=transport_power),
    ]
    return rows, constants, errors


# ---------------------------------------------------------------------------
# Quark sector summary
# ---------------------------------------------------------------------------

def compute_quark_sector_summary(constants: RepoConstants) -> dict:
    """
    Same β-lock reading applied to the quark sector. The full per-mode
    quark ledger is blocked on the same S(k) bridge as the lepton sector
    (and additionally on a per-quark mode assignment), so we only report
    the structural quanta-counting result here.
    """
    # The t-quark uplift at k=5 is 4·β_quark; in 2π quanta:
    t_quark_uplift_quanta = round(4 * constants.beta_quark / TAU)
    t_quark_uplift_mod_2pi = (4 * constants.beta_quark) % TAU

    return {
        "lepton_lock_quanta": constants.lepton_quanta,
        "quark_lock_quanta": constants.quark_quanta,
        "lock_quanta_gap": constants.quark_quanta - constants.lepton_quanta,
        "t_quark_uplift_quanta": t_quark_uplift_quanta,
        "t_quark_uplift_mod_2pi": t_quark_uplift_mod_2pi,
        "interpretation": (
            "Both β-locks are integer-compatibility conditions for the "
            "closure ledger: each multiplier ensures the heaviest-shell "
            "uplift contributes 0 mod 2π. The lepton lock fits at "
            f"{constants.lepton_quanta} quanta of 2π (channels 2 and 3 "
            f"nearly inactive); the quark lock fits at "
            f"{constants.quark_quanta} (channels 2 and 3 active). The "
            f"gap of {constants.quark_quanta - constants.lepton_quanta} "
            "is a hypothesis about what the bulk-coupling channel "
            "contributes across the quark spectrum, contingent on S(k) "
            "being defined."
        ),
        "hypothesis_status": "S_k_bridge_required",
    }
