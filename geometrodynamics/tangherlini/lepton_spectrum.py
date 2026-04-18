"""
Charged-lepton mass spectrum from 5D Tangherlini eigenvalues.

Each charged lepton is modelled as a non-orientable wormhole defect on S³
— a particle/antiparticle pair linked through the embedding dimension
(Wheeler's "mass without mass, charge without charge").  The Tangherlini
radial/S³ eigenproblem then supplies a tower of stable defect frequencies
ω_{l,n} that we identify with the lepton mass ladder via

    m_{l,n} = ω_{l,n} · (ℏc / R_throat)

where R_throat is the single dimensionful length in the geometry.  Because
exactly one length scale controls the whole tower, one lepton mass is used
as the calibration anchor; the remaining two masses are genuine predictions
of the bare geometry.

The purpose of this module is *calibration*, not a claim of agreement.
It quantifies how much of the lepton mass ratios the simplest Tangherlini
defect can explain, and — just as importantly — how much residual structure
(Möbius/non-orientable corrections, Hopf phase sectors, multi-pass throat
knots, condensate dressing) must come from elsewhere to close the gap.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

from geometrodynamics.tangherlini.radial import solve_radial_modes


# ── PDG 2024 charged-lepton masses (MeV/c²) ──────────────────────────────────
PDG_LEPTON_MASSES_MEV: dict[str, float] = {
    "electron": 0.51099895000,
    "muon": 105.6583755,
    "tau": 1776.86,
}

# ℏc in MeV·fm — converts ω/R (inverse length) into an energy.
HBAR_C_MEV_FM: float = 197.3269804


# ── Assignments: which (l, n) is identified with which generation ────────────

@dataclass(frozen=True)
class LeptonAssignment:
    """Identify a charged-lepton generation with a Tangherlini mode."""

    name: str
    l: int
    n: int


# Radial-ladder hypothesis: all three generations share the same l=1
# angular shape (same charge, same spin-½ holonomy) and differ only by
# radial overtone n = 0, 1, 2.  Physically: same defect geometry, more
# internal nodes for higher generations.
DEFAULT_ASSIGNMENT_RADIAL: tuple[LeptonAssignment, ...] = (
    LeptonAssignment("electron", l=1, n=0),
    LeptonAssignment("muon", l=1, n=1),
    LeptonAssignment("tau", l=1, n=2),
)

# Angular-ladder hypothesis: all three generations are throat-ground-states
# (n=0) and differ by S³ angular momentum l = 1, 3, 5 (odd l preserved by
# the orientation-reversing non-orientable throat map).
DEFAULT_ASSIGNMENT_ANGULAR: tuple[LeptonAssignment, ...] = (
    LeptonAssignment("electron", l=1, n=0),
    LeptonAssignment("muon", l=3, n=0),
    LeptonAssignment("tau", l=5, n=0),
)


# ── Report dataclasses ───────────────────────────────────────────────────────

@dataclass(frozen=True)
class LeptonMode:
    name: str
    l: int
    n: int
    omega: float
    mass_pred_mev: float
    mass_pdg_mev: float
    rel_error: float


@dataclass(frozen=True)
class LeptonSpectrumReport:
    assignment: tuple[LeptonAssignment, ...]
    calibration_lepton: str
    calibration_length_fm: float
    modes: tuple[LeptonMode, ...]
    ratios_pred: dict[str, float]
    ratios_pdg: dict[str, float]
    ratios_rel_err: dict[str, float]
    worst_mass_rel_err: float


# ── Solver wrapper ───────────────────────────────────────────────────────────

def _collect_omegas(
    ls: Iterable[int],
    n_max: int,
    N: int,
) -> dict[tuple[int, int], float]:
    """Run solve_radial_modes for each l and collect ω[(l, n)] for n ≤ n_max."""
    table: dict[tuple[int, int], float] = {}
    n_req = n_max + 1
    for l in ls:
        oms, _, _ = solve_radial_modes(l, N=N, n_modes=n_req + 2)
        for n, w in enumerate(oms):
            if n <= n_max:
                table[(int(l), int(n))] = float(w)
    return table


# ── Public API ───────────────────────────────────────────────────────────────

def compute_lepton_spectrum(
    assignment: Sequence[LeptonAssignment] = DEFAULT_ASSIGNMENT_RADIAL,
    calibrate_to: str = "electron",
    N: int = 140,
) -> LeptonSpectrumReport:
    """Map Tangherlini eigenfrequencies onto the charged-lepton ladder.

    The anchor lepton sets the throat length R_throat via

        R_throat [fm] = ω_anchor · ℏc / m_anchor[MeV]

    and every other mode's predicted mass is ω_{l,n} · ℏc / R_throat.
    """
    names = {a.name for a in assignment}
    unknown = names - set(PDG_LEPTON_MASSES_MEV)
    if unknown:
        raise ValueError(f"unknown lepton name(s): {sorted(unknown)}")
    if calibrate_to not in names:
        raise ValueError(
            f"calibrate_to={calibrate_to!r} is not in assignment names {sorted(names)}"
        )

    ls_needed = sorted({a.l for a in assignment})
    n_max = max(a.n for a in assignment)
    omegas = _collect_omegas(ls_needed, n_max=n_max, N=N)
    for a in assignment:
        if (a.l, a.n) not in omegas:
            raise RuntimeError(
                f"solver returned too few eigenmodes for l={a.l}, n={a.n}; "
                f"increase N."
            )

    anchor = next(a for a in assignment if a.name == calibrate_to)
    m_anchor = PDG_LEPTON_MASSES_MEV[anchor.name]
    w_anchor = omegas[(anchor.l, anchor.n)]
    R_throat_fm = w_anchor * HBAR_C_MEV_FM / m_anchor
    mass_scale = HBAR_C_MEV_FM / R_throat_fm  # MeV per unit ω

    modes: list[LeptonMode] = []
    for a in assignment:
        w = omegas[(a.l, a.n)]
        m_pred = w * mass_scale
        m_pdg = PDG_LEPTON_MASSES_MEV[a.name]
        rel = (m_pred - m_pdg) / m_pdg
        modes.append(LeptonMode(a.name, a.l, a.n, w, m_pred, m_pdg, rel))

    ratios_pred: dict[str, float] = {}
    ratios_pdg: dict[str, float] = {}
    ratios_err: dict[str, float] = {}
    for m in modes:
        if m.name == anchor.name:
            continue
        key = f"{m.name}/{anchor.name}"
        ratios_pred[key] = m.omega / w_anchor
        ratios_pdg[key] = m.mass_pdg_mev / m_anchor
        ratios_err[key] = (ratios_pred[key] - ratios_pdg[key]) / ratios_pdg[key]

    worst = max((abs(m.rel_error) for m in modes), default=0.0)

    return LeptonSpectrumReport(
        assignment=tuple(assignment),
        calibration_lepton=anchor.name,
        calibration_length_fm=R_throat_fm,
        modes=tuple(modes),
        ratios_pred=ratios_pred,
        ratios_pdg=ratios_pdg,
        ratios_rel_err=ratios_err,
        worst_mass_rel_err=worst,
    )


def format_report(report: LeptonSpectrumReport) -> str:
    """Human-readable summary of a LeptonSpectrumReport."""
    hdr = (
        f"Lepton spectrum (calibrated to {report.calibration_lepton}):\n"
        f"  Derived throat length R_throat = {report.calibration_length_fm:.4f} fm\n"
        f"{'  name':<12}{'l':>4}{'n':>4}{'omega':>12}"
        f"{'m_pred [MeV]':>16}{'m_PDG [MeV]':>16}{'rel err':>12}"
    )
    rows = [
        f"  {m.name:<10}{m.l:>4d}{m.n:>4d}{m.omega:>12.6f}"
        f"{m.mass_pred_mev:>16.6f}{m.mass_pdg_mev:>16.6f}{m.rel_error:>+12.2%}"
        for m in report.modes
    ]
    ratio_hdr = "  ratios (predicted / PDG / rel err):"
    ratio_rows = [
        f"    {k:<16}{report.ratios_pred[k]:>12.4f}"
        f"{report.ratios_pdg[k]:>12.4f}{report.ratios_rel_err[k]:>+12.2%}"
        for k in report.ratios_pred
    ]
    return "\n".join([hdr, *rows, ratio_hdr, *ratio_rows])
