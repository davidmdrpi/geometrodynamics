"""
Closed-history framework — the unifying principle.

This module encodes the core axiom of the geometrodynamic program:

    **Only globally closed, conservation-respecting histories become events.**

A history is a set of Events connected by Worldlines on S³, subject to:
  1. Topological closure: every worldline endpoint is an event
  2. Conservation: charges, energy, momentum balance at every vertex
  3. Phase closure: total phase around every loop ≡ 0 or π (mod 2π)
  4. Stationarity: the history has extremal action

This framework unifies:
  - Bell correlations (two-event detection history)
  - QED transactions (offer/confirm handshake history)
  - BH thermodynamics (many-throat condensate history)

through a single closure principle applied to different topologies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Dict, List, Optional, Tuple

import numpy as np


class EventType(Enum):
    """Types of events in a history."""
    CREATION = auto()        # pair creation from vacuum / throat nucleation
    DETECTION = auto()       # measurement / absorption at a detector
    SCATTERING = auto()      # interaction vertex (exchange, deflection)
    ANNIHILATION = auto()    # pair annihilation / throat closure
    EMISSION = auto()        # photon / graviton emission from mouth
    ABSORPTION = auto()      # photon / graviton absorption at mouth


@dataclass
class Event:
    """A single event in a history.

    An event is a point on S³ × time where something happens:
    creation, detection, scattering, emission, absorption, annihilation.
    """

    event_id: int
    event_type: EventType
    p4: np.ndarray               # position on S³ (unit 4-vector)
    t: float                     # coordinate time
    orientation: int = 0         # +1 particle, −1 antiparticle, 0 neutral
    phase: float = 0.0           # accumulated phase at this event
    detector_setting: Optional[float] = None  # θ for detection events
    outcome: Optional[int] = None             # ±1 for detection events
    charges: Dict[str, float] = field(default_factory=dict)  # conserved charges


@dataclass
class Worldline:
    """A worldline connecting two events.

    Carries the topological and dynamical data for transport between
    events: orientation, phase accumulation, and the throat/cavity
    data if the worldline passes through a wormhole.
    """

    start_event_id: int
    end_event_id: int
    orientation: int = +1        # +1 forward, −1 backward (time-reversed)
    phase_accumulated: float = 0.0
    is_throat_transport: bool = False  # passes through non-orientable throat
    geodesic_angle: float = 0.0       # angular separation on S³


@dataclass
class ClosureResult:
    """Result of checking whether a history satisfies global closure."""

    is_closed: bool
    phase_mismatch: float          # total phase error (0 = perfect closure)
    branch: int                    # 0 or 1 (π-branch)
    conservation_error: float      # max charge/energy violation
    weight: float                  # closure weight (Gaussian in mismatch)
    n_events: int
    n_worldlines: int


@dataclass
class History:
    """A complete history: events + worldlines + closure check.

    A history is a candidate physical process.  It becomes actual
    only if it satisfies the closure conditions: topological closure,
    conservation, and phase quantization.
    """

    events: Dict[int, Event] = field(default_factory=dict)
    worldlines: List[Worldline] = field(default_factory=list)

    def add_event(self, event: Event) -> None:
        self.events[event.event_id] = event

    def add_worldline(self, wl: Worldline) -> None:
        self.worldlines.append(wl)

    def total_phase(self) -> float:
        """Sum of all phase contributions around the history loop."""
        phase = 0.0
        for ev in self.events.values():
            phase += ev.phase
        for wl in self.worldlines:
            phase += wl.phase_accumulated
            if wl.is_throat_transport:
                # Non-orientable transport contributes π/2 to the SU(2) phase
                phase += np.pi / 2.0
        return phase

    def check_conservation(self) -> float:
        """Check charge conservation across the history.

        Each event carries unsigned charge magnitudes.  The orientation
        sign determines whether the charge is positive or negative.
        Conservation requires:  Σ (orientation_i × charge_i) = 0
        for each charge type across all events.
        """
        if not self.events:
            return 0.0

        all_charges = set()
        for ev in self.events.values():
            all_charges.update(ev.charges.keys())

        max_err = 0.0
        for charge_type in all_charges:
            total = sum(
                ev.orientation * ev.charges.get(charge_type, 0.0)
                for ev in self.events.values()
            )
            max_err = max(max_err, abs(total))
        return max_err

    def check_closure(self, sigma: float = 0.6) -> ClosureResult:
        """Check all closure conditions for this history.

        Returns a ClosureResult with the closure status, phase
        mismatch, branch assignment, and closure weight.
        """
        total_phase = self.total_phase()
        cons_err = self.check_conservation()

        # Phase closure: total ≡ 0 or π (mod 2π)
        wrapped = ((total_phase + np.pi) % (2.0 * np.pi)) - np.pi
        mm0 = abs(wrapped)
        mmpi = abs(((wrapped - np.pi) + np.pi) % (2.0 * np.pi) - np.pi)

        if mm0 <= mmpi:
            mismatch = mm0
            branch = 0
        else:
            mismatch = mmpi
            branch = 1

        weight = float(np.exp(-mismatch ** 2 / (2.0 * sigma ** 2)))

        # Conservation must also hold
        if cons_err > 0.01:
            weight *= np.exp(-cons_err ** 2 / 0.01)

        is_closed = weight > 0.01

        return ClosureResult(
            is_closed=is_closed,
            phase_mismatch=mismatch,
            branch=branch,
            conservation_error=cons_err,
            weight=weight,
            n_events=len(self.events),
            n_worldlines=len(self.worldlines),
        )


# ── History constructors for specific physics ────────────────────────────────

def make_bell_history(
    theta_a: float,
    theta_b: float,
    outcome_a: int,
    outcome_b: int,
    geodesic_angle: float = np.pi,
) -> History:
    """Construct a Bell measurement history.

    Two-event history: creation → detection_A + detection_B.
    The pair is created at the throat (implicit) and detected at
    two antipodal points with specified detector settings and outcomes.

    The phase contributions are:
      - Geodesic transport: θ_AB / 2 (spin-½)
      - Throat transport: π/2 (from non-orientable wrap)
      - Detector projection: outcome_a × θ_a/2 + outcome_b × θ_b/2
    """
    from geometrodynamics.transaction.s3_geometry import hsp, antipode4
    from geometrodynamics.bell.hopf_phases import (
        geodesic_spin_phase,
        detector_holonomy_phase,
    )

    p4_a = hsp(0.5, 0.3, 0.0)
    p4_b = antipode4(p4_a)

    # Detection events
    ev_a = Event(
        event_id=0,
        event_type=EventType.DETECTION,
        p4=p4_a,
        t=0.0,
        orientation=+1,
        detector_setting=theta_a,
        outcome=outcome_a,
        phase=0.0,  # detector physics is in the SU(2) amplitude, not here
        charges={"electric": 1.0},
    )

    ev_b = Event(
        event_id=1,
        event_type=EventType.DETECTION,
        p4=p4_b,
        t=0.0,
        orientation=-1,
        detector_setting=theta_b,
        outcome=outcome_b,
        phase=0.0,  # detector physics is in the SU(2) amplitude, not here
        charges={"electric": 1.0},
    )

    # Worldline through the throat
    wl = Worldline(
        start_event_id=0,
        end_event_id=1,
        orientation=+1,
        phase_accumulated=geodesic_spin_phase(geodesic_angle),
        is_throat_transport=True,
        geodesic_angle=geodesic_angle,
    )

    h = History()
    h.add_event(ev_a)
    h.add_event(ev_b)
    h.add_worldline(wl)
    return h


def make_transaction_history(
    src_p4: np.ndarray,
    dst_p4: np.ndarray,
    q_src: float,
    q_dst: float,
    t_emit: float,
    t_absorb: float,
) -> History:
    """Construct a Wheeler-Feynman transaction history.

    Four-event history: emission → propagation → absorption → advanced_return.
    The closed loop requires phase closure for the full round trip.
    """
    from geometrodynamics.transaction.s3_geometry import geo4

    theta = geo4(src_p4, dst_p4)

    ev_emit = Event(
        event_id=0,
        event_type=EventType.EMISSION,
        p4=src_p4,
        t=t_emit,
        orientation=+1,
        charges={"electric": abs(q_src)},
    )

    ev_absorb = Event(
        event_id=1,
        event_type=EventType.ABSORPTION,
        p4=dst_p4,
        t=t_absorb,
        orientation=-1,
        charges={"electric": abs(q_dst)},
    )

    # Retarded worldline
    wl_ret = Worldline(
        start_event_id=0,
        end_event_id=1,
        orientation=+1,
        phase_accumulated=theta / 2.0,  # spin transport
        geodesic_angle=theta,
    )

    # Advanced worldline (time-reversed return)
    wl_adv = Worldline(
        start_event_id=1,
        end_event_id=0,
        orientation=-1,
        phase_accumulated=-theta / 2.0,
        geodesic_angle=theta,
    )

    h = History()
    h.add_event(ev_emit)
    h.add_event(ev_absorb)
    h.add_worldline(wl_ret)
    h.add_worldline(wl_adv)
    return h


# ── Branch enumeration for Bell ──────────────────────────────────────────────

@dataclass
class BranchResult:
    """One branch of a Bell measurement."""
    outcome_a: int
    outcome_b: int
    closure: ClosureResult
    probability: float


def enumerate_bell_branches(
    theta_a: float,
    theta_b: float,
    geodesic_angle: float = np.pi,
) -> List[BranchResult]:
    """Enumerate all four outcome branches for a Bell measurement.

    For each (outcome_a, outcome_b) ∈ {±1}²:
    1. Compute the SU(2) amplitude from throat transport + detectors
    2. Construct the history and check phase closure
    3. Branch weight = |amplitude|² × closure_weight

    The SU(2) amplitude comes from the Bell analyzer module (which
    uses the Hopf-derived throat transport T = iσ_y).  The closure
    weight comes from the history framework.  Together they produce
    correlations that are both geometrically derived AND globally
    closed.

    Returns sorted by probability (highest first).
    """
    from geometrodynamics.bell.analyzers import singlet_amplitude

    branches = []
    total_weight = 0.0

    for sa in [+1, -1]:
        for sb in [+1, -1]:
            # SU(2) amplitude from throat transport + detector projection
            amp = singlet_amplitude(theta_a, theta_b, sa, sb)
            amp_sq = float(abs(amp) ** 2)

            # History closure check
            h = make_bell_history(theta_a, theta_b, sa, sb, geodesic_angle)
            cl = h.check_closure()

            # Combined weight: amplitude × closure
            weight = amp_sq * cl.weight

            branches.append(BranchResult(
                outcome_a=sa,
                outcome_b=sb,
                closure=cl,
                probability=weight,
            ))
            total_weight += weight

    # Normalise
    if total_weight > 1e-30:
        for b in branches:
            b.probability /= total_weight

    branches.sort(key=lambda b: b.probability, reverse=True)
    return branches
