"""
Wormhole-network traversal: the advanced confirmation as a globally
causal, everywhere-future-directed path (PR #216).

``handshake.advanced_confirm_amplitude`` models the Wheeler–Feynman
confirmation *phenomenologically* — a conjugated ("advanced") phase
assigned by fiat.  This module supplies the explicit mechanism: waves
only ever propagate forward in time, but a wormhole network can carry a
retarded wave into the global past.  A mouth absorbs the wave at the
future antipode, the throat transmits it (the PR #215 greybody),
traversal runs *forward* in mouth-clock time, and the far mouth —
desynchronized by differential aging (a deep gravitational well; the
Morris–Thorne–Yurtsever mechanism) — re-emits it at an earlier global
time.  Propagating forward again on S³, the wave intersects the
original emission event.

The central identity (derived, and machine-checked by the tests and by
``experiments/closure_ledger/wormhole_network_confirm_probe.py``):
expressing the locally-evolved phase in global exterior time forces the
throat's global-frame transfer factor to be

    t_AB(w) * U_BA(w) ,   U_BA(w) = e^{i w Delta_BA} * decorations ,

so the projected two-point kernel between absorption and return is the
analytic continuation of retarded propagation to a *negative* exterior
interval — the advanced kernel, with a greybody form factor.  The
clock-offset mouth IS the advanced half of the transaction.

Conventions: exterior S³ of unit radius, unit wave speed, global
monochromatic form ``e^{-i w t}``; mouth-local clocks
``tau_i = clock_rate_i * (t - clock_offset_i)``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional

import numpy as np


# ── Network elements ─────────────────────────────────────────────────────────

@dataclass
class NetworkMouth:
    """One mouth of a network throat, embedded in the S³ exterior.

    Parameters
    ----------
    mouth_id : str
        Label.
    psi : float
        Angular position on the S³ exterior.
    link_id : str
        Identifier of the throat this mouth terminates.
    clock_rate : float
        d(tau)/d(t): local proper-time rate against global exterior
        time (gravitational time dilation of the mouth's environment).
    clock_offset : float
        o_i in ``tau_i = clock_rate * (t - o_i)``.  Differential aging
        of the two mouths of a link (one parked in a deep well) shows
        up as unequal offsets — the time-machine structure of the
        frozen network.
    orientation : int
        Deck orientation (+1/-1) carried through the mouth (the
        #48/#63 Z2 label).
    transfer_phase : float
        Fixed geometric phase acquired crossing this mouth.
    """

    mouth_id: str
    psi: float
    link_id: str
    clock_rate: float = 1.0
    clock_offset: float = 0.0
    orientation: int = +1
    transfer_phase: float = 0.0

    def local_time(self, t: float) -> float:
        return self.clock_rate * (t - self.clock_offset)

    def global_time(self, tau: float) -> float:
        return tau / self.clock_rate + self.clock_offset


@dataclass
class NetworkThroat:
    """A throat linking two mouths — a TWO-PORT scattering device.

    Each mouth carries its own barrier (its ``MouthPort``): a wave
    entering through A's interface partially transmits, bounces between
    the two interior barrier faces (repeated loops), and leaks out
    through B's own interface as an echo train.  The composite
    transmission ``t_AB(w)`` is *derived* from the two ports by summing
    the interior loops (the Fabry–Pérot geometric series), not supplied.

    Parameters
    ----------
    mouth_A, mouth_B : NetworkMouth
        Entry and exit mouths.
    tau_th : float
        Interior transit time between the two ports' reference planes
        (mouth-clock time), > 0: each pass is future-directed in the
        throat's own clock.  A loop-k emergence takes (2k+1) tau_th.
    port_A, port_B : MouthPort
        One-sided scattering data of each mouth's barrier (PR #215
        greybody: t, r_out, r_in per interface).
    """

    mouth_A: NetworkMouth
    mouth_B: NetworkMouth
    tau_th: float
    port_A: "MouthPort"
    port_B: "MouthPort"

    @property
    def delta_BA(self) -> float:
        """Clock desynchronization o_B - o_A (< 0: B deep in the past)."""
        return self.mouth_B.clock_offset - self.mouth_A.clock_offset

    def U_BA(self, w: float) -> complex:
        """Global-frame transfer factor forced by clock continuity.

        A wave entering A at global t and leaving B at global
        t' = t + tau_th + Delta_BA carries local phase e^{-i w tau_th};
        re-expressed against the global form e^{-i w t'} this is an
        extra factor e^{+i w Delta_BA} — the clock offset *is* the
        transfer phase (times the fixed mouth decorations).
        """
        deco = (self.mouth_A.orientation * self.mouth_B.orientation
                * np.exp(1j * (self.mouth_A.transfer_phase
                               + self.mouth_B.transfer_phase)))
        return complex(np.exp(1j * w * self.delta_BA) * deco)

    # ── the derived two-port scattering matrix ───────────────────────

    def _loop_factor(self, w: float) -> complex:
        """One interior round trip: reflect off B's inner face, off A's
        inner face, with the round-trip propagation phase."""
        return complex(self.port_A.r_in(w) * self.port_B.r_in(w)
                       * np.exp(2j * w * self.tau_th))

    def t_AB(self, w: float) -> complex:
        """Composite excess transmission A -> B with ALL interior loops
        summed: t_A t_B / (1 - r_inA r_inB e^{2 i w tau_th}).  This is
        the factor in excess of free interior propagation (transparent
        ports give exactly 1); the free transit phase e^{-i w tau_th}
        is carried by the traversal leg."""
        return complex(self.port_A.t(w) * self.port_B.t(w)
                       / (1.0 - self._loop_factor(w)))

    def r_AA(self, w: float) -> complex:
        """Composite reflection back into A's exterior: the direct
        bounce plus every path that enters, loops k times, and exits
        back through A (reciprocity: interior->exterior transmission of
        a port equals its t)."""
        lf = self._loop_factor(w)
        return complex(self.port_A.r_out(w)
                       + self.port_A.t(w) ** 2 * self.port_B.r_in(w)
                       * np.exp(2j * w * self.tau_th) / (1.0 - lf))

    def loop_expansion(self, w: float, kmax: int) -> list:
        """The echo amplitudes: loop k transmits with
        t_A t_B (r_inA r_inB e^{2 i w tau_th})^k after a local transit
        of (2k+1) tau_th.  Partial sums converge to t_AB(w)."""
        t0 = self.port_A.t(w) * self.port_B.t(w)
        lf = self._loop_factor(w)
        return [complex(t0 * lf ** k) for k in range(kmax + 1)]


@dataclass
class MouthPort:
    """One-sided scattering data of a mouth's barrier interface.

    ``t(w)``      transmission amplitude (either direction — 1D
                  reciprocity makes exterior->interior and
                  interior->exterior transmissions equal);
    ``r_out(w)``  reflection for exterior incidence;
    ``r_in(w)``   reflection for interior incidence (the face the
                  repeated loops bounce off).  Unitarity ties them:
                  |r_in| = |r_out| and r_in = -conj(r_out) t/conj(t).
    """

    t: Callable[[float], complex]
    r_out: Callable[[float], complex]
    r_in: Callable[[float], complex]


def transparent_port() -> MouthPort:
    """A fully open interface: t = 1, no reflections (the single-
    interface limit of PR #216 v1 — no repeated loops)."""
    one = lambda w: 1.0 + 0.0j
    zero = lambda w: 0.0 + 0.0j
    return MouthPort(t=one, r_out=zero, r_in=zero)


# ── Traversal bookkeeping ────────────────────────────────────────────────────

@dataclass
class TraversalLeg:
    """One leg of a network path, with its local-clock ledger."""

    name: str
    t_start: float           # global exterior time at leg start
    t_end: float             # global exterior time at leg end
    local_duration: float    # elapsed time on the leg's own clock (> 0)
    factor: complex          # amplitude/phase factor of the leg


@dataclass
class NetworkTraversal:
    """Result of a full source -> antipode -> throat -> past -> source
    confirmation path."""

    legs: List[TraversalLeg] = field(default_factory=list)
    amp: complex = 1.0 + 0.0j
    t_emit: float = 0.0
    t_absorb: float = 0.0
    t_emerge: float = 0.0
    t_return: float = 0.0

    @property
    def locally_forward(self) -> bool:
        """Every leg future-directed in its own clock."""
        return all(leg.local_duration > 0 for leg in self.legs)

    @property
    def globally_advanced(self) -> bool:
        """Emergence precedes absorption in global time — the network
        reaches the past without any locally-backward propagation."""
        return self.t_emerge < self.t_absorb


def s3_leg(w: float, dist: float, t_start: float,
           name: str = "exterior") -> TraversalLeg:
    """Retarded exterior propagation over angular distance ``dist``."""
    return TraversalLeg(
        name=name,
        t_start=t_start,
        t_end=t_start + dist,
        local_duration=dist,
        factor=complex(np.exp(-1j * w * dist)),
    )


def traverse_throat(throat: NetworkThroat, w: float, t_entry: float,
                    loop: Optional[int] = None) -> TraversalLeg:
    """Forward traversal A -> B.

    ``loop=None`` (default): the composite factor with all interior
    loops coherently summed — the monochromatic steady state; the leg
    is labelled with the primary (k = 0) transit times.

    ``loop=k``: the k-th echo alone — local duration (2k+1) tau_th
    (one crossing plus k interior round trips, every one of them
    future-directed in the throat clock), global exit
    t_entry + (2k+1) tau_th + Delta_BA, amplitude
    t_A t_B (r_inA r_inB e^{2 i w tau_th})^k.
    """
    deco = (throat.mouth_A.orientation * throat.mouth_B.orientation
            * np.exp(1j * (throat.mouth_A.transfer_phase
                           + throat.mouth_B.transfer_phase)))
    if loop is None:
        local = throat.tau_th
        amp = throat.t_AB(w)
        name = "throat (loops summed)"
    else:
        local = (2 * loop + 1) * throat.tau_th
        amp = throat.loop_expansion(w, loop)[loop]
        name = f"throat (echo k={loop})"
    t_exit = t_entry + local + throat.delta_BA
    factor = amp * np.exp(-1j * w * local)
    return TraversalLeg(
        name=name,
        t_start=t_entry,
        t_end=t_exit,
        local_duration=local,
        factor=complex(factor * deco),
    )


def emergence_train(throat: NetworkThroat, w: float, t_entry: float,
                    kmax: int) -> List[TraversalLeg]:
    """The full echo train of one absorption event: emergences at
    t_entry + (2k+1) tau_th + Delta_BA, geometrically damped by the
    interior loop factor, each locally future-directed."""
    return [traverse_throat(throat, w, t_entry, loop=k)
            for k in range(kmax + 1)]


def emergent_frequency(throat: NetworkThroat, w: float) -> float:
    """Global frequency at emergence: w * rate_A / rate_B (the throat is
    stationary in mouth-clock time, so the *local* frequency is
    conserved through traversal)."""
    return w * throat.mouth_A.clock_rate / throat.mouth_B.clock_rate


def closure_offset(d_A: float, d_B: float, tau_th: float) -> float:
    """The clock offset Delta_BA for which the network path returns
    exactly to the emission event (time closure, group level up to the
    throat's Wigner delay): Delta_BA = -(d_A + d_B + tau_th)."""
    return -(d_A + d_B + tau_th)


def network_confirmation(throat: NetworkThroat, w: float,
                         t_emit: float, d_A: float,
                         d_B: float) -> NetworkTraversal:
    """The explicit mechanism replacing
    ``handshake.advanced_confirm_amplitude``: one retarded C-wave,
    emitted at ``t_emit``, absorbed at the future mouth after exterior
    distance ``d_A``, transmitted through the throat (greybody),
    traversed forward, re-emitted through the clock-offset mouth in the
    global past, and propagated forward again over ``d_B`` back to the
    source point.  Every leg is future-directed in its own clock; the
    return can precede the emission only through the mouths' frozen
    clock offset."""
    out = NetworkTraversal(t_emit=t_emit)
    leg1 = s3_leg(w, d_A, t_emit, name="offer: source -> mouth A")
    out.legs.append(leg1)
    out.t_absorb = leg1.t_end

    leg2 = traverse_throat(throat, w, leg1.t_end)
    out.legs.append(leg2)
    out.t_emerge = leg2.t_end

    leg3 = s3_leg(w, d_B, leg2.t_end, name="return: mouth B -> source")
    out.legs.append(leg3)
    out.t_return = leg3.t_end

    out.amp = leg1.factor * leg2.factor * leg3.factor
    return out


def projected_kernel(throat: NetworkThroat, w: float, d_B: float) -> dict:
    """The absorption -> return segment expressed in exterior labels.

    Global displacement:  D' = tau_th + Delta_BA + d_B   (< 0 when the
    offset outweighs the forward legs).  The accumulated factor is

        t_AB(w) * e^{-i w (tau_th + d_B)} * decorations
      = [t_AB(w) * U_BA-decorations] * e^{-i w D'} ,

    i.e. exactly the retarded phase rule analytically continued to the
    negative interval D' — the ADVANCED kernel — times the greybody
    weight.  Returns both sides for machine comparison."""
    d_prime = throat.tau_th + throat.delta_BA + d_B
    kernel = (traverse_throat(throat, w, 0.0).factor
              * np.exp(-1j * w * d_B))
    advanced_phase = complex(np.exp(-1j * w * d_prime))
    weight = kernel / advanced_phase          # the confirmation weight
    return {
        'interval': d_prime,
        'kernel': complex(kernel),
        'advanced_phase': advanced_phase,
        'confirmation_weight': complex(weight),
        'greybody_magnitude': abs(throat.t_AB(w)),
    }
