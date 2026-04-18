"""
S³ antipodal worldline-crossing visualisation (v10).

Continuous paths on S³, discrete transaction events at isolated
worldline ↔ antipodal-trace crossings.  Each momentum packet fires only
when a *triple gate* aligns:

  (1) cavity amplitude ``|bₙ| > b_min``
  (2) phase closure ``ΔΦₙ mod 2π`` near ``0`` or ``π``
  (3) worldline crossing ``χᵢⱼ(t,t') = arccos(xⱼ(t)·x̄ᵢ(t')) < ε_cross``

The discreteness is *geometric*: two smooth 1D worldlines intersect the
antipodal map of each other only at isolated instants (codimension-3 in
(t,t')-space).  No ``ℏ`` is inserted — quantisation drops out of the S³
topology.

The :class:`AntipodalCrossingSim` class runs headlessly and is directly
testable.  :func:`run_animation` drives a matplotlib :class:`FuncAnimation`
view with the main S³ scene plus χ/energy side panels.

Refactored from the JSX ``AntipodalCrossing`` prototype (v10).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from geometrodynamics.transaction.s3_geometry import antipode4, geo4, hsp


# ── Physical / gate constants ────────────────────────────────────────────────
C: float = 1.0                  # ring propagation speed on S³
GAMMA: float = 0.08             # cavity decay per bounce
AMP_MIN: float = 0.06           # cavity amplitude floor
AMP_EMIT: float = 0.90          # initial ring amplitude

MATCH_TOL: float = 0.50         # geodesic tolerance for antipodal proximity
EPS_CROSS: float = 0.38         # worldline-crossing tolerance (radians)
EPS_PHASE: float = 0.55         # phase-closure tolerance (radians)

EMIT_CD: float = np.pi * 0.75   # emission cooldown
TX_CD: float = np.pi * 0.65     # post-TX absorber cooldown
N_BOUNCE: int = 9               # max bounces before a ring is killed
TRAIL_LEN: int = 55             # worldline history length (samples)
TRAIL_DT: float = 0.10          # worldline sampling cadence

DEFAULT_DT: float = 0.018


# ── Default particle set: two near-antipodal pairs ───────────────────────────
_DEFAULT_PARTICLES: Tuple[dict, ...] = (
    dict(pid=0, chi=1.20, th=1.08, ph=0.32,
         dchi= 0.0026, dth= 0.0030, dph= 0.0050,
         color="#38d8ff", label="matter A"),
    dict(pid=1, chi=1.94, th=2.06, ph=3.46,
         dchi=-0.0022, dth=-0.0026, dph=-0.0044,
         color="#ff6040", label="anti-A"),
    dict(pid=2, chi=1.52, th=1.62, ph=1.28,
         dchi= 0.0018, dth= 0.0022, dph= 0.0054,
         color="#70ff38", label="matter B"),
    dict(pid=3, chi=1.62, th=1.52, ph=4.42,
         dchi=-0.0016, dth=-0.0020, dph=-0.0048,
         color="#e050ff", label="anti-B"),
)


def w_cross(chi: float) -> float:
    """Gaussian crossing weight ``exp(-χ²/2ε²)``."""
    return float(np.exp(-(chi * chi) / (2.0 * EPS_CROSS * EPS_CROSS)))


def w_phase(phi: float) -> float:
    """Phase-closure weight: Gaussian peak at the 0/π branch."""
    p0 = phi % (2 * np.pi)
    pm = min(p0, abs(np.pi - p0), abs(2 * np.pi - p0))
    return float(np.exp(-(pm * pm) / (2.0 * EPS_PHASE * EPS_PHASE)))


@dataclass
class _Particle:
    pid: int
    chi: float
    th: float
    ph: float
    dchi: float
    dth: float
    dph: float
    color: str
    label: str
    p4: np.ndarray = field(default_factory=lambda: np.zeros(4))
    emit_cd: float = 0.0
    tx_cd: float = 0.0
    energy: float = 1.0
    trail: List[Tuple[float, np.ndarray]] = field(default_factory=list)
    last_trail_t: float = -1.0e9

    def refresh_p4(self) -> None:
        self.p4 = hsp(self.chi, self.th, self.ph)


@dataclass
class _Ring:
    rid: int
    src_id: int
    p_front: np.ndarray
    geo_r: float = 0.02
    bounce: int = 0
    amp: float = AMP_EMIT
    omega: float = 1.0
    phase: float = 0.0
    hit_set: set = field(default_factory=set)
    done: bool = False
    t_done: Optional[float] = None


@dataclass
class Crossing:
    """A recorded candidate worldline/antipodal-trace crossing."""

    t: float
    xj: np.ndarray
    xi_bar: np.ndarray
    src_id: int
    dst_id: int
    chi: float
    weight: float
    color: str


@dataclass
class Absorption:
    """A fired discrete transaction event."""

    t: float
    p4: np.ndarray
    src_id: int
    dst_id: int
    color: str
    bounce: int
    amp: float
    chi: float
    weight: float


# ═════════════════════════════════════════════════════════════════════════════
# Headless simulation
# ═════════════════════════════════════════════════════════════════════════════


class AntipodalCrossingSim:
    """Pure-Python simulation of the v10 triple-gate transaction model.

    Use :meth:`step` to advance by ``dt``.  Inspect :attr:`particles`,
    :attr:`rings`, :attr:`crossings`, :attr:`absorptions` and
    :attr:`energy_history` between steps for rendering or diagnostics.

    Parameters
    ----------
    seed
        Seed for the internal PRNG (controls the random kicks on TX).
    particles
        Sequence of keyword dicts matching :class:`_Particle`.  Defaults to
        the two near-antipodal pairs from the reference prototype.
    """

    def __init__(
        self,
        seed: int = 0,
        particles: Optional[Sequence[dict]] = None,
    ) -> None:
        self._rng = np.random.default_rng(seed)
        specs = particles if particles is not None else _DEFAULT_PARTICLES
        self.particles: List[_Particle] = []
        for spec in specs:
            p = _Particle(**spec)
            p.refresh_p4()
            p.emit_cd = float(self._rng.random() * 1.5)
            self.particles.append(p)
        self.t: float = 0.0
        self.rings: List[_Ring] = []
        self.crossings: List[Crossing] = []
        self.absorptions: List[Absorption] = []
        self.energy_history: List[Tuple[float, float, float]] = []
        self._ring_counter: int = 0

    # ── Stepping ─────────────────────────────────────────────────────────
    def step(self, dt: float = DEFAULT_DT) -> None:
        self.t += dt
        self._move_particles(dt)
        self._emit_rings()
        self._propagate_rings(dt)
        self._record_energy()

    def run(self, n_steps: int, dt: float = DEFAULT_DT) -> None:
        for _ in range(n_steps):
            self.step(dt)

    # ── Diagnostics ──────────────────────────────────────────────────────
    def tx_count(self) -> int:
        return len(self.absorptions)

    def best_chi_now(self) -> float:
        """Smallest current-instant crossing angle across ordered pairs."""
        best = np.pi
        for j in self.particles:
            for i in self.particles:
                if i.pid == j.pid or not i.trail:
                    continue
                for _tp, xi_tp in i.trail:
                    chi = geo4(j.p4, antipode4(xi_tp))
                    if chi < best:
                        best = chi
        return float(best)

    # ── Internals ────────────────────────────────────────────────────────
    def _move_particles(self, dt: float) -> None:
        for p in self.particles:
            p.chi = float(np.clip(p.chi + p.dchi * dt, 0.10, np.pi - 0.10))
            p.th = float(np.clip(p.th + p.dth * dt, 0.06, np.pi - 0.06))
            p.ph = float((p.ph + p.dph * dt) % (2 * np.pi))
            if p.chi < 0.12 or p.chi > np.pi - 0.12:
                p.dchi = -p.dchi
            if p.th < 0.08 or p.th > np.pi - 0.08:
                p.dth = -p.dth
            p.refresh_p4()
            p.emit_cd = max(0.0, p.emit_cd - dt)
            p.tx_cd = max(0.0, p.tx_cd - dt)
            if self.t - p.last_trail_t >= TRAIL_DT:
                p.trail.append((self.t, p.p4.copy()))
                if len(p.trail) > TRAIL_LEN:
                    p.trail.pop(0)
                p.last_trail_t = self.t

    def _emit_rings(self) -> None:
        for A in self.particles:
            if A.emit_cd > 0:
                continue
            near_anti = any(
                B.pid != A.pid and geo4(B.p4, antipode4(A.p4)) < MATCH_TOL
                for B in self.particles
            )
            if not near_anti:
                continue
            self.rings.append(
                _Ring(
                    rid=self._ring_counter,
                    src_id=A.pid,
                    p_front=A.p4.copy(),
                    omega=1.0 + (self._ring_counter % 3),
                )
            )
            self._ring_counter += 1
            A.emit_cd = EMIT_CD

    def _propagate_rings(self, dt: float) -> None:
        for ring in self.rings:
            if ring.done:
                continue
            ring.geo_r += C * dt
            ring.phase += ring.omega * dt
            if ring.geo_r < np.pi:
                continue

            anti_pt = antipode4(ring.p_front)
            absorber = self._find_absorber(ring, anti_pt)

            amp_ok = ring.amp > AMP_MIN * 2
            phase_ok = w_phase(ring.phase) > np.exp(-0.5 * EPS_PHASE ** 2 * 0.25)

            cross_data = None
            if absorber is not None and amp_ok:
                src = self.particles[ring.src_id]
                cross_data = self._best_crossing(absorber, src, ring)
                if cross_data is not None and cross_data["chi"] < EPS_CROSS * 2.5:
                    self.crossings.append(
                        Crossing(
                            t=self.t,
                            xj=absorber.p4.copy(),
                            xi_bar=cross_data["xi_bar"].copy(),
                            src_id=ring.src_id,
                            dst_id=absorber.pid,
                            chi=cross_data["chi"],
                            weight=cross_data["w"],
                            color=absorber.color,
                        )
                    )
            self.crossings = [c for c in self.crossings if self.t - c.t < 1.5]

            fired = (
                absorber is not None
                and amp_ok
                and phase_ok
                and cross_data is not None
                and cross_data["chi"] < EPS_CROSS
            )
            if fired:
                self._fire_tx(ring, absorber, cross_data)  # type: ignore[arg-type]
            else:
                self._bounce_ring(ring, absorber)

        self.rings = [
            r for r in self.rings
            if not r.done or self.t - (r.t_done or self.t) < 0.7
        ]

    def _find_absorber(self, ring: _Ring, anti_pt: np.ndarray) -> Optional[_Particle]:
        for P in self.particles:
            if P.pid == ring.src_id or P.tx_cd > 0 or P.pid in ring.hit_set:
                continue
            if geo4(P.p4, anti_pt) < MATCH_TOL:
                return P
        return None

    def _best_crossing(
        self, pj: _Particle, pi: _Particle, ring: _Ring,
    ) -> Optional[dict]:
        if not pi.trail:
            return None
        best: Optional[dict] = None
        best_w = 0.0
        xj_now = pj.p4
        wp = w_phase(ring.phase)
        for _t, xi_tp in pi.trail:
            xi_bar = antipode4(xi_tp)
            chi = geo4(xj_now, xi_bar)
            if chi >= EPS_CROSS * 2.5:
                continue
            wc = w_cross(chi)
            w = wc * wp * ring.amp
            if w > best_w:
                best_w = w
                best = dict(chi=chi, w_cross=wc, w_phase=wp, xi_bar=xi_bar, w=w)
        return best if best is not None and best_w > 0.05 else None

    def _fire_tx(self, ring: _Ring, absorber: _Particle, cd: dict) -> None:
        ring.done = True
        ring.t_done = self.t
        W = cd["w_cross"] * cd["w_phase"]
        kick = ring.amp * W * 0.24
        rng = self._rng
        absorber.dph += (rng.random() - 0.5) * kick * 2
        absorber.dth += (rng.random() - 0.5) * kick
        absorber.dchi += (rng.random() - 0.5) * kick * 0.5
        absorber.energy = min(2.5, absorber.energy + ring.amp * W * 0.18)
        src = self.particles[ring.src_id]
        src.energy = max(0.3, src.energy - ring.amp * W * 0.14)
        absorber.tx_cd = TX_CD
        self.absorptions.insert(
            0,
            Absorption(
                t=self.t,
                p4=absorber.p4.copy(),
                src_id=ring.src_id,
                dst_id=absorber.pid,
                color=absorber.color,
                bounce=ring.bounce,
                amp=ring.amp,
                chi=cd["chi"],
                weight=W,
            ),
        )
        self.absorptions = self.absorptions[:32]

    def _bounce_ring(self, ring: _Ring, absorber: Optional[_Particle]) -> None:
        if absorber is not None:
            ring.hit_set.add(absorber.pid)
        ring.p_front = antipode4(ring.p_front)
        ring.geo_r = 0.02
        ring.bounce += 1
        ring.amp *= 1.0 - GAMMA
        ring.hit_set.clear()
        if ring.bounce > N_BOUNCE or ring.amp < AMP_MIN:
            ring.done = True
            ring.t_done = self.t

    def _record_energy(self) -> None:
        e_cav = float(
            sum(
                0.5 * r.amp ** 2 * (1 + r.omega ** 2)
                for r in self.rings
                if not r.done
            )
        )
        e_kin = float(sum(0.5 * p.energy for p in self.particles))
        self.energy_history.append((self.t, e_cav, e_kin))
        if len(self.energy_history) > 600:
            self.energy_history.pop(0)


# ═════════════════════════════════════════════════════════════════════════════
# Matplotlib renderer
# ═════════════════════════════════════════════════════════════════════════════


def _stereo(q: np.ndarray, clip: float = 3.6) -> Optional[np.ndarray]:
    """Stereographic projection S³ → R³ from the north pole (0,0,0,1)."""
    d = 1.0 - q[3]
    if abs(d) < 0.08:
        return None
    p = np.array([q[0] / d, q[1] / d, q[2] / d])
    if np.dot(p, p) > clip * clip:
        return None
    return p


def run_animation(
    sim: Optional[AntipodalCrossingSim] = None,
    n_steps: int = 4000,
    interval_ms: int = 16,
    show: bool = True,
):
    """Animate the simulation with matplotlib.

    Parameters
    ----------
    sim
        An existing :class:`AntipodalCrossingSim`.  A fresh one is built if
        ``None`` is passed.
    n_steps
        Number of animation frames to schedule.
    interval_ms
        Frame interval in milliseconds.
    show
        If ``True``, call :func:`plt.show`.  Set to ``False`` when driving
        the animation programmatically (e.g. saving to a file).

    Returns
    -------
    (matplotlib.animation.FuncAnimation, AntipodalCrossingSim)
    """
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    sim = sim if sim is not None else AntipodalCrossingSim()

    fig = plt.figure(figsize=(12, 7), facecolor="#010106")
    gs = fig.add_gridspec(
        3, 2, width_ratios=(2.2, 1.0), height_ratios=(3.0, 1.0, 1.2),
        hspace=0.35, wspace=0.20,
        left=0.04, right=0.98, top=0.94, bottom=0.06,
    )
    ax3d = fig.add_subplot(gs[:, 0], projection="3d")
    ax_chi = fig.add_subplot(gs[0, 1])
    ax_e = fig.add_subplot(gs[1, 1])
    ax_log = fig.add_subplot(gs[2, 1])

    for ax in (ax_chi, ax_e, ax_log):
        ax.set_facecolor("#02020a")
        for spine in ax.spines.values():
            spine.set_color("#1a2838")
        ax.tick_params(colors="#6a8aad", labelsize=7)

    ax3d.set_facecolor("#02020a")
    ax3d.xaxis.pane.set_edgecolor("#0a1422")
    ax3d.yaxis.pane.set_edgecolor("#0a1422")
    ax3d.zaxis.pane.set_edgecolor("#0a1422")
    ax3d.xaxis.pane.fill = False
    ax3d.yaxis.pane.fill = False
    ax3d.zaxis.pane.fill = False
    ax3d.grid(False)
    for axis in (ax3d.xaxis, ax3d.yaxis, ax3d.zaxis):
        axis.line.set_color("#0a1422")
    ax3d.tick_params(colors="#23405f", labelsize=6)

    # Reference S³ equator (stereographically the unit sphere at origin)
    u = np.linspace(0, 2 * np.pi, 60)
    v = np.linspace(0, np.pi, 30)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))
    ax3d.plot_wireframe(
        sx, sy, sz, color="#1a3a75", linewidth=0.4, alpha=0.18,
        rcount=10, ccount=16,
    )

    title = fig.suptitle(
        "S³ antipodal worldline crossing  —  v10   "
        "χᵢⱼ(t,t') = arccos(xⱼ(t)·x̄ᵢ(t'))",
        color="#c8b018", fontsize=11, y=0.98,
    )
    title.set_alpha(0.92)

    def _draw() -> None:
        # ── 3D scene ─────────────────────────────────────────────────
        ax3d.cla()
        ax3d.set_facecolor("#02020a")
        ax3d.set_xlim(-2.8, 2.8)
        ax3d.set_ylim(-2.8, 2.8)
        ax3d.set_zlim(-2.8, 2.8)
        ax3d.set_xticks([])
        ax3d.set_yticks([])
        ax3d.set_zticks([])
        ax3d.set_box_aspect((1, 1, 1))
        ax3d.grid(False)
        for axis in (ax3d.xaxis, ax3d.yaxis, ax3d.zaxis):
            axis.pane.set_edgecolor("#0a1422")
            axis.pane.fill = False
            axis.line.set_color("#0a1422")
        ax3d.plot_wireframe(
            sx, sy, sz, color="#1a3a75", linewidth=0.4, alpha=0.18,
            rcount=10, ccount=16,
        )

        for p in sim.particles:
            if len(p.trail) >= 2:
                pts = [_stereo(q) for _t, q in p.trail]
                # Worldline (solid)
                seg_x, seg_y, seg_z = [], [], []
                for pt in pts:
                    if pt is None:
                        seg_x.append(np.nan)
                        seg_y.append(np.nan)
                        seg_z.append(np.nan)
                    else:
                        seg_x.append(pt[0]); seg_y.append(pt[1]); seg_z.append(pt[2])
                ax3d.plot(seg_x, seg_y, seg_z, color=p.color, alpha=0.55, linewidth=1.2)
                # Antipodal trace (dashed)
                apts = [_stereo(antipode4(q)) for _t, q in p.trail]
                ax_, ay_, az_ = [], [], []
                for pt in apts:
                    if pt is None:
                        ax_.append(np.nan); ay_.append(np.nan); az_.append(np.nan)
                    else:
                        ax_.append(pt[0]); ay_.append(pt[1]); az_.append(pt[2])
                ax3d.plot(
                    ax_, ay_, az_,
                    color=p.color, alpha=0.25, linewidth=0.9, linestyle="--",
                )

            cur = _stereo(p.p4)
            if cur is not None:
                ax3d.scatter(
                    [cur[0]], [cur[1]], [cur[2]],
                    color=p.color, s=45 + p.energy * 18,
                    edgecolors="white", linewidths=0.4,
                    depthshade=False,
                )
            anti = _stereo(antipode4(p.p4))
            if anti is not None:
                ax3d.scatter(
                    [anti[0]], [anti[1]], [anti[2]],
                    facecolors="none", edgecolors=p.color,
                    s=30, alpha=0.35, linewidths=0.6,
                    depthshade=False,
                )

        # Crossings
        for cr in sim.crossings:
            pt = _stereo(cr.xj)
            if pt is None:
                continue
            age = max(0.0, 1.0 - (sim.t - cr.t) / 1.5)
            proximity = max(0.0, 1.0 - cr.chi / EPS_CROSS)
            alpha = min(1.0, 0.25 + 0.75 * proximity * age)
            ax3d.scatter(
                [pt[0]], [pt[1]], [pt[2]],
                marker="o", s=80 + 240 * proximity,
                facecolors="none",
                edgecolors=(1.0, 1.0, 0.82, alpha),
                linewidths=1.4 + 1.4 * proximity,
                depthshade=False,
            )

        # Ring front + antipodal target
        for r in sim.rings:
            if r.done:
                continue
            fr = _stereo(r.p_front)
            if fr is not None:
                ax3d.scatter(
                    [fr[0]], [fr[1]], [fr[2]],
                    marker="^", s=22 + 40 * r.amp,
                    color="#c8a028", alpha=min(0.85, 0.3 + r.amp * 0.6),
                    depthshade=False,
                )
            anti = _stereo(antipode4(r.p_front))
            if anti is not None:
                ax3d.scatter(
                    [anti[0]], [anti[1]], [anti[2]],
                    marker="v", s=22 + 40 * r.amp,
                    color="#ffd040",
                    alpha=min(0.85, 0.25 + r.amp * 0.6),
                    depthshade=False,
                )

        # TX bursts (recent)
        for ab in sim.absorptions:
            age = sim.t - ab.t
            if age > 1.1:
                continue
            pt = _stereo(ab.p4)
            if pt is None:
                continue
            alive = 1.0 - age / 1.1
            ax3d.scatter(
                [pt[0]], [pt[1]], [pt[2]],
                marker="*", s=140 + 200 * alive,
                color="#ffe037", alpha=alive,
                edgecolors="#ff8020", linewidths=0.8,
                depthshade=False,
            )

        # ── χ panel ─────────────────────────────────────────────────
        ax_chi.cla()
        ax_chi.set_facecolor("#02020a")
        ax_chi.set_xlim(0, 1)
        ax_chi.set_ylim(0, np.pi)
        ax_chi.set_ylabel("χᵢⱼ", color="#6a8aad", fontsize=8)
        ax_chi.tick_params(colors="#6a8aad", labelsize=7)
        ax_chi.axhline(EPS_CROSS, color="#ffff70", alpha=0.55, linestyle="--", linewidth=0.8)
        pair_cols = ("#38d8ff", "#70ff38", "#ff8040", "#e050ff")
        pairs = ((0, 1), (2, 3), (0, 2), (1, 3))
        for (ai, bi), col in zip(pairs, pair_cols):
            if ai >= len(sim.particles) or bi >= len(sim.particles):
                continue
            A = sim.particles[ai]
            B = sim.particles[bi]
            if len(B.trail) < 2:
                continue
            xi_bar_now = antipode4(A.p4)
            ys = [geo4(bp4, xi_bar_now) for _t, bp4 in B.trail]
            xs = np.linspace(0, 1, len(ys))
            ax_chi.plot(xs, ys, color=col, linewidth=0.9, alpha=0.75)
        ax_chi.set_title(
            "χᵢⱼ(t, t')  crossing angle",
            color="#6a8aad", fontsize=8, loc="left", pad=4,
        )

        # ── Energy panel ────────────────────────────────────────────
        ax_e.cla()
        ax_e.set_facecolor("#02020a")
        ax_e.tick_params(colors="#6a8aad", labelsize=7)
        if sim.energy_history:
            arr = np.asarray(sim.energy_history)
            ts = arr[:, 0]
            ecav = arr[:, 1]
            ekin = arr[:, 2]
            ax_e.fill_between(ts, 0, ecav, color="#c8a020", alpha=0.35, linewidth=0)
            ax_e.plot(ts, ecav, color="#c8a020", linewidth=1.0, label="E_cav")
            ax_e.plot(ts, ekin, color="#22c0d0", linewidth=1.0, label="E_kin")
            ax_e.plot(ts, ecav + ekin, color="#aaaaff",
                      linewidth=0.7, linestyle="--", alpha=0.6, label="E_tot")
            for ab in sim.absorptions:
                if ts[0] <= ab.t <= ts[-1]:
                    ax_e.axvline(ab.t, color="#ffe037", alpha=0.45, linewidth=0.8)
            ax_e.legend(
                loc="upper left", fontsize=6, frameon=False,
                labelcolor="#6a8aad", ncol=3,
            )
            ax_e.set_xlim(ts[0], ts[-1] + 1e-6)
        ax_e.set_title(
            "E_cav ↔ E_kin   │ = TX",
            color="#6a8aad", fontsize=8, loc="left", pad=4,
        )

        # ── Log panel ───────────────────────────────────────────────
        ax_log.cla()
        ax_log.set_facecolor("#02020a")
        ax_log.axis("off")
        lines = [
            f"t = {sim.t:6.2f}   TX = {sim.tx_count():>3d}   "
            f"χ_min = {sim.best_chi_now():.3f}"
        ]
        for ab in sim.absorptions[:6]:
            lines.append(
                f"t={ab.t:5.2f}  P{ab.src_id}→P{ab.dst_id}  "
                f"b={ab.bounce}  χ={ab.chi:.3f}  W={ab.weight:.3f}"
            )
        ax_log.text(
            0.01, 0.98, "\n".join(lines),
            color="#cdb015", fontsize=8, family="monospace",
            va="top", ha="left",
            transform=ax_log.transAxes,
        )

    def _update(_frame: int):
        sim.step()
        _draw()
        return ()

    _draw()
    ani = FuncAnimation(
        fig, _update, frames=n_steps, interval=interval_ms,
        blit=False, repeat=False, cache_frame_data=False,
    )
    if show:
        import matplotlib.pyplot as _plt
        _plt.show()
    return ani, sim


__all__ = [
    "AntipodalCrossingSim",
    "Crossing",
    "Absorption",
    "run_animation",
    "w_cross",
    "w_phase",
    "AMP_EMIT",
    "AMP_MIN",
    "EPS_CROSS",
    "EPS_PHASE",
    "MATCH_TOL",
    "GAMMA",
    "C",
    "EMIT_CD",
    "TX_CD",
    "N_BOUNCE",
    "TRAIL_LEN",
    "TRAIL_DT",
    "DEFAULT_DT",
]
