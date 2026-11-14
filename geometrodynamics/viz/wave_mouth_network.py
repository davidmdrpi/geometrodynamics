"""A wave that moves mouths, collapses a throat and is booked in a ledger.

A comparator for the visualization of PR #340, whose drawing is moved by a
solved rotor model but propagates no wave along the handle, has no collapse
and keeps no ledger. Everything drawn here is integrated from ONE Hamiltonian:

    H = sum_i dx [ p_i^2/2 + (psi_{i+1}-psi_i)^2/(2 dx^2) + m^2 psi_i^2/2 ]
        + sum_k [ sqrt(M^2 + P_k^2) - M            relativistic mouth
                  + pq_k^2/(2 mu) + S(q_k)         throat opening and support
                  + kappa/2 (q_k - q0_k - <psi>_k)^2 ]   reciprocal coupling

    <psi>_k = sum_i dx g(x_i - X_k) psi_i,    g a normalised Gaussian.

``psi`` is a real wave on a ring, ``X_k`` a mouth position, ``q_k`` its
throat opening. The coupling has the mismatch form of
``waves/moving_mouth_pairs.py``. ``S`` is a support potential with a
metastable open state at ``q = 1`` holding energy ``E_s`` above the collapsed
state ``q = 0``.

What is solved: the wave's force on each mouth, the trigger of a collapse,
the wave the collapse emits, and every term of the energy ledger.
What is assumed: this Hamiltonian. The metastable support in particular is
put in by hand; no GR support dynamics, no Einstein solution, no S^3.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

COMPONENTS = ("wave", "mouth_kinetic", "throat", "coupling")
POSITIONS = (25.0, 45.0, 62.0)
AMPLITUDE = 0.35            # the documented scenario
SUBTHRESHOLD = 0.2          # control: pushes the mouths, no collapse
DURATION = 70.0
STEP = 0.005


@dataclass(frozen=True)
class Network:
    sites: int = 800
    length: float = 80.0
    field_mass: float = 0.2
    mouth_mass: float = 6.0
    shape_mass: float = 8.0
    kappa: float = 3.0
    width: float = 0.4
    stored: tuple = (4.0, 0.0, 0.0)          # E_s per mouth; 0 = stable throat
    barrier: float = 0.3                     # h in S; must exceed 3/16
    stiffness: float = 4.0                   # stable throats: K/2 (q-1)^2

    def __post_init__(self):
        if self.barrier <= 3 / 16:
            raise ValueError("barrier must exceed 3/16 for a metastable state")
        if min(self.mouth_mass, self.shape_mass, self.kappa, self.width) <= 0:
            raise ValueError("positive masses, coupling and width required")

    @property
    def dx(self) -> float:
        return self.length / self.sites

    @property
    def x(self) -> np.ndarray:
        return self.dx * np.arange(self.sites)

    @property
    def mouths(self) -> int:
        return len(self.stored)

    @property
    def dim(self) -> int:
        return 2 * self.sites + 4 * self.mouths

    # ---- state -----------------------------------------------------------

    def unpack(self, z):
        n, k = self.sites, self.mouths
        return (z[:n], z[n:2 * n], z[2 * n:2 * n + k],
                z[2 * n + k:2 * n + 2 * k], z[2 * n + 2 * k:2 * n + 3 * k],
                z[2 * n + 3 * k:])

    def pack(self, psi, p, pos, mom, q, pq):
        return np.concatenate([psi, p, pos, mom, q, pq])

    def start(self, positions=POSITIONS, amplitude=0.0, centre=10.0,
              wavenumber=1.3, spread=3.0, direction=1.0):
        """Mouths at rest, throats open, and a travelling packet."""
        x = self.x
        d = (x - centre + self.length / 2) % self.length - self.length / 2
        envelope = amplitude * np.exp(-d ** 2 / (2 * spread ** 2))
        omega = np.sqrt(wavenumber ** 2 + self.field_mass ** 2)
        psi = envelope * np.cos(wavenumber * d)
        p = direction * omega * envelope * np.sin(wavenumber * d)
        k = self.mouths
        return self.pack(psi, p, np.asarray(positions, float), np.zeros(k),
                         np.ones(k), np.zeros(k))

    # ---- pieces of H -----------------------------------------------------

    def kernel(self, pos):
        d = (self.x[None, :] - pos[:, None] + self.length / 2) % self.length \
            - self.length / 2
        g = np.exp(-d ** 2 / (2 * self.width ** 2)) / (
            self.width * np.sqrt(2 * np.pi))
        return g, d

    def support(self, q):
        """Metastable double well where energy is stored; a stable well otherwise."""
        e, h = np.asarray(self.stored), self.barrier
        well = e * (16 * h * q ** 2 * (1 - q) ** 2 + q ** 2 * (3 - 2 * q))
        return np.where(e > 0, well, 0.5 * self.stiffness * (q - 1) ** 2)

    def support_force(self, q):
        e, h = np.asarray(self.stored), self.barrier
        well = -e * (32 * h * q * (1 - q) * (1 - 2 * q) + 6 * q * (1 - q))
        return np.where(e > 0, well, -self.stiffness * (q - 1))

    def mismatch(self, psi, pos, q):
        g, d = self.kernel(pos)
        return q - 1.0 - self.dx * g @ psi, g, d

    def components(self, z) -> dict:
        psi, p, pos, mom, q, pq = self.unpack(z)
        grad = (np.roll(psi, -1) - psi) / self.dx
        mis, _, _ = self.mismatch(psi, pos, q)
        return {
            "wave": float(self.dx * np.sum(
                0.5 * p ** 2 + 0.5 * grad ** 2
                + 0.5 * self.field_mass ** 2 * psi ** 2)),
            "mouth_kinetic": np.sqrt(self.mouth_mass ** 2 + mom ** 2)
            - self.mouth_mass,
            "throat": pq ** 2 / (2 * self.shape_mass) + self.support(q),
            "coupling": 0.5 * self.kappa * mis ** 2,
        }

    def energy(self, z) -> float:
        return float(sum(np.sum(v) for v in self.components(z).values()))

    def momentum(self, z) -> float:
        """Continuum Noether momentum: wave plus mouths. Not exact on a lattice."""
        psi, p, _, mom, _, _ = self.unpack(z)
        grad = (np.roll(psi, -1) - np.roll(psi, 1)) / (2 * self.dx)
        return float(-self.dx * np.sum(p * grad) + np.sum(mom))

    def forces(self, z):
        """Minus the gradient of the potential part of H."""
        psi, _, pos, _, q, _ = self.unpack(z)
        mis, g, d = self.mismatch(psi, pos, q)
        lap = (np.roll(psi, -1) - 2 * psi + np.roll(psi, 1)) / self.dx ** 2
        # H is a sum over sites of dx * density, so dp/dt = -(1/dx) dH/dpsi.
        f_psi = lap - self.field_mass ** 2 * psi + self.kappa * (mis @ g)
        dg = g * d / self.width ** 2                 # d g / d X
        f_pos = self.kappa * mis * (self.dx * dg @ psi)
        f_q = self.support_force(q) - self.kappa * mis
        return f_psi, f_pos, f_q

    def velocity(self, mom):
        return mom / np.sqrt(self.mouth_mass ** 2 + mom ** 2)

    def rhs(self, z):
        psi, p, pos, mom, q, pq = self.unpack(z)
        f_psi, f_pos, f_q = self.forces(z)
        return self.pack(p, f_psi, self.velocity(mom), f_pos,
                         pq / self.shape_mass, f_q)

    # ---- integrator ------------------------------------------------------

    def step(self, z, h):
        """Velocity Verlet; H is kinetic(p, P, p_q) plus potential(psi, X, q)."""
        psi, p, pos, mom, q, pq = (v.copy() for v in self.unpack(z))
        f_psi, f_pos, f_q = self.forces(z)
        p += 0.5 * h * f_psi
        mom += 0.5 * h * f_pos
        pq += 0.5 * h * f_q
        psi += h * p
        pos += h * self.velocity(mom)
        q += h * pq / self.shape_mass
        f_psi, f_pos, f_q = self.forces(self.pack(psi, p, pos, mom, q, pq))
        p += 0.5 * h * f_psi
        mom += 0.5 * h * f_pos
        pq += 0.5 * h * f_q
        return self.pack(psi, p, pos, mom, q, pq)


@dataclass
class History:
    model: Network
    time: np.ndarray
    states: np.ndarray
    ledger: dict = field(default_factory=dict)

    @property
    def total(self) -> np.ndarray:
        return sum(np.sum(np.atleast_2d(v.T).T, axis=1)
                   for v in self.ledger.values())

    def part(self, name: str) -> np.ndarray:
        n, k = self.model.sites, self.model.mouths
        offsets = {"psi": (0, n), "p": (n, 2 * n),
                   "X": (2 * n, 2 * n + k), "P": (2 * n + k, 2 * n + 2 * k),
                   "q": (2 * n + 2 * k, 2 * n + 3 * k),
                   "pq": (2 * n + 3 * k, 2 * n + 4 * k)}
        lo, hi = offsets[name]
        return self.states[:, lo:hi]


def evolve(model: Network, z0, duration=DURATION, h=STEP, sample=20) -> History:
    steps = int(round(duration / h))
    if steps % sample:
        raise ValueError("duration must hold a whole number of samples")
    z = np.asarray(z0, float).copy()
    times, states, rows = [0.0], [z.copy()], [model.components(z)]
    for j in range(steps):
        z = model.step(z, h)
        if (j + 1) % sample == 0:
            times.append((j + 1) * h)
            states.append(z.copy())
            rows.append(model.components(z))
    ledger = {name: np.array([np.atleast_1d(r[name]) for r in rows])
              for name in COMPONENTS}
    return History(model, np.array(times), np.array(states), ledger)


def summary(hist: History) -> dict:
    """Numbers every caption is allowed to quote."""
    m = hist.model
    total = hist.total
    q, pos, mom = hist.part("q"), hist.part("X"), hist.part("P")
    wave = hist.ledger["wave"][:, 0]
    collapsed = q < 0.5
    events = [float(hist.time[c.argmax()]) if c.any() else None
              for c in collapsed.T]
    momentum = np.array([m.momentum(z) for z in hist.states])
    return {
        "energy_initial": float(total[0]),
        "energy_max_relative_error": float(np.abs(total - total[0]).max()
                                           / abs(total[0])),
        "momentum_initial": float(momentum[0]),
        "momentum_max_drift": float(np.abs(momentum - momentum[0]).max()),
        "wave_energy_initial": float(wave[0]),
        "wave_energy_final": float(wave[-1]),
        "wave_energy_max": float(wave.max()),
        "wave_gain": float(wave[-1] / wave[0]) if wave[0] > 0 else None,
        "collapse_times": events,
        "final_openings": q[-1].tolist(),
        "displacements": (pos[-1] - pos[0]).tolist(),
        "final_velocities": m.velocity(mom[-1]).tolist(),
        "max_speed": float(np.abs(m.velocity(mom)).max()),
        "stored_released": float(np.sum(m.support(q[0]) - m.support(q[-1]))),
    }


def render(hist: History, path, title=None):
    """Three rows: spacetime diagram, throat openings, energy ledger."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m, t = hist.model, hist.time
    psi, pos, q = hist.part("psi"), hist.part("X"), hist.part("q")
    s = summary(hist)
    colours = ["#b3261e", "#1f6f8b", "#6a4c93", "#2e7d32", "#8d6e63"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=False,
                             gridspec_kw={"height_ratios": [3.2, 1.2, 1.8]})

    ax = axes[0]
    scale = np.percentile(np.abs(psi), 99.5) or 1.0
    ax.imshow(psi, origin="lower", aspect="auto", cmap="RdBu_r",
              vmin=-scale, vmax=scale,
              extent=[0, m.length, t[0], t[-1]], interpolation="nearest")
    for k in range(m.mouths):
        wrapped = pos[:, k] % m.length
        ax.plot(wrapped, t, color=colours[k % 5], lw=2.2,
                label=f"mouth {k + 1} worldline")
        if s["collapse_times"][k] is not None:
            tc = s["collapse_times"][k]
            ax.plot(np.interp(tc, t, wrapped), tc, marker="*", ms=16,
                    color="k", zorder=5)
    ax.set_xlabel("position on the ring")
    ax.set_ylabel("time")
    ax.set_title("Solved wave (colour) and mouth worldlines; star = throat "
                 "collapse", fontsize=11)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    ax = axes[1]
    for k in range(m.mouths):
        ax.plot(t, q[:, k], color=colours[k % 5], lw=1.6,
                label=f"throat {k + 1} (stored {m.stored[k]:g})")
    ax.axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax.axhline(0.0, color="0.6", lw=0.8, ls=":")
    ax.set_ylabel("throat opening q")
    ax.set_xlabel("time")
    ax.legend(loc="center right", fontsize=8)

    ax = axes[2]
    wave = hist.ledger["wave"][:, 0]
    kinetic = hist.ledger["mouth_kinetic"].sum(axis=1)
    throat = hist.ledger["throat"].sum(axis=1)
    coupling = hist.ledger["coupling"].sum(axis=1)
    ax.stackplot(t, wave, throat, coupling, kinetic,
                 labels=["wave", "throat support and motion", "coupling",
                         "mouth kinetic"],
                 colors=["#4f83cc", "#c98f2b", "#9e9e9e", "#b3261e"],
                 alpha=0.9)
    ax.plot(t, hist.total, color="k", lw=1.2, label="total")
    ax.set_ylabel("energy")
    ax.set_xlabel("time")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=8,
              ncol=5, frameon=False)
    ax.set_title(f"Ledger: total conserved to "
                 f"{s['energy_max_relative_error']:.1e} relative", fontsize=11)

    fig.suptitle(title or "One Hamiltonian: a wave moves mouths, a throat "
                 "collapses into a larger wave, the ledger closes",
                 fontsize=12, y=0.995)
    fig.text(0.5, 0.005, "Comparator model on a ring. The metastable support "
             "is assumed, not derived; this is not a GR solution.",
             ha="center", fontsize=8, color="0.3")
    fig.tight_layout(rect=(0, 0.015, 1, 0.985))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return s


def scenario(amplitude=AMPLITUDE, stored=(4.0, 0.0, 0.0), duration=DURATION,
             h=STEP, sample=20) -> History:
    """The documented run: one packet, three mouths, one loaded throat."""
    model = Network(stored=stored)
    return evolve(model, model.start(amplitude=amplitude), duration, h, sample)
