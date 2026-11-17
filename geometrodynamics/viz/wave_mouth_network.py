"""A wave that moves mouths, tips a loaded throat and is booked in a ledger.

A comparator for the visualization of PR #340, whose drawing is moved by a
solved rotor model but propagates no wave along the handle, has no collapse
and keeps no ledger. Everything drawn here is integrated from ONE Hamiltonian:

    H = sum_i dx [ p_i^2/2 + (psi_{i+1}-psi_i)^2/(2 dx^2) + m^2 psi_i^2/2 ]
        + sum_k [ sqrt(M^2 + P_k^2) - M            relativistic mouth
                  + pq_k^2/(2 mu) + S(q_k)         throat coordinate, support 
                  + kappa/2 (q_k - q0_k - <psi>_k)^2 ]   reciprocal coupling

    <psi>_k = sum_i dx g(x_i - X_k) psi_i,    g a normalised Gaussian.

``psi`` is a real wave on a ring, ``X_k`` a mouth position, ``q_k`` a signed
reaction coordinate of its throat (``1`` open, ``0`` collapsed; it overshoots
below zero, so it is not a nonnegative opening). The coupling has the mismatch form of
``waves/moving_mouth_pairs.py``. ``S`` is a support potential with a
metastable open state at ``q = 1`` holding energy ``E_s`` above the collapsed
state ``q = 0``.

What is solved: the wave's force on each mouth, the trigger of a collapse,
the field energy the transition releases, and every term of the energy ledger.
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


def barrier_top(model: Network) -> float:
    """Location of the support barrier between the open and collapsed wells."""
    grid = np.linspace(0.0, 1.0, 2001)
    h = model.barrier
    well = 16 * h * grid ** 2 * (1 - grid) ** 2 + grid ** 2 * (3 - 2 * grid)
    inner = well[1:-1]
    peak = (inner > well[:-2]) & (inner > well[2:])
    return float(grid[1:-1][peak][0]) if peak.any() else 1.0


def transitions(hist: History) -> list:
    """Time each LOADED throat enters its collapsed well and stays there.

    A transition needs a collapsed well to exist (stored energy > 0) and the
    coordinate to remain below the barrier top until the end of the run. A
    stable throat swinging below a threshold is not one; ``summary`` reports
    those separately as ``threshold_crossings``.
    """
    m, q = hist.model, hist.part("q")
    top = barrier_top(m)
    out = []
    for k in range(m.mouths):
        below = q[:, k] < top
        if m.stored[k] <= 0 or not below[-1]:
            out.append(None)
            continue
        above = np.flatnonzero(~below)
        start = above[-1] + 1 if above.size else 0
        out.append(float(hist.time[start]))
    return out


def far_field_energy(hist: History, index: int, radius: float | None = None
                     ) -> float:
    """Scalar-field energy at sites farther than ``radius`` from every mouth.

    Not an outgoing flux: on a ring it includes waves that have gone round.
    It separates field energy away from the mouths from energy held in the
    field around them.
    """
    m = hist.model
    radius = 5 * m.width if radius is None else radius
    psi, p = hist.part("psi")[index], hist.part("p")[index]
    pos = hist.part("X")[index]
    d = np.abs((m.x[None, :] - pos[:, None] + m.length / 2) % m.length
               - m.length / 2)
    mask = d.min(axis=0) > radius
    grad = (np.roll(psi, -1) - psi) / m.dx
    density = 0.5 * p ** 2 + 0.5 * grad ** 2 + 0.5 * m.field_mass ** 2 * psi ** 2
    return float(m.dx * np.sum(density[mask]))


def summary(hist: History) -> dict:
    """Numbers every caption is allowed to quote."""
    m = hist.model
    total = hist.total
    q, pos, mom = hist.part("q"), hist.part("X"), hist.part("P")
    field_energy = hist.ledger["wave"][:, 0]
    crossed = q < 0.5
    crossings = [float(hist.time[c.argmax()]) if c.any() else None
                 for c in crossed.T]
    momentum = np.array([m.momentum(z) for z in hist.states])
    return {
        "energy_initial": float(total[0]),
        "energy_max_relative_error": float(np.abs(total - total[0]).max()
                                           / abs(total[0])),
        "momentum_initial": float(momentum[0]),
        "momentum_max_drift": float(np.abs(momentum - momentum[0]).max()),
        # Total scalar-field energy over the whole ring, mouths included.
        "field_energy_initial": float(field_energy[0]),
        "field_energy_final": float(field_energy[-1]),
        "field_energy_max": float(field_energy.max()),
        "field_energy_gain": float(field_energy[-1] / field_energy[0])
        if field_energy[0] > 0 else None,
        "far_field_energy_initial": far_field_energy(hist, 0),
        "far_field_energy_final": far_field_energy(hist, -1),
        "barrier_top": barrier_top(m),
        "transition_times": transitions(hist),
        "threshold_crossings": crossings,
        "min_coordinate": q.min(axis=0).tolist(),
        "final_coordinates": q[-1].tolist(),
        "displacements": (pos[-1] - pos[0]).tolist(),
        "final_velocities": m.velocity(mom[-1]).tolist(),
        "max_speed": float(np.abs(m.velocity(mom)).max()),
        "stored_released": float(np.sum(m.support(q[0]) - m.support(q[-1]))),
    }


def validate(duration=DURATION, steps=(STEP, STEP / 2)) -> dict:
    """The documented run against an independent integrator, full window."""
    from scipy.integrate import solve_ivp

    m = Network()
    z0 = m.start(amplitude=AMPLITUDE)
    ref = solve_ivp(lambda t, z: m.rhs(z), (0.0, duration), z0,
                    method="DOP853", rtol=1e-10, atol=1e-12).y[:, -1]
    errors = []
    for h in steps:
        z = evolve(m, z0, duration, h=h, sample=int(round(duration / h))
                   ).states[-1]
        errors.append(float(np.linalg.norm(z - ref) / np.linalg.norm(ref)))
    e0 = m.energy(z0)
    return {"duration": duration, "steps": list(steps),
            "relative_final_state_error": errors,
            "improvement": errors[0] / errors[1],
            "dop853_energy_error": float(abs(m.energy(ref) - e0) / e0)}


BOX = (20.0, 40.0)          # fixed boundaries enclosing mouth 1 all run long


def bond_flux(model: Network, z, bonds) -> np.ndarray:
    """Rightward field-energy flux through the bonds (i-1, i), exact on the lattice.

    With site density ``e_i = p_i^2/2 + (psi_{i+1}-psi_i)^2/(2 dx^2) + m^2
    psi_i^2/2`` the free field obeys ``d(dx e_i)/dt = J_i - J_{i+1}`` with
    ``J_i = -p_i (psi_i - psi_{i-1}) / dx``. Sources act only where a mouth
    kernel is nonzero.
    """
    psi, p = model.unpack(z)[:2]
    bonds = np.asarray(bonds)
    return -p[bonds] * (psi[bonds] - psi[bonds - 1]) / model.dx


def box_energy(model: Network, z, bonds, inside=(0,)) -> tuple:
    """(field energy on sites a..b-1, energy of the mouths inside the box)."""
    psi, p = model.unpack(z)[:2]
    a, b = bonds
    grad = (np.roll(psi, -1) - psi) / model.dx
    density = 0.5 * p ** 2 + 0.5 * grad ** 2 + 0.5 * model.field_mass ** 2 * psi ** 2
    c = model.components(z)
    local = sum(float(c[n][k]) for n in ("mouth_kinetic", "throat", "coupling")
                for k in inside)
    return float(model.dx * density[a:b].sum()), local


def flux_ledger(amplitude=AMPLITUDE, stored=(4.0, 0.0, 0.0), box=BOX,
                duration=DURATION, h=STEP, sample=20, model=None) -> dict:
    """Time-integrated energy flux through two fixed boundaries around mouth 1.

    The flux is sampled at EVERY step and integrated by the trapezoid rule.
    ``net_outward`` is what left the box minus what entered; ``gross_*``
    split each boundary's instantaneous flux by sign. The balance

        net_outward(t) = E_inside(0) - E_inside(t)

    with ``E_inside`` the box field energy plus mouth 1's kinetic, throat
    and coupling energy is the check that the flux is the right one.
    """
    m = model or Network(stored=tuple(stored))
    a, b = (int(round(v / m.dx)) for v in box)
    bonds = np.array([a, b])
    z = m.start(amplitude=amplitude)
    steps = int(round(duration / h))
    out_sign = np.array([-1.0, 1.0])               # outward normal at a, b
    j = out_sign * bond_flux(m, z, bonds)
    cum = np.zeros(2)
    gross_out, gross_in = np.zeros(2), np.zeros(2)
    field0, local0 = box_energy(m, z, (a, b))
    rows = [(0.0, 0.0, 0.0, field0, local0, float(m.unpack(z)[2][0]))]
    for n in range(steps):
        z = m.step(z, h)
        j_new = out_sign * bond_flux(m, z, bonds)
        inc = 0.5 * h * (j + j_new)
        cum += inc
        gross_out += np.clip(inc, 0, None)
        gross_in += np.clip(-inc, 0, None)
        j = j_new
        if (n + 1) % sample == 0:
            f, loc = box_energy(m, z, (a, b))
            rows.append(((n + 1) * h, cum[0], cum[1], f, loc,
                         float(m.unpack(z)[2][0])))
    r = np.array(rows)
    inside = r[:, 3] + r[:, 4]
    net = r[:, 1] + r[:, 2]
    residual = net - (inside[0] - inside)
    return {
        "box": list(box), "step": h, "duration": duration,
        "time": r[:, 0], "outward_left": r[:, 1], "outward_right": r[:, 2],
        "net_outward": net, "box_field_energy": r[:, 3],
        "mouth_energy_inside": r[:, 4], "mouth_position": r[:, 5],
        "net_outward_final": float(net[-1]),
        "net_outward_left_final": float(cum[0]),
        "net_outward_right_final": float(cum[1]),
        "gross_outward": gross_out.tolist(), "gross_inward": gross_in.tolist(),
        "inside_initial": float(inside[0]), "inside_final": float(inside[-1]),
        "box_field_energy_final": float(r[-1, 3]),
        "balance_max_residual": float(np.abs(residual).max()),
        "mouth_margin": float(min((r[:, 5] - box[0]).min(),
                                  (box[1] - r[:, 5]).min())),
    }


def render(hist: History, path, title=None):
    """Three rows: spacetime diagram, throat coordinates, energy ledger."""
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
        if s["transition_times"][k] is not None:
            tc = s["transition_times"][k]
            ax.plot(np.interp(tc, t, wrapped), tc, marker="*", ms=16,
                    color="k", zorder=5)
    ax.set_xlabel("position on the ring")
    ax.set_ylabel("time")
    ax.set_title("Solved wave (colour) and mouth worldlines; star = loaded "
                 "throat enters its collapsed well", fontsize=11)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)

    ax = axes[1]
    for k in range(m.mouths):
        ax.plot(t, q[:, k], color=colours[k % 5], lw=1.6,
                label=f"throat {k + 1} (stored {m.stored[k]:g})")
    ax.axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax.axhline(0.0, color="0.6", lw=0.8, ls=":")
    ax.axhline(s["barrier_top"], color="0.4", lw=0.8, ls="--")
    ax.set_ylabel("throat coordinate q (signed)")
    ax.set_xlabel("time")
    ax.legend(loc="center right", fontsize=8)

    ax = axes[2]
    wave = hist.ledger["wave"][:, 0]
    kinetic = hist.ledger["mouth_kinetic"].sum(axis=1)
    throat = hist.ledger["throat"].sum(axis=1)
    coupling = hist.ledger["coupling"].sum(axis=1)
    ax.stackplot(t, wave, throat, coupling, kinetic,
                 labels=["scalar field (whole ring)", "throat support and motion", "coupling",
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

    fig.suptitle(title or "One Hamiltonian: a wave moves mouths and releases a "
                 "throat's stored energy into the field; the ledger closes",
                 fontsize=12, y=0.995)
    fig.text(0.5, 0.005, "Comparator model on a ring. The metastable support "
             "is assumed, not derived; this is not a GR solution.",
             ha="center", fontsize=8, color="0.3")
    fig.tight_layout(rect=(0, 0.015, 1, 0.985))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    return s


def render_flux(flux: dict, path, transition=None):
    """Cumulative outward energy flux through the two fixed boundaries."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = flux["time"]
    ink, muted = "#0b0b0b", "#52514e"
    series = [("outward_left", "through x = %g" % flux["box"][0], "#2a78d6", 2.0),
              ("outward_right", "through x = %g" % flux["box"][1], "#eb6834", 2.0),
              ("net_outward", "net, both boundaries", ink, 2.4)]
    fig, ax = plt.subplots(figsize=(9.5, 5.2), facecolor="#fcfcfb")
    ax.set_facecolor("#fcfcfb")
    ax.axhline(0.0, color=muted, lw=0.8)
    for key, label, colour, lw in series:
        ax.plot(t, flux[key], color=colour, lw=lw, label=label)
        ax.annotate(f"{flux[key][-1]:+.2f}", (t[-1], flux[key][-1]),
                    xytext=(6, 0), textcoords="offset points", va="center",
                    fontsize=9, color=ink)
    if transition is not None:
        ax.axvline(transition, color=muted, lw=0.8, ls="--")
        ax.annotate("throat 1 enters\nits collapsed well",
                    (transition, 0.82 * flux["net_outward"].max()),
                    xytext=(6, 0), textcoords="offset points", fontsize=9,
                    color=muted, va="center")
    ax.set_xlabel("time")
    ax.set_ylabel("cumulative outward energy flux")
    ax.set_xlim(0, t[-1] * 1.07)
    ax.grid(axis="y", color="#e6e5e0", lw=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    ax.set_title("Energy leaving the box around mouth 1, integrated from the "
                 "lattice flux", fontsize=11, loc="left", color=ink)
    fig.text(0.01, 0.01, "Negative = energy entering (the trigger packet). "
             "Balance against the energy inside the box closes to "
             f"{flux['balance_max_residual']:.1e}. Comparator model; ring, so "
             "late inflow includes reflections.", fontsize=8, color=muted)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def scenario(amplitude=AMPLITUDE, stored=(4.0, 0.0, 0.0), duration=DURATION,
             h=STEP, sample=20) -> History:
    """The documented run: one packet, three mouths, one loaded throat."""
    model = Network(stored=stored)
    return evolve(model, model.start(amplitude=amplitude), duration, h, sample)
