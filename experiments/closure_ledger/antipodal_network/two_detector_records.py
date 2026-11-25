"""Two nonlinear string terminations driven by the two outputs of the elliptic mixer."""
import json

import numpy as np
from scipy.integrate import quad
from scipy.stats import qmc

from .elliptic_transport import THETAS, split_weights

GAMMA = 2.
DURATION = 2.5
READ_TIME = 24.
SADDLE = .75
LABELS = ("neither", "right_only", "left_only", "both", "uncertified")


def potential(q):
    return q*q*(6*q*q-14*q+9)


def slope(q):
    return 24*q*(q-SADDLE)*(q-1)


def incoming(t):
    if t <= 0 or t >= DURATION:
        return 0.
    x = t/DURATION
    q = x**3*(10-15*x+6*x*x)
    v = 30*x*x*(1-x)**2/DURATION
    a = 60*x*(1-x)*(1-2*x)/DURATION**2
    return (a+GAMMA*v+slope(q))/(2*GAMMA)


PULSE_ENERGY = quad(lambda t: GAMMA*incoming(t)**2, 0., DURATION)[0]


def preparation(count, seed, e):
    if e == 0:
        return np.zeros((1, 2)), np.zeros((1, 2))
    power = int(np.log2(count))
    if count != 2**power or count < 2:
        raise ValueError("use a power-of-two number of launches")
    ranks = qmc.Sobol(4, scramble=True, seed=seed).random_base2(power-1)
    radius = np.sqrt(2*e*ranks[:, ::2])
    phase = 2*np.pi*ranks[:, 1::2]
    q = radius*np.cos(phase)/np.sqrt(18.)
    p = radius*np.sin(phase)
    return np.vstack((q, q[:, ::-1])), np.vstack((p, p[:, ::-1]))


def inverse_mass(p, mu):
    return (p + mu*p.sum(axis=-1, keepdims=True))/(1+2*mu)


def mechanical(q, v, mu):
    return (potential(q)+v*v/2).sum(axis=-1) + mu*(v[..., 0]-v[..., 1])**2/2


def certify(q, energy):
    written = q > SADDLE
    # At a first escape the crossing port reaches B. The other port still
    # has at least the minimum potential of its current basin.
    escape = potential(SADDLE) + np.min(written, axis=-1)*potential(1.)
    mask = written[..., 0].astype(int)+2*written[..., 1]
    return np.where(energy < escape, mask, 4)


def evolve(fractions, q0, p0, mu, dt):
    """fractions has one (right, left) incident-energy pair per measurement input."""
    amplitudes = np.sqrt(np.asarray(fractions, dtype=float))[:, None, :]
    shape = (len(amplitudes), len(q0), 2)
    q = np.broadcast_to(q0, shape).copy()
    v = np.broadcast_to(inverse_mass(p0, mu), shape).copy()
    initial = mechanical(q, v, mu)
    net_work = np.zeros_like(initial)
    incident = np.zeros((len(amplitudes), 1))
    max_balance = 0.
    steps = round(READ_TIME/dt)
    if not np.isclose(steps*dt, READ_TIME) or not np.isclose(round(DURATION/dt)*dt, DURATION):
        raise ValueError("grid must include pulse cutoff and read time")

    def stage(q, v, f):
        a = inverse_mass(2*GAMMA*f-GAMMA*v-slope(q), mu)
        # Incoming minus outgoing power for the string solution at the endpoint.
        power = GAMMA*(f*f-(f-v)**2).sum(axis=-1)
        return a, power

    for j in range(steps):
        t = dt*j
        f0 = amplitudes*incoming(t)
        fm = amplitudes*incoming(t+dt/2)
        f1 = amplitudes*incoming(t+dt)
        a1, w1 = stage(q, v, f0)
        v2 = v+dt*a1/2
        a2, w2 = stage(q+dt*v/2, v2, fm)
        v3 = v+dt*a2/2
        a3, w3 = stage(q+dt*v2/2, v3, fm)
        v4 = v+dt*a3
        a4, w4 = stage(q+dt*v3, v4, f1)
        q += dt*(v+2*v2+2*v3+v4)/6
        v += dt*(a1+2*a2+2*a3+a4)/6
        net_work += dt*(w1+2*w2+2*w3+w4)/6
        incident += dt*GAMMA*(f0*f0+4*fm*fm+f1*f1).sum(axis=-1)/6
        energy = mechanical(q, v, mu)
        max_balance = max(max_balance, np.max(abs(energy-initial-net_work)))
    if not np.all(np.isfinite(q)) or not np.all(np.isfinite(energy)):
        raise FloatingPointError("nonfinite terminal state")
    labels = certify(q, energy)
    rows = []
    for i, codes in enumerate(labels):
        counts = np.bincount(codes, minlength=len(LABELS))
        rows.append(dict(counts=dict(zip(LABELS, counts.tolist())),
                         initial_energy_max=float(initial[i].max()),
                         final_energy_range=[float(energy[i].min()), float(energy[i].max())],
                         outgoing_energy_mean=float((incident[i]-net_work[i]).mean()),
                         mean_position=q[i].mean(axis=0).tolist()))
    return rows, labels, float(max_balance/PULSE_ENERGY)


def experiment(weights, mu, e, count=512, seed=363, dt=.025):
    q, p = preparation(count, seed, e)
    rows, labels, balance = evolve(np.column_stack((weights, 1-weights)), q, p, mu, dt)
    right = (labels == 1).mean(axis=1)
    for row, theta, weight, prob in zip(rows, THETAS, weights, right):
        row.update(theta=float(theta), right_incident_fraction=float(weight),
                   right_only_probability=float(prob), target=float(np.cos(theta)**2))
    born_error = np.max(abs(right-np.cos(THETAS)**2))
    report = dict(mu=mu, ready_energy_parameter=e, count_per_angle=len(q), seed=seed, dt=dt,
                  rows=rows, max_budget_residual=balance, max_born_error=float(born_error),
                  max_error_from_mixer=float(np.max(abs(right-weights))),
                  all_exclusive=bool(np.all((labels == 1) | (labels == 2))),
                  born_window=bool(born_error <= .01))
    return report, labels


def probe():
    weights = np.array([split_weights(theta, .02)[0] for theta in THETAS])
    runs = {}
    primary, coarse = experiment(weights, .25, .02)
    fine, refined = experiment(weights, .25, .02, dt=.0125)
    runs["primary"], runs["refined"] = primary, fine
    for name, mu, e in (("uncoupled", 0., .02), ("broader", .25, .06), ("cold", .25, 0.)):
        runs[name], _ = experiment(weights, mu, e)
    runs["second_ensemble"], _ = experiment(weights, .25, .02, count=1024, seed=364)
    q, p = preparation(1, 363, 0.)
    controls, _, residual = evolve([(0., 0.), (1., 0.), (1., 1.)], q, p, 0., .025)
    counts_equal = all(a["counts"] == b["counts"] for a, b in zip(primary["rows"], fine["rows"]))
    return dict(pulse_energy=PULSE_ENERGY, read_time=READ_TIME, barrier=potential(SADDLE),
                written_minimum=potential(1.), runs=runs,
                controls=dict(zip(("dark", "one_full", "two_full"), controls)),
                control_budget_residual=residual,
                checks=dict(refined_counts_equal=counts_equal,
                            individual_changes=int(np.count_nonzero(coarse != refined)),
                            budget_pass=max([r["max_budget_residual"] for r in runs.values()]
                                            + [residual]) < 1e-5,
                            calibration_pass=controls[1]["counts"]["right_only"] == 1))


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2, allow_nan=False))
