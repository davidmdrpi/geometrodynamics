"""Lossless, integer-delay ring. The string waves are stored at cell midpoints."""
import json

import numpy as np

from . import two_detector_records as base

LABELS = (*base.LABELS[:4], "changing_or_unsettled")
START, END = 12., 36.


def counts(codes):
    return dict(zip(LABELS, np.bincount(codes, minlength=5).tolist()))


def basin_pair(q):
    return (q[0] > base.SADDLE).astype(int) + 2*(q[1] > base.SADDLE)


def persistent(low, high):
    ready = high < .65
    written = low > .85
    return np.where((ready | written).all(axis=0),
                    written[0].astype(int)+2*written[1], 4)


def evolve(pair, q0, p0, delay=4., dt=.025, closed=True, mu=.25):
    """One incident energy pair, an ensemble of receiver preparations, and one trace."""
    steps, slots = round(END/dt), round(delay/dt)
    if not all(np.isclose(round(t/dt)*dt, t) for t in (START, END, delay, base.DURATION)):
        raise ValueError("time grid must contain the pulse, delay and observation boundaries")
    if delay < base.DURATION or min(pair) < 0:
        raise ValueError("initial pulses must fit on the two arms")
    ampl = np.sqrt(pair)
    q = q0.T.copy()
    v = base.inverse_mass(p0, mu).T.copy()
    waveform = np.array([base.incoming((j+.5)*dt) for j in range(slots)])
    initial = waveform[:, None, None]*ampl[None, :, None]
    # Each entry contains both identical arms' incoming wave for a given port.
    line = np.broadcast_to(initial, (slots, 2, len(q0))).copy()
    incident = base.GAMMA*dt*np.sum(initial**2)
    total0 = base.mechanical(q.T, v.T, mu)+incident
    escaped = np.zeros(len(q0))
    trace = np.empty((steps+1, 2))
    trace[0] = q[:, 0]
    low = np.full_like(q, np.inf)
    high = -low
    integral_net = np.zeros_like(q)
    sample_steps = {round(t/dt) for t in (base.DURATION, delay, 2*delay, 3*delay, START, 24., END)}
    snapshots = []
    max_residual = 0.
    max_iteration_error = 0.
    a = 1 + dt*base.GAMMA/2
    # (M + dt*gamma*I/2)^-1 acting on symmetric and antisymmetric parts.
    diag = (a+mu)/(a*(a+2*mu))
    cross = mu/(a*(a+2*mu))
    for j in range(steps):
        slot = j % slots
        f = line[slot]
        drive = v + mu*(v-v[::-1]) + dt*base.GAMMA*f
        mid_v = v
        for iteration in range(5):
            jump = dt*mid_v
            mid_q = q+jump/2
            gradient = base.slope(mid_q) + (6*mid_q-3.5)*jump**2
            rhs = drive - dt*gradient/2
            trial = diag*rhs + cross*rhs[::-1]
            if iteration == 4:
                max_iteration_error = max(max_iteration_error, np.max(abs(trial-mid_v)))
            mid_v = trial
        outgoing = mid_v-f
        integral_net += base.GAMMA*dt*(f*f-outgoing*outgoing)
        if closed:
            line[slot] = outgoing[::-1]
        else:
            escaped += base.GAMMA*dt*np.sum(outgoing**2, axis=0)
            line[slot] = 0.
        q += dt*mid_v
        v = 2*mid_v-v
        trace[j+1] = q[:, 0]
        if j+1 in sample_steps:
            snapshots.append(dict(t=(j+1)*dt, counts=counts(basin_pair(q))))
        if j+1 >= round(START/dt):
            np.minimum(low, q, out=low)
            np.maximum(high, q, out=high)
        if (j+1) % slots == 0 or j+1 == steps:
            stored = base.GAMMA*dt*np.einsum("jkn,jkn->n", line, line)
            balance = base.mechanical(q.T, v.T, mu)+stored+escaped-total0
            max_residual = max(max_residual, np.max(abs(balance))/base.PULSE_ENERGY)
    if not np.all(np.isfinite(q)) or not np.isfinite(max_residual):
        raise ArithmeticError("nonfinite integration")
    stable = persistent(low, high)
    snap = basin_pair(q)
    # A supplementary readout margin, using |qdot_i| <= sqrt(2E_initial)
    # and half a sample interval to the nearest measured endpoint.
    extra = dt*np.sqrt(2*total0)/2
    robust = persistent(low-extra, high+extra)
    both_single = (stable == 1) | (stable == 2)
    row = dict(delay=delay, closed=closed, mu=mu, dt=dt, launches=len(q0),
               incident_fractions=list(pair), sampled_incident_energy=float(incident),
               persistent_counts=counts(stable), endpoint_counts=counts(snap),
               buffered_counts=counts(robust), exclusive=float(both_single.mean()),
               right_only=float((stable == 1).mean()), max_energy_residual=float(max_residual),
               mean_net_work_at_ports=integral_net.mean(axis=1).tolist(),
               basin_snapshots=snapshots,
               max_iteration_change=float(max_iteration_error))
    return row, trace, stable


def run_angles(weights, count, seed, dt, indices=None, closed=True):
    q, p = base.preparation(count, seed, .02)
    rows = []
    for i in range(len(weights)) if indices is None else indices:
        row, _, _ = evolve([weights[i], 1-weights[i]], q, p, dt=dt, closed=closed)
        row.update(index=int(i), theta=float(base.THETAS[i]), born=float(np.cos(base.THETAS[i])**2))
        rows.append(row)
    return rows


def difference(a, b):
    return max(abs(a[k]-b[k]) for k in ("exclusive", "right_only"))


def gates(rows, open_balanced):
    balanced = next(r for r in rows if r["index"] == 4)
    error = max(abs(r["right_only"]-r["born"]) for r in rows)
    return dict(above_half_and_open=balanced["exclusive"] > max(.5, open_balanced),
                near_unit_all_angles=min(r["exclusive"] for r in rows) >= .99,
                born_window=error <= .01, max_born_error=error)


def probe():
    weights = [base.split_weights(t, .02)[0] for t in base.THETAS]
    primary = run_angles(weights, 512, 363, .025)
    fine = run_angles(weights, 512, 363, .0125)
    repeat = run_angles(weights, 1024, 364, .025)
    open_rows = run_angles(weights, 512, 363, .025, closed=False)
    repeat_open = run_angles(weights, 1024, 364, .025, indices=[4], closed=False)[0]
    differences = [difference(a, b) for a, b in zip(primary, fine)]
    offending = [i for i, d in enumerate(differences) if d > .02]
    finer = run_angles(weights, 512, 363, .00625, indices=offending)
    q, p = base.preparation(512, 363, .02)
    controls = []
    for delay in (3., 4., 6.):
        closed, track, _ = evolve([.5, .5], q, p, delay=delay)
        opened, reference, _ = evolve([.5, .5], q, p, delay=delay, closed=False)
        err = abs(track-reference).max(axis=1)
        hits = np.flatnonzero(err > 1e-8)
        slots = round(delay/.025)
        samples = np.arange(0, len(track), 10)
        controls.append(dict(closed=closed, open=opened, pre_echo_max=float(max(err[:slots+1])),
                             first_deviation=None if not len(hits) else float(hits[0]*.025),
                             trace_t=(samples*.025).tolist(), loop_trace=track[samples].tolist(),
                             open_trace=reference[samples].tolist()))
    cold = np.zeros((1, 2))
    balanced, track, _ = evolve([.5, .5], cold, cold)
    dark, _, _ = evolve([0., 0.], cold, cold)
    unlinked, _, _ = evolve([.5, .5], q, p, mu=0.)
    all_rows = (primary+fine+repeat+open_rows+finer+[repeat_open, dark, balanced, unlinked]
                +[row[k] for row in controls for k in ("closed", "open")])
    return dict(primary=primary, fine=fine, repeat=repeat, open=open_rows,
                repeat_open_balanced=repeat_open, extra_refinement=finer,
                length_controls=controls, cold_balanced=balanced, cold_dark=dark,
                balanced_no_mutual_inertia=unlinked,
                gates={name:gates(rows, op) for name, rows, op in
                       (("primary", primary, open_rows[4]["exclusive"]),
                        ("fine", fine, open_rows[4]["exclusive"]),
                        ("repeat", repeat, repeat_open["exclusive"]))},
                checks=dict(max_energy_residual=max(r["max_energy_residual"] for r in all_rows),
                            max_step_fraction_difference=max(differences),
                            large_difference_indices=offending,
                            pre_echo_pass=all(c["pre_echo_max"] < 1e-10 for c in controls),
                            cold_symmetry_max=float(max(abs(track[:, 0]-track[:, 1])))))


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2, allow_nan=False))
