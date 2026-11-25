"""Two proper-time receivers on #364's ring with prescribed expansion."""
import json
import math
from pathlib import Path

import numpy as np

from . import closed_loop_records as ring

base = ring.base
D = 4.
RATES = (0., 1e-6, .1, .5, 1., 3.)


def scale(eta, rate):
    return np.sqrt(1 + 2*rate*np.maximum(eta-base.DURATION, 0.)/D)


def proper_time(eta, rate):
    # Stable integral of sqrt(1+c*x), including its zero-rate limit.
    x = np.maximum(eta-base.DURATION, 0.)
    a = scale(eta, rate)
    return np.minimum(eta, base.DURATION) + 2*x*(a*a+a+1)/(3*(a+1))


def advance(q, v, f, dt):
    """#364's discrete-gradient update, with local proper step dt."""
    mu = .25
    a = 1 + dt*base.GAMMA/2
    diag, cross = (a+mu)/(a*(a+2*mu)), mu/(a*(a+2*mu))
    drive = v + mu*(v-v[::-1]) + dt*base.GAMMA*f
    mid = v
    for _ in range(5):
        jump = dt*mid
        center = q+jump/2
        rhs = drive - dt*(base.slope(center)+(6*center-3.5)*jump**2)/2
        mid = diag*rhs + cross*rhs[::-1]
    return q+dt*mid, 2*mid-v, mid


def evolve(pair, rate, mode="redshift", limit=.025, count=512):
    if mode not in ("redshift", "timing", "open") or rate < 0 or limit <= 0:
        raise ValueError("invalid transport or grid")
    slots = 8*math.ceil(D*scale(ring.END, rate)/(8*limit))
    step = D/slots
    steps = round(ring.END/step)
    eta = np.arange(steps+1)*step
    time = proper_time(eta, rate)
    delta = np.diff(time)
    if delta.max() > limit*(1+1e-9):
        raise ArithmeticError("proper time resolution exceeded")
    a = scale(eta, rate)
    # Exact cell-average scale factor; matches the proper-time quadrature.
    amid = delta/step
    exponent = .5 if mode == "timing" else 1.
    wave_unit = a**(1-2*exponent)
    mid_unit = amid**(1-2*exponent)
    gain = amid**exponent
    q0, p0 = base.preparation(count, 363, .02)
    q, v = q0.T.copy(), base.inverse_mass(p0, .25).T.copy()
    pulse = np.array([base.incoming((j+.5)*step) for j in range(slots)])
    initial = pulse[:, None, None]*np.sqrt(pair)[None, :, None]
    line = np.broadcast_to(initial, (slots, 2, count)).copy()
    packet_sum = np.full(count, base.GAMMA*step*np.sum(initial**2))
    start = base.mechanical(q.T, v.T, .25)+packet_sum
    loss, escape = np.zeros(count), np.zeros(count)
    low, high = np.full_like(q, np.inf), np.full_like(q, -np.inf)
    trace = np.empty((2*slots+1, 2))
    trace[0] = q[:, 0]
    residual = 0.
    for j, dt in enumerate(delta):
        slot = j % slots
        stored = line[slot]
        f = stored/gain[j]
        q, v, mid = advance(q, v, f, dt)
        outgoing = mid-f
        loss += packet_sum*(wave_unit[j]-mid_unit[j])
        incoming_norm = np.sum(stored*stored, axis=0)
        if mode == "open":
            escape += base.GAMMA*dt*np.sum(outgoing*outgoing, axis=0)
            stored.fill(0.)
        else:
            stored[:] = gain[j]*outgoing[::-1]
        packet_sum += base.GAMMA*step*(np.sum(stored*stored, axis=0)-incoming_norm)
        loss += packet_sum*(mid_unit[j]-wave_unit[j+1])
        if j+1 <= 2*slots:
            trace[j+1] = q[:, 0]
        if j+1 == slots:
            first = ring.basin_pair(q)
            alive = (first == 1) | (first == 2)
        elif j+1 > slots:
            alive &= ring.basin_pair(q) == first
        if j+1 >= 3*slots:
            np.minimum(low, q, out=low)
            np.maximum(high, q, out=high)
        if (j+1) % slots == 0:
            inventory = base.GAMMA*step*np.einsum("jkn,jkn->n", line, line)
            defect = base.mechanical(q.T,v.T,.25)+wave_unit[j+1]*inventory+loss+escape-start
            residual = max(residual, float(np.max(np.abs(defect))/base.PULSE_ENERGY))
    stable = ring.persistent(low, high)
    same = alive & (stable == first)
    if not all(np.isfinite(x).all() for x in (q,v,loss,low,high)) or residual > 1e-8:
        raise ArithmeticError(f"failed energy/finite check: {residual}")
    return dict(rate=rate, mode=mode, limit=limit, slots=slots,
                first_singles=int(np.isin(first, [1,2]).sum()),
                unbroken_singles=int(alive.sum()), unbroken_and_durable=int(same.sum()),
                persistent=ring.counts(stable), final=ring.counts(ring.basin_pair(q)),
                lost_fraction_mean=float(np.mean(loss/start)), max_budget_residual=residual,
                first_echo_time=float(time[slots]), end_time=float(time[-1])), trace


def probe():
    angles = base.THETAS
    pairs = {str(i): [w, 1-w] for i,t in enumerate(angles)
             for w in [base.split_weights(t, .02)[0]]}
    pairs["bright"] = [.875,.875]
    rows = []
    pre_open = {}
    controls = []
    def run(key, rate, mode="redshift", limit=.025):
        row, trace = evolve(pairs[key], rate, mode, limit)
        row["input"] = key
        rows.append(row)
        if rate == 1. and limit == .025 and key in ("0", "4", "bright"):
            if mode == "redshift":
                pre_open[key] = trace[:row["slots"]+1]
            elif mode == "open":
                controls.append(float(np.max(np.abs(trace[:row["slots"]+1]-pre_open[key]))))
        return row
    for rate in RATES:
        keys = list(pairs) if rate in (0.,1.) else ["0","4","bright"]
        for key in keys:
            run(key, rate)
    for key in ("0","4","bright"):
        run(key, 1., "open")
        run(key, 1., "timing")
    coarse = {r["input"]: r for r in rows if r["rate"]==1 and r["mode"]=="redshift"}
    refined = {}
    def fractions(row):
        c = row["persistent"]
        return np.array([c["right_only"], c["right_only"]+c["left_only"],
                         row["unbroken_and_durable"]])/512
    flags = []
    for key in pairs:
        refined[key] = run(key, 1., limit=.0125)
        if max(abs(fractions(coarse[key])-fractions(refined[key]))) > .02:
            flags.append(key)
            run(key, 1., limit=.00625)
    archive = Path(ring.__file__).resolve().parents[1]/"runs/20261125_closed_loop_records/probe.json"
    previous = json.loads(archive.read_text())["primary"]
    for i,r in enumerate(previous):
        row = next(v for v in rows if v["rate"]==0 and v["input"]==str(i))
        assert row["persistent"] == r["persistent_counts"]
        assert row["final"] == r["endpoint_counts"]
    assert len(controls)==3 and max(controls) < 1e-10
    def gates(angular):
        stats = np.array([fractions(angular[str(i)]) for i in range(9)])
        error = max(abs(stats[:,0]-np.cos(angles)**2))
        opened = next(r for r in rows if r["input"]=="4" and r["mode"]=="open")
        return dict(balanced_above_half_and_open=bool(stats[4,1]>max(.5, fractions(opened)[1])),
                    near_unit_every_angle=bool(min(stats[:,1])>=.99),
                    born_window=bool(error<=.01), max_born_error=float(error))
    return dict(rows=rows, pre_echo_differences=controls, static_archive_match=True,
                refinement_flags=flags, coarse_gates=gates(coarse), fine_gates=gates(refined))


def strong_rate_followup():
    """Post-primary checks of the surviving high-rate cell and its controls."""
    rows = []
    for mode, limit in (("redshift",.025),("redshift",.0125),("timing",.025),
                        ("timing",.0125),("open",.025)):
        for key in [str(i) for i in range(9)] + ["bright"]:
            if mode != "redshift" and key not in ("0","4","bright"):
                continue
            if key == "bright":
                pair = [.875,.875]
            else:
                weight = base.split_weights(base.THETAS[int(key)], .02)[0]
                pair = [weight,1-weight]
            row, _ = evolve(pair,3.,mode,limit)
            rows.append(dict(row,input=key))
    return rows


if __name__ == "__main__":
    print(json.dumps(probe(), allow_nan=False))
