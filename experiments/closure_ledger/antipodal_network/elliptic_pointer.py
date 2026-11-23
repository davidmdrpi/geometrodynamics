"""Classical winding doublet coupled to a switched mechanical pointer (PR #241 mixer)."""
import json

import numpy as np
from scipy.integrate import solve_ivp

DAMPING = 0.4
GAIN = 0.1 * (2 - np.sqrt(2))
FREQUENCY = np.sqrt(1 - (DAMPING / 2)**2)
WRITE_TIME = np.pi / FREQUENCY
DECAY = np.exp(-DAMPING * WRITE_TIME / 2)
MATCHED_WIDTH = (1 + DECAY) * GAIN / DECAY
INPUTS = np.linspace(-1, 1, 9)


def channels(winding, angle):
    """Amplitudes in the two eigenchannels of the unit elliptic mixer."""
    a, b = winding
    return np.array([a + np.exp(1j * angle) * b,
                     a - np.exp(1j * angle) * b]) / np.sqrt(2)


def integrate(rhs, duration, rows):
    shape = rows.shape

    def flat_rhs(t, y):
        return np.asarray(rhs(y.reshape(shape))).ravel()

    result = solve_ivp(flat_rhs, (0, duration), rows.ravel(), method="DOP853",
                       t_eval=[duration], rtol=2e-11, atol=2e-13, max_step=0.5)
    if not result.success:
        raise RuntimeError(result.message)
    return result.y[:, -1].reshape(shape)


def latch(q, v, heat):
    def rhs(y):
        x, speed, _ = y
        return speed, x*(1-x*x) - DAMPING*speed, DAMPING*speed*speed

    x, speed, total_heat = integrate(rhs, 80, np.stack([q, v, heat]))
    energy = speed**2/2 + (x*x-1)**2/4
    locked = (energy < 0.25) & (abs(abs(x)-1) < 0.001) & (abs(speed) < 0.001)
    return x, speed, total_heat, energy, np.where(locked, np.sign(x), 0).astype(int)


def evolve(q0, plus, minus, duration=WRITE_TIME):
    """Integrate the reciprocal amplitude and pointer equations for every launch."""
    q0, plus, minus = np.broadcast_arrays(q0, plus, minus)
    zeros = np.zeros_like(q0)
    z0 = abs(plus)**2 - abs(minus)**2
    energy0 = q0**2 / 2 - GAIN*q0*z0

    def write(y):
        q, v, a, b, c, d, _ = y
        phase_rate = GAIN*q
        force = GAIN*(a*a + b*b - c*c - d*d)
        return (v, force - q - DAMPING*v, -phase_rate*b, phase_rate*a,
                phase_rate*d, -phase_rate*c, DAMPING*v*v)

    q, v, a, b, c, d, heat = integrate(write, duration, np.stack(
        [q0, zeros, plus.real, plus.imag, minus.real, minus.imag, zeros]))
    plus_end, minus_end = a + 1j*b, c + 1j*d
    z_end = abs(plus_end)**2 - abs(minus_end)**2
    before = (v*v + q*q)/2 - GAIN*q*z_end
    after = v*v/2 + (q*q-1)**2/4
    work = after - before

    x, speed, total_heat, final_energy, records = latch(q, v, heat)
    return dict(records=records,
                q_switch=q, v_switch=v, q_final=x, v_final=speed, plus_end=plus_end,
                minus_end=minus_end, work=work, heat=total_heat,
                ledger_error=np.maximum(abs(before+heat-energy0),
                                        abs(final_energy+total_heat-energy0-work)),
                norm_error=abs(abs(plus_end)**2 + abs(minus_end)**2 - 1),
                population_error=abs(z_end-z0))


def ensemble(count, width, density, time_factor=1):
    if density not in ("uniform", "circle_projection"):
        raise ValueError(f"Unknown pointer preparation: {density}")
    u = (np.arange(count)+0.5)/count
    initial = width*(2*u-1 if density == "uniform" else -np.cos(np.pi*u))
    initial, z = np.broadcast_arrays(initial, INPUTS[:, None])
    result = evolve(initial, np.sqrt((1+z)/2), np.sqrt((1-z)/2), time_factor*WRITE_TIME)
    records = result["records"]
    positive = np.sum(records == 1, axis=1)
    p = positive/count
    summary = dict(count=count, width=float(width), density=density,
                   write_time_factor=time_factor, positive=positive.tolist(),
                   negative=np.sum(records == -1, axis=1).tolist(),
                   no_record=np.sum(records == 0, axis=1).tolist(),
                   max_born_error=float(max(abs(p-(1+INPUTS)/2))),
                   max_energy_residual=float(result["ledger_error"].max()),
                   max_norm_error=float(result["norm_error"].max()),
                   max_population_drift=float(result["population_error"].max()))
    if time_factor == 1:
        predicted = (1+DECAY)*GAIN*z - DECAY*initial
        threshold = np.clip(MATCHED_WIDTH*INPUTS/width, -1, 1)
        limiting_p = ((1+threshold)/2 if density == "uniform"
                      else 0.5 + np.arcsin(threshold)/np.pi)
        summary.update(max_switch_position_error=float(abs(result["q_switch"]-predicted).max()),
                       max_switch_speed=float(abs(result["v_switch"]).max()),
                       max_abs_switch_position=float(abs(result["q_switch"]).max()),
                       basin_disagreements=int(np.sum(records != np.sign(predicted))),
                       continuum_positive=limiting_p.tolist())
    return summary


def probe():
    centered = evolve(np.array([0.]), np.sqrt([0.5]), np.sqrt([0.5]))
    controls = []
    for phase in (0, np.pi/2, np.pi):
        amplitudes = channels(np.array([1, np.exp(1j*phase)])/np.sqrt(2), 0)
        controls.append(dict(relative_phase=float(phase),
                             channel_difference=float(abs(amplitudes[0])**2-abs(amplitudes[1])**2)))
    seeds = np.array([-0.06, 0.0, 0.06])
    balanced = evolve(seeds, np.sqrt(0.5), np.sqrt(0.5))
    examples = [dict(q0=float(seed), q_switch=float(balanced["q_switch"][i]),
                     q_final=float(balanced["q_final"][i]), record=int(balanced["records"][i]),
                     switch_work=float(balanced["work"][i]), radiated=float(balanced["heat"][i]))
                for i, seed in enumerate(seeds)]
    return dict(balanced_examples=examples, inputs=INPUTS.tolist(), gain=float(GAIN), write_time=float(WRITE_TIME),
                decay=float(DECAY), matched_width=float(MATCHED_WIDTH),
                balanced_center=dict(q_final=centered["q_final"].item(),
                                     record=centered["records"].item()),
                equal_winding_controls=controls,
                primary=[ensemble(n, 0.12, "uniform") for n in (256, 512)],
                matched_uniform=ensemble(512, MATCHED_WIDTH, "uniform"),
                matched_circle=ensemble(512, MATCHED_WIDTH, "circle_projection"),
                timing_controls=[ensemble(256, 0.12, "uniform", f) for f in (0.9, 1.1)])


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2, allow_nan=False))
