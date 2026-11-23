"""Independently locate a nonlinear KG core and integrate the two candidate flows."""
import json

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

PLUS_WEIGHT = 0.8535130337631147
WINDINGS = np.array([[1.], [-1.]])
HALF_BOX = 32.


def profile(x):
    """Rest profile and its derivative for inverse decay length one."""
    denominator = 2 + np.sqrt(2)*np.cosh(2*x)
    amplitude = np.sqrt(2/denominator)
    derivative = -amplitude*np.sqrt(2)*np.sinh(2*x)/denominator
    return amplitude, derivative


def interp(values, x, h):
    """Local cubic interpolation on the periodic uniform grid."""
    index = (x+HALF_BOX)/h
    cell = int(np.floor(index))
    z = index-cell
    weights = np.array([-z*(z-1)*(z-2), 3*(z+1)*(z-1)*(z-2),
                        -3*(z+1)*z*(z-2), (z+1)*z*(z-1)])/6
    return values[..., (cell+np.arange(-1, 3)) % values.shape[-1]] @ weights


def core_location(x, density):
    i = np.argmax(density)
    if i < 4 or i+4 >= x.size:
        raise ValueError("core too close to edge")
    band = slice(i-4, i+5)
    spline = CubicSpline(x[band], density[band])
    peak = brentq(spline.derivative(), x[i-1], x[i+1])
    local = (density > np.roll(density, 1)) & (density >= np.roll(density, -1))
    local[i] = False
    secondary = np.max(density[local], initial=0.)/spline(peak)
    return peak, spline(peak).item(), spline(peak, 2).item(), secondary


def run(c=8., n=1024, dt=.001, stride=10, force=1., offset=0.,
        weight=PLUS_WEIGHT, boost=0.):
    """c=None runs the matched nonlinear Schrödinger limit."""
    x = np.linspace(-HALF_BOX, HALF_BOX, n, endpoint=False)
    h = x[1]-x[0]
    q = 2*np.pi*np.fft.fftfreq(n, h)
    # The taper reaches zero before the periodic seam.
    taper = np.cos(np.pi/2*np.clip((abs(x)-16)/12, 0, 1))**2
    gradient_shape = -force*x*taper
    gamma = 1. if c is None else 1/np.sqrt(1-(boost/c)**2)
    omega = 0. if c is None else c*np.sqrt(c*c-1)
    amplitude, derivative = profile(gamma*x)
    phase_gradient = boost if c is None else omega*gamma*boost/c**2
    basis = np.sqrt([[weight], [1-weight]])*np.exp(1j*phase_gradient*x)
    u = basis*amplitude
    z = basis*(-gamma*boost*derivative + 1j*omega*(1-gamma)*amplitude)
    state = np.stack([u, z])
    nstep = round(2/dt)
    if not np.isclose(nstep*dt, 2) or nstep % stride:
        raise ValueError("observation steps must include the end of the drive")

    def rhs(t, state):
        u, z = state
        potential = offset + gradient_shape*np.sin(np.pi*t/2)**2
        s = np.sum(abs(u)**2, axis=0)
        lap = np.fft.ifft(-q*q*np.fft.fft(u))
        spatial = lap + (-1 + 4*s - 1.5*s*s)*u
        rotation = -1j*WINDINGS*potential
        if c is None:
            return np.stack([.5j*spatial + rotation*u, np.zeros_like(z)])
        return np.stack([z+rotation*u, (2j*omega+rotation)*z+c*c*spatial])

    tracers = np.zeros(2)
    history = []
    previous_velocity = None
    work = 0.
    previous_power = 0.
    for step in range(nstep+1):
        t = step*dt
        if step % stride == 0:
            u, z = state
            s = np.sum(abs(u)**2, axis=0)
            peak, height, curvature, secondary = core_location(x, s)
            grad = np.fft.ifft(1j*q*np.fft.fft(u))
            j = np.imag(u.conj()*grad).sum(axis=0)
            if c is None:
                numerators, denominators = (j, j), (s, s)
                total_charge, total_energy, power = h*s.sum(), None, 0.
            else:
                pi = z-1j*omega*u
                potential = offset + gradient_shape*np.sin(np.pi*t/2)**2
                number = -np.imag(u.conj()*pi)/c**2
                signed = (WINDINGS*number).sum(axis=0)
                energy = (.5*((abs(pi)**2/c**2 + abs(grad)**2).sum(axis=0)
                              + c*c*s - 2*s*s + .5*s**3) + potential*signed)
                flux = -(grad.conj()*(pi-1j*WINDINGS*potential*u)).real.sum(axis=0)
                numerators, denominators = (flux, j), (energy, number.sum(axis=0))
                total_charge, total_energy = h*number.sum(), h*energy.sum()
                power = h*np.sum(signed*gradient_shape)*np.pi/2*np.sin(np.pi*t)

            def speeds(positions):
                result = []
                for num, den, pos in zip(numerators, denominators, positions):
                    norm = interp(den, pos, h)
                    if norm <= 0:
                        raise ArithmeticError("nonpositive tracer density")
                    result.append(interp(num, pos, h)/norm)
                return np.array(result)

            observation_dt = stride*dt
            if previous_velocity is None:
                tracers[:] = peak
            else:
                trial = tracers + observation_dt*previous_velocity
                tracers += observation_dt*(previous_velocity+speeds(trial))/2
                work += observation_dt*(previous_power+power)/2
            # A conserved number median independently tests the number tracer.
            # Periodic tails are negligible here; the median starts at the origin.
            charge_density = denominators[1]
            nodes = np.append(x, HALF_BOX)
            spline = CubicSpline(nodes, np.append(charge_density, charge_density[0]),
                                 bc_type="periodic").antiderivative()
            half_charge = (spline(-HALF_BOX)+spline(HALF_BOX))/2
            median = brentq(lambda pos: spline(pos)-half_charge, -8., 8.)

            # Independent trigonometric peak and instantaneous velocity from s_x=0.
            # The exact band-limited interpolant uses grid index 0 at -HALF_BOX.
            spectrum = np.fft.fft(s)/n
            spectral_peak = peak
            for _ in range(3):
                phase = np.exp(1j*q*(spectral_peak+HALF_BOX))
                slope = np.sum(1j*q*spectrum*phase).real
                second = np.sum(-q*q*spectrum*phase).real
                spectral_peak -= slope/second
            if abs(spectral_peak-peak) > h/2:
                raise ArithmeticError("peak locators disagree on the selected core")
            if c is None:
                density_rate = -np.fft.ifft(1j*q*np.fft.fft(j)).real
            else:
                density_rate = 2*(u.conj()*z).real.sum(axis=0)
            phase = np.exp(1j*q*(spectral_peak+HALF_BOX))
            mixed_derivative = np.sum(1j*q*np.fft.fft(density_rate)*phase).real/n
            curvature_spectral = np.sum(-q*q*spectrum*phase).real
            core_speed = -mixed_derivative/curvature_spectral
            local_energy_speed, local_number_speed = speeds([spectral_peak]*2)

            previous_velocity = speeds(tracers)
            previous_power = power
            history.append(dict(t=t, core=peak, energy_path=tracers[0], number_path=tracers[1],
                                height=height, curvature=curvature, secondary_ratio=secondary,
                                spectral_peak=spectral_peak, core_speed=core_speed,
                                energy_speed_at_peak=local_energy_speed,
                                number_speed_at_peak=local_number_speed, number_median=median,
                                total_charge=total_charge, total_energy=total_energy, work=work))
        if step < nstep:
            k1 = rhs(t, state)
            k2 = rhs(t+dt/2, state+dt*k1/2)
            k3 = rhs(t+dt/2, state+dt*k2/2)
            k4 = rhs(t+dt, state+dt*k3)
            state += dt*(k1+2*k2+2*k3+k4)/6
    return history


def summarize(rows):
    def difference(a, b):
        return max(abs(r[a]-r[b]) for r in rows)

    q0 = rows[0]["total_charge"]
    e0 = rows[0]["total_energy"]
    return dict(core_energy_error=difference("core", "energy_path"),
                core_number_error=difference("core", "number_path"),
                predictor_separation=difference("energy_path", "number_path"),
                median_number_path_error=difference("number_median", "number_path"),
                peak_locator_disagreement=difference("core", "spectral_peak"),
                max_secondary=max(r["secondary_ratio"] for r in rows),
                min_height_ratio=min(r["height"]/rows[0]["height"] for r in rows),
                min_curvature_magnitude=min(-r["curvature"] for r in rows),
                max_charge_drift=max(abs(r["total_charge"]/q0-1) for r in rows),
                max_work_residual=None if e0 is None else
                max(abs(r["total_energy"]-e0-r["work"])/e0 for r in rows),
                endpoint=rows[-1])


def probe():
    settings = dict(primary={}, fine=dict(n=2048, dt=.0005), c12=dict(c=12.),
                    c16=dict(c=16.), nonrelativistic=dict(c=None), drive_off=dict(force=0.),
                    reversed=dict(force=-1.), pure=dict(weight=1.), offset=dict(offset=.1),
                    boost=dict(boost=.5, force=0.))
    # After the registered comparison, isolate the one-winding relativistic residual.
    settings.update(pure_fine=dict(weight=1., n=2048, dt=.0005),
                    pure12=dict(weight=1., c=12.), pure16=dict(weight=1., c=16.),
                    pure_nr=dict(weight=1., c=None))
    histories = {name: run(**params) for name, params in settings.items()}
    primary, fine = histories["primary"], histories["fine"]
    change = max(abs(a[key]-b[key]) for a, b in zip(primary, fine[::2])
                 for key in ("core", "energy_path", "number_path"))
    tol = max(1e-5, 5*change)
    results = {k: summarize(v) for k, v in histories.items()}
    ref = results["fine"]
    resolved = ref["predictor_separation"] > 2*tol and ref["max_secondary"] < .5
    matches = [law for law in ("number", "energy") if ref["core_"+law+"_error"] <= tol]
    decision = "unresolved"
    if resolved:
        decision = "+".join(matches) if matches else "neither"
    trace_keys = ("t", "core", "energy_path", "number_path", "number_median",
                  "spectral_peak", "core_speed", "energy_speed_at_peak", "number_speed_at_peak")
    traces = {}
    for name in ("primary", "fine", "nonrelativistic"):
        rows = histories[name]
        step = round(.05/(rows[1]["t"]-rows[0]["t"]))
        traces[name] = [{k: row[k] for k in trace_keys} for row in rows[::step]]
    return dict(settings=settings, defaults=dict(c=8., n=1024, dt=.001, stride=10, weight=PLUS_WEIGHT),
                results=results, max_refinement_change=change, tolerance=tol, verdict=decision,
                # Summary maxima use every observation; save trajectories every 0.05.
                traces=traces)


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2, allow_nan=False))
