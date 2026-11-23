"""Hilbert energy current of the two winding channels in a minimally coupled 5D KK lift."""
import json

import numpy as np

from .elliptic_transport import END, FORCE, SIGMA, TAU, packet, winding_weights

K = np.array([[1.], [-1.]])
TIMES = (0., 0.5, 1., 2., 5., 10., END)


def connection(t, x, offset=0.):
    """Linear SG force where the packets are driven; taper before the periodic seam."""
    taper = np.cos(np.pi/2*np.clip((abs(x)-36)/18, 0, 1))**2
    amplitude = 2*FORCE*np.sin(np.pi*t/TAU)**2 if t < TAU else 0.
    return offset - amplitude*x*taper


def currents(u, pi, dxu, c, a):
    """pi is the carrier-restored D_t derivative (fast common phase factored out)."""
    dt = pi - 1j*K*a*u
    norm = np.sum(abs(u)**2, axis=0)
    norm_flux = np.sum(np.imag(u.conj()*dxu), axis=0)
    lag_trace = np.sum(-abs(pi)**2/c**2+abs(dxu)**2+c**2*abs(K*u)**2, axis=0)
    energy = np.sum((pi.conj()*dt).real, axis=0)/c**2 + lag_trace/2
    flux = -np.sum((dxu.conj()*dt).real, axis=0)
    return norm, norm_flux, energy, flux


def direct_fiber_contraction(u, pi, dxu, c, a):
    """Average T^A_t over chi using the 5x5 metric, before taking the ratio."""
    radius = 1/c
    metric = np.diag([-c*c, 1., 1., 1., radius**2])
    metric[0,0] += (radius*a)**2
    metric[0,4] = metric[4,0] = -radius**2*a
    inverse = np.linalg.inv(metric)
    phase = np.exp(1j*K*(np.arange(16)*2*np.pi/16))
    grad = np.zeros((16,5), complex)
    grad[:,0] = (pi-1j*K[:,0]*a*u) @ phase
    grad[:,1] = dxu @ phase
    grad[:,4] = (1j*K[:,0]*u) @ phase
    product = np.einsum("na,nb->nab",grad.conj(),grad).real
    trace = np.einsum("ab,nab->n",inverse,product)
    mixed = inverse @ (product.mean(axis=0)-metric*trace.mean()/2)
    return -mixed[0,0], -mixed[1,0]


def exact_control(c=8., p=1., offset=0.1):
    """Exact stationary local KG data at x=0, including a possible nonzero A'(0)."""
    intrinsic = np.sqrt(c**4+c*c*p*p)
    u = np.full((2,1),1/np.sqrt(2))
    pi = -1j*intrinsic*u
    dxu = 1j*K*p*u
    rho,j,e,s = currents(u,pi,dxu,c,offset)
    return dict(c=c, momentum=p, channel_potential_gap=2*offset,
                norm_velocity=float((j/rho)[0]), energy_velocity=float((s/e)[0]),
                analytic=c*c*p*offset/intrinsic**2, energy_density=float(e[0]))


def archived_field_check(t, c, weight, n=1024):
    """Evaluate stress on #358's actual Schrödinger envelopes as an NR lift.
    The Schrödinger equation fixes their full time derivatives, including phase.
    """
    x = np.linspace(-56.,56.,n,endpoint=False)
    d, speed, width, dilation = packet(t)
    other = x-K*d
    u = np.sqrt(np.array([[weight],[1-weight]])/(np.sqrt(2*np.pi)*width)) * np.exp(
        -other**2/(4*width**2)+1j*(K*speed*x+dilation*other**2/2))
    log_x = -other/(2*width**2)+1j*(K*speed+dilation*other)
    dxu = u*log_x
    # D_t removes the k A potential in this lift. The common rest phase is retained.
    pi = -1j*c*c*u + .5j*u*(log_x**2-1/(2*width**2)+1j*dilation)
    a = connection(t,x)
    rho,j,e,s = currents(u,pi,dxu,c,a)
    mask = rho > rho.max()*1e-6
    error = s[mask]/e[mask]-j[mask]/rho[mask]
    return dict(t=t, rms_energy_minus_358=float(np.sqrt(np.average(error**2,weights=rho[mask]))),
                midpoint_energy=float(e[n//2]))


def run(c, dt, n, weight, offset=0.):
    """Strang split exact free KG and KK shift. Positive-frequency initialization."""
    x = np.linspace(-56.,56.,n,endpoint=False)
    spacing = x[1]-x[0]
    q = 2*np.pi*np.fft.fftfreq(n,spacing)
    omega = np.sqrt(c**4+c*c*q*q)
    g = np.exp(-x*x/(4*SIGMA**2))/(2*np.pi*SIGMA**2)**0.25
    u = np.sqrt([[weight],[1-weight]])*g
    pi = np.fft.ifft(-1j*omega*np.fft.fft(u))
    carrier = np.exp(1j*c*c*dt)
    co = carrier*np.cos(omega*dt)
    si = carrier*np.sin(omega*dt)
    snapshots = []
    checkpoints = {round(t/dt):t for t in TIMES}
    mid = n//2
    one = np.argmin(abs(x-1))
    stationary_energies = []
    work = 0.
    drive_energies = []
    for step in range(round(END/dt)+1):
        if step*dt <= TAU:
            # Explicit metric work: derivative of the Hamiltonian at fixed (u,pi).
            charge = (pi.conj()*(-1j*K*u)).real.sum(axis=0)/c**2
            a_dot = connection(TAU/2,x)*np.pi/TAU*np.sin(2*np.pi*step*dt/TAU)
            work += dt*spacing*np.sum(a_dot*charge)  # endpoint powers are zero
        if step in checkpoints:
            t = checkpoints[step]
            a = connection(t,x,offset)
            dxu = np.fft.ifft(1j*q*np.fft.fft(u))
            rho,j,e,s = currents(u,pi,dxu,c,a)
            # Avoid numerical vacuum in all reported velocity errors.
            mask = rho > 1e-6*rho.max()
            v, ve = j[mask]/rho[mask], s[mask]/e[mask]
            d, speed, width, dilation = packet(t)
            dens = np.array([weight,1-weight])[:,None]*np.exp(-(x[mask]-K*d)**2/(2*width**2))
            reference = (dens*(K*speed+dilation*(x[mask]-K*d))).sum(axis=0)/dens.sum(axis=0)
            def rms(error):
                return float(np.sqrt(np.average(error**2,weights=rho[mask])))
            ed, fl = direct_fiber_contraction(u[:,one],pi[:,one],dxu[:,one],c,a[one])
            snapshots.append(dict(t=t, rms_energy_minus_norm=rms(ve-v),
                rms_energy_minus_358=rms(ve-reference), rms_norm_minus_358=rms(v-reference),
                midpoint_energy=float(e[mid]), midpoint_energy_v=float(s[mid]/e[mid]),
                midpoint_norm_v=float(j[mid]/rho[mid]), x_sample=float(x[one]),
                potential_gap_at_sample=float(2*a[one]),
                sample_energy_v=float(s[one]/e[one]), sample_norm_v=float(j[one]/rho[one]),
                direct_tensor_residual=float(max(abs(ed-e[one]),abs(fl-s[one]))),
                minimum_energy_over_density=float(min(e[mask]/rho[mask]))))
            if t in (0.,TAU):
                drive_energies.append(float(spacing*e.sum()))
            if t >= TAU:
                stationary_energies.append(float(spacing*e.sum()))
        if step == round(END/dt):
            break
        half = np.exp(-.5j*K*connection((step+.5)*dt,x,offset)*dt)
        uk, pk = np.fft.fft(half*u), np.fft.fft(half*pi)
        u = half*np.fft.ifft(co*uk+si*pk/omega)
        pi = half*np.fft.ifft(co*pk-si*omega*uk)
    return dict(c=c,dt=dt,n=n,plus_weight=float(weight),offset=offset,
                metric_work=float(work), drive_energy_change=drive_energies[1]-drive_energies[0],
                drive_budget_relative=abs(drive_energies[1]-drive_energies[0]-work)/drive_energies[0],
                free_energy_relative_drift=float(np.ptp(stationary_energies)/stationary_energies[0]),
                snapshots=snapshots)


def probe():
    weights = winding_weights(np.pi/8,.02)
    w = weights[1]/(weights[1]+weights[-1])
    return dict(exact_equal=exact_control(offset=0.), exact_split=exact_control(),
                archived_field=[archived_field_check(t,8.,w) for t in TIMES],
                discarded_population=float(1-weights[1]-weights[-1]),
                balanced=run(8.,.01,1024,.5), weak_offset=run(8.,.01,1024,.5,.1),
                mixed=[run(c,.01,1024,w) for c in (8.,12.,16.)],
                refinement=run(8.,.005,2048,w))


if __name__ == "__main__":
    print(json.dumps(probe(),indent=2,allow_nan=False))
