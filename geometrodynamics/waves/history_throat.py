"""Reciprocal history-coupled scalar/neck action; see its frozen protocol.

Exterior time t and neck time sigma are distinct. Clock energies are canonical
variables that change the mouth inertia. A reversed clock map requires a global
history solve, not a forward scheduled impulse. This is an effective model.
"""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial import legendre
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator

from .coupled_mouth_throat import CoupledMouthThroat, Parameters as OldParameters, real_basis


@dataclass(frozen=True)
class Parameters:
    lmax: int = 8
    nodes: int = 17
    length: float = .565
    mass: float = 1.
    charge: float = .5
    mouth_angle: float = .3
    footprint: float = .20
    packet_width: float = .35
    packet_energy: float = .02
    coupling: float = .5
    spatial_parity: int = -1
    time_orientation: int = 1
    horizon: float = 8.
    neck_horizon: float = 8.
    spacing: float = .02
    birth_a: float = .30
    birth_b: float = -.30
    endpoint_load: float = 0.
    blocked: bool = False
    rtol: float = 2e-9
    fixed_point_tolerance: float = 2e-7
    max_iterations: int = 80
    relaxation: float = .5


class EndpointHistory:
    def __init__(self, times, values, velocities):
        self.times = np.asarray(times)
        self.curve = CubicHermiteSpline(times, values, velocities, axis=0, extrapolate=False)

    def sample(self, clocks):
        clocks = np.asarray(clocks)
        if np.min(clocks) < self.times[0]-1e-10 or np.max(clocks) > self.times[-1]+1e-10:
            raise ValueError('Requested neck event lies outside the computed history')
        # Roundoff only; physical out-of-range events were rejected above.
        clocks = np.clip(clocks, self.times[0], self.times[-1])
        return np.diag(self.curve(clocks)), np.diag(self.curve(clocks, 1))


class HistoryThroat:
    def __init__(self, parameters=Parameters()):
        p = self.par = parameters
        if (p.lmax < 2 or p.nodes < 3 or min(p.mass, p.length, p.footprint,
                p.packet_width, p.horizon, p.spacing) <= 0 or p.neck_horizon < p.horizon
                or min(p.packet_energy, p.coupling, p.endpoint_load) < 0
                or p.spatial_parity not in (-1, 1) or p.time_orientation not in (-1, 1)
                or max(abs(p.birth_a), abs(p.birth_b)) >= 1
                or not 0 < p.relaxation <= 1 or p.max_iterations < 1):
            raise ValueError('Invalid scale, resolution, parity or history domain')
        basis_model = CoupledMouthThroat(OldParameters(lmax=p.lmax, offset=0., nodes=3,
                                                     footprint=p.footprint))
        self.n = basis_model.ns
        self.ell, self.d = basis_model.ell, basis_model.d
        self.omega = self.ell+.5
        self.smooth = np.exp(-.5*p.footprint**2*self.ell*(self.ell+1))
        # Exact finite Fourier representation of each harmonic on this great circle.
        count = 2*p.lmax+1
        angles = np.arange(count)*2*np.pi/count
        values = real_basis(p.lmax, np.arccos(np.cos(angles)),
                            np.where(np.sin(angles) >= 0, 0., np.pi))*self.smooth
        coef = np.fft.fft(values, axis=0)/count
        self.freq = np.fft.fftfreq(count, 1/count)
        self.coef = coef
        self.eta = np.array([1., p.spatial_parity])
        self.epsilon = np.array([1., p.time_orientation])
        self.anchor = np.array([0., p.neck_horizon if p.time_orientation == -1 else 0.])
        h = p.length/(p.nodes-1)
        self.mn = np.full(p.nodes, h); self.mn[[0, -1]] *= .5
        self.kn = np.zeros((p.nodes, p.nodes))
        for j in range(p.nodes-1):
            self.kn[j:j+2, j:j+2] += np.array([[1., -1.], [-1., 1.]])/h
        ls = np.arange(p.lmax+1)
        self.coulomb_coef = np.zeros(p.lmax+1)
        self.coulomb_coef[1:] = (np.exp(-p.footprint**2*ls[1:]*(ls[1:]+1))
            *(2*ls[1:]+1)/(4*np.pi*ls[1:]*(ls[1:]+1)))
        self.coulomb_derivative = legendre.legder(self.coulomb_coef)
        self.sigma = np.linspace(0., p.neck_horizon, int(np.ceil(p.neck_horizon/p.spacing))+1)

    def footprint(self, q):
        e = np.exp(1j*np.asarray(q)[..., None]*self.freq)
        return (e@self.coef).real, ((e*(1j*self.freq))@self.coef).real

    def coulomb(self, q):
        p = self.par; delta = q[0]-q[1]; co = np.cos(delta)
        value = p.charge**2*(legendre.legval(1., self.coulomb_coef)
                            -legendre.legval(co, self.coulomb_coef))
        derivative = p.charge**2*legendre.legval(co, self.coulomb_derivative)*np.sin(delta)
        return value, np.array([derivative, -derivative])

    def unpack(self, state):
        n = self.n
        return state[:n], state[n:2*n], *state[2*n:].reshape(4, 2)

    def initial(self):
        p = self.par
        state = np.zeros(2*self.n+8)
        lam = self.ell*(self.ell+1)
        u = lam*np.exp(-.5*p.packet_width**2*lam)*real_basis(p.lmax, 0., 0.)
        u *= np.sqrt(p.packet_energy/(.5*np.dot(self.omega*u, self.omega*u)))
        state[:self.n] = u
        state[2*self.n:2*self.n+2] = [p.mouth_angle, -p.mouth_angle]
        speeds = np.array([p.birth_a, p.birth_b])
        state[2*self.n+2:2*self.n+4] = p.mass*speeds/np.sqrt(1-speeds**2)
        return state

    def quantities(self, state, neck):
        u, momentum, q, p, clocks, clock_energy = self.unpack(state)
        f, fq = self.footprint(q)
        values, derivatives = neck.sample(self.epsilon*clocks+self.anchor)
        gap = values-self.eta*(f@u)
        potential = .5*self.par.coupling*gap**2
        mu = self.par.mass+clock_energy+potential
        if np.min(mu) <= 0:
            raise ValueError('Nonpositive effective mouth mass')
        energy = np.sqrt(p*p+mu*mu)
        return f, fq, gap, derivatives, potential, mu, energy

    def hamiltonian(self, state, neck):
        u, momentum, q, _, _, _ = self.unpack(state)
        *_, energy = self.quantities(state, neck)
        return float(np.sum(energy)+self.coulomb(q)[0]
                     +.5*(momentum@momentum+np.dot(self.omega*u, self.omega*u)))

    def angular(self, state):
        u, momentum, _, p, _, _ = self.unpack(state)
        return float(np.sum(p)-momentum@self.d@u)

    def rhs(self, t, state, neck):
        u, momentum, q, p, _, _ = self.unpack(state)
        f, fq, gap, zdot, _, mu, energy = self.quantities(state, neck)
        rate = mu/energy
        g = self.par.coupling
        source = rate*g*self.eta*gap
        pdot = -self.coulomb(q)[1]+source*(fq@u)
        kdot = -rate*g*gap*self.epsilon*zdot
        return np.r_[momentum, -self.omega**2*u+source@f, p/energy, pdot, rate, kdot]

    def exterior(self, neck):
        def contact(t, state):
            q = state[2*self.n:2*self.n+2]
            return np.arccos(np.cos(q[0]-q[1]))-2*self.par.mouth_angle
        contact.terminal = True; contact.direction = -1
        sol = solve_ivp(lambda t, y: self.rhs(t, y, neck), (0., self.par.horizon),
            self.initial(), method='DOP853', rtol=self.par.rtol, atol=self.par.rtol*.01,
            max_step=min(.04, 2*self.par.spacing), dense_output=True, events=contact)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    def neck(self, exterior):
        p = self.par; n = p.nodes
        times = np.linspace(0., exterior.t[-1], int(np.ceil(exterior.t[-1]/p.spacing))+1)
        samples = exterior.sol(times).T
        clocks = samples[:, 2*self.n+4:2*self.n+6]
        inverse = [PchipInterpolator(clocks[:, i], times, extrapolate=False) for i in range(2)]
        domains = [sorted([self.anchor[i], self.epsilon[i]*clocks[-1, i]+self.anchor[i]])
                   for i in range(2)]

        def rhs(sigma, state):
            z, v = state[:n], state[n:2*n]
            force = -self.kn@z
            for i, endpoint in enumerate([0, n-1]):
                if domains[i][0] <= sigma <= domains[i][1]:
                    proper = (sigma-self.anchor[i])/self.epsilon[i]
                    t = float(inverse[i](np.clip(proper, 0., clocks[-1, i])))
                    y = exterior.sol(t)
                    u, _, q, _, _, _ = self.unpack(y)
                    f, _ = self.footprint(q[i:i+1])
                    force[endpoint] -= p.coupling*(z[endpoint]-self.eta[i]*float(f[0]@u))
            force[[0, -1]] -= p.endpoint_load*v[[0, -1]]
            return np.r_[v, force/self.mn, p.endpoint_load*np.sum(v[[0, -1]]**2)]

        # Split at coverage edges, where a finite worldtube segment ends.
        cuts = np.unique(np.r_[0., p.neck_horizon, np.ravel(domains)])
        cuts = cuts[(cuts >= 0) & (cuts <= p.neck_horizon)]
        y = np.zeros(2*n+1); pieces = []
        for a, b in zip(cuts[:-1], cuts[1:]):
            if b-a < 1e-12:
                continue
            sol = solve_ivp(rhs, (a, b), y, method='DOP853', rtol=p.rtol,
                atol=p.rtol*.01, max_step=min(.02, p.spacing), dense_output=True)
            if not sol.success:
                raise RuntimeError(sol.message)
            pieces.append(sol); y = sol.y[:, -1]
        out = np.empty((len(self.sigma), 2*n+1))
        for sol in pieces:
            mask = (self.sigma >= sol.t[0]-1e-12) & (self.sigma <= sol.t[-1]+1e-12)
            out[mask] = sol.sol(np.clip(self.sigma[mask], sol.t[0], sol.t[-1])).T
        return out

    def solve(self, verbose=False):
        p = self.par
        values = np.zeros((len(self.sigma), 2)); velocities = values.copy()
        residuals = []; converged = False
        for iteration in range(p.max_iterations):
            history = EndpointHistory(self.sigma, values, velocities)
            exterior = self.exterior(history)
            if p.blocked or p.coupling == 0:
                states = np.zeros((len(self.sigma), 2*p.nodes+1))
            else:
                states = self.neck(exterior)
            new_values = states[:, [0, p.nodes-1]]
            new_velocities = states[:, [p.nodes, 2*p.nodes-1]]
            error = max(float(np.max(abs(new_values-values))),
                        float(np.max(abs(new_velocities-velocities))))
            residuals.append(error)
            if verbose:
                print(f'iteration={iteration+1} residual={error:.6g}', flush=True)
            if error < p.fixed_point_tolerance:
                converged = True
                break
            values = (1-p.relaxation)*values+p.relaxation*new_values
            velocities = (1-p.relaxation)*velocities+p.relaxation*new_velocities
        times = np.linspace(0., exterior.t[-1], int(np.ceil(exterior.t[-1]/p.spacing))+1)
        ys = exterior.sol(times).T
        h = np.array([self.hamiltonian(y, history) for y in ys])
        j = np.array([self.angular(y) for y in ys])
        n = p.nodes
        neck_energy = .5*np.sum(states[:, n:2*n]**2*self.mn, axis=1)
        neck_energy += .5*np.einsum('ti,ij,tj->t', states[:, :n], self.kn, states[:, :n])
        k = ys[-1, 2*self.n+6:2*self.n+8]
        work_error = float(neck_energy[-1]+states[-1, -1]-self.epsilon@k)
        energy_error = float(np.max(abs(h-h[0]))/max(abs(h[0]), 1e-30))
        angular_error = float(np.max(abs(j-j[0])))
        summary = dict(converged=converged, iterations=len(residuals), residuals=residuals,
            exterior_final_time=float(times[-1]), stop='contact' if len(exterior.t_events[0]) else 'horizon',
            initial_exterior_energy=float(h[0]), energy_relative_residual=energy_error,
            angular_absolute_residual=angular_error, neck_clock_work_residual=work_error,
            neck_final_energy=float(neck_energy[-1]), endpoint_load_energy=float(states[-1, -1]),
            final_clock_energies=k.tolist(),
            ledger_pass=bool(energy_error < 1e-5 and angular_error < 1e-6
                            and abs(work_error) < 1e-5*max(p.packet_energy, 1e-12)))
        return dict(model=self, exterior=exterior, endpoint_history=history, t=times,
            states=ys, sigma=self.sigma, neck_states=states, neck_energy=neck_energy,
            exterior_energy=h, angular=j, summary=summary)

    def transfer_event(self, result, emission, source, destination, path_length=None):
        """Kinematic ray prediction on the solved history, with no extrapolation.

        Finite-element packets are dispersive: these are geometric event
        predictions, not a measurement of a tagged transmitted packet.
        """
        t = result['t']; states = result['states']
        if not t[0] <= emission <= t[-1]:
            return None
        clocks = states[:, 2*self.n+4:2*self.n+6]
        source_clock = float(PchipInterpolator(t, clocks[:, source])(emission))
        sigma_in = self.epsilon[source]*source_clock+self.anchor[source]
        sigma_out = sigma_in+(self.par.length if path_length is None else path_length)
        target = (sigma_out-self.anchor[destination])/self.epsilon[destination]
        if not clocks[0, destination] <= target <= clocks[-1, destination]:
            return None
        arrival = float(PchipInterpolator(clocks[:, destination], t)(target))
        return dict(emission=float(emission), arrival=arrival, source=int(source),
            destination=int(destination), neck_time_in=float(sigma_in), neck_time_out=float(sigma_out),
            elapsed_exterior=arrival-float(emission),
            time_map_orientation=int(self.epsilon[source]*self.epsilon[destination]))
