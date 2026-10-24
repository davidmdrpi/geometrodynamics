"""Positive-energy deformable-defect benchmark, NOT the BAM interaction law.

L = 1/2 sum_a (u_a,t^2-u_a,x^2) - V,
V = (phi^2-1)^2/2 + m2*(psi1^2+psi2^2)/2
    + g*phi*(psi1^2-psi2^2)/2 + lam*(psi1^2+psi2^2)^2/4.

The field phi deforms and translates by its own PDE; there are no assigned
particle forces. The spatial deck map (phi,psi1,psi2)->(-phi,psi2,psi1)
preserves the action. It is a two-channel analogue, not a GR Alice metric.
The local continuous coupling is intentionally isolated from antipodal routing.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Parameters:
    length: float = 200.
    dx: float = .1
    mass_squared: float = 4.
    coupling: float = 1.
    quartic: float = 1.
    twisted: bool = False


class DeformableDefect:
    def __init__(self, parameters=Parameters()):
        self.par = parameters
        p = parameters
        if (not np.isfinite([p.length,p.dx,p.mass_squared,p.coupling,p.quartic]).all()
                or min(p.length,p.dx,p.quartic) <= 0 or p.mass_squared <= abs(p.coupling)):
            raise ValueError('Positive scales and propagating stable vacuum masses required')
        self.n = int(round(p.length/p.dx))
        if self.n < 32 or abs(self.n*p.dx-p.length) > 1e-10:
            raise ValueError('Length must be an integer multiple of dx, with >=32 cells')
        self.x = (np.arange(self.n)-self.n/2)*p.dx

    @staticmethod
    def deck(u):
        return np.array([-u[0],u[2],u[1]])

    def neighbors(self, u):
        left, right = np.roll(u,1,axis=-1), np.roll(u,-1,axis=-1)
        if self.par.twisted:
            left[:,0] = self.deck(u[:,-1])
            right[:,-1] = self.deck(u[:,0])
        return left, right

    def potential(self, u):
        phi, a, b = u
        r2 = a*a+b*b
        p = self.par
        return .5*(phi*phi-1)**2 + .5*p.mass_squared*r2 + .5*p.coupling*phi*(a*a-b*b) + .25*p.quartic*r2*r2

    def acceleration(self, u):
        left, right = self.neighbors(u)
        phi, a, b = u
        p = self.par
        r2 = a*a+b*b
        dv = np.array([2*phi*(phi*phi-1)+.5*p.coupling*(a*a-b*b),
                       (p.mass_squared+p.coupling*phi+p.quartic*r2)*a,
                       (p.mass_squared-p.coupling*phi+p.quartic*r2)*b])
        return (left-2*u+right)/p.dx**2-dv

    def step(self, u, v, dt):
        # Velocity Verlet; exact time reversal of the discrete update.
        vhalf = v+.5*dt*self.acceleration(u)
        unew = u+dt*vhalf
        return unew, vhalf+.5*dt*self.acceleration(unew)

    def energy(self, u, v):
        _, right = self.neighbors(u)
        return float(self.par.dx*np.sum(.5*np.sum(v*v+(right-u)**2/self.par.dx**2,axis=0)+self.potential(u)))

    def momenta(self, u, v):
        left, right = self.neighbors(u)
        return -np.sum(v*(right-left)/2,axis=1)

    def centers(self, u):
        # Restrict to the central half so boundary chart changes are not defects.
        phi = u[0]
        ids = np.where(phi[:-1]*phi[1:] < 0)[0]
        exact = np.where(phi==0)[0]
        roots = list(self.x[exact])
        roots += list(self.x[ids]-self.par.dx*phi[ids]/(phi[ids+1]-phi[ids]))
        return sorted(x for x in roots if abs(x)<self.par.length/4)

    def initial(self, pair=False, separation=14., speed=0.):
        x = self.x
        if pair:
            if self.par.twisted:
                raise ValueError('A kink-antikink pair uses the untwisted sector')
            gamma = 1/np.sqrt(1-speed*speed)
            uleft, uright = gamma*(x+separation/2),gamma*(x-separation/2)
            phi = np.tanh(uleft)-np.tanh(uright)-1
            velocity = speed*gamma*(1/np.cosh(np.clip(uleft,-350,350))**2
                                    +1/np.cosh(np.clip(uright,-350,350))**2)
        else:
            if not self.par.twisted:
                raise ValueError('One kink on the circle requires the twisted sector')
            phi = np.tanh(x)
            velocity = np.zeros_like(x)
        u = np.array([phi,np.zeros_like(x),np.zeros_like(x)])
        v = np.array([velocity,np.zeros_like(x),np.zeros_like(x)])
        return u, v

    def add_packet(self, u, v, center=25., width=5., omega=2.8, amplitude=.04,
                   direction=-1, channel=1):
        """Finite real packet; spectral velocity in its asymptotic mass branch.

        Uses the discrete spatial dispersion for initial data. The packet must
        be well separated from defects and the chart seam at preparation.
        """
        if channel not in (1,2) or direction not in (-1,1) or width<=0:
            raise ValueError('Invalid packet channel, direction or width')
        vacuum = float(u[0,np.argmin(abs(self.x-center))])
        mass2 = self.par.mass_squared+self.par.coupling*vacuum*(1 if channel==1 else -1)
        k2 = omega*omega-mass2
        if k2<=0 or .5*self.par.dx*np.sqrt(k2)>=1:
            raise ValueError('Packet frequency outside propagating grid band')
        k = 2/self.par.dx*np.arcsin(.5*self.par.dx*np.sqrt(k2))
        packet = amplitude*np.exp(-.5*((self.x-center)/width)**2)*np.cos(k*(self.x-center))
        wavenumber = 2*np.pi*np.fft.fftfreq(self.n,d=self.par.dx)
        frequency = np.sqrt(mass2+(2*np.sin(wavenumber*self.par.dx/2)/self.par.dx)**2)
        velocity = np.fft.ifft(-1j*direction*np.sign(wavenumber)*frequency*np.fft.fft(packet)).real
        u[channel] += packet
        v[channel] += velocity


def step_scattering(omega, mass_squared=4., coupling=1., extent=18.):
    """Independent ODE scattering on phi=tanh(x), incident from the right.

    Unit outgoing amplitude at the left; match both value and derivative on
    the right. Return energy fractions and +x impulse per incident energy.
    The analytic tanh-step reflection is a separate validation target.
    """
    if omega*omega <= mass_squared+abs(coupling):
        raise ValueError('Both asymptotic channels must propagate')
    kl = np.sqrt(omega*omega-mass_squared+coupling)
    kr = np.sqrt(omega*omega-mass_squared-coupling)
    u0 = np.exp(1j*kl*extent)
    result = solve_ivp(lambda x,y:[y[1],-(omega*omega-mass_squared-coupling*np.tanh(x))*y[0]],
                       [-extent,extent],[u0,-1j*kl*u0],rtol=2e-11,atol=2e-12)
    value, derivative = result.y[:,-1]
    incident = .5*(value+1j*derivative/kr)*np.exp(1j*kr*extent)
    reflected = .5*(value-1j*derivative/kr)*np.exp(-1j*kr*extent)
    reflection = abs(reflected/incident)**2
    transmission = kl/kr/abs(incident)**2
    analytic = (np.sinh(np.pi*(kl-kr)/2)/np.sinh(np.pi*(kl+kr)/2))**2
    return dict(reflection=float(reflection),transmission=float(transmission),
                flux_residual=float(reflection+transmission-1),
                analytic_reflection=float(analytic),
                impulse_per_energy=float((transmission*kl-(1+reflection)*kr)/omega),
                k_left=float(kl),k_right=float(kr))


def evolve(model, u, v, end, dt, sample_dt=.25):
    if dt<=0 or dt>model.par.dx*.5 or end<0:
        raise ValueError('Positive dt <= dx/2 and nonnegative duration required')
    count=int(round(end/dt)); stride=max(1,int(round(sample_dt/dt)))
    if abs(count*dt-end)>1e-9:
        raise ValueError('Duration must be an integer number of steps')
    e0=model.energy(u,v); p0=model.momenta(u,v)
    rows=[]
    for step in range(count+1):
        if step%stride==0 or step==count:
            momentum=model.momenta(u,v)
            rows.append(dict(time=step*dt,energy=model.energy(u,v),energy_error=model.energy(u,v)-e0,
                             momenta=momentum.tolist(),momentum_error=float(np.sum(momentum-p0)),
                             centers=model.centers(u)))
        if step<count:
            u,v=model.step(u,v,dt)
    return u,v,rows
