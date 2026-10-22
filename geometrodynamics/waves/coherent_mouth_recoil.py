"""Conservative local recoil and a fixed antipodal two-port path experiment.

This is a mechanism discriminator, not a moving-mouth GR solution. The two
ports are energy-normalized scalar travelling modes. In the centred m=0
reduction their resultant is neck-normal: the tangential force is ZERO by
azimuthal symmetry, not merely uncomputed. In particular the optional Klein
compression sign does not implement the tube's angular-momentum map.
"""
from dataclasses import dataclass

import numpy as np


def size_impedance(neck_over_mouth, spatial_dimension):
    """Equal-speed one-mode admittance: rim length in 2D, area in 3D."""
    ratio = float(neck_over_mouth)
    if not np.isfinite(ratio) or ratio <= 0 or spatial_dimension not in (2, 3):
        raise ValueError('Positive radius ratio and dimension 2 or 3 required')
    return ratio**(spatial_dimension-1)


def electric_stress_cross(static_field, wave_field):
    """Electric Maxwell-stress cross term in units epsilon=1.

    T(E0+e)-T(E0)-T(e). This local field identity does not supply the global
    static field or its dynamical exchange with the BAM throat.
    """
    e0, e = np.asarray(static_field, float), np.asarray(wave_field, float)
    if e0.shape != e.shape or e.shape not in ((2,), (3,)) or not np.isfinite([e0,e]).all():
        raise ValueError('Finite matching two- or three-component fields required')
    return np.outer(e0,e)+np.outer(e,e0)-np.eye(len(e))*np.dot(e0,e)


def junction(impedance_ratio):
    """Lossless reciprocal junction from continuity of value and flux.

    Z = throat / surface characteristic admittance; equal wave speed and a
    single matched transverse mode are assumptions. Z=(neck/mouth radius)^2
    only in the three-spatial-dimensional area model used in the probe.
    """
    z = float(impedance_ratio)
    if not np.isfinite(z) or z <= 0:
        raise ValueError('Positive finite impedance ratio required')
    r = (1-z)/(1+z)
    t = 2*np.sqrt(z)/(1+z)
    return np.array([[r, t], [t, -r]])


def amplitudes(value):
    a = np.asarray(value, dtype=complex)
    if a.shape != (2,) or not np.all(np.isfinite(a)):
        raise ValueError('Two finite incident amplitudes required')
    return a


def scatter(incident, impedance_ratio):
    """Fixed-mouth stress impulse, in c=1 units, for one matched packet.

    Incident directions are (+1,-1); outgoing directions are (-1,+1).
    Thus dp = (Es-Et)+(Es_out-Et_out), with NO signed-amplitude force.
    """
    a = amplitudes(incident)
    b = junction(impedance_ratio) @ a
    ein = np.abs(a)**2
    eout = np.abs(b)**2
    return dict(outgoing=b, energy=float(ein.sum()),
                impulse=float(ein[0]-ein[1]+eout[0]-eout[1]),
                energy_residual=float(eout.sum()-ein.sum()))


def recoil_operator(impedance_ratio):
    s = junction(impedance_ratio)
    j = np.diag([1., -1.])
    return j+s.T@j@s


def translated_junction(position, omega, impedance_ratio):
    """Translation phases: S(x)=exp(i*k*J*x) S(0) exp(i*k*J*x).

    With c=1, (-i/omega) S^dagger dS/dx at x=0 equals recoil_operator.
    This is local longitudinal virtual work, not a global mouth separation.
    """
    if not np.isfinite([position, omega]).all():
        raise ValueError('Finite position and frequency required')
    d = np.diag(np.exp(1j*omega*position*np.array([1., -1.])))
    return d@junction(impedance_ratio)@d


def finite_mass_scatter(incident, impedance_ratio, mass=1., momentum=0.):
    """One encounter with exact 1+1-dimensional energy/momentum closure.

    Boost the inputs to the initial mouth rest frame, scatter with the stated
    junction, and rescale all outgoing packet energies by one positive factor
    to put the recoiling mouth back on its mass shell. Then boost back.
    The common rescaling is an explicit constitutive completion, not a derived
    moving-boundary spectrum. It cannot establish inter-event phase coherence.
    """
    a = amplitudes(incident)
    if not np.isfinite([mass, momentum]).all() or mass <= 0:
        raise ValueError('Positive finite mass and finite momentum required')
    em = float(np.hypot(mass, momentum))
    gamma, v = em/mass, momentum/em
    din = np.array([1., -1.]); dout = -din
    ar = a*np.sqrt(gamma*(1-v*din))
    br = junction(impedance_ratio)@ar
    energy = float(np.vdot(ar, ar).real)
    pin = float(np.dot(din, np.abs(ar)**2))
    pout = float(np.dot(dout, np.abs(br)**2))
    if energy == 0:
        scale = 1.
    else:
        # (M+E-kE)^2 - (pin-k*pout)^2 = M^2.
        aa = (mass+energy)*energy-pin*pout
        bb = max(0., energy*energy-pout*pout)
        cc = 2*mass*energy+energy*energy-pin*pin
        disc = aa*aa-bb*cc
        if disc < -1e-12*aa*aa:
            raise ValueError('No elastic mass-shell completion')
        scale = cc/(aa+np.sqrt(max(0., disc)))
    pr = pin-scale*pout
    er = float(np.hypot(mass, pr))
    final_p = gamma*(pr+v*er)
    final_e = gamma*(er+v*pr)
    bout = br*np.sqrt(scale*gamma*(1+v*dout))
    eout = np.abs(bout)**2
    ein = np.abs(a)**2
    return dict(outgoing=bout, momentum=float(final_p), mouth_energy=float(final_e),
                impulse=float(final_p-momentum), rest_energy_scale=float(scale),
                energy_residual=float(final_e+eout.sum()-em-ein.sum()),
                momentum_residual=float(final_p+np.dot(dout, eout)-momentum-np.dot(din, ein)))


def focus_overlap(radius_over_width):
    """Captured energy of a normalized circular Gaussian focus (2D section).

    This is an aperture overlap model, not a diffraction law or an S3 tensor
    focusing calculation. Unmatched energy remains in the bulk.
    """
    x = np.asarray(radius_over_width, dtype=float)
    if np.any(~np.isfinite(x)) or np.any(x < 0):
        raise ValueError('Nonnegative finite focus ratio required')
    return -np.expm1(-x*x/2)


@dataclass(frozen=True)
class PathParameters:
    surface_bins: int = 24
    throat_bins: int = 6
    impedance_ratio: float = (np.sqrt(2)-1)**2
    handle: str = 'klein'
    capture_phase: str = 'after_caustic'
    seed_phase: float = np.pi
    throat_seed_fraction: float = 0.
    throat_seed_phase: float = np.pi
    pulse_bins: int = 1
    periods: int = 48


def closed_paths(parameters=PathParameters(), coherent=True, sample=False):
    """Finite seed, closed passive delay network; no pump or repeated source.

    Two fixed antipodal nodes have a surface leg of duration pi and a throat
    leg of L=pi*Nt/Ns. Each local scattering is unitary. A surface traversal
    carries -i for post-caustic capture (1 for pre-caustic capture). A throat
    traversal optionally carries the ASSUMED compression holonomy -1 (Klein).
    Empty-focus geometry and moving aperture hits are handled by #347, not by
    this aligned-node reduction. Delays are exact integer multiples of dt.

    Bins are complex envelopes of a positive-frequency carrier, NOT real
    broadband flashes. For a real pulse the caustic is a Hilbert transform;
    its reshaping and the resulting overlap are absent from this envelope map.
    The recorded impulses are neck-normal stresses on fixed centred m=0 nodes;
    their tangential resultant is zero. Supports absorb their normal impulse.
    The incoherent control transports intensities with |S|^2, also conserving
    the initial energy. It is an ensemble dephasing control, not absorption.
    """
    p = parameters
    if (p.handle not in ('klein', 'torus')
            or p.capture_phase not in ('after_caustic', 'before_caustic')
            or any(not isinstance(v, (int, np.integer)) or v < 1
                   for v in (p.surface_bins, p.throat_bins, p.pulse_bins, p.periods))
            or p.pulse_bins > min(p.surface_bins, p.throat_bins)
            or not np.isfinite([p.seed_phase, p.throat_seed_phase, p.throat_seed_fraction]).all()
            or not 0 <= p.throat_seed_fraction <= 1):
        raise ValueError('Invalid path parameters')
    s = junction(p.impedance_ratio)
    buffers = [np.zeros((n, 2), dtype=complex if coherent else float)
               for n in (p.surface_bins, p.throat_bins)]
    # Finite initial inventory, optionally split between both return channels.
    envelope = np.exp(-.5*((np.arange(p.pulse_bins)-(p.pulse_bins-1)/2)
                           /max(1., p.pulse_bins/3))**2)
    envelope /= np.linalg.norm(envelope)*np.sqrt(2)
    seed = envelope[:, None]*np.array([1., np.exp(1j*p.seed_phase)])
    surface_seed = seed*np.sqrt(1-p.throat_seed_fraction)
    throat_seed = seed*np.sqrt(p.throat_seed_fraction)*np.exp(1j*p.throat_seed_phase)
    for buf, value in zip(buffers, (surface_seed, throat_seed)):
        buf[:p.pulse_bins] = value if coherent else np.abs(value)**2
    phases = np.array([-1j if p.capture_phase == 'after_caustic' else 1.,
                       -1. if p.handle == 'klein' else 1.])
    impulse = np.zeros(2); exposure = np.zeros(2); interference = np.zeros(2)
    max_energy_error = 0.; rows = []
    for step in range(p.periods*p.surface_bins):
        idx = [step % len(b) for b in buffers]
        incoming = np.stack([b[i].copy() for b, i in zip(buffers, idx)], axis=1)
        if coherent:
            outgoing = incoming@s.T
            ein, eout = np.abs(incoming)**2, np.abs(outgoing)**2
            cross = 4*s[0, 0]*s[0, 1]*np.real(incoming[:, 0].conj()*incoming[:, 1])
        else:
            outgoing = incoming@(s*s).T
            ein, eout = incoming, outgoing
            cross = np.zeros(2)
        dp = ein[:, 0]-ein[:, 1]+eout[:, 0]-eout[:, 1]
        impulse += dp; interference += cross; exposure += ein.sum(axis=1)
        for k, (b, i) in enumerate(zip(buffers, idx)):
            # Opposite node is reached only after the full edge delay.
            b[i] = outgoing[::-1, k]*(phases[k] if coherent else 1.)
        inventory = sum(float(np.sum(np.abs(b)**2) if coherent else np.sum(b)) for b in buffers)
        max_energy_error = max(max_energy_error, abs(inventory-1.))
        if sample and (np.any(ein) or step == 0):
            rows.append([float(step*np.pi/p.surface_bins), *impulse.tolist(), inventory])
    return dict(impulse=impulse.tolist(), exposure=exposure.tolist(),
                interference_impulse=interference.tolist(),
                impulse_per_exposure=(impulse/np.maximum(exposure, 1e-300)).tolist(),
                initial_energy=1., final_energy=inventory,
                max_energy_error=max_energy_error, duration=p.periods*np.pi,
                sample=rows)


def path_sum(omega, impedance_ratio, length, survival=.9, handle='klein',
             capture_phase='after_caustic', source=(1., 0., 0., 0.), terms=256,
             clock_offset=0.):
    """Frequency-domain sum of ALL walks; externally fed diagnostic only.

    a = b + survival*P*S*a. survival<1 accounts for escape into untracked modes,
    whose power is reported. It is NOT a closed steady-state free-energy source.
    The finite-seed closed_paths experiment is the no-drive check. A unitary
    closed loop at survival=1 needs a homogeneous eigenmode/boundary solution.
    The default source is one-sided. The former antisymmetric surface-only
    source makes impulse/exposure blind to the surface phase at zero offset;
    it is retained in the probe as a degeneracy control, not the main source.
    sigma_A=t, sigma_B=t+offset give directed throat delays L-offset and
    L+offset. Negative delays are permitted ONLY in this frequency-domain
    boundary diagnostic, never inserted as events in closed_paths. Reversing
    the orientation of time itself would mix frequency sectors and is absent.
    """
    if (not 0 <= survival < 1 or handle not in ('klein', 'torus')
            or capture_phase not in ('before_caustic', 'after_caustic')
            or not np.isfinite([omega, length, clock_offset]).all() or length <= 0 or terms < 1):
        raise ValueError('Invalid path sum; closed poles require a boundary solution')
    s = np.kron(np.eye(2), junction(impedance_ratio))
    phase = np.array([(-1j if capture_phase == 'after_caustic' else 1)*np.exp(1j*omega*np.pi),
                      (-1 if handle == 'klein' else 1)*np.exp(1j*omega*length)])
    prop = np.zeros((4, 4), dtype=complex)
    prop[0, 2] = prop[2, 0] = phase[0]
    prop[1, 3] = phase[1]*np.exp(1j*omega*clock_offset)  # B -> A
    prop[3, 1] = phase[1]*np.exp(-1j*omega*clock_offset) # A -> B
    u = prop@s
    b = np.asarray(source, dtype=complex)
    if b.shape != (4,) or not np.isfinite(b).all() or np.linalg.norm(b) == 0:
        raise ValueError('Nonzero finite four-port source required')
    a = np.linalg.solve(np.eye(4)-survival*u, b)
    total = b.copy(); term = b.copy()
    for _ in range(1, terms):
        term = survival*u@term; total += term
    outgoing = s@a
    energy = float(np.vdot(a, a).real)
    dp = [scatter(a[i:i+2], impedance_ratio)['impulse'] for i in (0, 2)]
    return dict(incident=a, impulse=dp, exposure=energy,
                source_energy=float(np.vdot(b,b).real),
                impulse_per_source_energy=float(sum(dp)/np.vdot(b,b).real),
                impulse_per_exposure=float(sum(dp)/energy),
                walk_sum_relative_error=float(np.linalg.norm(total-a)/np.linalg.norm(a)),
                escape_power=float((1-survival**2)*np.vdot(outgoing, outgoing).real),
                throat_delays=[float(length-clock_offset), float(length+clock_offset)],
                advanced_boundary=bool(abs(clock_offset)>length),
                unitary_error=float(np.linalg.norm(u.conj().T@u-np.eye(4))))
