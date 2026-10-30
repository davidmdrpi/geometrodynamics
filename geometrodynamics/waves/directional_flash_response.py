"""Directional diagnostics for the existing transverse Maxwell birth flash.

These functions do not evolve mouths or insert continuous forces into BAM.
The round-sphere Maxwell propagator is a reference field, not the strict
antipodal/throat boundary solution. A sampled force or finite-window integral
must not be applied as an event kick without that boundary solution.
"""
import numpy as np

from .antipodal_throat_transport import AntipodalTransport, inertial_mouth, unit


def tangent_frame(center, first=None):
    center = unit(center)
    if first is None:
        first = np.eye(3)[np.argmin(abs(center))]
    first = np.asarray(first, float)
    first = unit(first - (first @ center)*center)
    return np.stack([first, np.cross(center, first)])


def electric(model, momentum, points):
    """Receiver-footprint averaged E_T; sum source copies before calling."""
    return -model._rot(model.footprints(np.atleast_2d(points)), momentum/model.omega)


def electric_integral(model, q0, p0, points, start, stop):
    """Exact integral at FIXED receiver points: int E dt = A(start)-A(stop).

    The interval is specified diagnostic data, not an interaction duration.
    This identity does not hold along a moving receiver without convective terms.
    """
    qa, _ = model.free_wave(q0, p0, start)
    qb, _ = model.free_wave(q0, p0, stop)
    return model._rot(model.footprints(np.atleast_2d(points)), (qa-qb)/model.omega)


def rim_geometry(center, radius, samples=128, first=None):
    """Geodesic ring and its parallel-transported tangent basis on unit S2."""
    if not 0 < radius < np.pi/2 or samples < 8:
        raise ValueError('Require 0 < radius < pi/2 and at least 8 samples')
    center = unit(center)
    frame = tangent_frame(center, first)
    phi = 2*np.pi*np.arange(samples)/samples
    radial = np.cos(phi)[:, None]*frame[0] + np.sin(phi)[:, None]*frame[1]
    points = np.cos(radius)*center + np.sin(radius)*radial
    er = -np.sin(radius)*center + np.cos(radius)*radial
    ep = -np.sin(phi)[:, None]*frame[0] + np.cos(phi)[:, None]*frame[1]
    return points, er, ep, phi, frame


def rim_modes(vectors, er, ep, phi):
    """m=0 and real m=1 components in the local radial/azimuthal basis.

    Return the mean vector transported back to the ring centre, in its frame.
    A uniform transported vector has radial AND azimuthal m=1, and zero m=0.
    This is a field projection, not a derived throat scattering matrix.
    """
    radial = np.einsum('ij,ij->i', vectors, er)
    azimuthal = np.einsum('ij,ij->i', vectors, ep)
    basis = np.stack([np.cos(phi), np.sin(phi)], axis=1)
    r1 = 2*np.mean(radial[:, None]*basis, axis=0)
    t1 = 2*np.mean(azimuthal[:, None]*basis, axis=0)
    mean = .5*np.array([r1[0]-t1[1], r1[1]+t1[0]])
    return dict(radial_m0=float(np.mean(radial)), azimuthal_m0=float(np.mean(azimuthal)),
                radial_m1=r1.tolist(), azimuthal_m1=t1.tolist(),
                transported_mean=mean.tolist())


def power_audit(model, field_momentum, points, velocities):
    """Same footprint on force and current gives exact reciprocal electric power.

    This checks the existing Maxwell coupling at an instant. It does not switch
    that coupling on at isolated events or claim a finite recoil solution.
    """
    force = model.charges[:, None]*electric(model, field_momentum, points)
    mechanical = float(np.sum(force*velocities))
    current = model.source(points, velocities).sum(axis=0)
    field = float(np.vdot(field_momentum, current).real)
    return dict(mechanical_power=mechanical, field_power=field, residual=mechanical+field)


def coasting_focus_audit(angle, speed):
    """First unperturbed flight only; stop each source at its first open capture.

    Mouths start at +/-a and coast outwards until the opposite rims touch at
    t=(pi-2a)/v. Empty foci alternate between +/- the emission point. No force
    or throat continuation is inferred. Contact itself is excluded.
    """
    if not 0 < angle < np.pi/4 or not 0 < speed < 1:
        raise ValueError('Require 0 < a < pi/4 and 0 < v < 1')
    end = (np.pi-2*angle)/speed
    mouths = [inertial_mouth(angle, speed, horizon=end, radius=angle),
              inertial_mouth(-angle, -speed, horizon=end, radius=angle)]
    transport = AntipodalTransport(mouths)
    events = []
    for source in range(2):
        origin, _ = mouths[source].frame(0.)
        for n in range(1, int(np.ceil(end/np.pi))):
            time = n*np.pi
            focus = (-1)**n*origin
            target, local = transport.mouth_at_focus(time, focus)
            events.append(dict(source=source, return_number=n, time=float(time),
                               target=target, event='empty_focus' if target is None else 'open_capture',
                               offset=None if target is None else float(np.linalg.norm(local))))
            if target is not None:
                break
    return dict(angle=angle, speed=speed, contact_time=float(end), events=events)
