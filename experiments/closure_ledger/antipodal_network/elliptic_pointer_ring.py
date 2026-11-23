"""Post-freeze control using all eight sites and the nonlinear metric hopping of PR #241."""
import json

import numpy as np

from .elliptic_pointer import DAMPING, INPUTS, WRITE_TIME, integrate, latch

SITES = 8
CHI = 2*np.pi*np.arange(SITES)/SITES
QUADRUPOLE = np.cos(2*(CHI + np.pi/SITES))
WINDING = np.exp(1j*CHI)/np.sqrt(SITES)


def ring_terms(q, field):
    """Return H(q)ψ, the field energy, and its negative q derivative."""
    difference = np.roll(field, -1, axis=0) - field
    radius = 1 - 0.1*q*QUADRUPOLE[:, None]
    current = difference/radius**2
    h_field = np.roll(current, 1, axis=0) - current
    density = abs(difference)**2
    energy = np.sum(density/radius**2, axis=0)
    force = -np.sum(0.2*QUADRUPOLE[:, None]*density/radius**3, axis=0)
    return h_field, energy, force


def run(count=256):
    q0 = 0.12*(2*(np.arange(count)+0.5)/count-1)
    initial, z = np.broadcast_arrays(q0, INPUTS[:, None])
    q, z = initial.ravel(), z.ravel()
    plus, minus = np.sqrt((1+z)/2), np.sqrt((1-z)/2)
    field = (WINDING[:, None]*(plus + minus)
             + WINDING.conj()[:, None]*(plus - minus))/np.sqrt(2)
    _, field_energy0, _ = ring_terms(q, field)
    zeros = np.zeros_like(q)

    def rhs(y):
        position, speed = y[0].real, y[1].real
        h_field, _, force = ring_terms(position, y[3:])
        return np.vstack((speed, force-position-DAMPING*speed,
                          DAMPING*speed**2, -1j*h_field))

    end = integrate(rhs, WRITE_TIME, np.vstack((q, zeros, zeros, field)))
    position, speed, heat = end[:3].real
    wave = end[3:]
    _, field_energy, _ = ring_terms(position, wave)
    write_energy = (position**2+speed**2)/2 + field_energy
    residual = write_energy + heat - q**2/2 - field_energy0
    # At disconnection the free ring retains its q=0 energy. Its evolution
    # decouples, so the free pointer can be tested without integrating it again.
    _, _, heat_end, final_energy, records = latch(position, speed, heat)
    latch_residual = final_energy + heat_end - heat - speed**2/2 - (position**2-1)**2/4
    overlap = abs(WINDING.conj() @ wave)**2 + abs(WINDING @ wave)**2
    norm = np.sum(abs(wave)**2, axis=0)
    records = records.reshape(initial.shape)
    positive = (records == 1).sum(axis=1)
    return dict(count=count, inputs=INPUTS.tolist(), positive=positive.tolist(),
                negative=(records == -1).sum(axis=1).tolist(),
                no_record=(records == 0).sum(axis=1).tolist(),
                max_born_error=float(max(abs(positive/count-(1+INPUTS)/2))),
                max_leakage_at_switch=float(max(norm-overlap)),
                max_norm_error=float(max(abs(norm-1))),
                max_energy_residual=float(max(abs(residual).max(), abs(latch_residual).max())))


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
