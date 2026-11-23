"""Fixed-action launches from oscillator ensembles; no incoming radiation after launch."""
import json

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from .elliptic_pointer import DECAY, GAIN, INPUTS, MATCHED_WIDTH, latch


def measure(q0, v0):
    # Exact half-period propagator of the write oscillator, including momentum.
    q = (1+DECAY)*GAIN*INPUTS[:, None] - DECAY*q0
    v = np.broadcast_to(-DECAY*v0, q.shape)
    _, _, heat, energy, record = latch(q, v, np.zeros_like(q))
    positive = (record == 1).sum(axis=1)
    initial_energy = v*v/2 + (q*q-1)**2/4
    return dict(positive=positive.tolist(), no_record=(record == 0).sum(axis=1).tolist(),
                max_error=float(max(abs(positive/q.shape[1]-(1+INPUTS)/2))),
                max_latch_energy_residual=float(abs(heat+energy-initial_energy).max()))


def run(count):
    # Symmetrized scrambled Sobol points; independent coordinates under the
    # target product measure. Antithetic pairs exactly balance the symmetric input.
    if count < 2 or count & (count-1):
        raise ValueError("thermal quadrature requires a power of two")
    unit = qmc.Sobol(2, scramble=True, seed=357).random_base2(int(np.log2(count))-1)
    normal = ndtri(unit)
    normal = np.concatenate((normal, -normal), axis=0)
    phase = 2*np.pi*(np.arange(count)+0.5)/count
    rows = {}
    for factor in (0.5, 1., 2.):
        q, v = (factor*MATCHED_WIDTH*normal/np.sqrt(3)).T
        rows[f"thermal_rms_factor_{factor:g}"] = measure(q,v)
    q = MATCHED_WIDTH*np.cos(phase)
    rows["phase_at_rest"] = measure(q,np.zeros(count))
    rows["oscillator_phase"] = measure(q,-MATCHED_WIDTH*np.sin(phase))
    # Haar S² has dΩ/(4π)=d(cos θ)dφ/(4π); only this coordinate enters.
    q = MATCHED_WIDTH*(2*(np.arange(count)+0.5)/count-1)
    rows["sphere_projection_at_rest"] = measure(q,np.zeros(count))
    return rows


if __name__ == "__main__":
    print(json.dumps({n: run(n) for n in (1024, 4096)}, indent=2))
