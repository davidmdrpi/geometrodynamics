"""Review diagnostics; leave the frozen solver and nine-case archives unchanged.

Default: postprocess archived CSVs only. --ports-off also runs free-field
controls to distinguish the coordinate angular balance from a rotation charge.
"""
import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.linalg import cho_solve

from geometrodynamics.waves.coupled_mouth_throat import (
    CoupledMouthThroat, Parameters, real_basis,
)

ROOT = Path('experiments/closure_ledger/runs/20261020_coupled_mouth_throat_replay')


def window_mean(t, y, start, end):
    """Integrate the piecewise-linear archive with interpolated endpoints."""
    t, y = np.asarray(t), np.asarray(y)
    if len(t) != len(y) or np.any(np.diff(t) <= 0):
        raise ValueError('Require matching samples and strictly increasing times')
    if not t[0] <= start < end <= t[-1]:
        raise ValueError('Require a complete window within the archive')
    x = np.r_[start, t[(t > start) & (t < end)], end]
    return float(np.trapezoid(np.interp(x, t, y), x)/(end-start))


def radius_diagnostics(t, ratio, period):
    """Full radial-period averages, sampled at each eligible archive time.

    This is a stated post-hoc smoothing diagnostic, not a new success gate.
    The final possible window is included; no partial-window normalization.
    """
    if period <= 0 or period > t[-1]-t[0]:
        raise ValueError('No complete positive-period window')
    starts = np.unique(np.r_[t[t <= t[-1]-period], t[-1]-period])
    means = [window_mean(t, ratio, s, s+period) for s in starts]
    windows = []
    # CSV timestamps are rounded to 12 significant digits.
    for j in range(int(np.floor((t[-1]-t[0]+1e-9)/(2*np.pi)))):
        a, b = t[0]+2*np.pi*j, min(t[-1], t[0]+2*np.pi*(j+1))
        samples = np.r_[np.interp(a, t, ratio), ratio[(t > a) & (t < b)],
                        np.interp(b, t, ratio)]
        windows.append(dict(start=float(a), end=float(b),
            mean=window_mean(t, ratio, a, b), half_range=float(np.ptp(samples)/2)))
    return dict(instantaneous_minimum=float(np.min(ratio)), radial_period=float(period),
        minimum_full_radial_period_mean=float(min(means)),
        minimum_mean_window_start=float(starts[np.argmin(means)]),
        return_windows=windows)


class FreeRotationCharge:
    """Physical sphere rotation flux on a tilted slice, with model J's sign.

    Only the free sphere is diagnosed. This is not a derived total physical
    charge for the constitutive moving-port model.
    """
    def __init__(self, model):
        self.model = model
        lmax = model.par.lmax
        z, wt = np.polynomial.legendre.leggauss(2*lmax+5)
        phi = np.arange(4*lmax+8)*2*np.pi/(4*lmax+8)
        theta, ph = np.broadcast_arrays(np.arccos(z)[:, None], phi[None, :])
        self.weights = np.broadcast_to(wt[:, None]*2*np.pi/len(phi), theta.shape).ravel()
        self.y, self.yt, self.yp = [a.reshape(-1, model.ns) for a in
                                  real_basis(lmax, theta, ph, True)]
        self.sin = np.sin(theta).ravel()
        self.cos = np.cos(theta).ravel()
        self.cp, self.sp = np.cos(ph).ravel(), np.sin(ph).ravel()

    def __call__(self, state):
        m = self.model
        u, p, _, _ = m.unpack(state)
        u, p = u[:m.ns], p[:m.ns]
        v = self.y @ cho_solve(m.factor, p-m.bs@u)
        value = self.y@u
        gt = self.yt@u + m.kappa*self.sin*v
        gp = (self.yp@u)/self.sin
        xi_grad = self.cp*gt-self.cos*self.sp*gp
        xi_tilt = -m.kappa*self.sin*self.cp
        lagrangian = .5*(v*v-gt*gt-gp*gp-.25*value*value)
        density = (v-m.kappa*self.sin*gt)*xi_grad+xi_tilt*lagrangian
        return float(-self.weights@density)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ports-off', action='store_true')
    parser.add_argument('--out', type=Path, default=Path('docs/reviews/pr346/numerical_review.json'))
    args = parser.parse_args()
    result = dict(source_commit='2073e49741d9c1cb5c1d2c3da650f56138ba3f2a',
        scope='Post-hoc diagnostics; original acceptance gates and archives unchanged.',
        radius_mean_definition='Piecewise-linear CSV integral over 2*pi/sqrt(k/mu); full windows only, starting at eligible sample times plus the last possible start.',
        archive_sha256={}, cases={})
    for path in sorted(ROOT.glob('*.csv')):
        config_path = path.with_suffix('.json')
        config = json.loads(config_path.read_text())['parameters']
        d = np.genfromtxt(path, names=True, delimiter=',')
        period = 2*np.pi/np.sqrt(config['stiffness']/config['radial_mass'])
        result['cases'][path.stem] = radius_diagnostics(d['tau'], d['radius']/config['radius0'], period)
        for source in [path, config_path]:
            result['archive_sha256'][str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
    if args.ports_off:
        result['ports_off'] = []
        for offset in [-.3, -1.2]:
            for lmax in [10, 14, 18]:
                m = CoupledMouthThroat(replace(Parameters(), coupling=0., offset=offset, lmax=lmax))
                t, states, stop = m.evolve(horizon=2*np.pi)
                charge = FreeRotationCharge(m)
                q = np.array([charge(s) for s in states])
                j = np.array([m.angular(s) for s in states])
                row = dict(offset=offset, lmax=lmax, horizon=float(t[-1]), stop=stop,
                    initial_physical_charge=float(q[0]), physical_charge_range=float(np.ptp(q)),
                    canonical_charge_range=float(np.ptp(j)),
                    coordinate_reaction_range=float(np.ptp(states[:, -1])),
                    canonical_corrected_range=float(np.ptp(j+states[:, -1])))
                result['ports_off'].append(row)
                print(json.dumps(row), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
