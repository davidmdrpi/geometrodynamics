"""Mass balance on a reflection center in the charged dynamical loop."""
import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .charged_characteristics import propagate


OUTPUT = Path(__file__).resolve().parents[1]/'runs/20261126_cold_loop_mass/probe.json'
GAPS = (0.5, 0.2, 0.05, 0.01)
SAMPLES = (0.5, 0., -0.5, -1.)
PRIMARY_RADII = (.198, .197, .1965, .19625, .196, .1959)


def cold_family(epsilon):
    if not 0 < epsilon < 1:
        raise ValueError("the three reference horizons must be distinct and positive")
    b = 0.4
    a = b*(1-epsilon)
    c = 2-a-b
    roots = (a, b, c)
    surface_gravities = [math.prod(abs(s-t) for t in roots if t != s)/(2*s**1.5)
                         for s in roots]
    return dict(squared_horizons=roots, kappa=surface_gravities,
                charge=math.sqrt(a*b*c/2), mu=(a*b+a*c+b*c)/2)


def simulate(case):
    gap, n = case
    if gap is None:
        run = propagate(n, targets=PRIMARY_RADII, remap=True, audit_mass=True)
        reference = None
    else:
        reference = cold_family(gap)
        a, b, _ = reference['squared_horizons']
        targets = [math.sqrt(max(0, a+fraction*(b-a))) for fraction in SAMPLES]
        run = propagate(n, charge=reference['charge'], radius=math.sqrt(b),
                        amplitude=.02*gap, targets=targets, time_limit=12., audit_mass=True)
    return dict(epsilon=gap, reference=reference, **run)


def spatial_control():
    """Sample invariant masses across the periodic slice on uniform null grids."""
    family = cold_family(.01)
    rows = []
    for n in (512, 1024, 2048):
        run = propagate(n, charge=family['charge'], radius=math.sqrt(.4),
                        amplitude=.0002, targets=(), time_limit=2*math.pi+.02,
                        profile_times=(0., math.pi/2, math.pi, 2*math.pi))
        rows.append(dict(n=n, epsilon=.01, profiles=run['mass_profiles']))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    cases = [(None,n) for n in (8192,16384,32768)]
    cases += [(gap,n) for gap in GAPS for n in (512,1024,2048)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(simulate,cases))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(dict(runs=rows, spatial=spatial_control()),
                                 indent=1, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
