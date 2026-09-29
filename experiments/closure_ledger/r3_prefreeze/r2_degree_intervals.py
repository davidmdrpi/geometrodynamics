"""Post-review static-map diagnostic; not an Einstein evolution or R2 verdict.

Original pre-freeze bytes live at 2e984ac. Root search is not exhaustive;
degree quadrature needs refinement and is never rounded to an integer here.
"""
import argparse
import numpy as np
from experiments.closure_ledger.r3_prefreeze.r2_degree_scratch import (
    grid, degree, dphi, hopf, tangential_zeros,
)
from experiments.closure_ledger.r3_review_controls import exact_event_time


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--grid', type=int, default=72)
    ap.add_argument('--seeds', type=int, default=4000)
    args = ap.parse_args()
    if args.grid < 4 or args.seeds < 1:
        ap.error('grid >= 4 and seeds >= 1 required')
    X, h = grid(args.grid)
    Z = tangential_zeros(dphi, seeds=args.seeds)
    paired = all(any(np.linalg.norm(x+y) < 1e-6 for y, _ in Z) for x, _ in Z)
    print('KINEMATIC_ONLY; root-search completeness NOT_ESTABLISHED')
    print('roots found:', len(Z), 'antipodal pairing of found roots:', paired)
    for eps in (.001, .01, .05):
        times = sorted(float(exact_event_time(eps, lam)) for _, lam in Z)
        print('eps:', eps, 'exact central-branch tau/eps:', np.array(times)/eps)
    eps = .01
    # Group approximate antipodal duplicates only for choosing diagnostic times.
    # The tolerance is not an event-count or integer-degree certification.
    events = np.unique(np.round([exact_event_time(eps, lam) for _, lam in Z], 10))
    events = events[(events > -.05) & (events < .05)]
    edges = np.r_[-.05, events, .05]
    D = dphi(X)
    print('oriented identity degree control:', -degree(X, h))
    for lo, hi in zip(edges, edges[1:]):
        t = (lo+hi)/2
        R = -np.sqrt(3)/2*np.sin(2*t)
        F = R*X+eps*D
        print('tau:', t, 'degree estimate:', -degree(F, h),
              'sampled min|phi|:', np.linalg.norm(F, axis=-1).min())
    for t in (-.02, 0., .02):
        R = -np.sqrt(3)/2*np.sin(2*t)
        F = R*X+eps*hopf(X)
        print('Hopf tau:', t, 'degree estimate:', -degree(F, h),
              'analytic min|phi|:', np.hypot(R, eps))


if __name__ == '__main__':
    main()
