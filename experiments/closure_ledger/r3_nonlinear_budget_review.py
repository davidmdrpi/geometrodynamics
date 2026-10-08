"""Post-hoc archive diagnostics for the review of PR #321; no new evolution.

Keep registered labels intact. Independently validate the completed return
inside each failed two-return step, which the original scorer skips.
"""
import json
import numpy as np
from scipy.interpolate import PPoly
from experiments.closure_ledger import r3_nonlinear_budget_probe as p
from experiments.closure_ledger import r3_nonlinear_budget_replay as replay


def scale_range(family):
    """Global range of A on the reference spline, including interior extrema."""
    curve = PPoly(family.curve.c[:, :, 0], family.curve.x, extrapolate=False)
    roots = curve.derivative().roots(extrapolate=False)
    values = curve(np.r_[curve.x, roots[np.isfinite(roots)]])
    return float(min(values)), float(max(values))


def partial_return(family, run, d0):
    if run['arm'] != 'unforced' or run['terminal'] != 'NUMERICALLY_UNRESOLVED':
        raise ValueError('not an unresolved unforced run')
    if len(run['steps']) != 2:
        raise ValueError('wrong step count')
    step = run['steps'][1]
    if step['step'] != 2 or 'error' not in step or len(step['returns']) != 1:
        raise ValueError('expected one completed return in failed step')
    start = p.arr(step['start'], (29,))
    p.close(start, p.arr(run['steps'][0]['after'], (29,)))
    history = step['returns'][0]
    times = p.arr(history['times'], (33,))
    states = p.arr(history['states'], (33, 29))
    if times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError('invalid return times')
    p.close(states[0], start)
    p.close(states[:, 2], start[2]+times)
    residual = 0.
    for row in states:
        p.d.ingredients(row)
        residual = max(residual, float(np.max(abs(p.d.constraints(row)['residual']))))
    if residual > 1e-8:
        raise ValueError('partial-return constraint accuracy')
    end = states[-1]
    if abs(end[3]) > 1e-9 or end[7] >= 0:
        raise ValueError('not a descending section')
    low, high = scale_range(family)
    bound = max(low-end[0], end[0]-high, 0.)/d0
    return dict(method=run['method'], A=float(end[0]),
                scale_distance_lower_bound_over_d0=float(bound),
                max_constraint=residual, endpoint=end.tolist())


def diagnostics(directory=p.RUN):
    registered = replay.replay(directory)  # Authenticate bytes, sources and decisions.
    family = p.family()
    cases = []
    for k, summary in enumerate(registered['cases']):
        case = p.read(directory/f'case_{k}.json.gz.b64')
        d0 = case['initial']['d0']
        partial = [partial_return(family, r, d0) for r in case['runs'] if r['arm']=='unforced']
        difference = float(np.max(abs(np.array(partial[0]['endpoint'])-partial[1]['endpoint'])))
        if difference > 1e-4*d0:
            raise ValueError('partial-return integrator disagreement')
        costs = [r['steps'][0]['proposed_cost'] for r in case['runs'] if r['arm']=='controlled']
        cases.append(dict(detuning=case['detuning'],
                          cost_over_d0_squared=[c/d0**2 for c in costs],
                          partial_returns=partial, endpoint_max_difference=difference,
                          partial_section_outside_tube=bool(all(r['scale_distance_lower_bound_over_d0']>2 for r in partial)),
                          registered_unforced=summary['unforced']))
    return dict(scope='post-hoc diagnostics; registered labels unchanged',
                registered_label=registered['label'], reference_A_range=scale_range(family),
                minimum_detuning_from_horizon_cap=1/(2*abs(family.twist)*256), cases=cases)


if __name__ == '__main__':
    print(json.dumps(diagnostics(), indent=2, allow_nan=False))
