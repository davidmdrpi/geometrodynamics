"""Equation-level preflight of the submitted three-law null-horizon protocol.

This is NOT a spacetime packet evolution. Trace quadratures test coordinate
covariance and orientation; no observer waveform or backreaction is computed.
The original submitted protocol is preserved byte-for-byte in docs/protocols.
See docs/null_horizon_continuation_audit.md for derivations and conventions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import quad


def tortoise(r, b=1.0):
    if b <= 0 or r <= b:
        raise ValueError("tortoise requires b > 0 and r > b")
    return r + b / 2 * math.log((r - b) / (r + b))


def kruskal_product(r, b=1.0):
    return -b*b * (r-b)/(r+b) * math.exp(2*r/b)


def metric_factor(r, b=1.0):
    """Dimensionless F for dimensionful U=-b exp(-u/b), V=b exp(v/b)."""
    return ((r+b)/r)**2 * math.exp(-2*r/b)


def wave_coefficient(r, ell, b=1.0):
    """Q in u_UV + Q u = 0, evaluated without division by f at the horizon."""
    f = 1 - b*b/(r*r)
    bracket = ell*(ell+2)/r**2 + 3*f/(4*r**2) + 3*b*b/r**4
    return metric_factor(r, b)*bracket/4


def junction(V, g, gp, gpp=0.0, b=1.0):
    """H+ -> exterior H- in common parameter V; return formal BI quantities.

Scaled stress means 8*pi*G5 times mu or p. A negative slope cannot join the
specified retained exterior sides: its formal jumps are NOT physical stress.
"""
    if V <= 0 or g >= 0 or gp == 0 or b <= 0:
        raise ValueError("require V>0, g<0, gprime!=0, b>0")
    minus = -V/b
    plus = -g/(b*gp)
    return dict(nr_minus=minus, nr_plus=plus, jump_nr=plus-minus,
                formal_mu_8piG=-3*(plus-minus)/b,
                formal_p_8piG=-gpp/gp,
                future_generator_preserved=bool(gp > 0),
                retained_exterior_sides_compatible=bool(gp > 0),
                continuous_metric_transverse_slope=1/gp,
                physical_stress_interpretation=bool(gp > 0))


def packet(V, b=1.0, amplitude=0.01):
    """Submitted Gaussian: finite energy, NOT compact support."""
    z = (V-2*b)/(.25*b)
    a = amplitude*math.exp(-z*z/2)
    return a, -z/(.25*b)*a


def trace_control(kind, ell=0, boost=1.0, b=1.0, alpha=1.3, c=None):
    """Direct quadrature in transformed coordinates; no wave propagation.

    Horizon angular harmonics have unit integral Y^2. Integrals cover V>0.
    A second real variation q=(V-2b)*a/b gives a nonzero symplectic pairing.
    Passive boosts transform both the incoming data and the law parameters.
    """
    if boost <= 0 or b <= 0 or ell < 0 or int(ell) != ell:
        raise ValueError("require boost,b>0 and nonnegative integer ell")
    c = b*b if c is None else c
    if alpha <= 0 or c <= 0 or kind not in ("linear", "reciprocal"):
        raise ValueError("unknown map or nonpositive parameter")
    parity = (-1)**ell
    alpha_tilde = alpha/boost**2

    def quantities(w):
        V = w/boost
        a, ap = packet(V, b)
        aw = ap/boost
        q = (V-2*b)*a/b
        qw = (a/b+(V-2*b)*ap/b)/boost
        if kind == "linear":
            U, gp = -alpha_tilde*w, -alpha_tilde
        else:
            U, gp = -c/w, c/w**2
        out, outp = parity*a, parity*aw/gp
        qout, qoutp = parity*q, parity*qw/gp
        jac = abs(gp)  # increasing U orientation, including reversed maps
        return (aw*aw, outp*outp*jac,
                w/b*aw*aw, -U/b*outp*outp*jac,
                a*qw-q*aw, (out*qoutp-qout*outp)*jac)

    # Split around the narrow Gaussian so adaptive integration cannot miss it.
    edges = [0., boost*b, 2*boost*b, 3*boost*b, 4*boost*b, math.inf]
    vals = [sum(quad(lambda w: quantities(w)[i], lo, hi,
                     epsabs=1e-13, epsrel=1e-11)[0]
                for lo, hi in zip(edges[:-1], edges[1:])) for i in range(6)]
    return dict(kind=kind, ell=ell, boost=boost,
                alpha_in_this_chart=alpha_tilde if kind == "linear" else None,
                raw_in=vals[0], raw_out=vals[1], killing_in=vals[2],
                killing_out=vals[3], symplectic_in=vals[4], symplectic_out=vals[5])


def run():
    root = Path(__file__).resolve().parents[2]
    geometry = []
    for b in (0.5, 1., 2.):
        r = 1.7*b
        f = 1-b*b/r**2
        F_from_jacobian = -f*b*b/kruskal_product(r, b)
        geometry.append(dict(b=b, F_horizon=metric_factor(b, b),
            metric_jacobian_residual=abs(metric_factor(r,b)-F_from_jacobian),
            submitted_first_F_ratio=4., submitted_second_F_ratio=4*b*b,
            horizon_Q=[wave_coefficient(b, l, b) for l in (0,1)],
            horizon_K=72/b**4))
    linear = [dict(V=V, **junction(V, -1.3*V, -1.3)) for V in (1.,2.,3.)]
    # Exploratory control, explicitly NOT substituted for frozen Law B.
    reciprocal = [dict(V=V, **junction(V, -1/V, 1/V**2, -2/V**3))
                  for V in (1.,2.,3.)]
    traces = [trace_control(kind,l,boost) for kind in ("linear","reciprocal")
              for l in (0,1) for boost in (.5,1.,2.,4.)]
    delays = [dict(epsilon=eps,
                   earliest_round_trip_Killing_time=2*(tortoise(2.)-tortoise(1.+eps)),
                   earliest_round_trip_observer_proper_time=math.sqrt(.75)*2*(tortoise(2.)-tortoise(1.+eps)))
              for eps in (.1,.01,.001,.0001)]
    boosts = [dict(boost=L,
                   coordinate_time_origin_shift=math.log(L),
                   linear_map_relative_error_if_alpha_not_transformed=abs(L*L-1))
              for L in (.5,1.,2.,4.)]
    errors = []
    for kind in ("linear","reciprocal"):
        group = [x for x in traces if x['kind']==kind]
        base = next(x for x in group if x['ell']==0 and x['boost']==1.)
        errors += [abs(x['killing_out']/base['killing_out']-1) for x in group]
    checks = dict(
        metric_jacobian=max(x['metric_jacobian_residual'] for x in geometry)<1e-14,
        horizon_wave_coefficient_nonzero=all(min(x['horizon_Q'])>0 for x in geometry),
        linear_formal_jump_zero=all(abs(x['jump_nr'])<1e-14 for x in linear),
        linear_sides_incompatible=all(not x['retained_exterior_sides_compatible'] for x in linear),
        reciprocal_sides_compatible=all(x['retained_exterior_sides_compatible'] for x in reciprocal),
        reciprocal_negative_mu=all(x['formal_mu_8piG']<0 for x in reciprocal),
        trace_killing_energy_covariant=max(errors)<1e-10,
        trace_killing_energy_balance=max(abs(x['killing_out']/x['killing_in']-1) for x in traces)<1e-10,
        trace_symplectic_orientation=max(abs(x['symplectic_out']/x['symplectic_in']-
            (-1 if x['kind']=='linear' else 1)) for x in traces)<1e-10,
        cutoff_delay_increasing=all(y['earliest_round_trip_Killing_time']>x['earliest_round_trip_Killing_time']
                                   for x,y in zip(delays[:-1],delays[1:])))
    paths = ['docs/protocols/null_horizon_continuation_submitted.txt',
             'experiments/throat_dynamics/null_horizon_continuation_probe.py',
             'tests/test_null_horizon_continuation.py',
             'docs/null_horizon_continuation_audit.md']
    return dict(label='SUBMITTED_PROTOCOL_REQUIRES_CORRECTION' if all(checks.values())
                else 'AUDIT_NUMERICALLY_UNRESOLVED',
        experiment_kind='analytic preflight and characteristic-trace quadrature; no PDE evolution',
        original_protocol_preserved=True, checks=checks, geometry=geometry,
        laws=dict(
            A=dict(name='kruskal_analytic_quotient',
                G1='FAIL_RETURN_REQUIREMENT: local lift enters interior; quotient reverses time orientation',
                G2='LOCAL_SMOOTH_QUOTIENT: no junction shell; r=0 excluded, not regularized',
                G3='NOT_RUN: no retarded return law'),
            B=dict(name='direct_null_bounce',
                G1='FAIL: gprime<0 reverses horizon generator ordering',
                G2='INVALID_JUNCTION: retained exterior sides incompatible; formal zero jump is not admission',
                G3='NOT_RUN: only formal trace covariance and orientation checked'),
            C=dict(name='stretched_horizon_cutoff',
                G1='EVEN_BC_NOT_STATIC_NORMAL: corrected wall has finite-epsilon return with divergent delay',
                G2='UNSPECIFIED: test-field wall supplies no gravitational matching or material action',
                G3='NOT_RUN: original Neumann vector not boost covariant')),
        linear_junction=linear, reciprocal_diagnostic_control=reciprocal,
        trace_quadratures=traces, cutoff_null_ray_bounds=delays, boost_controls=boosts,
        uncomputed=['finite-packet PDE evolution','observer waveform','Misner-Sharp mass response',
                    'global causal consistency','surface matter action and stability',
                    'finite-thickness regularity','quantum event law'],
        sources={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths},
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path,
        default=Path('experiments/throat_dynamics/runs/20261014_null_horizon_preflight/result.json'))
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('label','checks','laws','cutoff_null_ray_bounds')}, indent=2))
