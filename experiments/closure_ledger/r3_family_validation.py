"""Off-diagonal validation follow-up (added 2026-10-01, #319 review).

The frozen S3 check probed real parts of the off-diagonal eigenvectors. At the
defective unit multipliers those collapse onto the rotation eigenspace and
miss the generalised (angular-momentum) directions. This supplement probes,
at the first registered sample:
1. six independent off-diagonal coordinate directions, by central differences
   of P∘P on the 12-coordinate map;
2. the strongest secular-shear direction (top right singular vector of
   N = M2_off - I) on the full system, with the first-order gravitational
   angular momentum cancelled by quartet components q_1..q_3, so that every
   constraint holds.
It also records ||N||, ||N^2|| and ||M2_off^100|| (secular, linear growth).
The frozen producer, archives and labels are unchanged.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp
from geometrodynamics.waves import r3_family as rf
from geometrodynamics.waves import nonlinear_supported_tt as d
from experiments.closure_ledger import r3_family_probe as probe

EPS = 1e-6
OUT = probe.RUN_DIR/'validation_followup.json'


def P2(z):
    return rf.P(rf.P(z))


def compensated_state(z):
    """Full state for section point z with the momentum constraint satisfied by q_1..q_3."""
    y = rf.to_state(z)
    base = d.constraints(y)['residual'][1:4]                  # gravity - matter, matter = 0 here
    J = np.zeros((3, 3))
    for k in range(3):
        yk = y.copy()
        yk[4+k] = 1.
        J[:, k] = d.constraints(yk)['residual'][1:4]-base     # linear in q_i at fixed q'
    y[4:7] = np.linalg.solve(J, -base)
    E = d.constraints(y)['residual'][0]
    y[7] = -np.sqrt(y[7]**2-2*E)                               # re-solve the Hamiltonian constraint
    return y


def full_two_returns(y0):
    def clock(t, y):
        return y[3]
    clock.direction = -1
    y = y0
    for _ in range(2):
        s = solve_ivp(lambda t, u: d.conformal_rhs(u), (0., np.pi+.8), y, method='DOP853', rtol=1e-12, atol=1e-14,
                      events=[clock], max_step=.05)
        y = s.y_events[0][[k for k, t in enumerate(s.t_events[0]) if t > 1.][0]]
    return y


def main():
    S = json.loads((probe.RUN_DIR/'stage_S.json').read_text())
    smp = S['samples'][0]
    M2 = np.array(smp['M2'])
    z = np.r_[np.array(smp['v'][:6]), np.zeros(6)]
    B = M2[6:, 6:]
    N = B-np.eye(6)
    coord = []
    for k in range(6):
        e = np.zeros(12)
        e[6+k] = 1.
        lin = (P2(z+EPS*e)-P2(z-EPS*e))/(2*EPS)
        coord.append(float(np.linalg.norm(lin-M2 @ e)/np.linalg.norm(M2 @ e)))
    v = np.r_[np.zeros(6), np.linalg.svd(N)[2][0]]
    yp, ym = compensated_state(z+EPS*v), compensated_state(z-EPS*v)
    cons = max(np.abs(d.constraints(y)['residual']).max() for y in (yp, ym))
    fp, fm = full_two_returns(yp), full_two_returns(ym)
    lin = (rf.to_section(fp)-rf.to_section(fm))/(2*EPS)
    shear = float(np.linalg.norm(lin-M2 @ v)/np.linalg.norm(M2 @ v))
    final_cons = max(np.abs(d.constraints(y)['residual']).max() for y in (fp, fm))
    rec = dict(sample_index=smp['index'], coordinate_directions_relative_error=coord,
               compensated_shear_relative_error=shear, compensated_initial_constraint=float(cons),
               compensated_final_constraint=float(final_cons),
               norm_N=float(np.linalg.norm(N, 2)), norm_N2=float(np.linalg.norm(N @ N, 2)),
               norm_M2off_pow100=float(np.linalg.norm(np.linalg.matrix_power(B, 100), 2)))
    if OUT.exists():
        raise FileExistsError('append-only: '+str(OUT))
    OUT.write_text(json.dumps(rec, indent=1)+'\n')
    print(json.dumps(rec, indent=1))


if __name__ == '__main__':
    main()
