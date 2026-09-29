"""PRE-FREEZE checks for the R3 extension: equations, parametrisation, linear map.
No cubic ESU normal form and no new ESU circle is computed here."""
import time
import numpy as np
from geometrodynamics.waves import r3_extension as rx, nonlinear_supported_tt as d, esu_floquet as fl
from geometrodynamics.waves.jets import variables, linear_part, Jet
from scipy.linalg import expm

rng = np.random.default_rng(5)
worst = 0.
for _ in range(10):
    A, Ap, q, qp = 1+.1*rng.normal(), .1*rng.normal(), .3*rng.normal(), rng.normal()
    b = sum(.1*rng.normal()*e for e in rx.E)
    M = expm(2*b)
    L = rng.normal(size=(3, 3))*.2
    y = d.pack(A, Ap, np.r_[q, 0, 0, 0], np.r_[qp, 0, 0, 0], M, L)
    f = d.conformal_rhs(y)
    _, _, _, _, Mp, Lp = d.unpack(f)
    r = rx.rhs([A, Ap, q, qp, M.tolist(), L.tolist()])
    worst = max(worst, abs(r[1]-f[1]), abs(r[3]-f[7]), np.abs(np.array(r[4])-Mp).max(), np.abs(np.array(r[5])-Lp).max(),
                abs(rx.constraint([A, Ap, q, qp, M.tolist(), L.tolist()])-d.constraints(y)['residual'][0]))
print('matrix rhs/constraint vs conformal_rhs: %.1e' % worst)
z = variables(np.r_[1., np.zeros(11)], 3)
back = rx.state_to_section(rx.section_to_state(z))
print('section round trip (order 3): %.1e' % max(np.abs((b-a).c).max() for a, b in zip(z, back)))
print('constraint on section data: %.1e' % np.abs(rx.constraint(rx.section_to_state(z)).c).max())
t = time.time()
out = rx.jet_return_map(1, 2048)
A1 = linear_part(out['P'])
tr310 = np.trace(fl.monodromy('T', 2))
blocks = [np.trace(A1[2+2*k:4+2*k, 2+2*k:4+2*k]) - tr310 for k in range(5)]
mask = np.zeros_like(A1, bool)
for bl in [(0, 1)]+[(2+2*k, 3+2*k) for k in range(5)]:
    mask[np.ix_(bl, bl)] = True
print('order-1 12D map: block traces - #310:', np.round(blocks, 12), 'off-block %.1e' % np.abs(A1[~mask]).max(),
      'hyperbolic', np.round(np.linalg.eigvals(A1[:2, :2]), 6), '%.0fs' % (time.time()-t))
