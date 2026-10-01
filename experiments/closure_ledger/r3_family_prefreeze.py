"""PRE-FREEZE method checks for the family/transverse-stability test.
No computation at points of the two-return family is made here."""
import time
import numpy as np
from geometrodynamics.waves import r3_family as rf, esu_floquet as fl

rng = np.random.default_rng(17)
z = np.r_[1.002, .003, .1*rng.normal(size=10)]
print('round trip (12D, anisotropic): %.1e' % np.abs(rf.to_section(rf.to_state(z))-z).max())
Z0 = np.r_[1., np.zeros(11)]
t = time.time()
p, J = rf.DP(Z0, 12)
ev = np.linalg.eigvals(J)
tr310 = np.trace(fl.monodromy('T', 2))
print('DP(z*): fixed %.1e, multipliers %.6f / %.6f, max block trace - #310 %.1e, %.0fs' % (
    np.abs(p-Z0).max(), abs(ev).max(), abs(ev).min(),
    max(abs(J[2+2*k, 2+2*k]+J[3+2*k, 3+2*k]-tr310) for k in range(5)), time.time()-t))
zt = np.r_[1.0004, -.002, .05, .03, -.04, .1, .02, -.05, .01, .02, -.03, .04]   # arbitrary point, off the family
for dims in (6, 12):
    p, Jj = rf.DP(zt[:dims], dims)
    Jf = rf._fd(lambda u: rf.P(np.r_[u, np.zeros(12-dims)])[:dims], zt[:dims], 1e-6)
    print('dims %d: |P_jet - P_event| %.1e, |DP_jet - DP_fd|/|DP| %.1e' % (
        dims, np.abs(p-rf.P(np.r_[zt[:dims], np.zeros(12-dims)])[:dims]).max(), np.abs(Jj-Jf).max()/np.abs(Jj).max()))
