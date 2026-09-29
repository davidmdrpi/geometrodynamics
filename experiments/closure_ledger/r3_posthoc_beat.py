"""POST-HOC diagnostic (not registered): why N3 failed in the R3 run.

2*rho0 is close to an integer, so the per-period phase increments carry a
beat of period 1/(1-frac(2 rho0)) ~ 33 clock periods, longer than K/2 = 24.
Prints the beat period and, per run, the raw mean shift and half-window
differences. No verdict depends on this script.
"""
import json
import numpy as np
from geometrodynamics.waves import r3_resonance as r3
from experiments.closure_ledger.r3_resonance_probe import RUN_DIR

rec = json.loads((RUN_DIR/'r3_resonance.json').read_text())
rho0 = rec['result']['rho0']
print('2*rho0 mod 1 =', 2*rho0 % 1, ' beat period =', 1/(1-(2*rho0 % 1)))
for run in rec['raw']['runs']:
    if run['integrator'] != 'primary':
        continue
    inc = run['increments']
    print(run['eps'], run['pol'],
          'raw mean - rho0 = %+.2e' % (np.mean(inc)/(2*np.pi)-rho0),
          'B48-B(first24) = %+.2e' % (r3.birkhoff(inc)-r3.birkhoff(inc[:24])),
          'B(last24)-B48 = %+.2e' % (r3.birkhoff(inc[24:])-r3.birkhoff(inc)))
