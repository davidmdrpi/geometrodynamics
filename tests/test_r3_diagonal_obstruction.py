import numpy as np
import gmpy2 as gp
from gmpy2 import mpfr as mp
from geometrodynamics.waves import r3_diagonal_taylor as dt
from geometrodynamics.waves import r3_diagonal_symmetry as ds
from geometrodynamics.waves import r3_high_precision_shooting as sh
from experiments.closure_ledger import r3_diagonal_obstruction_controls as controls
from experiments.closure_ledger import r3_diagonal_obstruction_probe as probe


def test_both_formulations_match_original_rhs_off_section():
    with gp.context(gp.get_context(),precision=224):
        z=['1.0004','-.002','.05','.03','-.04','.1']
        for form in ('momentum','velocity'):
            y=dt.to_state(z,form);y[2]=mp('.17');h=y[0]**2-y[2]**2/6
            velocity=[y[6]/h,y[7]/h] if form=='momentum' else y[6:]
            state=np.array(y[:6]+velocity,float);expected=ds.reduced_rhs(0,state)
            c=dt.coefficients(y,2,form);actual=np.array([a[1] for a in c],float)
            if form=='momentum':actual[6:]=(actual[6:]-(2*float(y[0]*y[1])-float(y[2]*y[3])/3)*np.array(velocity,float))/float(h)
            np.testing.assert_allclose(actual,expected,atol=5e-15,rtol=2e-14)


def test_independent_high_precision_maps_agree():
    r=controls.off_family();assert r['ok']


def test_known_twelfth_harmonic_and_zero_controls():
    assert controls.synthetic()['ok'];assert controls.synthetic('0')['ok']


def test_high_precision_inverse_and_isolated_root_control():
    with gp.context(gp.get_context(),precision=336):
        e=mp('1e-24')
        def P(z):
            th,I=z;p=I+e*gp.sin(12*th);return [th+p,p]
        result=probe.isolated_diagnostics(P,P,[[mp(0),mp(0)]])
        assert result['ok'];assert mp(result['inverse_norm'])>mp('1e22')


def test_returned_residual_is_recomputed_after_iteration_cap():
    with gp.context(gp.get_context(),precision=224):
        def P(z):return [z[0]+z[1],z[1]+mp('.125')]
        pt=dict(nodes=[['0','0']],lam='0',g=[0,-1],c0=[0,0],t0=[1,0],J=[[0,1,0],[0,0,1],[1,0,0]])
        result=sh.solve(P,pt,224,'1e-50',maxit=1)
        assert result['ok'] and mp(result['residual'])<mp('1e-50')


def test_failed_phase_cannot_score_as_suppression():
    rows=[dict(j=j,error='deliberate convergence failure') for j in range(probe.N)]
    result=probe.scan_score(rows,dict(ok=True))
    assert result['label']=='NUMERICALLY_UNRESOLVED'
    assert len(result['failed_phases'])==probe.N
