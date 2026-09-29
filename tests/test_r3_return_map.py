import numpy as np
import pytest
from geometrodynamics.waves import r3_return_map as rm
from geometrodynamics.waves import nonlinear_supported_tt as d
from geometrodynamics.waves import esu_floquet as fl
from geometrodynamics.waves.jets import Jet, variables, linear_part

LAM, A0, BETA = 3., 2.5, .7
C = dict(a=.3, b=-.4, c=.2, e=.25, f=-.15)


def _trig(v, name):
    return getattr(v, name)() if isinstance(v, Jet) else getattr(np, name)(v)


def toy(z, alpha0=A0):
    """Hyperbolic x exact-twist map conjugated by nonlinear symplectic shears."""
    def SP(z, s):
        Q1, P1, Q2, P2 = z
        return [Q1, P1+s*(C['b']*Q2*Q2+2*C['c']*Q1*Q2), Q2, P2+s*(3*C['a']*Q2*Q2+2*C['b']*Q1*Q2+C['c']*Q1*Q1)]

    def SQ(z, s):
        Q1, P1, Q2, P2 = z
        return [Q1+s*C['f']*P2*P2, P1, Q2+s*(3*C['e']*P2*P2+2*C['f']*P1*P2), P2]

    def R(z):
        Q1, P1, Q2, P2 = z
        al = alpha0+2*BETA*(Q2*Q2+P2*P2)/2
        return [LAM*Q1, P1/LAM, Q2*_trig(al, 'cos')+P2*_trig(al, 'sin'), -Q2*_trig(al, 'sin')+P2*_trig(al, 'cos')]
    return SQ(SP(R(SP(SQ(z, -1), -1)), 1), 1)


def test_jet_arithmetic_matches_functions():
    x, y = variables([.3, -.2], 4)
    f = ((x*y+1).exp()*(x+2).sqrt()/(1+y*y)).sin()
    g = lambda u, v: np.sin(np.exp(u*v+1)*np.sqrt(u+2)/(1+v*v))
    for dx, dy in ((1e-4, 2e-4), (-2e-4, 1e-4)):   # order-4 truncation: error ~ d^5
        taylor = f(np.array([dx, dy]))
        assert abs(taylor-g(.3+dx, -.2+dy)) < 1e-15


@pytest.mark.parametrize('alpha0', [A0, 3.0])
def test_normal_form_recovers_exact_twist_and_action(alpha0):
    z = variables(np.zeros(4), 3)
    nf = rm.normal_form(toy(z, alpha0), point=np.zeros(4))
    assert abs(nf['nu']-BETA/np.pi) < 1e-10 and abs(nf['dissipative']) < 1e-10
    assert abs(nf['theta']-alpha0) < 1e-12


def test_invariant_circles_recover_exact_rotation_and_action():
    P = lambda v: np.array(toy(list(v)), dtype=float)
    A = linear_part(toy(variables(np.zeros(4), 1)))
    K, om = rm.linear_circle(A, .05, 31, point=np.zeros(4))
    c = rm.invariant_circle(P, .05, K, om, 31)
    assert c['residual'] < 1e-12
    assert abs(c['omega']-A0-2*BETA*abs(c['action'])) < 1e-12


def test_reduced_equations_match_full_system():
    rng = np.random.default_rng(3)
    for _ in range(10):
        A, Ap, q, qp, x, xp = 1+.1*rng.normal(), .1*rng.normal(), .3*rng.normal(), rng.normal(), .1*rng.normal(), .3*rng.normal()
        y = d.pack(A, Ap, np.r_[q, 0, 0, 0], np.r_[qp, 0, 0, 0], np.diag(np.exp(2*x*np.diag(rm.B0))), xp*rm.B0)
        f = d.conformal_rhs(y)
        _, _, _, _, _, Lp = d.unpack(f)
        red = rm.rhs([A, Ap, q, qp, x, xp])
        assert np.allclose(red, [f[0], f[1], f[3], f[7], xp, np.trace(Lp @ rm.B0)], atol=1e-13)
        assert abs(rm.constraint([A, Ap, q, qp, x, xp])-d.constraints(y)['residual'][0]) < 1e-13


def test_section_data_satisfy_constraint():
    z = variables(rm.ZSTAR, 2)
    y = rm.section_to_state([z[0]+.01, z[1]-.02, z[2]+.03, z[3]+.01])
    assert abs(rm.constraint(y).c).max() < 1e-13


def test_linear_return_map_matches_310_and_is_symplectic():
    A = linear_part(rm.jet_return_map(1, 2048)['P'])
    assert abs(A[2, 2]+A[3, 3]-np.trace(fl.monodromy('T', 2))) < 1e-10
    assert np.abs(A.T @ rm.OMEGA @ A-rm.OMEGA).max() < 1e-11
    assert np.abs(A[:2, 2:]).max() == 0 and np.abs(A[2:, :2]).max() == 0


ARCHIVE = __import__('pathlib').Path(__file__).resolve().parents[1]/'experiments/closure_ledger/runs/20260929_r3_return_map/return_map.json'


@pytest.mark.skipif(not ARCHIVE.exists(), reason='archive not generated')
def test_return_map_archive_rescores_and_binds_sources():
    import json
    from experiments.closure_ledger import r3_return_map_probe as probe
    from experiments.closure_ledger.esu_floquet_probe import close
    rec = json.loads(ARCHIVE.read_text())
    assert rec['sources'] == probe.sources()
    again = json.loads(json.dumps(probe.score(rec['m1'], rec['m2'], rec['result']['theta0'])))
    assert again['label'] == rec['result']['label'] and again['checks'] == rec['result']['checks']
    assert close(again, rec['result'], 1e-9)
    bad = json.loads(ARCHIVE.read_text())
    bad['m2']['circles'][0]['omega'] += 1e-6
    assert not close(json.loads(json.dumps(probe.score(bad['m1'], bad['m2'], bad['result']['theta0']))), rec['result'], 1e-9)
