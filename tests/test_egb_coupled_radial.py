"""Coupled field equations, independent action Hessian, charge, and gauge checks."""
import math
import numpy as np
import pytest
from geometrodynamics.bulk.egb_coupled_radial import CoupledRadial,egb_tensor
from geometrodynamics.bulk.egb_shell_dynamics import junction_tensor,RoundEGBShell
from geometrodynamics.bulk.egb_shell_action import ScalarShellAction

BRACKETS=[(-6.5,-6.2),(-2.6,-2.4),(1.8,2.1),(6.2,6.6),(30.,32.)]

@pytest.fixture(scope='module')
def modes():
    m=CoupledRadial()
    return m,[m.root(b) for b in BRACKETS]


def test_bulk_equations_by_orthonormal_riemann_contraction(modes):
    m,roots=modes
    for value in roots[:3]:
        profile,_=m.mode(value)
        for u in [.4,.8,1.3]:
            f,fp=profile(u);fpp=-3/math.tan(u)*fp-(value+4+m.e)*f
            A=(-fpp-1/math.tan(u)*fp-2*f)/m.R**2
            B=(-2*f-2/math.tan(u)*fp)/m.R**2
            def tensor(eps):
                sec={(0,1):eps*(value*f-2*fpp/m.C)/m.R**2}
                for i in range(2,5):
                    sec[0,i]=eps*(value*f-2/math.tan(u)*fp/m.C)/m.R**2
                    sec[1,i]=1/m.R**2+eps*A
                    for j in range(i+1,5):sec[i,j]=1/m.R**2+eps*B
                return egb_tensor(sec,m.shell.alpha)
            derivative=(tensor(1e-5)-tensor(-1e-5))/2e-5
            density=m.action.q2*(value+m.e)*f
            np.testing.assert_allclose(derivative,np.eye(5)*density,rtol=3e-7,atol=1e-8)


def test_momentum_constraint_by_mixed_riemann_contraction(modes):
    m,roots=modes
    for value in roots[:3]:
        profile,_=m.mode(value)
        fp=float(profile(.8)[1])
        sectionals={(i,j):1/m.R**2 for i in range(1,5) for j in range(i+1,5)}
        mixed={(0,i,1,i):-fp/m.R for i in range(2,5)}
        tensor=egb_tensor(sectionals,m.shell.alpha,mixed)
        assert tensor[0,1]==pytest.approx(-3*m.C*fp/m.R,rel=1e-11,abs=1e-11)


def test_junctions_by_full_davis_contraction(modes):
    m,roots=modes
    for value in roots[:3]:
        f,z=m.mode(value);fa,fp=f(m.a);xi=(z-fa)*math.tan(m.a)
        Q=1/(m.R*math.tan(m.a));dQ=(fp-fa/math.tan(m.a)-xi/math.sin(m.a)**2)/m.R
        k0=-value*xi/m.R-2*fp/(m.C*m.R)
        def stress(eps):
            return junction_tensor(eps*k0,Q+eps*dQ,-eps*value*z/m.R**2,(1-2*eps*z)/m.r**2,m.shell.alpha)
        derivative=(stress(1e-6)-stress(-1e-6))/2e-6
        h=value+m.e;expected=np.diag([m.W*h*fa/m.shell.eta]+[m.W*h*fa]*3)
        np.testing.assert_allclose(derivative,expected,rtol=1e-7,atol=1e-8)
        assert max(abs(v) for v in m.junction_residuals(value).values())<2e-8
        assert abs(m.charge_residual(value))<2e-9


def test_norm_from_spectral_derivative_and_two_fixed_shell_gauges(modes):
    m,roots=modes
    for value in roots:
        norm=m.physical_norm(value)
        assert m.spectral_norm(value)==pytest.approx(norm,rel=2e-8,abs=1e-9)
        for power in [1,3]:
            assert m.symplectic_norm(value,power)['total']==pytest.approx(norm,rel=2e-8,abs=1e-9)
    assert [np.sign(m.physical_norm(z)) for z in roots[:4]]==[-1,1,-1,1]
    assert roots[0]<roots[1]<0<roots[2]


def test_independent_finite_element_convergence_and_physical_inertia(modes):
    m,roots=modes
    coarse=m.finite_element_spectrum(80);fine=m.finite_element_spectrum(160)
    for row in (coarse,fine):
        assert row['maximum_imaginary_part']<1e-8
        assert row['physical_negative_directions']==2
        assert row['auxiliary_negative_directions']==1
    for i in range(4):
        err0=abs(coarse['low_eigenvalues'][i]-roots[i]);err1=abs(fine['low_eigenvalues'][i]-roots[i])
        assert err1<.27*err0
        assert err1<.002
    K,M=m.pencil(40)
    # Multiplying by the physical metric yields a symmetric action stiffness.
    V=K@np.linalg.solve(M,K)+m.e*K
    np.testing.assert_allclose(V,V.T,rtol=1e-10,atol=1e-8)


def test_structural_kinetic_sign_scaling_and_exceptional_pole(modes):
    m,roots=modes
    value=roots[2];k=2.1
    scaled=CoupledRadial(ScalarShellAction(RoundEGBShell(radius=k,alpha=-.06*k*k)))
    assert scaled.root(BRACKETS[2])==pytest.approx(value,rel=1e-10)
    assert scaled.physical_norm(value)==pytest.approx(k**4*m.physical_norm(value),rel=1e-9)
    for alpha in [-.04,-.06,-.1,-.2]:
        family=CoupledRadial(ScalarShellAction(RoundEGBShell(alpha=alpha)))
        assert family.D<0
        assert family.finite_element_spectrum(50)['physical_negative_directions']==2
    pole=-m.d/m.x
    assert np.isfinite(m.determinant(pole))
    assert abs(m.determinant(pole))>1e-5
    assert abs(m.determinant(-m.e))>1e-5 # no even gauge root


def test_odd_gauge_mode_is_null_and_physical_odd_modes_are_positive():
    m=CoupledRadial()
    # The exceptional cos(u) solution has h=0; with z=0 its integrated
    # physical form vanishes by the exact boundary identity (no quadrature).
    boundary=math.sin(m.a)**3*math.cos(m.a)*math.sin(m.a)
    subtract=math.sin(m.a)**3*math.tan(m.a)*math.cos(m.a)**2
    assert boundary==pytest.approx(subtract,abs=1e-15)
    # Dirichlet at mouth and equator for physical odd modes: eigenvalue > -e.
    assert m.e>0
    vals=np.linspace(-m.e+.01,80,180)
    pairs=[(lo,hi) for lo,hi in zip(vals[:-1],vals[1:]) if m.determinant(lo,'odd')*m.determinant(hi,'odd')<0]
    assert len(pairs)>=2
    for pair in pairs[:2]:
        root=m.root(pair,'odd')
        assert m.physical_norm(root,'odd')>0
        assert root>-m.e


def test_full_curvature_action_recovers_myers_kinetic_boundary_term():
    from experiments.closure_ledger.egb_coupled_action_check import check_action_hessian
    coarse=check_action_hessian(1e-3);fine=check_action_hessian(5e-4)
    assert coarse['relative_error']<1e-6
    assert fine['relative_error']<.3*coarse['relative_error']
    assert fine['relative_error']<2e-7


def test_result_provenance_and_separate_admission_failure():
    import hashlib
    import json
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    archive=root/'experiments/closure_ledger/runs/20261013_egb_coupled_radial/result.json'
    result=json.loads(archive.read_text())
    assert result['label']=='COUPLED_RADIAL_ADMISSION_FAILED'
    assert all(result['numerical_gates'].values())
    assert not result['auxiliary_mass_is_physical_norm']
    for path,digest in result['sources'].items():
        assert hashlib.sha256((root/path).read_bytes()).hexdigest()==digest
    shell_mode=result['even_modes'][2]
    assert shell_mode['omega2_R2']>0
    assert shell_mode['physical_norm_no_S3_over_kappa']<0
    assert result['odd_modes_in_scanned_range'][0]['omega2_R2']<0
