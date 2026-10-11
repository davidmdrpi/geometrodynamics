import numpy as np
import pytest
from geometrodynamics.waves import r3_diagonal_symmetry as ds
from geometrodynamics.waves import r3_family as rf
from geometrodynamics.waves import nonlinear_supported_tt as full

Z = np.array([1.0004, -.002, .05, .03, -.04, .1])


def test_spatial_representation_matches_matrix_conjugation():
    for R in (ds.CYCLE, ds.REFLECTION):
        y = rf.to_state(np.r_[Z, np.zeros(6)])
        A, Ap, q, qp, M, L = full.unpack(y)
        expected = rf.to_section(full.pack(A, Ap, q, qp, R@M@R.T, R@L@R.T))[:6]
        np.testing.assert_allclose(ds.spatial(Z, R), expected, atol=2e-15)
    np.testing.assert_allclose(ds.spatial(ds.spatial(ds.spatial(Z))), Z, atol=1e-15)
    np.testing.assert_allclose(ds.spatial(ds.spatial(Z, ds.REFLECTION), ds.REFLECTION), Z, atol=1e-15)
    assert ds.chirality(ds.spatial(Z, ds.REFLECTION)) == pytest.approx(-ds.chirality(Z))


def test_reduced_rhs_is_diagonal_restriction_of_matrix_equations():
    y = ds.reduced_initial(Z)
    y[2] = .17
    A, Ap, q, qp = y[:4]
    M = np.diag(np.exp(2*y[4:6]@ds.BASIS))
    L = np.diag(y[6:8]@ds.BASIS)
    dy = full.conformal_rhs(full.pack(A, Ap, np.r_[q,0,0,0], np.r_[qp,0,0,0], M,L))
    a,b,c,d,md,ld = full.unpack(dy)
    expected = np.r_[a,b,c[0],d[0],ds.BASIS@np.diag(np.linalg.solve(M,md)/2),ds.BASIS@np.diag(ld)]
    np.testing.assert_allclose(ds.reduced_rhs(0,y),expected,atol=3e-15)


def test_half_clock_square_and_independent_matrix_off_family():
    h = ds.clock_map(Z)['z']
    hh = ds.clock_map(h)['z']
    p = ds.clock_map(Z,half=False)['z']
    np.testing.assert_allclose(hh,p,atol=2e-11,rtol=0)
    np.testing.assert_allclose(h,ds.clock_map(Z,matrix=True)['z'],atol=2e-12,rtol=0)


def test_fourier_fit_does_not_impose_spatial_symmetry():
    th=np.linspace(0,2*np.pi,80,endpoint=False)
    nodes=np.zeros((80,6));nodes[:,0]=1+.1*np.sin(5*th)
    nodes[:,2]=np.cos(th);nodes[:,4]=np.sin(th)
    f,_=ds.fourier_fit(nodes,8)
    np.testing.assert_allclose(f(th),nodes,atol=2e-14)
    assert ds.first_harmonic(4,3)==12
    assert ds.first_harmonic(4,1)==4
    with pytest.raises(ValueError): ds.clock_map(np.zeros(12))
