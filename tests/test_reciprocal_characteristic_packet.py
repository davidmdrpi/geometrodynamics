"""Bulk PDE, Einstein mass response and fixed-null-map compatibility checks."""
import math
import numpy as np
import pytest
import sympy as s
from scipy.integrate import simpson
from geometrodynamics.tangherlini.characteristic_packet import (
    compact_packet,diamond,geometry,radius,seam_response,solve_pair,
)


def test_implicit_geometry_and_regular_horizon():
    for b in (.5,1.,2.):
        U=np.linspace(-8*b,0,41)[:,None]
        V=np.linspace(0,4*b,51)[None,:]
        r,F,Q=geometry(U,V,b=b)
        H=-b*b*(r-b)/(r+b)*np.exp(2*r/b)
        assert np.max(np.abs(H-U*V))/b**2<1e-12
        assert np.allclose(r[-1],b)
        assert np.allclose(F[-1],4/np.e**2)
        assert np.all(Q>0)


def test_diamond_against_exact_nonzero_potential_and_free_solution():
    errors=[]
    for n in (32,64,128):
        U=np.linspace(0,1,n+1)**1.2
        V=np.linspace(0,1,n+1)
        # exp(aU+bV) solves u_UV+Q u=0 for Q=-ab=2.
        exact=np.exp(U[:,None]-2*V[None,:])
        got=diamond(U,V,exact[0],exact[:,0],np.full((n,n),2.))
        errors.append(np.max(np.abs(got-exact)))
        free=U[:,None]**2+np.sin(V)[None,:]
        zero=diamond(U,V,free[0],free[:,0],np.zeros((n,n)))
        assert np.max(np.abs(zero-free))<2e-13
    assert all(3.7<a/b<4.3 for a,b in zip(errors[:-1],errors[1:]))
    with pytest.raises(ValueError,match='corner'):
        diamond(U,V,np.ones(len(V)),np.zeros(len(U)),np.zeros((n,n)))


def test_mass_law_from_spherical_einstein_constraints():
    # Differentiate the geometric definition and insert rr/angular equations.
    r,a,d,F,su,sv,pu,pv,k=s.symbols('r a d F su sv pu pv k', nonzero=True)
    mixed=-2/r*(a*d+F/4)
    ruu=2*su*a-k*r*pu**2/3
    rvv=2*sv*d-k*r*pv**2/3
    m=r*r*(1+4*a*d/F)
    mu=s.diff(m,r)*a+s.diff(m,a)*ruu+s.diff(m,d)*mixed+s.diff(m,F)*2*F*su
    mv=s.diff(m,r)*d+s.diff(m,a)*mixed+s.diff(m,d)*rvv+s.diff(m,F)*2*F*sv
    assert s.simplify(mu+4*k*r**3*d*pu**2/(3*F))==0
    assert s.simplify(mv+4*k*r**3*a*pv**2/(3*F))==0


def test_compact_support_and_seam_constraint_obstruction():
    V=np.linspace(.5,4,2049)
    a=compact_packet(V)
    assert np.all(a[(V<=1)|(V>=3)]==0)
    out=seam_response(V,a)
    ap=np.gradient(a,V,edge_order=2)
    # Independently integrate the positive-kernel mismatch at the final point.
    expected=-simpson((V[-1]-V)**2/V[-1]*ap**2/3,x=V)
    assert out['induced_radius_mismatch'][-1]==pytest.approx(expected,rel=1e-9)
    assert out['incoming_delta_r'][-1]<out['outgoing_delta_r'][-1]<0
    assert out['incoming_delta_m'][-1]>0
    assert out['mass_parameter_jump'][-1]<0
    assert np.all(seam_response(V,np.zeros_like(V))['induced_radius_mismatch']==0)


def test_bulk_refinement_mass_path_and_symplectic_conservation():
    runs=[solve_pair(n) for n in (64,128,256)]
    for key in ('field','outgoing'):
        errors=[np.linalg.norm(b[key][::2,::2]-a[key])/np.linalg.norm(b[key][::2,::2])
                for a,b in zip(runs[:-1],runs[1:])]
        assert 3.4<errors[0]/errors[1]<4.6
    for side in ('incoming_mass','outgoing_mass'):
        err=[r[side]['relative_path_discrepancy'] for r in runs]
        assert all(3.7<a/b<4.3 for a,b in zip(err[:-1],err[1:]))
        assert err[-1]<1e-4
    for side in ('incoming_symplectic','outgoing_symplectic'):
        assert runs[-1][side]['relative_residual']<1e-4
    # No vanishing mismatch as the grid refines.
    mismatch=[r['seam']['induced_radius_mismatch'][-1] for r in runs]
    assert all(x<-1e-4 for x in mismatch)
    assert abs(mismatch[-1]-mismatch[-2])<abs(mismatch[-2]-mismatch[-3])/3
    # Outgoing horizon mass is integrated in U, independently of the analytic
    # reciprocal pullback used by seam_response in V.
    expected=runs[-1]['seam']['outgoing_delta_m']
    measured=runs[-1]['outgoing_mass']['delta_m'][:,0]
    assert np.max(np.abs(measured-expected))/np.max(np.abs(expected))<1e-4


def test_covariance_amplitude_order_and_angular_scope():
    base=solve_pair(64)
    for boost in (.5,2.,4.):
        got=solve_pair(64,boost=boost)
        assert np.allclose(got['observer_t'],base['observer_t'],atol=1e-13)
        assert np.max(np.abs(got['observer_u']-base['observer_u']))<1e-13
        assert np.max(np.abs(got['outgoing_mass']['delta_m']-base['outgoing_mass']['delta_m']))<1e-13
    for scale in (.5,2.):
        got=solve_pair(64,amplitude=.01*scale)
        assert np.allclose(got['field'],scale*base['field'],atol=1e-14)
        assert np.allclose(got['incoming_mass']['delta_m'],scale**2*base['incoming_mass']['delta_m'],atol=1e-14)
    odd=solve_pair(64,ell=1)
    assert 'incoming_mass' not in odd
    assert np.array_equal(odd['outgoing'][:,0],-odd['field'][-1])
    vacuum=solve_pair(32,amplitude=0)
    assert not np.any(vacuum['field'])
    assert not np.any(vacuum['incoming_mass']['delta_m'])


def test_published_archive_provenance_and_scope():
    import hashlib
    import json
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    folder=root/'experiments/throat_dynamics/runs/20261014_reciprocal_packet'
    data=json.loads((folder/'result.json').read_text())
    assert data['label']=='RECIPROCAL_FIXED_NULL_MATCHING_FAILS_AT_ORDER_A2'
    assert all(data['numerical_gates'].values())
    assert data['mismatch_resolved']
    assert 'nonlinear coupled metric evolution' in data['scope']['uncomputed']
    for name,digest in data['sources'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
    for name,digest in data['artifacts'].items():
        assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest
    for ell in (0,1):
        with np.load(folder/f'fields_N64_l{ell}.npz') as fields:
            assert fields['field'].shape==(65,129,2)
            assert fields['outgoing'].shape==(129,129,2)
            assert np.all(np.isfinite(fields['outgoing']))
