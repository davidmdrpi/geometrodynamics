"""Post-freeze uniqueness extension to de55f3f; not a frozen prediction.

All smooth nontrivial single-real-scalar supports of the exact round ESU
are homogeneous. The proof, including zero slices, is in the companion doc.
"""

from functools import lru_cache

import numpy as np

from . import scalar_esu_support as support
from . import reciprocal_scalar_tt as rt


SEED = 2026090916
REQUIRED_CHECKS = (
    "coefficient_expansion", "rank_and_alignment", "regularity_contradiction",
    "independent_momentum_route", "degenerate_controls", "failure_paths",
)


def coefficient_conditions(A, Ad, Add, B, Bd, Bdd, radius=1.):
    """h^3 (Box-1/a^2)(1/h) = C + L.x + x^T M x on unit S3."""
    B, Bd, Bdd = map(np.asarray, (B, Bd, Bdd))
    if radius <= 0 or not np.isfinite(radius):
        raise ValueError("positive finite radius required")
    if any(v.shape != (4,) or not np.isfinite(v).all() for v in (B,Bd,Bdd)):
        raise ValueError("finite ambient four-vectors required")
    if not np.isfinite([A,Ad,Add]).all():
        raise ValueError("finite scalar jets required")
    C = A*Add-2*Ad**2+(2*(B@B)-A*A)/radius**2
    L = A*Bdd+Add*B-4*Ad*Bd+A*B/radius**2
    M = (np.outer(B,Bdd)+np.outer(Bdd,B))/2-2*np.outer(Bd,Bd)
    return C,L,M


@lru_cache(maxsize=1)
def exact_certificate():
    import sympy as s
    a = s.symbols("a", positive=True)
    A,Ad,Add = s.symbols("A Ad Add", real=True)
    B,V,W,x = [s.Matrix(s.symbols(prefix+"0:4", real=True)) for prefix in ("B","V","W","x")]
    h,ht,htt = A+B.dot(x), Ad+V.dot(x), Add+W.dot(x)
    M = (B*W.T+W*B.T)/2-2*V*V.T
    C = A*Add-2*Ad**2+(2*B.dot(B)-A*A)/a**2
    L = A*W+Add*B-4*Ad*V+A*B/a**2
    # Direct differentiated reciprocal numerator, using the unit-sphere
    # gradient norm and Laplacian of the ambient linear function B.x.
    numerator = h*htt-2*ht**2+(2*(B.dot(B)-B.dot(x)**2)+3*h*B.dot(x)-h*h)/a**2
    residuals={"three_conditions":s.expand(numerator-C-L.dot(x)-(x.T*M*x)[0])}
    H = B.row_join(V).row_join(W)
    K = s.Matrix([[0,0,s.Rational(1,2)],[0,-2,0],[s.Rational(1,2),0,0]])
    factor = M-H*K*H.T
    residuals.update({f"rank_factor_{i}{j}":s.expand(factor[i,j]) for i in range(4) for j in range(4)})
    # B != 0: choose an orthogonal frame B=(b,0,0,0), b>0.
    b = s.symbols("b", positive=True)
    aligned = s.Matrix([b,0,0,0])
    N = (aligned*W.T+W*aligned.T)/2-2*V*V.T
    for i in (1,2,3):
        residuals[f"perpendicular_square_{i}"] = N[i,i]+2*V[i]**2
        residuals[f"perpendicular_acceleration_{i}"] = N[0,i].subs({V[j]:0 for j in (1,2,3)})-b*W[i]/2
    residuals["parallel_acceleration"] = N[0,0]-(b*W[0]-2*V[0]**2)
    # Projected amplitude q=A/b; use bbdd=2bd^2 from M=0 on B != 0.
    d,q,qd,qdd = s.symbols("d q qd qdd", real=True)
    dd = 2*d*d/b
    AA,AA_d,AA_dd = b*q,d*q+b*qd,dd*q+2*d*qd+b*qdd
    cc = AA*AA_dd-2*AA_d**2+(2*b*b-AA*AA)/a**2
    ll = AA*dd+AA_dd*b-4*AA_d*d+AA*b/a**2
    residuals["regularity_contradiction"] = s.expand(cc/b**2-q*ll/b**2+2*(qd**2+(q*q-1)/a**2))
    # Zero-momentum route derived from the unexpanded improved stress.
    hh,hd,hx,hdx = s.symbols("h hd hx hdx", nonzero=True, real=True)
    phi,pt,px,ptx = 1/hh,-hd/hh**2,-hx/hh**2,2*hd*hx/hh**3-hdx/hh**2
    T0i = pt*px-s.Rational(1,6)*(2*pt*px+2*phi*ptx)
    residuals["momentum_reciprocal"] = s.simplify(T0i-hdx/(3*hh**3))
    t,amplitude,phase = s.symbols("t amplitude phase", real=True)
    AA = amplitude*s.cos(t/a+phase)
    constant = AA*s.diff(AA,t,2)-2*s.diff(AA,t)**2+(2*b*b-AA*AA)/a**2
    residuals["constant_B_wave_residual"] = s.trigsimp(constant-2*(b*b-amplitude**2)/a**2)
    reduced = {key:str(s.simplify(value)) for key,value in residuals.items()}
    return dict(residuals=reduced,all_zero=all(v=="0" for v in reduced.values()),
                sym_convention="sym(U V^T)=(U V^T+V U^T)/2",
                rank_bound="M=H K H^T with H of shape 4x3; rank(M)<=3",
                alignment_scope="B!=0: Bd=lambda B and Bdd=2 lambda^2 B; at B=0 only Bd=0 follows pointwise",
                regularity_identity="q^2+a^2 qdot^2=1 contradicts |q|>1",
                constant_B_residual="2 (b^2-A0^2)/a^2",
                global_proof="docs/scalar_esu_uniqueness.md",
                proof_status="WRITTEN_ANALYTIC_PROOF_WITH_SYMBOLIC_LOCAL_CHECKS")


def independent_pointwise_checks():
    rng = np.random.default_rng(SEED)
    eta = np.diag([-1.,1.,1.,1.])
    rows=[]
    for a in (.7,1.,2.):
        B = rng.normal(0,.2,4)
        Bd,Bdd = rng.normal(0,.3,(2,4))
        A,Ad,Add = 2.,-.3,.4
        C,L,M = coefficient_conditions(A,Ad,Add,B,Bd,Bdd,a)
        points = rng.normal(size=(40,4))
        points /= np.linalg.norm(points,axis=1)[:,None]
        frames = np.einsum("iab,pb->pia",rt.QUATERNION_DERIVATIVES,points)
        errors=[];momentum=[];wave_sizes=[]
        for x,E in zip(points,frames):
            h,hd,hdd = A+B@x,Ad+Bd@x,Add+Bdd@x
            dh,dhd = E@B/a,E@Bd/a
            Hh = -(B@x)*np.eye(3)/a**2
            phi = 1/h
            derivative = np.r_[-hd/h**2,-dh/h**2]
            hessian = np.zeros((4,4))
            hessian[0,0] = 2*hd*hd/h**3-hdd/h**2
            hessian[0,1:] = hessian[1:,0] = 2*hd*dh/h**3-dhd/h**2
            hessian[1:,1:] = 2*np.outer(dh,dh)/h**3-Hh/h**2
            G = np.diag([3.,-1.,-1.,-1.])/a**2
            T = support.improved_stress(eta,eta,G,phi,derivative,hessian)
            wave = np.sum(eta*hessian)-phi/a**2
            errors.append(abs(h**3*wave-(C+L@x+x@M@x)))
            momentum.append(np.linalg.norm(T[0,1:]-dhd/(3*h**3)))
            wave_sizes.append(abs(wave))
        rows.append(dict(radius=a,wave_expansion_absolute=float(max(errors)),
                         wave_scale=max(1.,abs(C),float(np.linalg.norm(L)),float(np.linalg.norm(M))),
                         momentum_absolute=float(max(momentum)),
                         non_solution_wave_size=float(max(wave_sizes))))
    return rows


def rank_controls():
    rng=np.random.default_rng(SEED+1)
    errors=[];ranks=[]
    for _ in range(40):
        B,V,W=rng.normal(size=(3,4))
        M=coefficient_conditions(1.,0.,0.,B,V,W)[2]
        ranks.append(int(np.linalg.matrix_rank(M,tol=1e-10)))
        u=rng.normal()
        balanced=coefficient_conditions(1.,0.,0.,B,u*B,2*u*u*B)[2]
        errors.append(np.linalg.norm(balanced)/max(1.,np.linalg.norm(B)**2*u*u))
    e=np.eye(4)
    instant=coefficient_conditions(1.,0.,0.,np.zeros(4),np.zeros(4),e[0])[2]
    transverse_velocity=coefficient_conditions(1.,0.,0.,e[0],e[1],np.zeros(4))[2]
    transverse_acceleration=coefficient_conditions(1.,0.,0.,e[0],np.zeros(4),e[1])[2]
    mu_control=coefficient_conditions(1.,0.,0.,e[0],e[1],e[2])[2]
    return dict(max_rank=max(ranks),min_generic_rank=min(ranks),alignment_max_relative=float(max(errors)),
                zero_B_nonzero_Bdd_matrix_norm=float(np.linalg.norm(instant)),
                zero_B_control_Bdd=e[0].tolist(),
                transverse_velocity_matrix_norm=float(np.linalg.norm(transverse_velocity)),
                transverse_acceleration_matrix_norm=float(np.linalg.norm(transverse_acceleration)),
                trace_matched_identity_residual=float(np.linalg.norm(mu_control-np.trace(mu_control)*np.eye(4)/4)),
                rank_step_is_analytic="samples are controls; rectangular factorization proves the bound")


def failed_checks(checks):
    return [name for name in REQUIRED_CHECKS if checks.get(name) is not True]


def verdict(checks):
    failed=failed_checks(checks)
    return dict(uniqueness="UNRESOLVED" if failed else "ONLY_HOMOGENEOUS_SUPPORT_IN_STATED_CLASS",
                failed_checks=failed,post_freeze=True,
                BAM_support_selection="NOT_DERIVED",Phi_selection="NOT_DERIVED",causality_gate="OPEN")
