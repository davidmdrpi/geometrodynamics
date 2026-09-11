"""Round 12: independent odd conformal components on the exact round ESU.

Fixed-background source responses are not coupled metric evolution. Exact
polynomial certificates precede numerical whitening and singular values.
"""
from dataclasses import dataclass
from functools import lru_cache
from fractions import Fraction
import math

import numpy as np

from . import reciprocal_scalar_tt as rt

PUBLIC_PREREG = "a12cb2f96d846bc3208931e02e8f7f6b59a05d21"
REVIEW_AMENDMENT = "9c33867b216ea7e76b436720c9163527314c4234"
SEED = 2026091012
PHYSICAL_VERDICTS = ("background_existence", "kinetic_regularity",
                     "component_count_bound", "diagonal_sensitivity", "full_preparation_kernel")
_B, _K, _M, _D, _P = PHYSICAL_VERDICTS
CHECK_TARGETS = {
    "harmonic_parity": (_B, _M, _D, _P),
    "addition_identities": PHYSICAL_VERDICTS,
    "full_stress_background": PHYSICAL_VERDICTS,
    "einstein_normalization": PHYSICAL_VERDICTS,
    "kinetic_matrix": (_K,),
    "component_bound": (_M,),
    "independent_field_controls": (_B, _M, _D, _P),
    "exact_preparation_kernel": (_P,),
    "normalized_sensitivity": (_D, _P),
    "gram_reconstruction": (_D, _P),
    "scope_and_order": PHYSICAL_VERDICTS,
    "failure_paths": PHYSICAL_VERDICTS,
}
REQUIRED_CHECKS = tuple(CHECK_TARGETS)
_EXPECTED_TARGETS = dict(CHECK_TARGETS)
SCOPE = dict(coupled_dynamical_stability="NOT_ESTABLISHED",
             field_content_preparation_selection="NOT_DERIVED",
             TT_frequency_transfer="NOT_ESTABLISHED", coupled_support_response="NOT_DERIVED",
             Phi_selection="NOT_DERIVED", causality_gate="OPEN")


def verdict(checks, kernel_results=None):
    """A complete gate set is necessary; affirmative kernel evidence is separate."""
    result = dict(zip(PHYSICAL_VERDICTS, (
        "ODD_MULTIPLET_EXACT_ESU_SUPPORT", "POSITIVE_FULL_KINETIC_MATRIX",
        "FOUR_MINIMAL_IN_COMMON_PHASE_SINGLE_DEGREE_CLASS",
        "FINITE_FIXED_BACKGROUND_SUSCEPTIBILITY", "UNRESOLVED")))
    failed = [key for key in REQUIRED_CHECKS if checks.get(key) is not True]
    valid_map = CHECK_TARGETS == _EXPECTED_TARGETS and tuple(CHECK_TARGETS) == REQUIRED_CHECKS
    if not valid_map or set(checks) != set(REQUIRED_CHECKS):
        result.update(dict.fromkeys(PHYSICAL_VERDICTS, "UNRESOLVED"))
        failed.append("gate_schema")
    if kernel_results is not None and set(kernel_results) == {1, 3, 5}:
        if all(r.get("certified") is True and isinstance(r.get("nullity"), int)
               and r["nullity"] >= 0 for r in kernel_results.values()):
            result[_P] = {str(k): ("EXACT_EQUAL_STRESS_FAMILY" if r["nullity"]
                                   else "ISOLATED_FIXED_TRACE_GRAM")
                          for k, r in kernel_results.items()}
    if result[_P] == "UNRESOLVED":
        failed.append("kernel_evidence")
    for key in failed:
        for target in CHECK_TARGETS.get(key, ()):
            result[target] = "UNRESOLVED"
    if not valid_map or "gate_schema" in failed:
        result.update(dict.fromkeys(PHYSICAL_VERDICTS, "UNRESOLVED"))
    return {**result, **SCOPE, "failed_checks": sorted(set(failed))}


@dataclass(frozen=True)
class MultipletSupport:
    degree: int = 1
    radius: float = 1.
    kappa: float = 1.

    def __post_init__(self):
        if self.degree not in (1, 3, 5):
            raise ValueError("odd support degrees are 1, 3, 5")
        if not np.isfinite([self.radius, self.kappa]).all() or min(self.radius, self.kappa) <= 0:
            raise ValueError("finite positive radius and kappa required")

    @property
    def size(self): return (self.degree+1)**2
    @property
    def omega(self): return (self.degree+1)/self.radius
    @property
    def amplitude2(self): return 3/(self.kappa*self.size**2)
    @property
    def density(self): return 3/(2*self.kappa*self.radius**2)
    @property
    def cosmological_constant(self): return 3/(2*self.radius**2)
    @property
    def volume(self): return 2*math.pi**2*self.radius**3

    def kinetic(self, phi):
        phi = np.asarray(phi, dtype=float)
        f = 1-self.kappa*np.sum(phi*phi, axis=-1)/6
        if np.any(f <= 0): raise ValueError("nonpositive gravitational coefficient")
        K = np.eye(phi.shape[-1])/f[..., None, None]
        K = K+self.kappa*phi[..., :, None]*phi[..., None, :]/(6*f[..., None, None]**2)
        return f/self.kappa, K


@lru_cache(None)
def exact_moment(exponents):
    if any(e % 2 for e in exponents): return Fraction(0)
    half = [e//2 for e in exponents]
    return Fraction(math.prod(math.prod(range(1, 2*k, 2)) for k in half),
                    math.prod(range(4, 4+2*sum(half), 2)))


def generator(degree, S):
    ex = rt.powers(degree)
    index = {e: i for i, e in enumerate(ex)}
    out = np.zeros((len(ex), len(ex)), dtype=np.int64)
    for j, e in enumerate(ex):
        for r, count in enumerate(e):
            if count:
                for s in range(4):
                    if S[r, s]:
                        target = list(e); target[r] -= 1; target[s] += 1
                        out[index[tuple(target)], j] += count*int(S[r, s])
    return out


@lru_cache(None)
def rational_basis(k):
    """Integer harmonic basis from an exact Euclidean Laplacian nullspace."""
    import sympy as sp
    ex = rt.powers(k); low = rt.powers(k-2) if k >= 2 else ()
    index = {e: i for i, e in enumerate(low)}
    lap = sp.zeros(len(low), len(ex))
    for j, e in enumerate(ex):
        for a, n in enumerate(e):
            if n >= 2:
                v = list(e); v[a] -= 2
                lap[index[tuple(v)], j] += n*(n-1)
    vectors = lap.nullspace()
    B = sp.Matrix.hstack(*(v*sp.ilcm(*[c.q for c in v]) for v in vectors))
    if B.cols != (k+1)**2 or lap*B != sp.zeros(lap.rows, B.cols):
        raise ArithmeticError("exact harmonic basis failed")
    return np.array(B, dtype=np.int64)


def product_matrix(A, B, degree):
    """Coefficients of sum_ij H_ij A_i B_j for independent upper entries H."""
    ex = rt.powers(degree); target = {e: i for i, e in enumerate(rt.powers(2*degree))}
    n = A.shape[1]; pairs = [(i,j) for i in range(n) for j in range(i,n)]
    out = np.zeros((len(target), len(pairs)), dtype=np.int64)
    sparse_a = [[(ex[r], int(A[r,i])) for r in np.flatnonzero(A[:,i])] for i in range(n)]
    sparse_b = [[(ex[r], int(B[r,i])) for r in np.flatnonzero(B[:,i])] for i in range(n)]
    for c,(i,j) in enumerate(pairs):
        for u,v in ((i,j),) if i==j else ((i,j),(j,i)):
            for x,a in sparse_a[u]:
                for y,b in sparse_b[v]:
                    out[target[tuple(p+q for p,q in zip(x,y))],c] += a*b
    return out


@lru_cache(None)
def polynomial_operators(k):
    """Z and 12*A in integer coefficients; homogeneous degree 2k is complete."""
    B = rational_basis(k)
    small = [generator(k,S) for S in rt.QUATERNION_DERIVATIVES]
    large = [generator(2*k,S) for S in rt.QUATERNION_DERIVATIVES]
    dB = [M@B for M in small]
    Z = product_matrix(B,B,k)
    tensor = {}
    for i in range(3):
        for j in range(i,3):
            # Symmetric spatial Hessian in the invariant orthonormal frame.
            U = product_matrix(dB[i],dB[j],k)
            tensor[i,j] = 12*U-(large[i]@large[j]+large[j]@large[i])@Z
    A = np.vstack([tensor[0,0]-tensor[2,2],tensor[1,1]-tensor[2,2],
                   tensor[0,1],tensor[0,2],tensor[1,2]])
    # Haar trace functional: mean Z. Scale to primitive integer coefficients.
    moment = [exact_moment(e) for e in rt.powers(2*k)]
    trace = [sum((moment[r]*int(Z[r,c]) for r in range(len(moment))),Fraction(0))
             for c in range(Z.shape[1])]
    den = math.lcm(*(v.denominator for v in trace))
    tr = np.array([int(v*den) for v in trace],dtype=np.int64)
    tr //= math.gcd(*map(int,tr))
    return dict(B=B,Z=Z,A=A,trace=tr,large=np.array(large),tensor=tensor)


def _exact_nullspace(matrix):
    """Fraction-free exact row reduction; returned witnesses are integer columns."""
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix
    dm = DomainMatrix.from_Matrix(sp.Matrix(matrix)).convert_to(sp.ZZ)
    reduced, denom, pivots = dm.rref_den(method="FF")
    null = reduced.nullspace_from_rref(pivots).to_Matrix().T
    columns = []
    for j in range(null.cols):
        v = null[:,j]
        gcd = math.gcd(*(int(x) for x in v))
        columns.append([int(x)//gcd for x in v])
    return len(pivots), columns, list(pivots)


def modular_rank(matrix, prime=1000003):
    """Exact finite-field lower bound for rational rank; never an upper bound."""
    M = np.asarray(matrix,dtype=np.int64).copy() % prime
    rank = 0
    for col in range(M.shape[1]):
        candidates = np.flatnonzero(M[rank:,col])
        if not len(candidates): continue
        pivot = rank+int(candidates[0]); M[[rank,pivot]] = M[[pivot,rank]]
        M[rank,col:] = M[rank,col:]*pow(int(M[rank,col]),-1,prime)%prime
        idx = np.flatnonzero(M[rank+1:,col])+rank+1
        # int64 safe: both factors below prime.
        M[idx,col:] = (M[idx,col:]-M[idx,col,None]*M[rank,col:])%prime
        rank += 1
        if rank==M.shape[0]: break
    return rank


def verify_certificate(matrix, certificate):
    """Rank lower bound + exact independent kernel => equality over Q."""
    n = matrix.shape[1]
    rank = certificate.get("rank")
    if type(rank) is not int or not 0 <= rank <= n:
        return False
    vectors = certificate.get("kernel_columns",[])
    if any(not isinstance(x,(int,np.integer)) or isinstance(x,(bool,np.bool_))
           for row in vectors for x in row):
        return False
    W = np.array(vectors,dtype=object).T if vectors else np.empty((n,0),dtype=object)
    if W.shape != (n,n-certificate.get("rank",-1)):
        return False
    if "nullity" in certificate and certificate["nullity"] != W.shape[1]:
        return False
    if "fixed_trace_rank" in certificate and certificate["fixed_trace_rank"] != certificate["rank"]-1:
        return False
    # A proved absolute bound permits fast int64 arithmetic without overflow;
    # larger external certificates retain arbitrary-precision verification.
    bound = (max((sum(abs(int(x)) for x in row) for row in matrix),default=0)
             *max((abs(int(x)) for x in W.flat),default=0))
    if bound <= np.iinfo(np.int64).max:
        exact_zero = not np.any(np.asarray(matrix,dtype=np.int64)@np.asarray(W,dtype=np.int64))
    else:
        exact_zero = not np.any(np.asarray(matrix,dtype=object)@W)
    independent = modular_rank(np.asarray(W % 1000003,dtype=np.int64))==W.shape[1]
    lower = modular_rank(matrix)
    return bool(exact_zero and independent and lower==certificate["rank"])


def exact_kernel_certificate(k):
    op = polynomial_operators(k)
    answer = {"degree":k,"input_dimension":op["Z"].shape[1]-1,"operators":{}}
    for name,M in (("A",op["A"]),("L",np.vstack([op["Z"],op["A"]]))):
        constrained = np.vstack([op["trace"],M])
        rank, columns, pivots = _exact_nullspace(constrained)
        cert = dict(rank=rank, fixed_trace_rank=rank-1, nullity=constrained.shape[1]-rank,
                    kernel_columns=columns, pivots=pivots, prime=1000003)
        cert["certified"] = verify_certificate(constrained,cert)
        if not cert["certified"]: raise ArithmeticError("exact rank certificate failed")
        answer["operators"][name] = cert
    return answer


@lru_cache(None)
def numerical_basis(k):
    B = rational_basis(k).astype(float)
    ex = rt.powers(k)
    M = np.array([[float(exact_moment(tuple(x+y for x,y in zip(a,b)))) for b in ex] for a in ex])
    moments = B.T@M@B
    val,vec = np.linalg.eigh(moments)
    W = (vec/np.sqrt(val))@vec.T
    return B@W,W,moments


def basis_jets(k, points):
    """Haar basis, invariant derivatives and symmetrized Hessian, unit radius."""
    B,_,_ = numerical_basis(k)
    mon = rt.monomials(points,rt.powers(k))
    operators = [generator(k,S) for S in rt.QUATERNION_DERIVATIVES]
    Y = mon@B
    grad = np.array([mon@M@B for M in operators]).transpose(1,2,0)
    hess = np.empty((len(points),B.shape[1],3,3))
    for i in range(3):
        for j in range(3):
            hess[:,:,i,j] = mon@((operators[i]@operators[j]+operators[j]@operators[i])/2)@B
    return Y,grad,hess


def component_stress(model, phase, points, C=None, jets=None):
    """Direct sum of off-shell improved stresses; arbitrary real component matrix.

    C maps Haar modes to component fields; default is the ESU preparation.
    Includes all physical-radius derivatives and xi*G*phi^2.
    """
    Y,dY,hY = basis_jets(model.degree,points) if jets is None else jets
    if C is None: C = np.eye(model.size)*math.sqrt(model.amplitude2)
    C = np.asarray(C,dtype=float)
    f = Y@C.T; df = np.einsum("pni,mn->pmi",dY,C)/model.radius
    hf = np.einsum("pnij,mn->pmij",hY,C)/model.radius**2
    c,s = math.cos(phase),math.sin(phase)
    phi = c*f; derivative = np.concatenate((-model.omega*s*f[:,:,None],c*df),axis=2)
    hessian = np.empty((*phi.shape,4,4))
    hessian[:,:,0,0] = -model.omega**2*phi
    hessian[:,:,0,1:] = -model.omega*s*df
    hessian[:,:,1:,0] = -model.omega*s*df
    hessian[:,:,1:,1:] = c*hf
    eta = np.diag([-1.,1.,1.,1.]); E = np.diag([3.,-1.,-1.,-1.])/model.radius**2
    # Vectorized version of improved_stress, checked against the scalar route.
    square = np.einsum("pmi,ij,pmj->pm",derivative,eta,derivative)
    box = np.einsum("ij,pmij->pm",eta,hessian)
    outer = derivative[:,:,:,None]*derivative[:,:,None,:]
    stresses = outer-eta*square[:,:,None,None]/2+(eta*2*(square+phi*box)[:,:,None,None]
        -2*(outer+phi[:,:,None,None]*hessian)+E*phi[:,:,None,None]**2)/6
    return stresses,dict(phi=phi,derivative=derivative,hessian=hessian,box=box,
                         einstein=E,metric=eta)


def gram_stress(model, phase, points, G, jets=None):
    """Independent polynomial reduction, with G in the Haar basis."""
    Y,dY,hY = basis_jets(model.degree,points) if jets is None else jets
    GY = Y@G
    GdY = np.einsum("nm,pmi->pni",G,dY,optimize=True)
    Z = np.einsum("pn,pn->p",Y,GY)
    grad = 2*np.einsum("pni,pn->pi",dY,GY)/model.radius
    U = np.einsum("pni,pnj->pij",dY,GdY)/model.radius**2
    hess = 2*(U+np.einsum("pnij,pn->pij",hY,GY)/model.radius**2)
    lap = np.trace(hess,axis1=1,axis2=2)
    c,s = math.cos(phase),math.sin(phase)
    rho = model.omega**2*Z/2+c*c*lap/12
    T = np.zeros((len(points),4,4)); T[:,0,0] = rho
    T[:,0,1:] = T[:,1:,0] = -model.omega*s*c*grad/6
    Q = U-hess/6
    Q -= np.trace(Q,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
    T[:,1:,1:] = rho[:,None,None]*np.eye(3)/3+c*c*Q
    return T


def factor_gram(G):
    G = np.asarray(G,dtype=float)
    if not np.isfinite(G).all() or np.linalg.norm(G-G.T)>1e-12*max(1.,np.linalg.norm(G)):
        raise ValueError("finite symmetric Gram required")
    val,vec = np.linalg.eigh(G)
    if val.min() < -1e-13*max(1.,np.linalg.norm(G)):
        raise ValueError("positive semidefinite Gram required")
    return np.sqrt(np.maximum(val,0))[:,None]*vec.T


def symmetric_directions(n, fixed_trace=True):
    """Frobenius-orthonormal symmetric basis; sqrt(n) realizes frozen norm."""
    basis = []
    if fixed_trace:
        for j in range(1,n):
            d = np.zeros(n); d[:j]=1; d[j]=-j
            basis.append(np.diag(d)/math.sqrt(j*(j+1)))
    else:
        basis.extend(np.diag(row) for row in np.eye(n))
    for i in range(n):
        for j in range(i+1,n):
            E = np.zeros((n,n)); E[i,j]=E[j,i]=1/math.sqrt(2)
            basis.append(E)
    return np.array(basis)*math.sqrt(n)


def sensitivity_operator(k):
    """Exact spatial moments and exact time averages; numerical input whitening."""
    op = polynomial_operators(k); n=(k+1)**2; _,W,_=numerical_basis(k)
    H = symmetric_directions(n)
    raw = np.einsum("ai,dij,bj->dab",W,H,W)
    pairs = [(i,j) for i in range(n) for j in range(i,n)]
    transform = np.array([raw[:,i,j] for i,j in pairs])
    Z = op["Z"]@transform
    D = op["large"].astype(float); lap = sum(M@M for M in D)
    # Divide by rho0 = n*s^2*omega^2/2; s^2 and a^-2 cancel.
    omega = k+1; scale = n*omega**2/2
    R0 = (omega**2*Z/2+lap@Z/24)/scale
    Rc = (lap@Z/24)/scale
    J = np.array([omega*M@Z/(12*scale) for M in D])
    # 12 times unprojected Q; explicit tensor trace subtraction.
    Q = np.array([op["tensor"][i,j]@transform/(12*scale)
                  for i,j in ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))])
    tr = sum(Q[:3])/3; Q[:3] -= tr
    ex=rt.powers(2*k)
    moments = np.array([[float(exact_moment(tuple(x+y for x,y in zip(a,b)))) for b in ex] for a in ex])
    def square(X): return X.T@moments@X
    anis = sum(square(q)*(1 if i<3 else 2) for i,q in enumerate(Q))*3/8
    full = anis+(square(R0)+square(Rc)/2)*4/3+sum(square(j) for j in J)
    result={"directions":H,"matrices":{"A":(anis+anis.T)/2,"L":(full+full.T)/2}}
    return result


def singular_analysis(k, certificate):
    operators=sensitivity_operator(k)
    result={}
    for name,M in operators["matrices"].items():
        vals,vecs=np.linalg.eigh(M)
        nullity=certificate["operators"][name]["nullity"]
        positive=vals[nullity:]
        if len(positive) and positive.min()<=0: raise ArithmeticError("positive spectrum unresolved")
        result[name]=dict(rank=len(positive),nullity=nullity,
            singular_values=[0.]*nullity+np.sqrt(positive).tolist(),
            zero_value_provenance="exact rational kernel; positive singular values from numerical moment Gram matrix",
            singular_reporting_cutoff=float(1e-10*np.sqrt(positive.max())),
            smallest_positive=float(np.sqrt(positive.min())),largest=float(np.sqrt(positive.max())),
            numerical_kernel_gram_residual=float(np.max(np.abs(vals[:nullity]))) if nullity else 0.,
            min_direction=np.einsum("d,dij->ij",vecs[:,nullity],operators["directions"]),
            max_direction=np.einsum("d,dij->ij",vecs[:,-1],operators["directions"]))
    return result,operators


def exact_addition_certificate(k):
    """Check addition identities over Q, without numerical whitening."""
    import sympy as sp
    op=polynomial_operators(k); B=sp.Matrix(op["B"]); n=B.cols
    ex=rt.powers(k)
    moment=sp.Matrix([[sp.Rational(exact_moment(tuple(x+y for x,y in zip(a,b))))
                       for b in ex] for a in ex])
    metric=B.T*moment*B
    inverse=metric.inv()
    pair=sp.Matrix([inverse[i,j] for i in range(n) for j in range(i,n)])
    expected=[]
    for e in rt.powers(2*k):
        expected.append(0 if any(v%2 for v in e) else
                        math.factorial(k)//math.prod(math.factorial(v//2) for v in e))
    radial=sp.Matrix(expected)
    z=sp.Matrix(op["Z"])*pair
    residuals={"sum_squares": z-n*radial == sp.zeros(len(radial),1),
               "odd_parity": k%2==1,
               "dimension":n==(k+1)**2}
    for i,S in enumerate(rt.QUATERNION_DERIVATIVES):
        residuals[f"value_gradient_{i}"]=sp.Matrix(op["large"][i])*z==sp.zeros(len(radial),1)
        for j,T in enumerate(rt.QUATERNION_DERIVATIVES):
            U=product_matrix(generator(k,S)@op["B"],generator(k,T)@op["B"],k)
            target=n*k*(k+2)*radial/3 if i==j else sp.zeros(len(radial),1)
            residuals[f"gradient_{i}_{j}"]=sp.Matrix(U)*pair==target
    return dict(degree=k, residuals={key:"0" if v else "NONZERO" for key,v in residuals.items()},
                all_zero=all(residuals.values()))


def exact_background_certificate():
    """All stress/Einstein and kinetic eigenvalue identities after addition."""
    import sympy as sp
    n,a,kap,c=sp.symbols("N a kappa c",positive=True)
    # n=(k+1)^2 and lambda=n-1. The addition identities fix all summed jets.
    amp=3/(kap*n*n); omega2=n/a**2
    phi2=n*amp*c*c; dt2=n*amp*omega2*(1-c*c)
    spatial=n*amp*(n-1)*c*c/(3*a*a)
    phi_dtt=-omega2*phi2
    phi_hess=-spatial
    # Hess(sum phi_i^2) spatial=0, but individual Hess terms do not vanish.
    square=-dt2+3*spatial
    phi_box=-phi_dtt+3*phi_hess
    rho=dt2+square/2+(-2*(square+phi_box)-2*(dt2+phi_dtt)+3*phi2/a**2)/6
    p=spatial-square/2+(2*(square+phi_box)-2*(spatial+phi_hess)-phi2/a**2)/6
    target=3/(2*kap*a*a); Lambda=3/(2*a*a)
    f=1-kap*phi2/6
    residuals=dict(rho=rho-target,pressure=p-target/3,
                   einstein_time=3/a**2-Lambda-kap*rho,
                   einstein_space=-1/a**2+Lambda-kap*p,
                   radial_kinetic=1/f+kap*phi2/(6*f*f)-1/(f*f),
                   f= f-(1-c*c/(2*n)))
    residuals={key:str(sp.factor(v)) for key,v in residuals.items()}
    return dict(residuals=residuals,all_zero=set(residuals.values())=={"0"},
                minimum_components=4,scope="common-phase single-degree exact pointwise support",
                isotropy_phase="cosine nonzero; constant-norm inference uses field-zero phase")


def raw_kernel_to_haar(k, column):
    n=(k+1)**2; G=np.zeros((n,n)); index=0
    for i in range(n):
        for j in range(i,n):
            G[i,j]=G[j,i]=float(column[index]); index+=1
    _,W,_=numerical_basis(k)
    inverse=np.linalg.inv(W)
    H=inverse@G@inverse.T
    return H/(np.linalg.norm(H)/math.sqrt(n))


def source_norm(T, weights, density=1., anisotropic=False):
    if anisotropic:
        Q=T[...,1:,1:].copy()
        Q-=np.trace(Q,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)/3
        square=np.sum(Q*Q,axis=(-2,-1))
    else:
        square=np.sum(T*T,axis=(-2,-1))
    return float(np.sqrt(max(0.,np.asarray(weights)@square))/density)


def exact_divergence_certificate(k):
    """Post-freeze explanation of ker(A)=ker(L), not a frozen prediction.

    div A = -grad(Delta Z + 12*(k+1)^2 Z)/36. Since Z has degrees <=2k,
    the operator in parentheses has no nonconstant kernel. Fixed mean Z=0
    makes A=0 imply Z=0. The converse is automatic from L=(Z,A).
    """
    op=polynomial_operators(k); D=op['large']; tensor=op['tensor']
    trace=sum(tensor[i,i] for i in range(3))
    lap=sum(M@M for M in D)
    residuals=[]
    for j in range(3):
        divergence=sum(D[i]@(3*tensor[min(i,j),max(i,j)]-(trace if i==j else 0)) for i in range(3))
        residual=divergence+D[j]@(lap+12*(k+1)**2*np.eye(lap.shape[0],dtype=np.int64))@op['Z']
        residuals.append(int(np.max(abs(residual))))
    gap=12*(k+1)**2-2*k*(2*k+2)
    return dict(degree=k,residuals=residuals,minimum_spectral_gap=gap,
                certified=all(r==0 for r in residuals) and gap>0,post_freeze=True)
