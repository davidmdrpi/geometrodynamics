"""Local Lorentz gluing and real-linear history closure diagnostics.

Spatial identification, time orientation and an optional internal charge
conjugation are independent data. These identities do not construct a global
Einstein metric or supply the stress required to support a traversable handle.
"""
import numpy as np

METRIC = np.diag([-1.,1.,1.,1.])


def boost(velocity):
    v=np.asarray(velocity,float)
    if v.shape!=(3,) or not np.isfinite(v).all() or v@v>=1:
        raise ValueError('Finite subluminal three-velocity required')
    vv=v@v; gamma=1/np.sqrt(1-vv)
    result=np.eye(4); result[0,0]=gamma
    result[0,1:]=result[1:,0]=-gamma*v
    if vv>0:
        result[1:,1:]+=(gamma-1)*np.outer(v,v)/vv
    return result


def gluing_jacobian(spatial, time_orientation=1, source_velocity=(0.,0.,0.),
                     destination_velocity=(0.,0.,0.)):
    """J maps destination coordinate vectors to source coordinate vectors.

    The spatial map includes any inward/outward normal reversal at the seam;
    it must not be identified with a 2D rim-chart map alone.
    """
    r=np.asarray(spatial,float)
    if (r.shape!=(3,3) or not np.isfinite(r).all() or
            not np.allclose(r.T@r,np.eye(3),atol=1e-12,rtol=0) or time_orientation not in (-1,1)):
        raise ValueError('Orthogonal spatial map and time orientation +/-1 required')
    d=np.zeros((4,4));d[0,0]=time_orientation;d[1:,1:]=r
    return boost(-np.asarray(source_velocity))@d@boost(destination_velocity)


def faraday(electric, magnetic):
    e,b=np.asarray(electric,float),np.asarray(magnetic,float)
    if e.shape!=(3,) or b.shape!=(3,) or not np.isfinite([e,b]).all():
        raise ValueError('Finite three-component fields required')
    f=np.zeros((4,4));f[0,1:]=e;f[1:,0]=-e
    bx,by,bz=b
    f[1:,1:]=[[0,-bz,by],[bz,0,-bx],[-by,bx,0]]
    return f


def fields(f):
    return f[0,1:].copy(),np.array([-f[2,3],f[1,3],-f[1,2]])


def stress(f):
    raised=METRIC@f@METRIC
    return f@METRIC@f.T-.25*METRIC*np.sum(f*raised)


def transport_faraday(f,jacobian,charge_conjugation=1):
    if charge_conjugation not in (-1,1):
        raise ValueError('Charge conjugation must be +/-1')
    return charge_conjugation*jacobian.T@f@jacobian


def frequency_map(spatial, omega, clock_offset=0., time_orientation=1, charge_conjugation=1):
    """Electric phasor pullback for t_source=epsilon*t_destination+offset.

    Real field convention E(t)=Re[a exp(-i omega t)]. A reversed clock maps
    positive frequency to its conjugate. Applicable to static seam frames.
    """
    j=gluing_jacobian(spatial,time_orientation)
    if omega<=0 or not np.isfinite([omega,clock_offset]).all() or charge_conjugation not in (-1,1):
        raise ValueError('Positive finite frequency and valid parity required')
    m=charge_conjugation*time_orientation*j[1:,1:].T
    zero=np.zeros((3,3),complex)
    if time_orientation==1:
        return m*np.exp(-1j*omega*clock_offset),zero
    return zero,m*np.exp(1j*omega*clock_offset)


def compose_maps(first, second):
    """Compose (U,V): a -> U a + V conjugate(a), first after second."""
    u1,v1=first;u2,v2=second
    return u1@u2+v1@v2.conj(),u1@v2+v1@u2.conj()


def solve_history(u,v,source,tolerance=1e-11):
    """Solve a=U a+V conjugate(a)+b by real SVD, retaining compatibility.

    Closed lossless loops can have no driven solution or undetermined modes.
    Least-squares residual is NEVER reported as an exact compatible history.
    The returned solution is the minimum-norm member when nonunique.
    """
    u,v,b=np.asarray(u,complex),np.asarray(v,complex),np.asarray(source,complex)
    n=len(b)
    if (u.shape!=(n,n) or v.shape!=(n,n) or b.shape!=(n,) or tolerance<=0
            or not all(np.isfinite(z).all() for z in [u,v,b])):
        raise ValueError('Finite matching linear maps and source required')
    linear=np.block([[u.real+v.real,-u.imag+v.imag],
                     [u.imag+v.imag,u.real-v.real]])
    matrix=np.eye(2*n)-linear
    rhs=np.r_[b.real,b.imag]
    left,singular,right=np.linalg.svd(matrix)
    threshold=tolerance*max(1.,singular[0])
    rank=int(np.sum(singular>threshold))
    coeff=left.T@rhs
    y=right[:rank].T@(coeff[:rank]/singular[:rank]) if rank else np.zeros(2*n)
    residual=float(np.linalg.norm(matrix@y-rhs))
    compatible=residual<=tolerance*max(1.,np.linalg.norm(rhs))
    solution=y[:n]+1j*y[n:]
    return dict(solution=solution,residual=residual,rank=rank,nullity=2*n-rank,
                singular_values=singular,
                status='incompatible' if not compatible else ('unique' if rank==2*n else 'nonunique'))


def time_deck(values):
    """Internal deck map plus a half-period shift on an even spatial cover.

    Quotient identification is (t,x,u) ~ (-t,x+L,D u), D(phi,a,b)=(-phi,b,a).
    This free isometry of flat 1+1 spacetime defines a time-nonorientable
    benchmark. It does not assert a supported four/five-dimensional GR handle.
    """
    u=np.asarray(values)
    if u.ndim!=2 or u.shape[0]!=3 or u.shape[1]%2:
        raise ValueError('Three fields on an even spatial cover required')
    shifted=np.roll(u,u.shape[1]//2,axis=1)
    return np.array([-shifted[0],shifted[2],shifted[1]])


def time_reversing_initial(cell_fields,cell_velocities):
    """Initial data satisfying u0=G u0 and v0=-G v0.

    Evolution by a deck-invariant reversible action then obeys u(-t)=G u(t)
    by uniqueness. Values at the end of the cell must satisfy the spatial
    deck matching; the localized packet construction enforces this to its tails.
    """
    u,v=np.asarray(cell_fields),np.asarray(cell_velocities)
    if u.shape!=v.shape or u.ndim!=2 or u.shape[0]!=3:
        raise ValueError('Matching three-field cell data required')
    d=lambda x:np.array([-x[0],x[2],x[1]])
    return np.concatenate([u,d(u)],axis=1),np.concatenate([v,-d(v)],axis=1)
