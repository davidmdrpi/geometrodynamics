"""Existing conformal quartet on a localized compact handle, freeze 66bce68.

Initial constraints only. The collar and sign line bundle are free data.
"""
from functools import lru_cache
import numpy as np
from scipy.integrate import solve_ivp, solve_bvp
from scipy.interpolate import CubicHermiteSpline, PPoly
from . import mouth_momentum as vacuum

PREREG='66bce68ce5c9e07f4369283b18cb1ac4d3e94f73'
LENGTHS=(3.5,4.5,5.5)
AMPLITUDES=(0.,.1,.3)
SCHEDULE=((129,1e-6),(257,1e-8),(513,1e-10))
Q=np.sqrt(3)/2
F=7/8
ALPHA=6/7
LAMBDA=1.5
POTENTIAL=LAMBDA/F**2


def sech(x):return 1/np.cosh(x)


def cutoff(s,L):
    z=2*(np.asarray(s)-L+1)
    # Stable evaluation of the C-infinity step; never divide by zero.
    a=np.zeros_like(z,dtype=float);b=np.zeros_like(z,dtype=float)
    with np.errstate(divide='ignore',over='ignore'):
        np.exp(-1/z,out=a,where=z>0)
        np.exp(-1/(1-z),out=b,where=z<1)
    return np.where(z<=0,1.,np.where(z>=1,0.,b/(a+b)))


def theta_prime(s,L):return sech(s)*cutoff(s,L)


def momentum(s,L,eta):return eta*sech(L)**3*np.sin(np.pi*np.asarray(s)/(2*L))


def serialize(poly):return dict(x=poly.x.tolist(),c=poly.c.tolist(),axis=poly.axis)
def restore(record):
    result=PPoly(np.asarray(record['c']),np.asarray(record['x']))
    result.axis=record['axis']
    return result


@lru_cache(None)
def profile(L):
    nodes=np.linspace(0,L,4097)
    t=solve_ivp(lambda s,y:[float(theta_prime(s,L))],(0,L),[0.],method='DOP853',
                rtol=1e-12,atol=1e-14,max_step=.01,dense_output=True)
    a=solve_ivp(lambda s,y:[-ALPHA*momentum(s,L,1.)*theta_prime(s,L)/2],(L,0),
                [sech(L)**3/8],method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01,dense_output=True)
    if not(t.success and a.success):raise ArithmeticError('profile integration failed')
    theta=CubicHermiteSpline(nodes,t.sol(nodes)[0],theta_prime(nodes,L))
    tensor=CubicHermiteSpline(nodes,a.sol(nodes)[0],-ALPHA*momentum(nodes,L,1.)*theta_prime(nodes,L)/2)
    return theta,tensor


def solve(L,eta,n,tolerance):
    theta,tensor=profile(L)
    def rhs(s,y):
        p=momentum(s,L,eta);a=eta*tensor(s);tp=theta_prime(s,L)
        linear=(2-ALPHA*(tp*tp+2*np.cos(theta(s))**2))/8
        source=(6*a*a+ALPHA*p*p)/8
        with np.errstate(over='ignore',divide='ignore',invalid='ignore'):
            dd=linear*y[0]-source*y[0]**-7-POTENTIAL*y[0]**5/4
        return np.array([y[1],dd])
    x=np.linspace(0,L,n)
    root=np.sqrt(sech(x));image=np.sqrt(sech(2*L-x))
    y=F**.25*np.array([root+image,-.5*np.tanh(x)*root+.5*np.tanh(2*L-x)*image])
    solution=solve_bvp(rhs,lambda ya,yb:np.array([ya[1],yb[1]]),x,y,tol=tolerance,max_nodes=40000)
    return dict(L=L,eta=eta,n_initial=n,tolerance=tolerance,status=int(solution.status),
                success=bool(solution.success),message=solution.message,nodes=len(solution.x),
                iterations=int(solution.niter),solution=serialize(solution.sol),
                rms_residuals=solution.rms_residuals.tolist())


class Data:
    def __init__(self,record,profiles=None):
        self.record=record
        self.L=record['L'];self.eta=record['eta']
        if profiles is None:self.theta,self.tensor=profile(self.L)
        else:self.theta,self.tensor=restore(profiles['theta']),restore(profiles['tensor'])
        self.sol=restore(record['solution'])

    def reduced(self,s):
        s=np.asarray(s);L=self.L
        winding=np.floor((s+L)/(2*L)).astype(int)
        v=(s+L)%(2*L)-L
        sign=np.where(winding%2, -1.,1.)
        theta=sign*np.sign(v)*self.theta(abs(v))
        tp=sign*theta_prime(abs(v),L)
        psi=self.sol(abs(v))[0]
        ps=np.sign(v)*self.sol(abs(v),1)[0]
        a=self.eta*self.tensor(abs(v))
        p=momentum(s,L,self.eta)
        return theta,tp,psi,ps,a,p

    def fields(self,x):
        s,t,phi=x
        theta,tp,psi,ps,a,p=self.reduced(s)
        n=np.array([np.sin(t)*np.cos(phi),np.sin(t)*np.sin(phi),np.cos(t)])
        fields=Q*np.r_[np.sin(theta),np.cos(theta)*n]
        tangent=Q*np.r_[np.cos(theta),-np.sin(theta)*n]
        PiE=psi**-6*p*tangent
        return fields,np.sqrt(F)*PiE

    def metric_tensor(self,x,wrong_weight=False,correction=1.):
        s,t,_=x
        theta,tp,psi,ps,a,p=self.reduced(s)
        g=psi**4/F*np.diag([1.,1.,np.sin(t)**2])
        K=psi**-2/np.sqrt(F)*np.diag([2*a,-a,-a*np.sin(t)**2])
        return g,K

    def residuals(self,s):
        s=np.asarray(s);theta=self.theta(s);tp=theta_prime(s,self.L)
        p=momentum(s,self.L,self.eta);a=self.eta*self.tensor(s)
        psi=self.sol(s)[0];dd=self.sol(s,2)[0]
        H=dd-(2-ALPHA*(tp*tp+2*np.cos(theta)**2))*psi/8
        H+=(6*a*a+ALPHA*p*p)*psi**-7/8+POTENTIAL*psi**5/4
        M=2*self.eta*self.tensor(s,1)+ALPHA*p*tp
        return H,M

    def diagnostics(self):
        s=np.linspace(0,1,501);psi=self.sol(s)[0]
        r=psi**2/np.sqrt(F);relative=r/sech(s)-1
        sample=np.array([0.,self.L-.1,self.L-.05,self.L])
        y=self.sol(sample);ps=self.sol(sample,1)[0]
        a=self.eta*self.tensor(sample)
        radii=y[0]**2/np.sqrt(F)
        H=4*np.sqrt(F)*ps/y[0]**3
        trace=-2*np.sqrt(F)*a*y[0]**-6
        return dict(bulk_s=s.tolist(),bulk_relative=relative.tolist(),
                    sample_s=sample.tolist(),radii=radii.tolist(),areas=(4*np.pi*radii*radii).tolist(),
                    mean_curvature=H.tolist(),theta_plus=(H-trace).tolist(),theta_minus=(-H-trace).tolist())


def physical_constraints(data,x,h,order=2):
    # Reuse the independent coordinate curvature engine, not its vacuum verdict.
    class Adapter:
        record={'Lambda':0.}
        metric_tensor=data.metric_tensor
    engine=vacuum.coordinate_constraints if order==2 else _coordinate_order4
    geometry=engine(Adapter(),x,h)
    g=np.asarray(geometry['metric']);inverse=np.linalg.inv(g)
    phi,Pi=data.fields(x);x=np.asarray(x)
    if order==2:
        gradient=np.column_stack([(data.fields(x+np.eye(3)[i]*h)[0]-data.fields(x-np.eye(3)[i]*h)[0])/(2*h) for i in range(3)])
    else:
        offsets,weights=(-2,-1,1,2),(1/12,-2/3,2/3,-1/12)
        gradient=np.column_stack([sum(w*data.fields(x+np.eye(3)[i]*h*o)[0] for o,w in zip(offsets,weights))/h for i in range(3)])
    kinetic=float(Pi@Pi);spatial=float(np.einsum('Ai,ij,Aj',gradient,inverse,gradient))
    scalar_current=Pi@gradient
    geometric=geometry['R']+geometry['trace_K']**2-geometry['K2']
    H=F*geometric-kinetic-spatial-2*LAMBDA
    M=F*np.asarray(geometry['momentum'])+scalar_current
    hscale=max(1.,F*(abs(geometry['R'])+geometry['trace_K']**2+geometry['K2'])+kinetic+spatial+2*LAMBDA)
    mscale=max(1.,F*sum(np.linalg.norm(t) for t in geometry['momentum_terms'])+np.linalg.norm(scalar_current))
    return dict(geometry=geometry,fields=phi.tolist(),Pi_J=Pi.tolist(),gradient=gradient.tolist(),
                kinetic=kinetic,spatial=spatial,current=scalar_current.tolist(),
                H=H,M=M.tolist(),H_normalized=abs(H)/hscale,M_normalized=float(np.linalg.norm(M)/mscale),
                wrong_f_H=abs(geometric-kinetic-spatial-2*LAMBDA)/hscale,
                wrong_f_M=float(np.linalg.norm(np.asarray(geometry['momentum'])+scalar_current)/mscale))


def seam(data,h=1e-5):
    J=np.diag([1.,-1.,1.]);L=data.L;rows=[]
    for t in (.43,.91,1.47,2.13):
        x=np.array([L,t,.37]);z=np.array([-L,np.pi-t,.37+np.pi])
        g,K=data.metric_tensor(x);gz,Kz=data.metric_tensor(z)
        ph,pi=data.fields(x);pz,piz=data.fields(z)
        def derivatives(point):
            plus=point+np.array([h,0,0]);minus=point-np.array([h,0,0])
            gp,kp=data.metric_tensor(plus);gm,km=data.metric_tensor(minus)
            fp,pp=data.fields(plus);fm,pm=data.fields(minus)
            return (gp-gm)/(2*h),(kp-km)/(2*h),(fp-fm)/(2*h),(pp-pm)/(2*h)
        dg,dK,df,dp=derivatives(x);dz,dKz,dfz,dpz=derivatives(z)
        rows.append(dict(t=t,metric=(g-J@gz@J).tolist(),K=(K-J@Kz@J).tolist(),
                         fields=(ph+pz).tolist(),Pi=(pi+piz).tolist(),
                         dmetric=(dg-J@dz@J).tolist(),dK=(dK-J@dKz@J).tolist(),
                         dfields=(df+dfz).tolist(),dPi=(dp+dpz).tolist(),
                         missing_sign=(ph-pz).tolist()))
    return rows


@lru_cache(None)
def symbolic():
    import sympy as s
    z=s.symbols('z',real=True);f=s.Rational(7,8);q2=s.Rational(3,4)
    alpha=q2/f;U=s.Rational(3,2)/f**2
    y=f**s.Rational(1,4)/s.sqrt(s.cosh(z))
    reference=s.diff(y,z,2)-(2-3*alpha/s.cosh(z)**2)*y/8+U*y**5/4
    # General nonminimal projections for constant norm and tangent momenta.
    RE,AE,PE,DE=s.symbols('RE AE PE DE')
    jordan=f*(f*RE-f*AE)-f*PE-f*DE-3
    einstein=f**2*(RE-AE-(PE+DE)/f-3/f**2)
    divE,jE=s.symbols('divE jE')
    mJ=f*s.sqrt(f)*divE+s.sqrt(f)*jE
    mE=f**s.Rational(3,2)*(divE+jE/f)
    radial=1/f+q2/(6*f*f)
    # ADM projections of the summed improved stress. For S=sum(phi^2),
    # g_nn Box S-H_nn=-tr_spatial H. Constant S and nS=0 make this zero;
    # mixed H_ni=D_i(nS)+K_i^j D_j S is zero as well.
    Snn,Hsp,Gnn,Gni,rho,kinni=s.symbols('Snn Hsp Gnn Gni rho kinni')
    box=-Snn+Hsp
    Tnn=rho+(q2*Gnn-box-Snn)/6
    Tni=kinni+q2*Gni/6
    projections=(s.expand((Gnn-Tnn).subs(Hsp,0)-(f*Gnn-rho)),
                 s.expand(Gni-Tni-(f*Gni-kinni)))
    angle=s.symbols('angle',real=True)
    norm=q2*(s.sin(angle)**2+s.cos(angle)**2)
    tangent=q2*(s.sin(angle)*s.cos(angle)-s.cos(angle)*s.sin(angle))
    return tuple(str(s.simplify(v)) for v in (reference,jordan-einstein,mJ-mE,radial-1/f**2,
                                              *projections,norm-q2,tangent))


def reconstruct(record,profiles):
    """Frozen quintic Hermite extension, never a new branch/trajectory solve."""
    d=Data(record,profiles);x=d.sol.x;y=d.sol(x)
    theta=d.theta(x);tp=theta_prime(x,d.L)
    a=d.eta*d.tensor(x);p=momentum(x,d.L,d.eta)
    dd=(2-ALPHA*(tp*tp+2*np.cos(theta)**2))*y[0]/8
    dd-=(6*a*a+ALPHA*p*p)*y[0]**-7/8+POTENTIAL*y[0]**5/4
    # Direct scaled Hermite coefficients avoid subtracting the large constant
    # Bernstein coefficient to recover tiny endpoint slopes. The mathematical
    # interpolant is the same registered unique quintic.
    h=np.diff(x).astype(np.longdouble)
    value=y[0].astype(np.longdouble);velocity=y[1].astype(np.longdouble)
    acceleration=dd.astype(np.longdouble)
    c0=value[:-1];c1=h*velocity[:-1];c2=h*h*acceleration[:-1]/2
    d0=value[1:]-c0-c1-c2
    d1=h*velocity[1:]-c1-2*c2
    d2=h*h*acceleration[1:]-2*c2
    coefficients=np.array([c0,c1,c2,10*d0-4*d1+d2/2,
                           -15*d0+7*d1-d2,6*d0-3*d1+d2/2])
    coefficients/=h[None,:]**np.arange(6)[:,None]
    scalar=PPoly(np.asarray(coefficients[::-1],dtype=float),x)
    derivative=scalar.derivative()
    c=np.zeros(scalar.c.shape+(2,));c[:,:,0]=scalar.c;c[1:,:,1]=derivative.c
    state=PPoly(c,x);state.axis=1
    out=dict(record);out['solution']=serialize(state)
    out['reconstruction']='quintic Hermite from original knots, velocities and ODE accelerations'
    return out


def _coordinate_order4(interpolant, x, h, wrong_weight=False, correction=1.):
    """Direct g/K differences -> Christoffel/Ricci/divergence, no conformal PDE."""
    x = np.asarray(x, dtype=float)
    g, K = interpolant.metric_tensor(x, wrong_weight, correction)
    inverse = np.linalg.inv(g)
    dg, dK = np.zeros((3, 3, 3)), np.zeros((3, 3, 3))
    ddg = np.zeros((3, 3, 3, 3))
    axes = np.eye(3)*h
    offsets=(-2,-1,1,2)
    first=(1/12,-2/3,2/3,-1/12)
    second=(-1/12,4/3,4/3,-1/12)
    for a in range(3):
        gs=[];ks=[]
        for o in offsets:
            gg,kk=interpolant.metric_tensor(x+axes[a]*o,wrong_weight,correction)
            gs.append(gg);ks.append(kk)
        dg[a]=sum(w*gg for w,gg in zip(first,gs))/h
        dK[a]=sum(w*kk for w,kk in zip(first,ks))/h
        ddg[a,a]=(sum(w*gg for w,gg in zip(second,gs))-2.5*g)/h**2
        for b in range(a):
            val=np.zeros((3,3))
            for oa,wa in zip(offsets,first):
                for ob,wb in zip(offsets,first):
                    gg=interpolant.metric_tensor(x+axes[a]*oa+axes[b]*ob,wrong_weight,correction)[0]
                    val+=wa*wb*gg
            ddg[a,b]=ddg[b,a]=val/h**2
    G = np.zeros((3, 3, 3))
    dG = np.zeros((3, 3, 3, 3))
    dinv = np.array([-inverse@dg[a]@inverse for a in range(3)])
    for a in range(3):
        for b in range(3):
            for c in range(3):
                terms = dg[b, :, c]+dg[c, :, b]-dg[:, b, c]
                G[a, b, c] = inverse[a]@terms/2
                for d in range(3):
                    dt = ddg[d, b, :, c]+ddg[d, c, :, b]-ddg[d, :, b, c]
                    dG[d, a, b, c] = (dinv[d, a]@terms+inverse[a]@dt)/2
    ricci = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            for a in range(3):
                ricci[i, j] += dG[a, a, i, j]-dG[j, a, i, a]
                for b in range(3):
                    ricci[i, j] += G[a, a, b]*G[b, i, j]-G[a, j, b]*G[b, i, a]
    R = float(np.sum(inverse*ricci))
    mixed = inverse@K
    K2 = float(np.trace(mixed@mixed))
    trK = float(np.trace(mixed))
    dmixed = np.array([dinv[a]@K+inverse@dK[a] for a in range(3)])
    terms = np.zeros((3, 3))
    for i in range(3):
        terms[0, i] = sum(dmixed[j, j, i] for j in range(3))-np.trace(dmixed[i])
        terms[1, i] = sum(G[j, j, k]*mixed[k, i] for j in range(3) for k in range(3))
        terms[2, i] = -sum(G[k, j, i]*mixed[j, k] for j in range(3) for k in range(3))
    momentum = terms.sum(axis=0)
    Lambda = interpolant.record['Lambda']
    H = R+trK*trK-K2-2*Lambda
    return dict(point=x.tolist(), h=h, R=R, K2=K2, trace_K=trK,
                metric=g.tolist(), extrinsic_curvature=K.tolist(),
                hamiltonian=H, momentum=momentum.tolist(),
                hamiltonian_normalized=abs(H)/max(1., abs(R)+trK*trK+K2+2*abs(Lambda)),
                momentum_normalized=float(np.linalg.norm(momentum)/max(1., sum(np.linalg.norm(t) for t in terms))),
                momentum_terms=terms.tolist())
