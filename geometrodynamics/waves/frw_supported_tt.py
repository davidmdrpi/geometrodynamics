"""Fully supported homogeneous TT response on an exact four-scalar FRW family.

Freeze c15acd4. Finite-interval classical transport is not Floquet stability,
an adiabatic invariant, a rotor, or a quantum mode prescription.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp, quad

from . import coupled_multiplet_response as cm
from . import multiplet_scalar_stability as ms
from .scalar_esu_support import curvature_from_jets, improved_stress

PUBLIC_PREREG='c15acd48e61d71378a1af2284767d52bfee94810'
BASELINE='8ac4fd3e100b2ab75f91c08b4fae9ef605885294'
SEED=2026091207
REQUIRED_CHECKS=('exact_background','action_derivation','full_stress_response',
                 'scalar_constraint_closure','independent_geometry','two_clocks',
                 'static_limit','canonical_transport','negative_controls','failure_paths')
IDENTITIES=('spatial_curvature','gradient_expansion','action','stress_action','normal_form',
            'proper_clock','lambda_volume','extrinsic_trace','static_mass','static_stiffness',
            'scalar_trace','scalar_divergence')
SCOPE=dict(general_tensor_harmonics='NOT_DERIVED',nonlinear_stability='NOT_ESTABLISHED',
           tensor_quadratic_scalar_vector_sources='NOT_DERIVED',coupled_rotor='NOT_DERIVED',
           adiabatic_invariant='NOT_TESTED',Bogoliubov_bases='NOT_SPECIFIED',
           preparation_selection='NOT_DERIVED',Phi_selection='NOT_DERIVED',causality_gate='OPEN')
J=np.array([[0.,1.],[-1.,0.]])


@dataclass(frozen=True)
class FRWSupport:
    radius: float=1.
    kappa: float=1.
    departure: float=.1
    sign: int=1
    phase: float=0.

    def __post_init__(self):
        if (not np.isfinite([self.radius,self.kappa,self.departure,self.sign,self.phase]).all()
            or min(self.radius,self.kappa)<=0 or self.departure<=-1 or self.sign not in (-1,1)):
            raise ValueError('positive radius/kappa/initial scale, finite parameters, sign +/-1 required')

    def jets(self,eta):
        if not np.isfinite(eta):raise ValueError('finite conformal time required')
        A,Ap,App=ms.exact_scale(eta,self.departure,self.radius,self.sign)
        q,qp,qpp=ms.field_jets(eta,self.phase)[:3]*self.radius/math.sqrt(self.kappa)
        return A,Ap,App,q,qp,qpp

    def coefficients(self,eta):
        A,Ap,App,q,qp,qpp=self.jets(eta);a=self.radius;kap=self.kappa
        m=(A*A-kap*q*q/6)/(a*a)
        if m<=0:raise ValueError('nonpositive supported tensor kinetic coefficient')
        mp=(2*A*Ap-kap*q*qp/3)/(a*a)
        mpp=(2*(Ap*Ap+A*App)-kap*(qp*qp+q*qpp)/3)/(a*a)
        k=(8*A*A+2*kap*q*q/3)/(a*a)
        return m,mp,mpp,k

    def interval_bound(self,start=0.,end=.8):
        # A is monotone on this exact branch; q^2 <= 3a^2/(4 kappa).
        low=min(self.jets(start)[0],self.jets(end)[0])/self.radius
        # Check the chart throughout via monotonic |w| (including a backwards interval).
        return low*low-1/8

    def generator(self,eta):
        m,_,_,k=self.coefficients(eta)
        return np.array([[0.,1/m],[-k,0.]])

    def normal_potential(self,eta):
        m,mp,mpp,k=self.coefficients(eta)
        return k/m-mpp/(2*m)+mp*mp/(4*m*m)

    def acceleration(self,eta,beta,velocity):
        m,mp,_,k=self.coefficients(eta)
        return -(mp*velocity+k*beta)/m


def exact_certificate():
    import sympy as s
    e,A,Ap,kap,q,qp,V,lam,a=s.symbols('epsilon A Ap kappa q qp V Lambda a',positive=True)
    b1,b2=s.symbols('b1 b2',real=True);eigen=(b1,b2,-b1-b2)
    B=s.diag(*eigen);norm=s.trace(B*B);speed=s.symbols('speed_squared')
    R3=2/A**2*sum(2*s.exp(-2*e*b)-s.exp(4*e*b) for b in eigen)
    r2=s.diff(R3,e,2).subs(e,0)/2
    grad2=s.diff(sum(s.exp(-2*e*b) for b in eigen),e,2).subs(e,0)/2
    M=A*A/kap-q*q/6;Mp=2*A*Ap/kap-q*qp/3;K=8*A*A/kap+s.Rational(2,3)*q*q
    L2=V*((A*A/(2*kap)-q*q/12)*(A*A*r2+speed)-q*q*grad2/2)
    Q=q*q/(A*A);Qp=2*q*qp/(A*A)-2*Ap*q*q/A**3
    b,bp,bpp=s.symbols('b bp bpp')
    G=(bpp+2*Ap*bp/A+8*b)/(A*A)
    T=-2*Q*b/(A*A)+Qp*bp/(6*A*A)+Q*G/6
    m,mp,mpp,k,y,yp,ypp=s.symbols('m mp mpp k y yp ypp',positive=True)
    bj=y/s.sqrt(m);bpj=yp/s.sqrt(m)-mp*y/(2*m**s.Rational(3,2))
    bppj=(ypp/s.sqrt(m)-mp*yp/m**s.Rational(3,2)
          +(3*mp*mp/(4*m**s.Rational(5,2))-mpp/(2*m**s.Rational(3,2)))*y)
    normal=(m*bppj+mp*bpj+k*bj)/s.sqrt(m)-ypp-(k/m-mpp/(2*m)+mp*mp/(4*m*m))*y
    bd,bdd=s.symbols('bd bdd')
    conformal=M*(A*A*bdd+Ap*bd)+Mp*A*bd+K*b
    proper=A*(A*M*bdd+(Ap*M/A+Mp)*bd+K*b/A)
    # Exact exponential volume, including time-dependent F: trace K depends
    # only on A. Jacobi's log-det formula extends this diagonal computation
    # to every time-dependent, noncommuting STF beta.
    volume=s.prod(s.exp(2*e*v) for v in eigen)**s.Rational(1,2)
    lambda2=s.diff(-lam*V*A**4*volume/kap,e,2).subs(e,0)/2
    linear_volume=s.sqrt(s.prod(1+2*e*v for v in eigen))
    bad_lambda2=s.diff(-lam*V*A**4*linear_volume/kap,e,2).subs(e,0)/2
    theta=s.symbols('theta',real=True)
    qstatic=a*s.sqrt(3/(4*kap))*s.cos(theta)
    # Covariant divergence in an orthonormal invariant frame has epsilon_ijk
    # connection coefficients; its contraction with symmetric beta vanishes.
    divergence=[]
    for j in range(3):
        divergence.append(sum(s.LeviCivita(i,j,l)*B[i,l] for i in range(3) for l in range(3)))
    values=dict(spatial_curvature=r2+8*norm/A**2,gradient_expansion=grad2-2*norm,
        action=L2-V*(M*speed-K*norm)/2,
        stress_action=G-kap*T-kap*(M*bpp+Mp*bp+K*b)/A**4,normal_form=normal,
        proper_clock=conformal-proper,lambda_volume=lambda2,
        extrinsic_trace=s.trace(B),static_mass=M.subs({A:a,q:qstatic})*kap/a**2-(1-s.cos(theta)**2/8),
        static_stiffness=K.subs({A:a,q:qstatic})*kap/a**2-(8+s.cos(theta)**2/2),
        scalar_trace=s.trace(B),scalar_divergence=sum(abs(v) for v in divergence))
    # General symmetric divergence, independent of diagonal action reduction.
    entries=s.symbols('b00 b01 b02 b11 b12 b22')
    general=s.Matrix([[entries[0],entries[1],entries[2]],[entries[1],entries[3],entries[4]],
                      [entries[2],entries[4],entries[5]]])
    values['scalar_divergence']=sum(abs(s.simplify(sum(s.LeviCivita(i,j,l)*general[i,l]
                                        for i in range(3) for l in range(3)))) for j in range(3))
    values={key:str(s.simplify(value)) for key,value in values.items()}
    return dict(residuals=values,all_zero=set(values.values())=={'0'},
                bad_linear_lambda_coefficient=str(s.factor(bad_lambda2)),
                bad_linear_lambda_identity=str(s.simplify(bad_lambda2-lam*V*A**4*norm/kap)),
                bad_linear_lambda_value=float(bad_lambda2.subs({A:1,lam:1,V:1,kap:1,b1:1,b2:0})))


def expected_response(model,eta,beta,velocity,acceleration):
    A,Ap,App,q,qp,qpp=model.jets(eta);kap=model.kappa
    Q=q*q/(A*A);Qp=2*q*qp/(A*A)-2*Ap*Q/A
    G=(acceleration+2*Ap*velocity/A+8*beta)/(A*A)
    T=-2*Q*beta/(A*A)+Qp*velocity/(6*A*A)+Q*G/6
    fullG=np.zeros((4,4));fullG[1:,1:]=G
    fullT=np.zeros((4,4));fullT[1:,1:]=T
    return dict(einstein=fullG,stress=fullT,residual=fullG-kap*fullT,KG=np.zeros(4),R=0.)


def full_geometry(model,eta,beta,velocity,acceleration,epsilon=0.,clock='conformal',angles=(.83,1.07,.61)):
    """Full metric jets and four off-shell improved stresses; no reduced ODE.

    Input tensor jets are conformal. The proper-clock route transforms every
    jet before forming curvature, including the lapse and scalar second jet.
    I+2 epsilon beta is used only as a first-variation continuation.
    """
    A,Ap,App,q,qp,qpp=model.jets(eta)
    P=q/A;Pp=qp/A-q*Ap/(A*A)
    Ppp=qpp/A-2*qp*Ap/(A*A)-q*App/(A*A)+2*q*Ap*Ap/A**3
    if clock=='conformal':N,Nd,Ndd=A,Ap,App;At,Att=Ap,App;Pt,Ptt=Pp,Ppp;bt,btt=velocity,acceleration
    elif clock=='proper':
        N,Nd,Ndd=1.,0.,0.;At=Ap/A;Att=App/(A*A)-Ap*Ap/A**3
        Pt=Pp/A;Ptt=Ppp/(A*A)-Pp*Ap/A**3
        bt=velocity/A;btt=acceleration/(A*A)-velocity*Ap/A**3
    else:raise ValueError('clock must be conformal or proper')
    B=np.eye(3)+2*epsilon*beta
    if np.linalg.eigvalsh(B).min()<=0:raise ValueError('nonpositive spatial continuation')
    Bt=2*epsilon*bt;Btt=2*epsilon*btt
    H=A*A*B;Ht=2*A*At*B+A*A*Bt
    Htt=2*(At*At+A*Att)*B+4*A*At*Bt+A*A*Btt
    x,dx,ddx,E,dE,ddE=cm.coframe_jets(tuple(angles))
    metric=np.zeros((4,4));metric[0,0]=-N*N;metric[1:,1:]=E.T@H@E
    dg=np.zeros((4,4,4));ddg=np.zeros((4,4,4,4))
    dg[0,0,0]=-2*N*Nd;ddg[0,0,0,0]=-2*(Nd*Nd+N*Ndd)
    dg[0,1:,1:]=E.T@Ht@E;ddg[0,0,1:,1:]=E.T@Htt@E
    for c in range(3):
        dg[c+1,1:,1:]=dE[c].T@H@E+E.T@H@dE[c]
        ddg[0,c+1,1:,1:]=ddg[c+1,0,1:,1:]=dE[c].T@Ht@E+E.T@Ht@dE[c]
        for d in range(3):
            ddg[c+1,d+1,1:,1:]=(ddE[c,d].T@H@E+dE[c].T@H@dE[d]
                                       +dE[d].T@H@dE[c]+E.T@H@ddE[c,d])
    inverse,connection,ricci,G,R=curvature_from_jets(metric,dg,ddg)
    T=np.zeros((4,4));kg=[]
    for i in range(4):
        derivative=np.r_[Pt*x[i],P*dx[i]]
        second=np.empty((4,4));second[0,0]=Ptt*x[i]
        second[0,1:]=second[1:,0]=Pt*dx[i];second[1:,1:]=P*ddx[i]
        Hess=second-np.einsum('kij,k->ij',connection,derivative)
        T+=improved_stress(metric,inverse,G,P*x[i],derivative,Hess)
        kg.append(np.sum(inverse*Hess)-R*P*x[i]/6)
    frame=np.zeros((4,4));frame[0,0]=N;frame[1:,1:]=A*E;fi=np.linalg.inv(frame)
    Tm=frame@inverse@T@fi;Gm=frame@inverse@G@fi
    return dict(stress=Tm,einstein=Gm,residual=Gm+3/(2*model.radius**2)*np.eye(4)-model.kappa*Tm,
                KG=np.array(kg),R=R,phi=P*x)


def evolve(model,start=0.,end=.8,rtol=1e-12,atol=1e-14,normal=False):
    if model.interval_bound(start,end)<=0:raise ValueError('interval has no certified positive kinetic bound')
    def transform(t):
        m,mp,_,_=model.coefficients(t)
        return np.array([[math.sqrt(m),0.],[mp/(2*math.sqrt(m)),1/math.sqrt(m)]])
    initial=transform(start) if normal else np.eye(2)
    def rhs(t,state):
        B=np.array([[0.,1.],[-model.normal_potential(t),0.]]) if normal else model.generator(t)
        return (B@state.reshape(2,2)).ravel()
    sol=solve_ivp(rhs,(start,end),initial.ravel(),method='DOP853',rtol=rtol,atol=atol,max_step=math.pi/100)
    if not sol.success or not np.isfinite(sol.y).all():raise ArithmeticError('tensor transport failed')
    result=sol.y[:,-1].reshape(2,2)
    return np.linalg.solve(transform(end),result) if normal else result


def proper_evolve(model,start=0.,end=.8,rtol=1e-12,atol=1e-14):
    """Independent proper-clock background IVP, normalized canonical pair."""
    bound=model.interval_bound(start,end)
    if bound<=0:raise ValueError('interval has no certified positive kinetic bound')
    elapsed=quad(lambda eta:model.jets(eta)[0],start,end,epsabs=1e-12,epsrel=1e-12)[0]
    a=model.radius;kap=model.kappa
    def rhs(t,y):
        A,eta=y[:2]
        q=a/math.sqrt(kap)*ms.field_jets(eta,model.phase)[0]
        m=(A*A-kap*q*q/6)/(a*a);k=(8*A*A+2*kap*q*q/3)/(a*a)
        if A<=0 or m<=0:raise ValueError('proper-time evolution left regular chart')
        generator=np.array([[0.,1/m],[-k,0.]])/A
        Ad=model.sign*(A*A-a*a)/(math.sqrt(2)*a*A)
        return np.r_[Ad,1/A,(generator@y[2:].reshape(2,2)).ravel()]
    sol=solve_ivp(rhs,(0,elapsed),np.r_[model.jets(start)[0],start,np.eye(2).ravel()],
                  method='DOP853',rtol=rtol,atol=atol,max_step=min(model.jets(start)[0],model.jets(end)[0])*math.pi/100)
    if not sol.success or not np.isfinite(sol.y).all():raise ArithmeticError('proper tensor/background evolution failed')
    last=sol.y[:,-1]
    return dict(map=last[2:].reshape(2,2),eta_error=abs(float(last[1]-end)),
                scale_error=abs(float(last[0]-model.jets(end)[0]))/a,elapsed=elapsed)


def valid_maps(evidence):
    try:
        if not isinstance(evidence,list) or not evidence:return False
        for row in evidence:
            maps=np.asarray(row['maps'],dtype=float)
            if maps.shape!=(4,2,2) or not np.isfinite(maps).all():return False
            scale=max(1.,np.linalg.norm(maps[0]))
            if max(np.linalg.norm(x-maps[0])/scale for x in maps[1:])>=1e-8:return False
            if max(np.linalg.norm(x.T@J@x-J) for x in maps)>=1e-8:return False
            if max(abs(np.linalg.det(x)-1) for x in maps)>=1e-8:return False
            endpoint=float(row['eta_error']);bound=float(row['kinetic_bound'])
            if not np.isfinite([endpoint,bound]).all() or not 0<=endpoint<1e-9 or bound<=.2:return False
        return True
    except (KeyError,ValueError,TypeError,OverflowError):return False


def verdict(checks,evidence=None):
    failed=[k for k in REQUIRED_CHECKS if checks.get(k) is not True]
    if set(checks)!=set(REQUIRED_CHECKS):failed.append('gate_schema')
    if not valid_maps(evidence):failed.append('map_evidence')
    names=('linear_tensor_sector','supported_operator','clock_agreement','finite_interval_transport')
    values=('CLOSED_FOR_ZERO_SCALAR_PERTURBATION_DATA','FULLY_SUPPORTED_FRW_EQUATION_VERIFIED',
            'CONFORMAL_AND_PROPER_AGREE','CANONICAL_TRANSPORT_VERIFIED')
    return {**dict(zip(names,values if not failed else ['UNRESOLVED']*4)),**SCOPE,'failed_checks':sorted(set(failed))}
