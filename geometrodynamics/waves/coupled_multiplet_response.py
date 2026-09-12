"""Linear homogeneous TT response of the four-component degree-1 ESU.

Freeze 64d7f9c. No general scalar/vector response or nonlinear stability is
inferred. The support's metric stress is varied, not assumed rigid.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import solve_ivp

from . import reciprocal_scalar_tt as rt
from .scalar_esu_support import curvature_from_jets, improved_stress

PUBLIC_PREREG = '64d7f9c8e41a9b8ef591eae7b04685f038ed49b4'
BASELINE = 'fe3d423d322fd1ea8760dc4aee09dff7fdf2acb9'
SEED = 2026091014
REQUIRED_CHECKS = ('background_and_parity','scalar_sector_closure','quadratic_action',
                  'independent_curvature_stress','constraint_completion','bare_frequency_control',
                  'canonical_propagation','period_map_convergence','scope_and_order','failure_paths')
PHYSICAL = ('linear_tensor_sector','coupled_tensor_equation','bare_frequency_transfer','period_map_type')
SCOPE = dict(general_scalar_vector_response='NOT_DERIVED',full_dynamical_stability='NOT_ESTABLISHED',
             nonlinear_persistence='NOT_ESTABLISHED',preparation_selection='NOT_DERIVED',
             Phi_selection='NOT_DERIVED',causality_gate='OPEN')


def classify_period(traces):
    traces=np.asarray(traces,dtype=float)
    if traces.shape!=(2,) or not np.isfinite(traces).all() or abs(traces[0]-traces[1])>=1e-8:
        return 'UNRESOLVED'
    if np.all(abs(traces)<2-1e-7): return 'NUMERICALLY_ELLIPTIC'
    if np.all(abs(traces)>2+1e-7): return 'NUMERICALLY_HYPERBOLIC'
    return 'UNRESOLVED'


def verdict(checks,period_evidence=None):
    failed=[k for k in REQUIRED_CHECKS if checks.get(k) is not True]
    if set(checks)!=set(REQUIRED_CHECKS): failed.append('gate_schema')
    kind='UNRESOLVED'
    if period_evidence is not None:
        try:
            maps=np.asarray(period_evidence['maps'],dtype=float)
            traces=np.asarray(period_evidence['traces'],dtype=float)
            if (maps.shape==(2,2,2) and traces.shape==(2,) and np.isfinite(maps).all()
                and np.max(abs(np.trace(maps,axis1=1,axis2=2)-traces))<1e-12
                and np.max(abs(np.linalg.det(maps)-1))<1e-8
                and np.linalg.norm(maps[0]-maps[1])<1e-8):
                kind=classify_period(traces)
        except (KeyError,TypeError,ValueError): pass
    if kind=='UNRESOLVED': failed.append('period_evidence')
    if failed: values=dict.fromkeys(PHYSICAL,'UNRESOLVED')
    else:
        values=dict(linear_tensor_sector='CLOSED_AT_LINEAR_ORDER_FOR_TENSOR_INITIAL_DATA',
                    coupled_tensor_equation='PERIODIC_HAMILTONIAN_RESPONSE_VERIFIED',
                    bare_frequency_transfer='BARE_EQUATION_FAILS_FOR_THIS_SUPPORT',period_map_type=kind)
    return {**values,**SCOPE,'failed_checks':sorted(set(failed))}


@dataclass(frozen=True)
class CoupledSupport:
    radius: float=1.
    kappa: float=1.
    phase: float=0.

    def __post_init__(self):
        if not np.isfinite([self.radius,self.kappa,self.phase]).all() or min(self.radius,self.kappa)<=0:
            raise ValueError('finite positive radius and kappa; finite phase required')

    @property
    def amplitude(self):return math.sqrt(3/(4*self.kappa))
    @property
    def volume(self):return 2*math.pi**2*self.radius**3
    @property
    def C(self):return self.volume/self.kappa
    @property
    def Lambda(self):return 3/(2*self.radius**2)

    def field_jets(self,time):
        theta=2*time/self.radius+self.phase
        P=self.amplitude*math.cos(theta)
        Pd=-2*self.amplitude*math.sin(theta)/self.radius
        return P,Pd,-4*P/self.radius**2

    def coefficients(self,tau,eta=1.):
        """eta!=1 is a mathematical operator control, not an ESU solution."""
        theta=2*np.asarray(tau)+self.phase
        c=np.cos(theta)
        return (1-eta*c*c/8,eta*np.sin(2*theta)/4,
                eta*np.cos(2*theta),8+eta*c*c/2)

    def generator(self,tau,eta=1.):
        f,fp,fpp,g=self.coefficients(tau,eta)
        return np.array([[0.,1/f],[-g,0.]])

    def residual(self,time,beta,velocity,acceleration):
        f,fp,_,g=self.coefficients(time/self.radius)
        return f*acceleration+fp*velocity/self.radius+g*beta/self.radius**2

    def predicted_stress(self,time,beta,velocity,acceleration):
        P,Pd,_=self.field_jets(time); Q=P*P; Qd=2*P*Pd
        G=acceleration+8*beta/self.radius**2
        return -2*Q*beta/self.radius**2+Qd*velocity/6+Q*G/6

    def quadratic_action(self,time,beta,velocity):
        f,_,_,g=self.coefficients(time/self.radius)
        return self.C*(f*np.sum(velocity*velocity)-g*np.sum(beta*beta)/self.radius**2)/2


def stf(M):
    return M-np.trace(M,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)/3


def exact_certificate():
    """Exact action, full-stress and scalar-closure identities."""
    import sympy as sp
    e,a,kap,Q,Qd,V=sp.symbols('e a kappa Q Qdot V',positive=True)
    b1,b2=sp.symbols('b1 b2',real=True)
    B=sp.diag(b1,b2,-b1-b2)
    R3=2/a**2*sum(2*sp.exp(-2*e*b)-sp.exp(4*e*b) for b in (b1,b2,-b1-b2))
    r2=sp.diff(R3,e,2).subs(e,0)/2
    grad2=sp.diff(sum(sp.exp(-2*e*b) for b in (b1,b2,-b1-b2)),e,2).subs(e,0)/2
    # Use exact symbolic norm rather than rationalizing numerical STF radicals.
    speed2=sp.symbols('speed2')
    f=1-kap*Q/6; g=8+sp.Rational(2,3)*kap*Q
    L2=V*((r2+speed2)/(2*kap)-Q*grad2/(2*a*a)-Q*(r2+speed2)/12)
    expected=V/(2*kap)*(f*speed2-g*sp.trace(B*B)/a**2)
    beta,bd,bdd=sp.symbols('beta bd bdd')
    Ein=bdd+8*beta/a**2
    stress=-2*Q*beta/a**2+Qd*bd/6+Q*Ein/6
    equation=f*bdd-kap*Qd*bd/6+g*beta/a**2
    # Degree-1 identities directly in the four-dimensional embedding.
    S=[sp.Matrix(M.astype(int)) for M in rt.QUATERNION_DERIVATIVES]
    hess=[S[i]*S[j]+S[j]*S[i]+2*(i==j)*sp.eye(4) for i in range(3) for j in range(3)]
    x=sp.Matrix(sp.symbols('x0:4'))
    gradients=sp.Matrix.hstack(*(M*x for M in S))
    norm=(x.T*x)[0]
    tensor=gradients.T*gradients-norm*sp.eye(3)
    theta=sp.symbols('theta',real=True)
    c=sp.cos(theta); physicalQ=3*c*c/(4*kap)
    kinetic=sp.simplify(f.subs(Q,physicalQ)-(1-c*c/8))
    stiffness=sp.simplify(g.subs(Q,physicalQ)-(8+c*c/2))
    traceB=sp.trace(B)
    residuals=dict(spatial_curvature=sp.expand(r2+8*sp.trace(B*B)/a**2),
                   gradient_expansion=sp.expand(grad2-2*sp.trace(B*B)),
                   quadratic_action=sp.expand(L2-expected),
                   stress_action=sp.expand(Ein-kap*stress-equation),
                   determinant_tangent=traceB,kinetic=kinetic,stiffness=stiffness,
                   bare_turning_point=sp.Rational(7,8)*(-8)+sp.Rational(17,2)-sp.Rational(3,2))
    for i,M in enumerate(hess):residuals['hessian_'+str(i)]=sum(abs(z) for z in M)
    for i,z in enumerate(tensor):residuals['gradient_'+str(i)]=sp.expand(z)
    # Contracting Hess x=-g*x with any STF beta and using delta R=0
    # yields zero KG source for every component; parity is polynomial degree 1.
    for i,xi in enumerate(x):residuals['KG_'+str(i)]=traceB*xi
    residuals={k:str(sp.simplify(value)) for k,value in residuals.items()}
    return dict(residuals=residuals,all_zero=set(residuals.values())=={'0'},
                parity_degree=1,components=4,minimum_f='7/8',
                curvature_kinetic_proof='tr(K)=0 for det(exp(2 epsilon beta))=1; tr(K^2)=epsilon^2 tr(beta_dot^2)+O(epsilon^3)')


@lru_cache(None)
def _embedding_functions():
    """Independent coordinate jets through order three; no tensor equation."""
    import sympy as sp
    u=sp.symbols('chi theta phi',real=True)
    ch,th,ph=u
    x=sp.Matrix([sp.cos(ch),sp.sin(ch)*sp.cos(th),sp.sin(ch)*sp.sin(th)*sp.cos(ph),
                 sp.sin(ch)*sp.sin(th)*sp.sin(ph)])
    dx=x.jacobian(u)
    dd=[[[sp.diff(x[i],u[j],u[k]) for k in range(3)] for j in range(3)] for i in range(4)]
    ddd=[[[[sp.diff(x[i],u[j],u[k],u[l]) for l in range(3)] for k in range(3)]
           for j in range(3)] for i in range(4)]
    return sp.lambdify(u,(x,dx,dd,ddd),'numpy',cse=True)


@lru_cache(None)
def coframe_jets(angles):
    x,dx,dd,ddd=_embedding_functions()(*angles)
    x=np.asarray(x).ravel(); dx=np.asarray(dx); dd=np.asarray(dd); ddd=np.asarray(ddd)
    S=rt.QUATERNION_DERIVATIVES
    E=np.einsum('iab,b,ac->ic',S,x,dx)
    dE=np.empty((3,3,3)); ddE=np.empty((3,3,3,3))
    for c in range(3):
        dE[c]=np.einsum('iab,b,ad->id',S,dx[:,c],dx)+np.einsum('iab,b,ad->id',S,x,dd[:,:,c])
        for d in range(3):
            ddE[c,d]=(np.einsum('iab,b,ae->ie',S,dd[:,c,d],dx)
                +np.einsum('iab,b,ae->ie',S,dx[:,c],dd[:,:,d])
                +np.einsum('iab,b,ae->ie',S,dx[:,d],dd[:,:,c])
                +np.einsum('iab,b,ae->ie',S,x,ddd[:,:,c,d]))
    return x,dx,dd,E,dE,ddE


def full_geometry(model,time,beta,velocity,acceleration,epsilon,angles=(.83,1.07,.61)):
    """Full coordinate curvature and all four component stresses, off shell.

    The linear metric continuation I+2 epsilon beta has the same tangent as
    exp(2 epsilon beta). Finite epsilon is not asserted to solve Einstein.
    """
    x,dx,ddx,E,dE,ddE=coframe_jets(tuple(angles))
    a=model.radius; M=np.eye(3)+2*epsilon*beta
    if np.linalg.eigvalsh(M).min()<=0:raise ValueError('spatial metric must be positive')
    Mt=2*epsilon*velocity; Mtt=2*epsilon*acceleration
    metric=np.zeros((4,4)); metric[0,0]=-1.; metric[1:,1:]=a*a*E.T@M@E
    dg=np.zeros((4,4,4)); ddg=np.zeros((4,4,4,4))
    dg[0,1:,1:]=a*a*E.T@Mt@E
    ddg[0,0,1:,1:]=a*a*E.T@Mtt@E
    for c in range(3):
        dg[c+1,1:,1:]=a*a*(dE[c].T@M@E+E.T@M@dE[c])
        ddg[0,c+1,1:,1:]=ddg[c+1,0,1:,1:]=a*a*(dE[c].T@Mt@E+E.T@Mt@dE[c])
        for d in range(3):
            ddg[c+1,d+1,1:,1:]=a*a*(ddE[c,d].T@M@E+dE[c].T@M@dE[d]
                +dE[d].T@M@dE[c]+E.T@M@ddE[c,d])
    inverse,gamma,ricci,einstein,R=curvature_from_jets(metric,dg,ddg)
    P,Pd,Pdd=model.field_jets(time)
    stress=np.zeros((4,4)); waves=[]
    for I in range(4):
        derivative=np.r_[Pd*x[I],P*dx[I]]
        second=np.empty((4,4));second[0,0]=Pdd*x[I]
        second[0,1:]=second[1:,0]=Pd*dx[I];second[1:,1:]=P*ddx[I]
        hessian=second-np.einsum('kij,k->ij',gamma,derivative)
        stress+=improved_stress(metric,inverse,einstein,P*x[I],derivative,hessian)
        waves.append(np.sum(inverse*hessian)-R*P*x[I]/6)
    # Mixed invariant components remove the background pressure/metric term.
    frame=np.zeros((4,4));frame[0,0]=1.;frame[1:,1:]=a*E
    frame_inv=np.linalg.inv(frame)
    T=frame@inverse@stress@frame_inv
    G=frame@inverse@einstein@frame_inv
    residual=G+model.Lambda*np.eye(4)-model.kappa*T
    return dict(stress=T,einstein=G,residual=residual,KG=np.array(waves),R=R,
                metric=metric,phi=P*x,coframe=E)


def finite_response(model,time,beta,velocity,acceleration,epsilon,angles):
    plus=full_geometry(model,time,beta,velocity,acceleration,epsilon,angles)
    minus=full_geometry(model,time,beta,velocity,acceleration,-epsilon,angles)
    return {key:(plus[key]-minus[key])/(2*epsilon) for key in ('stress','einstein','residual','KG','R')}


def evolve(model=CoupledSupport(),periods=1,rtol=1e-12,atol=1e-14,eta=1.,normal_form=False,max_step=np.pi/100):
    T=np.pi/2
    if normal_form:
        f,fp,_,_=model.coefficients(0.,eta)
        conversion=np.array([[math.sqrt(f),0.],[fp/(2*math.sqrt(f)),1/math.sqrt(f)]])
        initial=conversion
        def rhs(tau,state):
            f,fp,fpp,g=model.coefficients(tau,eta)
            q=g/f-fpp/(2*f)+fp*fp/(4*f*f)
            return (np.array([[0.,1.],[-q,0.]])@state.reshape(2,2)).ravel()
    else:
        initial=np.eye(2)
        def rhs(tau,state):return (model.generator(tau,eta)@state.reshape(2,2)).ravel()
    times=np.linspace(0,periods*T,100*periods+1)
    sol=solve_ivp(rhs,(0,periods*T),initial.ravel(),t_eval=times,method='DOP853',
                  rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success or not np.isfinite(sol.y).all():raise ArithmeticError('linear response integration failed')
    matrices=sol.y.T.reshape(-1,2,2)
    if normal_form:
        for i,tau in enumerate(times):
            f,fp,_,_=model.coefficients(tau,eta)
            back=np.array([[1/math.sqrt(f),0.],[-fp/(2*math.sqrt(f)),math.sqrt(f)]])
            matrices[i]=back@matrices[i]
    return times,matrices


def physical_evolve(model,time_end):
    """Independent physical-time (b,b_dot) evolution; proper time is t."""
    f0=model.coefficients(0.)[0]
    initial=np.diag([1.,1/(model.radius*f0)])
    def rhs(t,state):
        f,fp,_,g=model.coefficients(t/model.radius)
        A=np.array([[0.,1.],[-g/(model.radius**2*f),-fp/(model.radius*f)]])
        return (A@state.reshape(2,2)).ravel()
    sol=solve_ivp(rhs,(0,time_end),initial.ravel(),method='DOP853',rtol=1e-12,atol=1e-14,
                  max_step=model.radius*np.pi/100)
    if not sol.success:raise ArithmeticError('physical-time solve failed')
    f=model.coefficients(time_end/model.radius)[0]
    return np.diag([1.,model.radius*f])@sol.y[:,-1].reshape(2,2)
