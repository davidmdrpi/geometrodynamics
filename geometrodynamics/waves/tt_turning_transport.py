"""Classical turning certificate and two-column TT asymptotic transport.

Freeze 2b928d1. No quantization or absence of all classical invariants follows.
Primes denote conformal eta derivatives; x=eta_star-eta has opposite sign.
"""
from functools import lru_cache
import copy
import math
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PUBLIC_PREREG='2b928d1a550e24019130544e9a4ff326969df965'
BASELINE='6099a0baba1a14536ebbaad50323468c3281836e'
SEED=2026091218
J=np.array([[0.,1.],[-1.,0.]])
DEPENDENCIES=dict(normalization='BAFT',phase_reduction='BAFT',turning_location='B',
    adiabatic_diagnostic='A',instantaneous_actions='A',frobenius_basis='FT',past_basis='T',
    transport_convergence='T',symplectic_completion='T',negative_controls='BAFT',scope='BAFT',failure_paths='BAFT')
TARGETS=dict(B='all_phase_turning_bound',A='named_instantaneous_actions',F='future_basis',T='complete_transport')
SCOPE=dict(all_classical_invariants_excluded=False,quantization='NOT_DERIVED',
           preferred_complex_structure='NOT_SPECIFIED',coupled_rotor='NOT_DERIVED',
           general_support_independence='NOT_ESTABLISHED',Phi_selection='NOT_DERIVED',causality_gate='OPEN')


def radius(x):
    if not np.isfinite(x) or x<=0:raise ValueError('x must be finite and positive')
    return 1/math.tanh(x/math.sqrt(2))


def x_from_radius(R):
    if not np.isfinite(R) or R<=1:raise ValueError('expanding branch requires R>1')
    return math.log((R+1)/(R-1))/math.sqrt(2)


def coefficients(x,alpha,esu=False):
    R=1. if esu else radius(x)
    c=math.cos(alpha-2*x);s=math.sin(alpha-2*x)
    m=R*R-c*c/8
    mp=math.sqrt(2)*R*(R*R-1)+c*s/2
    mpp=3*R**4-4*R*R+1+c*c-s*s
    k=8*R*R+c*c/2
    kp=8*math.sqrt(2)*R*(R*R-1)-2*c*s
    return m,mp,mpp,k,kp


def potential(x,alpha):
    """W and analytic dW/deta, from its exact quadratic-form numerator."""
    R=radius(x);z=R*R;c=math.cos(alpha-2*x);s=math.sin(alpha-2*x)
    m,mp,_,_,_=coefficients(x,alpha)
    N=4*z*z*(9-z)+2*z+(.75*z*z-7*z+.25)*c*c+math.sqrt(2)*R*(z-1)*c*s
    Np=R*(z-1)/math.sqrt(2)*(24*z*(6-z)+3*(z-2)*c*c)+24*z*c*s
    return N/(4*m*m),Np/(4*m*m)-N*mp/(2*m**3)


def turning_point(alpha):
    lo=x_from_radius(3.03);hi=x_from_radius(2.97)
    x=brentq(lambda x:potential(x,alpha)[0],lo,hi,xtol=1e-13,rtol=1e-14)
    W,Wp=potential(x,alpha)
    return dict(alpha=float(alpha),x=x,R=radius(x),W=W,Wprime=Wp)


@lru_cache(None)
def _turning_certificate():
    """Exact rational/polynomial all-phase certificate; no phase scan."""
    import sympy as s
    R,c,t=s.symbols('R c s',real=True);z=s.symbols('z',real=True)
    m=R**2-c*c/8;mp=s.sqrt(2)*R*(R*R-1)+c*t/2
    mpp=3*R**4-4*R*R+1+c*c-t*t;k=8*R*R+c*c/2
    N=4*R**4*(9-R*R)+2*R*R+(s.Rational(3,4)*R**4-7*R*R+s.Rational(1,4))*c*c+s.sqrt(2)*R*(R*R-1)*c*t
    circle=s.Poly(s.expand(4*k*m-2*m*mpp+mp*mp-N),t,c).rem(s.Poly(t*t+c*c-1,t,c)).as_expr()
    Np=s.diff(N,R)*(R*R-1)/s.sqrt(2)+s.diff(N,c)*(-2*t)+s.diff(N,t)*2*c
    expected=R*(R*R-1)/s.sqrt(2)*(24*R*R*(6-R*R)+3*(R*R-2)*c*c)+24*R*R*c*t
    derivative=s.Poly(s.expand(Np-expected),t,c).rem(s.Poly(t*t+c*c-1,t,c)).as_expr()
    N0=4*z*z*(9-z);h22=N0+2*z;h11=N0+s.Rational(3,4)*z*z-5*z+s.Rational(1,4)
    P=16*z*z-147*z+12
    det=s.expand(h11*h22-z*(z-1)**2/2)
    zu=s.Rational(303,100)**2
    identities=dict(numerator=circle,eta_derivative=derivative,
        determinant=s.expand(det-z**3*(z-9)*P),
        at_three=s.Poly(s.expand(N.subs(R,3)-(4*c+3*s.sqrt(2)*t)**2),t,c).rem(s.Poly(t*t+c*c-1,t,c)).as_expr())
    facts=dict(P_at_1=P.subs(z,1),P_at_9=P.subs(z,9),P_second_derivative=s.diff(P,z,2),
        P_at_upper=P.subs(z,zu),P_derivative_at_upper=s.diff(P,z).subs(z,zu),
        negative_h22_factor_at_upper=(2*z*z-18*z-1).subs(z,zu),
        negative_h22_derivative_at_upper=(4*z-18).subs(z,zu),
        bracket_at_9=(24*z*(6-z)+3*(z-2)).subs(z,9),
        bracket_derivative_at_9=(147-48*z).subs(z,9),
        numerator_derivative_upper=-10032+12*zu,
        root_derivative_upper=-9900/(4*zu*zu))
    conditions=[facts['P_at_1']<0,facts['P_at_9']<0,facts['P_second_derivative']>0,
        facts['P_at_upper']>0,facts['P_derivative_at_upper']>0,
        facts['negative_h22_factor_at_upper']>0,facts['negative_h22_derivative_at_upper']>0,
        facts['bracket_at_9']==-627,facts['bracket_derivative_at_9']<0,
        facts['numerator_derivative_upper']<-9900,facts['root_derivative_upper']<-20]
    identities={key:str(s.simplify(v)) for key,v in identities.items()}
    return dict(identities=identities,facts={key:str(v) for key,v in facts.items()},
                inequalities=[bool(v) for v in conditions],verified=set(identities.values())=={'0'} and all(conditions),
                frozen_window=['297/100','303/100'],proof_method='quadratic_form_and_rational_sign_bounds')


def turning_certificate():return copy.deepcopy(_turning_certificate())


def verify_certificate(data):
    """Recompute exact identities/sign facts, rather than trust a verified flag."""
    try:return data==_turning_certificate()
    except (TypeError,ValueError):return False


@lru_cache(None)
def radius_series(order):
    import sympy as s
    x=s.symbols('x')
    expression=s.series(x*x*s.coth(x/s.sqrt(2))**2,x,0,order+1).removeO().expand()
    return tuple(float(expression.coeff(x,i)) for i in range(order+1))


def series_coefficients(alpha,order=12):
    """Columns b_c,b_d, using the equation for h=x^2 m, j=x^2 k."""
    if order<4:raise ValueError('order >=4 required')
    r=np.array(radius_series(order));h=r.copy();j=8*r.copy()
    for n in range(order-1):
        if n==0:value=(1+math.cos(2*alpha))/2
        else:
            tr=math.cos(2*alpha) if n%2==0 else math.sin(2*alpha)
            value=.5*4.**n/math.factorial(n)*(-1)**(n//2)*tr
        h[n+2]-=value/8;j[n+2]+=value/2
    b=np.zeros((order+1,2));b[0,0]=1.;b[3,1]=-1/6
    resonance=np.zeros(2)
    for n in range(1,order+1):
        source=sum((n-l)*(n-3)*h[l]*b[n-l] for l in range(1,n+1))
        if n>=2:source=source+sum(j[l]*b[n-l-2] for l in range(n-1))
        if n==3:resonance=source
        else:b[n]=-source/(2*n*(n-3))
    return b,resonance


def series_basis(x,alpha,order=12):
    b,resonance=series_coefficients(alpha,order)
    n=np.arange(len(b));values=(x**n)@b
    dx=(n[1:]*x**(n[1:]-1))@b[1:]
    dxx=(n[2:]*(n[2:]-1)*x**(n[2:]-2))@b[2:]
    m,mp,_,k,_=coefficients(x,alpha)
    residual=m*dxx-mp*dx+k*values
    scales=np.maximum(1.,abs(m*dxx)+abs(mp*dx)+abs(k*values))
    return dict(matrix=np.array([values,-m*dx]),resonance=resonance,
                residual=residual,normalized=abs(residual)/scales,coefficients=b)


def frobenius_certificate():
    import sympy as s
    c=s.symbols('c',real=True);n=s.symbols('n',integer=True)
    h2=s.Rational(2,3)-c*c/8;j2=s.Rational(16,3)+c*c/2
    b4=s.simplify(-(8*h2+64+j2)/8)
    # h1=j1=b1=0: the resonant n=3 source vanishes, leaving b3 free.
    return dict(indicial='2*n*(n-3)',indices=[0,3],b1='0',b2='4',resonance='0',
                cubic_is_free=True,constant_cubic_choice='0',decaying_cubic='-1/6',
                b4=str(b4),b4_identity=str(s.simplify(b4+s.Rational(28,3)-c*c/16)),
                eta_wronskian_limit='1',x_wronskian_limit='-1',
                action_growing_coefficient='sqrt(8)',action_decaying_coefficient='1/(4*sqrt(8))')


def generator(x,alpha,esu=False):
    m,_,_,k,_=coefficients(x,alpha,esu)
    # State is (beta,p_eta); changing independent variable reverses the flow.
    return np.array([[0.,-1/m],[k,0.]])


def transport(start,end,alpha,initial=None,esu=False,rtol=1e-12,atol=1e-14):
    if not esu and min(start,end)<=0:raise ValueError('cannot cross x=0')
    initial=np.eye(2) if initial is None else np.asarray(initial,dtype=float)
    if initial.shape!=(2,2):raise ValueError('both canonical columns required')
    step=math.pi/100 if esu else min(math.pi/100,min(start,end)/4)
    sol=solve_ivp(lambda x,y:(generator(x,alpha,esu)@y.reshape(2,2)).ravel(),(start,end),initial.ravel(),
                  method='DOP853',rtol=rtol,atol=atol,max_step=step)
    if not sol.success or not np.isfinite(sol.y).all():raise ArithmeticError('canonical propagation failed')
    return sol.y[:,-1].reshape(2,2)


def input_basis(L,alpha,match=1.,**kwargs):
    return transport(L,match,alpha,**kwargs)@transport(0.,L,alpha,esu=True,**kwargs)


def output_basis(alpha,x0=.01,order=12,match=1.,**kwargs):
    data=series_basis(x0,alpha,order)
    return transport(x0,match,alpha,initial=data['matrix'],**kwargs)


def actions(x,alpha,state):
    b,p=np.asarray(state);m,mp,mpp,k,kp=coefficients(x,alpha)
    omega=math.sqrt(k/m);c=m*omega;cp=c*(mp/m+kp/k)/2
    action=(p*p/c+c*b*b)/2
    derivative=cp*(b*b-p*p/(c*c))/2
    y=math.sqrt(m)*b;yp=mp*b/(2*math.sqrt(m))+p/math.sqrt(m)
    W,Wp=potential(x,alpha)
    return dict(J_beta=action,J_beta_prime=derivative,y=y,yprime=yp,
                J_y=(yp*yp+W*y*y)/(2*math.sqrt(W)) if W>0 else None,
                W=W,Wprime=Wp,WKB=abs(Wp)/(2*W**1.5) if W>0 else None)


def valid_transport(evidence):
    try:
        if not isinstance(evidence,list) or len(evidence)!=16:return False
        if not np.allclose([r['alpha'] for r in evidence],np.arange(16)*math.pi/16,rtol=0,atol=1e-14):return False
        json.dumps(evidence,allow_nan=False)
        for row in evidence:
            maps=np.asarray(row['comparison_maps'],dtype=float)
            if maps.shape!=(7,2,2) or not np.isfinite(maps).all():return False
            scale=max(1.,np.linalg.norm(maps[0]))
            if max(np.linalg.norm(M-maps[0])/scale for M in maps)>=1e-8:return False
            if max(abs(np.linalg.det(M)-1) for M in maps)>=1e-8:return False
            if max(np.linalg.norm(M.T@J@M-J) for M in maps)>=1e-8:return False
        return True
    except (TypeError,ValueError,KeyError,OverflowError):return False


def valid_future(evidence):
    try:
        if evidence['frobenius']!=frobenius_certificate():return False
        rows=evidence['series']
        if len(rows)!=144:return False
        json.dumps(rows,allow_nan=False)
        expected={(float(a),n,x) for a in np.arange(16)*math.pi/16 for n in (8,10,12) for x in (.04,.02,.01)}
        if {(r['alpha'],r['order'],r['x0']) for r in rows}!=expected:return False
        for row in rows:
            values=np.asarray(row['normalized']);resonance=np.asarray(row['resonance'])
            if values.shape!=(2,) or resonance.shape!=(2,):return False
            if np.any(values<0) or np.max(values)>=1e-8 or np.max(abs(resonance))!=0:return False
            b=np.asarray(row['coefficients'])
            if b.shape!=(row['order']+1,2) or not np.isfinite(b).all():return False
            if not np.array_equal(b[0],[1,0]) or not np.array_equal(b[3],[0,-1/6]):return False
        return True
    except (TypeError,ValueError,KeyError,OverflowError):return False


def valid_actions(evidence):
    try:
        if evidence['action_identity']!='0':return False
        diag=evidence['diagnostics'];tuned=evidence['tuned'];endpoint=evidence['endpoint'];actions=evidence['actions']
        if (len(diag),len(tuned),len(endpoint),len(actions))!=(80,32,64,112):return False
        json.dumps([diag,tuned,endpoint,actions],allow_nan=False)
        for row in diag:
            samples=row['samples']
            if len(samples)!=4 or abs(samples[-1]['scaled']/samples[-1]['expected']-1)>=1e-3:return False
            if not all(samples[j+1]['WKB']>samples[j]['WKB']>0 for j in range(3)):return False
        if not all(len(r['samples'])==4 and abs(r['samples'][-1]['scaled']-.5)<1e-3 for r in tuned):return False
        last=[r for r in endpoint if r['x']==.00125]
        if len(last)!=16 or any(len(r['ratios'])!=2 or max(abs(v-1) for v in r['ratios'])>=2e-4 for r in last):return False
        if any(len(r['J_beta'])!=16 for r in actions):return False
        return max(max(r['J_beta']) for r in actions)>2
    except (TypeError,ValueError,KeyError,OverflowError,ZeroDivisionError):return False


def verdict(checks,certificate=None,transport_evidence=None,action_evidence=None,future_evidence=None):
    failed={t:[] for t in TARGETS}
    for gate,targets in DEPENDENCIES.items():
        if checks.get(gate) is not True:
            for t in targets:failed[t].append(gate)
    if set(checks)-set(DEPENDENCIES):
        for t in TARGETS:failed[t].append('gate_schema')
    if not verify_certificate(certificate):failed['B'].append('global_certificate')
    if not valid_transport(transport_evidence):failed['T'].append('transport_evidence')
    if not valid_actions(action_evidence):failed['A'].append('action_evidence')
    if not valid_future(future_evidence):
        failed['F'].append('future_evidence');failed['T'].append('future_evidence')
    values=dict(B='CERTIFIED_UNIQUE_SIMPLE_ROOT_IN_FROZEN_WINDOW',
                A='NAMED_ACTION_PREDICTIONS_VERIFIED_NOT_A_GENERAL_NO_GO',
                F='TWO_COLUMN_FROBENIUS_BASIS_VERIFIED',T='REAL_SYMPLECTIC_ASYMPTOTIC_MAP_VERIFIED')
    return {**{TARGETS[t]:values[t] if not failed[t] else 'UNRESOLVED' for t in TARGETS},
            **SCOPE,'failed_checks':failed}
