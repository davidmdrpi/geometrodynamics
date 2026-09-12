"""Constrained low scalar modes of the four-field ESU; freeze 8eb33d5.

The cover dipole is a linear scalar-type block, not a nonlinear completion.
Independent coordinate curvature below does not use the reduced equations.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space, expm

from . import coupled_multiplet_response as cm
from .scalar_esu_support import curvature_from_jets, improved_stress

PUBLIC_PREREG = '8eb33d5348138cafd73c77f0a4ba66a285e99a73'
BASELINE = '9a0bfbdbc553e85a4bef5632ee5184d760db95dd'
SEED = 2026091106
REQUIRED_CHECKS = ('field_reduction', 'homogeneous_constraints', 'exact_continuation',
                  'clock_and_gauge', 'dipole_reduction', 'dipole_parity',
                  'independent_geometry', 'period_maps', 'negative_controls', 'failure_paths')
EXACT_IDENTITIES = ('density','pressure','KG','constraint_propagation',
                    'homogeneous_linear_constraint','homogeneous_linear_acceleration','dipole_H_propagation',
                    'dipole_J_propagation','dipole_trace_reduction','exact_branch')
HOMOGENEOUS = np.array([[0., 1.], [2., 0.]])
# State (u,u',v,v'), for dimensionless sqrt(kappa) delta phi.
DIPOLE = np.array([[0., 1., 0., 0.], [-7., 0., 6., 0.],
                   [0., 0., 0., 1.], [2., 0., -3., 0.]])
SCOPE = dict(vector_response='NOT_DERIVED', tensor_quadratic_sources='NOT_DERIVED',
             dipole_nonlinear_completion='NOT_ESTABLISHED', coupled_rotor='NOT_DERIVED',
             preparation_measure='NOT_SPECIFIED', Phi_selection='NOT_DERIVED', causality_gate='OPEN')


def field_jets(tau, phase=0.):
    p=math.sqrt(3/4)*math.cos(2*tau+phase)
    pd=-math.sqrt(3)*math.sin(2*tau+phase)
    return np.array([p,pd,-4*p,-4*pd])


def constraints(tau, phase=0.):
    p,pd,*_=field_jets(tau,phase)
    return np.array([[5*p,pd,-3*p,0.], [2*pd/3,-p/3,-pd,p]])


def dipole_coefficients(tau, u, v, lapse=(0.,0.,0.), phase=0.):
    """Off-shell Hamiltonian, covariant flux, pressure and KG coefficients.

    Metric is a^2[-(1+2 epsilon alpha Y)d tau^2+gamma]. Einstein
    residuals are (-rho, -J grad Y, ..., (-2 alpha-pressure)Y delta).
    u,v contain value/first/second derivative. Stress units 1/(kappa a^2).
    """
    p,pd,*_=field_jets(tau,phase); u,ud,udd=u;v,vd,vdd=v; al,ald,_=lapse
    rho=pd*ud+5*p*u-3*p*v-al*pd*pd
    flux=p*vd-p*ud/3+2*pd*u/3-pd*v+al*p*pd/3
    pressure=(pd*ud-p*udd-2*p*u-al*pd*pd+ald*p*pd)/3+p*v-3*al*p*p
    kg_xy=-(udd-vdd)-9*(u-v)-10*al*p+ald*pd
    kg_d=-vdd+2*u-3*v+al*p
    return dict(rho=rho,flux=flux,pressure=pressure,kg_xy=kg_xy,kg_d=kg_d,curvature=6*al)


def exact_certificate():
    """Symbolic constraints, propagation and off-shell FRW reduction."""
    import sympy as s
    A,Ad,Add,q,qd,qdd,kap,lam=s.symbols('A Ad Add q qd qdd kappa Lambda', nonzero=True)
    p=q/A; pt=qd/A-q*Ad/A**2
    ptt=qdd/A-2*qd*Ad/A**2-q*Add/A**2+2*q*Ad**2/A**3
    # Full improved stress before KG: canonical gradient has 3 p^2/A^2.
    Q=p*p; Qt=2*p*pt; Qtt=2*(pt*pt+p*ptt)
    h=Ad/A
    G00=3*(h*h+1)
    Gsp=-2*Add/A+h*h-1
    rho=(pt*pt+3*p*p)/(2*A*A)+(3*h*Qt+G00*Q)/(6*A*A)
    pressure=(pt*pt-p*p)/(2*A*A)+(-Qtt-h*Qt+Gsp*Q)/(6*A*A)
    R=6*(Add/A+1)/A**2
    kg=(-ptt-2*h*pt-3*p)/A**2-R*p/6
    E=(qd*qd+4*q*q)/2
    C=Ad*Ad+A*A-lam*A**4/3-kap*E/3
    Cprime=s.diff(C,A)*Ad+s.diff(C,Ad)*Add+s.diff(C,q)*qd+s.diff(C,qd)*qdd
    P,Pd,U,Ud,V,Vd,al,ald=s.symbols('P Pd U Ud V Vd alpha alphad')
    H=Pd*Ud+5*P*U-3*P*V
    J=P*Vd-P*Ud/3+2*Pd*U/3-Pd*V
    def evolve(expr):
        return sum(s.diff(expr,z)*dz for z,dz in zip((P,Pd,U,Ud,V,Vd),
                    (Pd,-4*P,Ud,-7*U+6*V,Vd,2*U-3*V)))
    # Lapse retained here; Einstein trace + KG + H force alpha=0.
    udd=-7*U+6*V-9*al*P+ald*Pd
    pr=(Pd*Ud-P*udd-2*P*U-al*Pd**2+ald*P*Pd)/3+P*V-3*al*P**2
    Hr=H-al*Pd**2
    a,epsilon,r,rd,rdd,z,zd=s.symbols('a epsilon r rd rdd z zd', nonzero=True)
    bgq,bgqd=s.symbols('bgq bgqd')
    clin=s.diff(C.subs({A:a*(1+epsilon*r),Ad:a*epsilon*rd,q:bgq+epsilon*z,qd:bgqd+epsilon*zd}),epsilon).subs(epsilon,0)
    alin=s.diff((Add+A-2*lam*A**3/3).subs({A:a*(1+epsilon*r),Add:a*epsilon*rdd}),epsilon).subs(epsilon,0)
    residuals=dict(
        density=s.factor(rho-E/A**4),
        pressure=s.factor(pressure.subs(qdd,-4*q)-E/(3*A**4)),
        KG=s.factor(kg+(qdd+4*q)/A**3),
        constraint_propagation=s.factor(Cprime.subs({qdd:-4*q,Add:-A+2*lam*A**3/3})),
        homogeneous_linear_constraint=s.factor(clin.subs(lam,3/(2*a*a))+kap*(bgqd*zd+4*bgq*z)/3),
        homogeneous_linear_acceleration=s.factor(alin.subs(lam,3/(2*a*a))-a*(rdd-2*r)),
        dipole_H_propagation=s.expand(evolve(H)+3*J),
        dipole_J_propagation=s.expand(evolve(J)-H/3),
        dipole_trace_reduction=s.expand(pr-Hr/3),
        exact_branch=s.factor(((A*A-a*a)**2/(2*a*a)+A*A-A**4/(2*a*a))-a*a/2))
    values={k:str(s.simplify(v)) for k,v in residuals.items()}
    return dict(residuals=values,all_zero=set(values.values())=={'0'})


def period_map(generator, period=math.pi, rtol=1e-12, atol=1e-14):
    n=len(generator)
    sol=solve_ivp(lambda t,y:(generator@y.reshape(n,n)).ravel(),(0,period),np.eye(n).ravel(),
                  method='DOP853',rtol=rtol,atol=atol,max_step=math.pi/100)
    if not sol.success or not np.isfinite(sol.y).all():raise ArithmeticError('period integration failed')
    return sol.y[:,-1].reshape(n,n)


def constrained_dipole_map(phase=0.,**kwargs):
    basis=null_space(constraints(0.,phase))
    M=period_map(DIPOLE,**kwargs)
    return basis.T@M@basis, float(np.linalg.norm(constraints(math.pi,phase)@M@basis))


def geometry_from_jets(metric,dg,ddg,phi,derivative,second,a,kap,E):
    inv,connection,ricci,G,R=curvature_from_jets(metric,dg,ddg)
    stress=np.zeros((4,4)); kg=[]
    for i in range(4):
        hessian=second[i]-np.einsum('kij,k->ij',connection,derivative[i])
        stress+=improved_stress(metric,inv,G,phi[i],derivative[i],hessian)
        kg.append(np.sum(inv*hessian)-R*phi[i]/6)
    frame=np.zeros((4,4));frame[0,0]=a;frame[1:,1:]=a*E
    fi=np.linalg.inv(frame)
    tm=frame@inv@stress@fi;gm=frame@inv@G@fi
    return dict(stress=tm,einstein=gm,residual=gm+3/(2*a*a)*np.eye(4)-kap*tm,
                KG=np.array(kg),R=R,phi=phi)


def spatial_jets(angles):
    x,dx,ddx,E,dE,ddE=cm.coframe_jets(tuple(angles))
    g=E.T@E
    dg=np.array([de.T@E+E.T@de for de in dE])
    ddg=np.array([[ddE[c,d].T@E+dE[c].T@dE[d]+dE[d].T@dE[c]+E.T@ddE[c,d]
                    for d in range(3)] for c in range(3)])
    return x,dx,ddx,E,g,dg,ddg


def frw_geometry(Ajets,qjets,a=1.,kap=1.,angles=(.83,1.07,.61)):
    """Off-shell FRW coordinate curvature and four independent field stresses."""
    A,Ad,Add=Ajets;q,qd,qdd=qjets
    if min(A,a,kap)<=0:raise ValueError('positive scale/radius/kappa required')
    x,dx,ddx,E,g,gs,gss=spatial_jets(angles)
    unit=np.zeros((4,4));unit[0,0]=-1;unit[1:,1:]=g
    metric=A*A*unit;dg=np.zeros((4,4,4));ddg=np.zeros((4,4,4,4))
    dg[0]=2*A*Ad*unit;ddg[0,0]=2*(Ad*Ad+A*Add)*unit
    for c in range(3):
        dg[c+1,1:,1:]=A*A*gs[c]
        ddg[0,c+1,1:,1:]=ddg[c+1,0,1:,1:]=2*A*Ad*gs[c]
        for d in range(3):ddg[c+1,d+1,1:,1:]=A*A*gss[c,d]
    P=q/A;Pd=qd/A-q*Ad/A**2
    Pdd=qdd/A-2*qd*Ad/A**2-q*Add/A**2+2*q*Ad*Ad/A**3
    derivative=np.column_stack((Pd*x,P*dx));second=np.zeros((4,4,4))
    second[:,0,0]=Pdd*x;second[:,0,1:]=second[:,1:,0]=Pd*dx;second[:,1:,1:]=P*ddx
    return geometry_from_jets(metric,dg,ddg,P*x,derivative,second,a,kap,E)


def expected_frw(Ajets,qjets,a=1.,kap=1.,angles=(.83,1.07,.61)):
    A,Ad,Add=Ajets;q,qd,qdd=qjets
    x=cm.coframe_jets(tuple(angles))[0]
    rho=(qd*qd+4*q*q)/(2*A**4)
    pressure=(qd*qd-2*q*qdd-4*q*q)/(6*A**4)
    T=np.diag([-rho,pressure,pressure,pressure])
    G=np.diag([-3*(Ad*Ad+A*A)/A**4]+[(Ad*Ad-2*A*Add-A*A)/A**4]*3)
    return dict(stress=T,einstein=G,residual=G+3/(2*a*a)*np.eye(4)-kap*T,
                KG=-(qdd+4*q)*x/A**3,R=6*(Add+A)/A**3)


def expected_homogeneous(tau,r,z,a=1.,kap=1.,phase=0.,angles=(.83,1.07,.61)):
    r,rd,rdd=r;z,zd,zdd=z;p,pd,*_=field_jets(tau,phase)
    delta_energy=pd*zd+4*p*z
    rho=delta_energy-6*r;pressure=(pd*zd-p*zdd)/3-2*r
    T=np.diag([-rho,pressure,pressure,pressure])/(kap*a*a)
    G=np.diag([6*r]+[2*r-2*rdd]*3)/(a*a)
    x=cm.coframe_jets(tuple(angles))[0]
    return dict(stress=T,einstein=G,residual=G-kap*T,
                KG=-(zdd+4*z)*x/(math.sqrt(kap)*a*a),R=6*(rdd-2*r)/(a*a))


def exact_scale(tau,epsilon,a=1.,sign=1.):
    """Local exact fixed-energy branch; raises before a pole or A<=0."""
    w=epsilon/(2+epsilon)*math.exp(sign*math.sqrt(2)*tau)
    if abs(w)>=1:raise ValueError('outside this local positive-scale chart')
    A=a*(1+w)/(1-w)
    return np.array([A,sign*(A*A-a*a)/(math.sqrt(2)*a),-A+A**3/a**2])


def gauge_jets(tau,phase=0.):
    """Minus Lie derivative by xi=T Y d_tau + L grad Y (not a solution ansatz)."""
    T=np.array([math.sin(tau),math.cos(tau),-math.sin(tau)])*.13
    L=np.array([math.cos(.7*tau),-.7*math.sin(.7*tau),-.49*math.cos(.7*tau),.343*math.sin(.7*tau)])*.17
    T3=-.13*math.cos(tau);p,pd,pdd,pddd=field_jets(tau,phase)
    u=-np.array([T[0]*pd,T[1]*pd+T[0]*pdd,T[2]*pd+2*T[1]*pdd+T[0]*pddd])
    v=-np.array([L[0]*p,L[1]*p+L[0]*pd,L[2]*p+2*L[1]*pd+L[0]*pdd])
    return dict(u=u,v=v,lapse=-np.array([T[1],T[2],T3]),spatial=L[:3],shift=T-L[1:])


def dipole_geometry(tau,u,v,lapse=(0.,0.,0.),spatial=(0.,0.,0.),shift=(0.,0.,0.),
                    epsilon=0.,a=1.,kap=1.,phase=0.,direction=(1.,0.,0.,0.),angles=(.83,1.07,.61)):
    """Full lapse/shift/trace metric jets and scalar-pair jets, before gauge fixing."""
    x,dx,ddx,E,g,gs,gss=spatial_jets(angles)
    d=np.asarray(direction,dtype=float);Y=d@x;dy=d@dx;ddy=np.einsum('i,iab->ab',d,ddx)
    dddx=np.asarray(cm._embedding_functions()(*angles)[3]);dddy=np.einsum('i,iabc->abc',d,dddx)
    metric=np.zeros((4,4));metric[0,0]=-a*a;metric[1:,1:]=a*a*g
    dg=np.zeros((4,4,4));ddg=np.zeros((4,4,4,4))
    for c in range(3):
        dg[c+1,1:,1:]=a*a*gs[c]
        for e in range(3):ddg[c+1,e+1,1:,1:]=a*a*gss[c,e]
    def add(i,j,f,K,Ks,Kss):
        metric[i,j]+=epsilon*a*a*f[0]*K
        dg[0,i,j]+=epsilon*a*a*f[1]*K
        dg[1:,i,j]+=epsilon*a*a*f[0]*Ks
        ddg[0,0,i,j]+=epsilon*a*a*f[2]*K
        ddg[0,1:,i,j]+=epsilon*a*a*f[1]*Ks
        ddg[1:,0,i,j]+=epsilon*a*a*f[1]*Ks
        ddg[1:,1:,i,j]+=epsilon*a*a*f[0]*Kss
    add(0,0,-2*np.asarray(lapse),Y,dy,ddy)
    for i in range(3):
        add(0,i+1,shift,dy[i],ddy[i],dddy[i]);add(i+1,0,shift,dy[i],ddy[i],dddy[i])
        for j in range(3):
            ks=gs[:,i,j]*Y+g[i,j]*dy
            kss=gss[:,:,i,j]*Y+np.outer(gs[:,i,j],dy)+np.outer(dy,gs[:,i,j])+g[i,j]*ddy
            add(i+1,j+1,2*np.asarray(spatial),g[i,j]*Y,ks,kss)
    p,pd,pdd,*_=field_jets(tau,phase);u=np.asarray(u);v=np.asarray(v);w=u-v
    xy=x*Y;dxy=dx*Y+x[:,None]*dy
    ddxy=ddx*Y+np.einsum('ia,b->iab',dx,dy)+np.einsum('ib,a->iab',dx,dy)+x[:,None,None]*ddy
    phi=(p*x+epsilon*(w[0]*xy+v[0]*d))/math.sqrt(kap)
    derivative=np.zeros((4,4));derivative[:,0]=pd*x+epsilon*(w[1]*xy+v[1]*d)
    derivative[:,1:]=p*dx+epsilon*w[0]*dxy
    second=np.zeros((4,4,4));second[:,0,0]=pdd*x+epsilon*(w[2]*xy+v[2]*d)
    second[:,0,1:]=second[:,1:,0]=pd*dx+epsilon*w[1]*dxy
    second[:,1:,1:]=p*ddx+epsilon*w[0]*ddxy
    return geometry_from_jets(metric,dg,ddg,phi,derivative/math.sqrt(kap),second/math.sqrt(kap),a,kap,E)


def expected_dipole(tau,u,v,lapse=(0.,0.,0.),a=1.,kap=1.,phase=0.,direction=(1.,0.,0.,0.),angles=(.83,1.07,.61)):
    c=dipole_coefficients(tau,u,v,lapse,phase)
    x,dx,_,E,*_=spatial_jets(angles);d=np.asarray(direction);Y=d@x
    grad=np.linalg.solve(E.T,d@dx)
    T=np.zeros((4,4));T[0,0]=-c['rho']*Y
    T[0,1:]=-c['flux']*grad;T[1:,0]=c['flux']*grad
    T[1:,1:]=c['pressure']*Y*np.eye(3)
    G=np.zeros((4,4));G[1:,1:]=-2*lapse[0]*Y*np.eye(3)
    return dict(stress=T/(kap*a*a),einstein=G/(a*a),residual=(G-T)/(a*a),
                KG=(c['kg_xy']*x*Y+c['kg_d']*d)/(math.sqrt(kap)*a*a),R=c['curvature']*Y/(a*a))


def valid_map_pair(values,expected):
    try:
        maps=np.asarray(values,dtype=float)
        if maps.shape!=(2,2,2) or not np.isfinite(maps).all():return False
        scale=max(1.,np.linalg.norm(expected))
        if np.max([np.linalg.norm(m-expected)/scale for m in maps])>=1e-8:return False
        if np.linalg.norm(maps[0]-maps[1])/scale>=1e-8:return False
        if np.max(abs(np.linalg.det(maps)-1))>=1e-7:return False
        return True
    except (TypeError,ValueError,OverflowError):return False


def valid_period_evidence(evidence):
    try:
        return (valid_map_pair(evidence['homogeneous'],expm(HOMOGENEOUS*math.pi))
                and valid_map_pair(evidence['dipole'],-np.eye(2)))
    except (KeyError,TypeError):return False


def verdict(checks,evidence=None):
    failed=[k for k in REQUIRED_CHECKS if checks.get(k) is not True]
    if set(checks)!=set(REQUIRED_CHECKS):failed.append('gate_schema')
    if not valid_period_evidence(evidence):failed.append('period_evidence')
    keys=('homogeneous_physical_block','homogeneous_constraint_completion','dipole_cover_block','dipole_restricted_admissibility')
    if failed:out=dict.fromkeys(keys,'UNRESOLVED')
    else:out=dict(zip(keys,('HYPERBOLIC_GROWING_MODE','EXACT_FRW_CONTINUATION',
                           'CONSTRAINED_LINEAR_NEUTRAL_SEMISIMPLE','EXCLUDED_UNDER_STATED_ANTIPODAL_RESTRICTIONS')))
    return {**out,**SCOPE,'failed_checks':sorted(set(failed))}
