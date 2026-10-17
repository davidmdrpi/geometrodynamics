"""Arbitrary-precision diagonal Einstein-quartet maps in two formulations.

Primary: canonical shear momenta p=h*x'. Check: shear velocities x'.
The coefficient recurrences use distinct state equations and the crossing
solvers are Newton and bracketed bisection respectively. Shared scalar series
arithmetic is explicit; agreement is a numerical check, not interval arithmetic.
"""
import gmpy2 as gp
from gmpy2 import mpfr as mp

PRIMARY = dict(bits=224, order=48, eps='1e-56', formulation='momentum')
CHECK = dict(bits=336, order=64, eps='1e-80', formulation='velocity')


def conv(a,b,k):
    return gp.fsum(a[j]*b[k-j] for j in range(k+1))


def horner(a,t):
    out=a[-1]
    for v in reversed(a[:-1]):out=out*t+v
    return out


def basis():
    s6,s2=gp.sqrt(mp(6)),gp.sqrt(mp(2))
    return [[1/s6,1/s6,-2/s6],[1/s2,-1/s2,mp(0)]]


def coefficients(y,order,formulation):
    if formulation not in ('momentum','velocity'):raise ValueError('unknown formulation')
    B=basis();A,V,q,W,x,u,a,b=([mp(v)] for v in y)
    xs=[x,u];s=[a,b]
    inv=[[gp.exp(-2*sum(B[j][i]*xs[j][0] for j in range(2)))] for i in range(3)]
    sq=[[gp.exp(4*sum(B[j][i]*xs[j][0] for j in range(2)))] for i in range(3)]
    aa,qq,h,hp,vv,r,tt,ell,total,potential,force=[],[],[],[],[[],[]],[],[],[],[],[],[[],[]]
    accel=[[],[]]
    for k in range(order):
        if k:
            for i in range(3):
                for target,scale in ((inv,-2),(sq,4)):
                    target[i].append(gp.fsum(j*scale*sum(B[d][i]*xs[d][j] for d in range(2))*target[i][k-j] for j in range(1,k+1))/k)
        aa.append(conv(A,A,k));qq.append(conv(q,q,k));h.append(aa[k]-qq[k]/6)
        hp.append(2*conv(A,V,k)-conv(q,W,k)/3)
        for d in range(2):
            vv[d].append((s[d][k]-gp.fsum(h[j]*vv[d][k-j] for j in range(1,k+1)))/h[0] if formulation=='momentum' else s[d][k])
        tt.append(sum(e[k] for e in inv));r.append(4*tt[k]-2*sum(e[k] for e in sq))
        ell.append(conv(vv[0],vv[0],k)+conv(vv[1],vv[1],k));total.append(r[k]+ell[k])
        potential.append(tt[k]+total[k]/6)
        for d in range(2):
            f=gp.fsum(B[d][i]*(conv(qq,inv[i],k)-4*conv(h,inv[i],k)-4*conv(h,sq[i],k)) for i in range(3))
            force[d].append(f)
            if formulation=='velocity':
                accel[d].append((f-conv(hp,vv[d],k)-gp.fsum(h[j]*accel[d][k-j] for j in range(1,k+1)))/h[0])
        A.append(V[k]/(k+1));V.append((-conv(A,total,k)/6+conv(aa,A,k))/(k+1))
        q.append(W[k]/(k+1));W.append(-conv(potential,q,k)/(k+1))
        for d in range(2):
            xs[d].append(vv[d][k]/(k+1))
            s[d].append((force[d][k] if formulation=='momentum' else accel[d][k])/(k+1))
    return [A,V,q,W,x,u,a,b]


def to_state(z,formulation):
    A,pA,x,p1,u,p2=map(mp,z);B=basis();h=A*A;V=-pA/6
    inv=[gp.exp(-2*(x*B[0][i]+u*B[1][i])) for i in range(3)]
    sq=[gp.exp(4*(x*B[0][i]+u*B[1][i])) for i in range(3)]
    r=4*sum(inv)-2*sum(sq);ell=(p1*p1+p2*p2)/(h*h)
    w2=6*V*V-h*ell+h*r-3*A**4
    if A<=0 or w2<=0:raise ValueError('invalid section chart')
    return [A,V,mp(0),-gp.sqrt(w2),x,u,p1 if formulation=='momentum' else p1/h,p2 if formulation=='momentum' else p2/h]


def constraint(y,formulation):
    A,V,q,W,x,u,a,b=y;B=basis();h=A*A-q*q/6
    inv=[gp.exp(-2*(x*B[0][i]+u*B[1][i])) for i in range(3)]
    sq=[gp.exp(4*(x*B[0][i]+u*B[1][i])) for i in range(3)]
    ell=(a*a+b*b)/(h*h) if formulation=='momentum' else a*a+b*b
    return -3*V*V+W*W/2+h*ell/2-h*(4*sum(inv)-2*sum(sq))/2+q*q*sum(inv)/2+3*A**4/2


def crossing(cs,h,formulation):
    if formulation=='momentum':
        s=h*cs[2][0]/(cs[2][0]-horner(cs[2],h))
        derivative=[j*v for j,v in enumerate(cs[2])][1:]
        for _ in range(100):
            delta=horner(cs[2],s)/horner(derivative,s);s-=delta
            if abs(delta)<mp(2)**(-gp.get_context().precision+12):break
        else:raise ArithmeticError('Newton crossing failed')
    else:
        lo,hi=mp(0),h;sign=horner(cs[2],lo)>0
        for _ in range(gp.get_context().precision+4):
            s=(lo+hi)/2
            if (horner(cs[2],s)>0)==sign:lo=s
            else:hi=s
        s=(lo+hi)/2
    if not 0<=s<=h:raise ArithmeticError('crossing outside step')
    return s


def clock_map(z,bits=224,order=48,eps='1e-56',formulation='momentum',half=False):
    with gp.context(gp.get_context(),precision=bits):
        y=to_state(z,formulation);time=mp(0);factor=mp(eps)**(mp(1)/order)*mp('.7');direction=1 if half else -1
        for step in range(10000):
            if time>mp('4.0'):raise ArithmeticError('no crossing by cap')
            if y[0]<=0 or y[0]**2-y[2]**2/6<=0:raise ValueError('left positive chart')
            cs=coefficients(y,order,formulation)
            radii=[abs(c[k])**(-mp(1)/k) for c in cs for k in (order-1,order) if c[k]]
            h=min(mp('.15'),min(radii)*factor)
            q0,q1=y[2],horner(cs[2],h)
            if time+h>mp('.5') and direction*q0<0 and direction*q1>=0:
                s=crossing(cs,h,formulation);y=[horner(c,s) for c in cs];time+=s
                A,V,q,W,x,u,a,b=y;hh=A*A-q*q/6
                out=[A,-6*V,x,a if formulation=='momentum' else hh*a,u,b if formulation=='momentum' else hh*b]
                return dict(z=out,constraint=abs(constraint(y,formulation)),time=time,steps=step+1)
            y=[horner(c,h) for c in cs];time+=h
        raise ArithmeticError('Taylor step cap')
