"""Mixed-precision augmented shooting and high-precision closure diagnostics."""
import numpy as np
import gmpy2 as gp
from gmpy2 import mpfr as mp
from . import r3_breaking as b


def residual(P,Z,lam,g,t0,c0):
    n=len(g);q=len(Z);out=[]
    for i in range(q):out.extend(a-c for a,c in zip(P(Z[i]),Z[(i+1)%q]))
    for k in range(n):out[(q-1)*n+k]-=lam*g[k]
    out.append(gp.fsum(t0[k]*(Z[0][k]-c0[k]) for k in range(n)))
    return out


def solve(P,point,bits,tol='1e-45',start=None,maxit=24):
    """Residual and state updates at high precision; fixed double Jacobian.

    Recompute the residual at the returned nodes, including iteration exhaustion.
    Retain exact string slices when supplied. A failed contraction is a failed
    solve, not permission to silently accept an old residual.
    """
    with gp.context(gp.get_context(),precision=bits):
        src=start or point;Z=[[mp(v) for v in z] for z in src['nodes']];lam=mp(src['lam'])
        g,t,c=([mp(v) for v in point[k]] for k in ('g','t0','c0'))
        J=np.array(point['J'],float);history=[];n=len(g)
        for _ in range(maxit):
            F=residual(P,Z,lam,g,t,c);r=max(map(abs,F));history.append(str(r))
            if r<=mp(tol):break
            correction=np.linalg.solve(J,-np.array(F,float))
            if not np.isfinite(correction).all():raise ArithmeticError('nonfinite correction')
            for i,z in enumerate(Z):
                for k in range(n):z[k]+=mp(correction[i*n+k])
            lam+=mp(correction[-1])
        F=residual(P,Z,lam,g,t,c);r=max(map(abs,F))
        return dict(nodes=[[str(v) for v in z] for z in Z],lam=str(lam),residual=str(r),
                    ok=bool(r<=mp(tol)),history=history,iterations=len(history),
                    inverse_row_norm=float(np.linalg.norm(np.linalg.inv(J)[-1],ord=1)))


def closure(P,Z):
    return [a-c for i,z in enumerate(Z) for a,c in zip(P(z),Z[(i+1)%len(Z)])]


def closure_jacobian(P,Z,step):
    n=len(Z[0]);q=len(Z);h=mp(step);out=[[mp(0) for _ in range(n*q)] for _ in range(n*q)]
    for i,z in enumerate(Z):
        for k in range(n):
            plus,minus=list(z),list(z);plus[k]+=h;minus[k]-=h
            col=[(a-c)/(2*h) for a,c in zip(P(plus),P(minus))]
            for row in range(n):out[i*n+row][i*n+k]=col[row]
        for k in range(n):out[i*n+k][((i+1)%q)*n+k]-=1
    return out


def inverse(A):
    """Partial-pivot Gauss-Jordan at the caller's precision."""
    n=len(A);a=[list(row)+[mp(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        pivot=max(range(k,n),key=lambda i:abs(a[i][k]))
        if not a[pivot][k]:raise ArithmeticError('singular closure Jacobian')
        a[k],a[pivot]=a[pivot],a[k];scale=a[k][k];a[k]=[v/scale for v in a[k]]
        for i in range(n):
            if i!=k:
                scale=a[i][k];a[i]=[x-scale*y for x,y in zip(a[i],a[k])]
    return [row[n:] for row in a]


def matrix_norm(A):return max(gp.fsum(map(abs,row)) for row in A)


def matvec(A,v):return [gp.fsum(x*y for x,y in zip(row,v)) for row in A]
