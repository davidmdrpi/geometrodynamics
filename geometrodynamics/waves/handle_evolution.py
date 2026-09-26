"""Short-time, coupled Einstein/quartet evolution; no imposed crossing kicks."""
import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicSpline
from . import localized_mouth as initial

L = 5.5
T = .01
PARITY = np.array([1, 1, 1, 1, -1, 1, -1, 1])


def shifted(y, offset):
    n = y.shape[-1]
    indices = np.arange(n)+offset
    return y[:, indices % n]*PARITY[:, None]**(indices//n % 2)


def derivatives(y, dx):
    mm, m, p, pp = [shifted(y, j) for j in (-2, -1, 1, 2)]
    return (mm-8*m+8*p-pp)/(12*dx), (-mm+16*m-30*y+16*p-pp)/(12*dx*dx)


def geometry(y, dy, ddy):
    A, B, k, l, u, v, P, Q = y
    As, Bs, _, ls, us, vs, _, _ = dy
    f = 1-(u*u+v*v)/6
    if not np.isfinite(y).all() or np.min(f) <= .1 or min(np.min(A), np.min(B)) <= .005:
        raise ArithmeticError('registered positivity/domain stop')
    K = k+2*l
    dot = u*P+v*Q
    radial = u*us+v*vs
    Z = (us*us+vs*vs)/f+radial*radial/(6*f*f)
    kinetic = (P*P+Q*Q)/f+dot*dot/(6*f*f)
    current = (P*us+Q*vs)/f+dot*radial/(6*f*f)
    angular = v*v/(f*B*B)
    U = 1.5/f**2
    D = ddy[1]-As/A*Bs
    Rs = -2*D/(A*A*B)
    Ro = (1-(Bs/A)**2-B*D/(A*A))/(B*B)
    rho = (kinetic+Z/(A*A)+2*angular)/2+U
    Lm = (kinetic-Z/(A*A)-2*angular)/2-U
    Ss, So = Z/(A*A)+Lm, angular+Lm
    H = Rs+2*Ro+K*K-k*k-2*l*l-2*rho
    M = -2*ls+2*Bs/B*(k-l)+current
    hs = np.maximum(1, abs(Rs)+2*abs(Ro)+K*K+k*k+2*l*l+2*abs(rho))
    ms = np.maximum(1, 2*abs(ls)+2*abs(Bs/B)*(abs(k)+abs(l))+abs(current))
    return dict(f=f, K=K, dot=dot, radial=radial, Z=Z, current=current,
                angular=angular, U=U, Rs=Rs, Ro=Ro, rho=rho, Ss=Ss, So=So,
                H=H, M=M, Hn=H/hs, Mn=M/ms)


def rhs_jets(y, dy, ddy):
    g = geometry(y, dy, ddy)
    A, B, k, l, u, v, P, Q = y
    As, Bs, _, _, us, vs, _, _ = dy
    f, K = g['f'], g['K']
    connection = 2*Bs/B-As/A
    reaction = -g['dot']/(3*f)
    spatial = g['radial']/(3*f*A*A)
    return np.array([-A*k, -B*l,
                     g['Rs']+K*k-g['Z']/(A*A)-g['U'],
                     g['Ro']+K*l-g['angular']-g['U'], P, Q,
                     K*P+(ddy[4]+connection*us)/(A*A)+reaction*P+spatial*us-u/f,
                     K*Q+(ddy[5]+connection*vs)/(A*A)-2*v/(B*B)+reaction*Q+spatial*vs-v/f]), g


def interpolate(values, positions):
    """Four-point cubic interpolation of periodic geometric observables."""
    n = values.shape[-1]
    z = (np.asarray(positions)+L)/(2*L)*n
    i = np.floor(z).astype(int)
    x = z-i
    weights = (-x*(x-1)*(x-2)/6, (x+1)*(x-1)*(x-2)/2,
               -(x+1)*x*(x-2)/2, (x+1)*x*(x-1)/6)
    return sum(w*values[..., (i+j) % n] for j, w in zip((-1, 0, 1, 2), weights))


def particle_rhs(y, dy, g, particles):
    x, p = particles
    A, As, f, fs = interpolate(np.array([y[0], dy[0], g['f'],
                                        -(y[4]*dy[4]+y[5]*dy[5])/3]), x)
    H = np.sqrt(1/f+p*p/(A*A))
    return np.array([p/(A*A*H), (fs/(2*f*f)+p*p*As/A**3)/H])


def tube(y, dy, g):
    points = np.linspace(L-.2, L, 129)
    A, B = y[:2]
    W = A*B*B
    inventory, flux, work = interpolate(np.array([
        -W*g['current'], W*g['Ss'],
        W*(dy[0]/A*g['Ss']+2*dy[1]/B*g['So'])]), points)
    P = 4*np.pi*simpson(inventory, x=points)
    flux_rate = -4*np.pi*(flux[-1]-flux[0])
    work_rate = 4*np.pi*simpson(work, x=points)
    return P, flux_rate+work_rate, flux_rate


def prepared(record, profiles, n):
    s = np.linspace(-L, L, n, endpoint=False)
    d = initial.Data(record, profiles)
    theta, _, psi, _, a, p = d.reduced(s)
    u, v = initial.Q*np.sin(theta), initial.Q*np.cos(theta)
    return s, np.array([psi**2, psi**2, 2*a*psi**-6, -a*psi**-6,
                        u, v, psi**-6*p*v, -psi**-6*p*u])


def worldline_summary(times, samples):
    samples = np.asarray(samples)
    results = []
    for j in range(6):
        sign = 1 if j < 3 else -1
        x = sign*samples[:, 0, j]
        if not x[0] < L < x[-1]:
            results.append(dict(crossed=False))
            continue
        # Monotone inverse interpolation locates the geometric surface.
        crossing = float(CubicSpline(x, times)(L))
        record = dict(crossed=True, time=crossing,
                      p_hat=float(CubicSpline(times, samples[:, 2, j])(crossing)), windows=[])
        for width in (T/16, T/32, T/64):
            if crossing-width < 0 or crossing+width > T:
                raise ArithmeticError('registered momentum window outside evolution')
            row = dict(half_width=width)
            for name, index in [('p_s', 1), ('p_hat', 2)]:
                curve = CubicSpline(times, samples[:, index, j])
                row[name] = float(curve(crossing+width)-curve(crossing-width))
            record['windows'].append(row)
        results.append(record)
    return results


def run(record, profiles, n):
    s, y = prepared(record, profiles, n)
    dx = 2*L/n
    steps = 200*n//512
    dt = T/steps
    x = np.array([L-.05]*3+[-L+.05]*3)
    speed = np.array([.25, .5, .75, -.25, -.5, -.75])
    f = 1-(y[4]**2+y[5]**2)/6
    p = interpolate(y[0]/np.sqrt(f), x)*speed/np.sqrt(1-speed*speed)
    particles = np.array([x, p])
    balance = np.zeros(2)
    rows, trajectories, snapshots = [], [], []
    def derivative(state, probes):
        dy, ddy = derivatives(state, dx)
        derivative, g = rhs_jets(state, dy, ddy)
        return derivative, particle_rhs(state, dy, g, probes), np.array(tube(state, dy, g)[1:])
    for step in range(steps+1):
        if step % (n//512) == 0:
            dy, ddy = derivatives(y, dx)
            g = geometry(y, dy, ddy)
            P, _, _ = tube(y, dy, g)
            jordan_r = y[1]/np.sqrt(g['f'])
            A, f = interpolate(np.array([y[0], g['f']]), particles[0])
            phat = particles[1]*np.sqrt(f)/A
            rows.append([step*dt, np.max(abs(g['Hn'])), np.max(abs(g['Mn'])),
                         np.max(abs(g['H'])), np.max(abs(g['M'])), np.min(y[0]), np.min(y[1]),
                         np.min(g['f']), jordan_r[0], jordan_r[n//2], jordan_r[0]/jordan_r[n//2],
                         P, balance[0], balance[1], min(jordan_r[1], jordan_r[-1])-jordan_r[0],
                         float(np.mean(np.sqrt(y[4, abs(s)<1]**2+y[5, abs(s)<1]**2)))])
            trajectories.append(np.array([particles[0], particles[1], phat, phat/np.sqrt(1+phat**2)]).tolist())
            # Common physical grid at every output time for evolution convergence.
            snapshots.append(y[:, ::n//512].copy())
        if step == steps:
            break
        a, b, c = derivative(y, particles)
        aa, bb, cc = derivative(y+dt*a/2, particles+dt*b/2)
        aaa, bbb, ccc = derivative(y+dt*aa/2, particles+dt*bb/2)
        aaaa, bbbb, cccc = derivative(y+dt*aaa, particles+dt*bbb)
        y += dt/6*(a+2*aa+2*aaa+aaaa)
        particles += dt/6*(b+2*bb+2*bbb+bbbb)
        balance += dt/6*(c+2*cc+2*ccc+cccc)
    return dict(N=n, eta=record['eta'], steps=steps, dt=dt,
                diagnostics=np.asarray(rows).tolist(), trajectories=trajectories,
                final_fields=y.tolist(), comparison_fields=np.asarray(snapshots),
                crossings=worldline_summary(np.asarray(rows)[:, 0], trajectories))


def homogeneous_controls():
    out = []
    for epsilon in (0., 1e-4, -1e-4):
        a = 1+epsilon
        y0 = [a, (a*a-1)/np.sqrt(2), initial.Q, 0.]
        def rhs(t, y):
            a, ad, q, qd = y
            f = 1-q*q/(6*a*a)
            return np.array([ad, -a+a**3, qd, -4*q])/(a*np.sqrt(f))
        times = np.linspace(0, T, 201)
        sol = solve_ivp(rhs, (0, T), y0, method='DOP853', t_eval=times, rtol=1e-11, atol=1e-13)
        if not sol.success:
            raise ArithmeticError('homogeneous control failed')
        a, ad, q, qd = sol.y
        constraint = ad*ad+a*a-.5*a**4-(qd*qd+4*q*q)/6
        out.append(dict(epsilon=epsilon, states=sol.y.tolist(), constraint_max=float(np.max(abs(constraint)))))
    return out


def validate_equations():
    """Independent target compatibility, exact breathing and coordinate Ricci."""
    import sympy as sp
    import mpmath as mp
    from .coordinate_budget import coordinate_check
    z = sp.symbols('z0:4')
    f = 1-sum(v*v for v in z)/6
    metric = sp.eye(4)/f+sp.Matrix(z)*sp.Matrix(z).T/(6*f*f)
    compatible = True
    for i in range(4):
        for j in range(i, 4):
            for k in range(4):
                predicted = sum(((int(m == k)*z[i]+int(m == i)*z[k])*metric[m,j]
                                 +(int(m == k)*z[j]+int(m == j)*z[k])*metric[i,m])/(6*f)
                                for m in range(4))
                compatible &= sp.simplify(sp.diff(metric[i,j], z[k])-predicted) == 0
    # Angular target connection vanishes because phi dot d_angle(phi)=0;
    # the sphere Laplacian sends n to -2n, so no omitted angular mode appears.
    th, az, u, v = sp.symbols('th az u v', real=True)
    n = sp.Matrix([sp.sin(th)*sp.cos(az), sp.sin(th)*sp.sin(az), sp.cos(th)])
    angular = sp.simplify((n.T*sp.diff(n, th))[0]) == 0
    angular &= sp.simplify((n.T*sp.diff(n, az))[0]) == 0
    lap = sp.diff(n, th, 2)+sp.cot(th)*sp.diff(n, th)+sp.diff(n, az, 2)/sp.sin(th)**2
    angular &= all(sp.trigsimp(a+2*b) == 0 for a, b in zip(lap, n))
    round_error = 0.
    geometric_error = 0.
    with mp.workdps(50):
        def round_state(tau, s):
            q = mp.sqrt(3)/2
            R = q*mp.cos(2*tau)
            Rtau = -2*q*mp.sin(2*tau)
            f = 1-R*R/6
            ftau = -R*Rtau/3
            k = -ftau/(2*f**mp.mpf('1.5'))
            return [mp.sqrt(f)/mp.cosh(s)]*2+[k,k,R*mp.tanh(s),R/mp.cosh(s),
                       Rtau/mp.sqrt(f)*mp.tanh(s),Rtau/mp.sqrt(f)/mp.cosh(s)]
        for tau in map(mp.mpf, ('0', '.003', '.01')):
            for s in map(mp.mpf, ('.2', '.7', '1.3', '3.1')):
                state = round_state(tau,s)
                first = [mp.diff(lambda x: round_state(tau,x)[i], s) for i in range(8)]
                second = [mp.diff(lambda x: round_state(tau,x)[i], s,2) for i in range(8)]
                f = 1-(state[4]**2+state[5]**2)/6
                exact = [mp.diff(lambda x: round_state(x,s)[i], tau)/mp.sqrt(f) for i in range(8)]
                measured, _ = rhs_jets(np.array(state,dtype=float),np.array(first,dtype=float),np.array(second,dtype=float))
                round_error = max(round_error,float(np.max(abs(measured-np.array(exact,dtype=float)))))
        s = mp.mpf('.7')
        def model(x):
            return [.8+.1*mp.sin(x), .6+.1*mp.cos(x), .2*mp.sin(x), .1*mp.cos(x),
                    mp.mpf('.1'),mp.mpf('.2'),mp.mpf('.03'),mp.mpf('.04')]
        y = model(s)
        first = [mp.diff(lambda x:model(x)[i],s) for i in range(8)]
        second = [mp.diff(lambda x:model(x)[i],s,2) for i in range(8)]
        def values(x,t,azimuth):
            A,B,k,l,u,v,P,Q = model(x)
            return [A*A,B*B,B*B*mp.sin(t)**2,A*A*k,B*B*l,B*B*l*mp.sin(t)**2]+[mp.mpf(0)]*8
        independent = coordinate_check(values,[s,mp.mpf('.91'),mp.mpf('.37')])
        g = geometry(np.array(y,dtype=float),np.array(first,dtype=float),np.array(second,dtype=float))
        expected_m = -2*first[3]+2*first[1]/y[1]*(y[2]-y[3])
        geometric_error = max(abs(float(independent['R'])-float(g['Rs']+2*g['Ro'])),
                              abs(float(independent['divergence'][0]-expected_m)))
    return dict(target_metric_compatible=bool(compatible), equivariant_closure=bool(angular),
                round_evolution_error=round_error, coordinate_geometry_error=geometric_error,
                initial_norm_acceleration=-4*initial.Q/initial.F)
