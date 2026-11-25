"""Symmetric intensity loads and a reciprocal overlap pointer on antipodal paths."""
import json

import numpy as np
from scipy.optimize import brentq, minimize

from .antipodal_two_time import arc_action, hopf, rank, unit

BARE = np.pi**2/2
LOADS = (-.3, -.1, 0., .1, .3)
SUSCEPTIBILITIES = (0., 2., 3.9, 4., 4.1, 6., 10.)


def direct_certificate(kappa, floating):
    """Action improvement relative to the entire complementary class, at u=.5
    in the pinned problem or against any complementary floating history."""
    delta = kappa/4
    excess = np.sin(delta) if floating else np.sin(2*delta)/2
    improvement = kappa*excess - 2*delta**2
    return dict(trial_displacement=delta, trial_gain=float(improvement),
                minimum_defect_bound=float(improvement/abs(kappa)) if kappa else 0.)


def pointer_prediction(chi):
    if chi < 0:
        raise ValueError("g²/k must be nonnegative")
    # Nonzero root obeys sinc(2*delta) = 4/chi with radians as argument.
    delta = brentq(lambda d: np.sinc(2*d/np.pi)-4/chi, 0., np.pi/2) if chi>4 else 0.
    overlap = np.sin(delta)
    return dict(overlap=float(overlap), pointer_abs=float(np.sqrt(chi)*overlap),
                action=float(BARE + 2*delta**2-chi*overlap**2/2),
                invariant_mean_defect=float(overlap/2), balanced_doubles=float(delta/np.pi),
                balanced_misses=float(delta/np.pi))


def arc_gradient(nodes):
    """Node gradient of #367's arc action at lambda=2 on unit endpoints."""
    n = len(nodes)-1
    x,y = nodes[:-1],nodes[1:]
    cos = np.clip((x*y).sum(axis=1),-1.,1.)
    d = np.arccos(cos)
    h = 1/np.sinc(d/np.pi)
    dh = np.full_like(d,-1/3.)
    regular = d>1e-4
    z = d[regular]
    dh[regular] = (z*np.cos(z)-np.sin(z))/np.sin(z)**3
    fiber = (hopf(x)*y).sum(axis=1)
    dc = n*(-h+fiber**2*h*dh)
    df = n*(fiber*h**2)
    grad = np.zeros_like(nodes)
    grad[:-1] += dc[:,None]*y-df[:,None]*hopf(y)
    grad[1:] += dc[:,None]*x+df[:,None]*hopf(x)
    return grad


def solve(kind, strength, floating, seed, n=8):
    """Unconstrained node coordinates with normalization and antipodal gluing.
    For the overlap load, include Q as an independent optimizer variable."""
    rng = np.random.default_rng(seed)
    x = np.array([1.,0.,1.,0.])/np.sqrt(2)
    v = rng.normal(size=4)
    v = unit(v-v.dot(x)*x)
    t = np.arange(n)[:,None]*np.pi/n
    guess = unit(np.cos(t)*x+np.sin(t)*v+.05*rng.normal(size=(n,4)))
    start = 0 if floating else 1
    initial = guess[start:].ravel()
    if kind=="pointer":
        initial = np.r_[initial,rng.normal()]
    elif kind=="ready_pair":
        initial = np.r_[initial,rng.normal(size=2)]
    def unpack(params):
        raw = params[:4*(n-start)].reshape((-1,4))
        nodes = unit(raw)
        if not floating:
            nodes = np.vstack((x,nodes))
        return raw, np.vstack((nodes,-nodes[0]))
    def objective(params):
        raw, nodes = unpack(params)
        action = arc_action(nodes,2.)
        grad = arc_gradient(nodes)
        a,b = nodes[0],nodes[n//2]
        if kind=="direct":
            action -= strength*(rank(a)+rank(b))
            grad[0,:2] -= 2*strength*a[:2]
            grad[n//2,:2] -= 2*strength*b[:2]
        elif kind=="ready_pair":
            # Unit mass, time and terminal spring. Linear pointer paths
            # exactly minimize the kinetic cost for any two terminal values.
            read = params[-2:]
            ranks = np.array([rank(a),rank(b)])
            mismatch = read-ranks
            action += (read.dot(read)+mismatch.dot(mismatch))/2
            grad[0,:2] -= 2*mismatch[0]*a[:2]
            grad[n//2,:2] -= 2*mismatch[1]*b[:2]
        elif kind=="pointer":
            q = params[-1]
            g = np.sqrt(strength)
            overlap = a.dot(b)
            action += q*q/2-g*q*overlap
            grad[0] -= g*q*b
            grad[n//2] -= g*q*a
        else:
            raise ValueError(f"unknown coupling {kind}")
        grad[0] -= grad[-1]
        grad = grad[start:n]
        pts = nodes[start:n]
        gradient = (grad-(grad*pts).sum(axis=1,keepdims=True)*pts)/np.linalg.norm(raw,axis=1,keepdims=True)
        flat = gradient.ravel()
        if kind=="pointer":
            flat = np.r_[flat,q-g*overlap]
        elif kind=="ready_pair":
            flat = np.r_[flat,read+mismatch]
        return action,flat
    result = minimize(objective,initial,jac=True,method="L-BFGS-B",
                      options=dict(ftol=1e-15,gtol=1e-9,maxiter=3000,maxcor=25))
    if not result.success:
        raise RuntimeError(result.message)
    _,nodes = unpack(result.x)
    a,b = nodes[0],nodes[n//2]
    return dict(action=float(result.fun), u=float(rank(a)), v=float(rank(b)),
                real_overlap=float(a.dot(b)), fiber_overlap=float(hopf(a).dot(b)),
                pointer=float(result.x[-1]) if kind=="pointer" else None,
                paired_readout=result.x[-2:].tolist() if kind=="ready_pair" else None,
                gradient_max=float(np.max(abs(result.jac))))


def probe():
    direct = []
    for kappa in LOADS:
        for floating in (False,True):
            direct.append(dict(kappa=kappa, floating=floating,
                               certificate=direct_certificate(kappa,floating),
                               fits=[solve("direct",kappa,floating,s) for s in (365,367,368)]))
    pointer = []
    for chi in SUSCEPTIBILITIES:
        pointer.append(dict(chi=chi, prediction=pointer_prediction(chi),
                            fits=[solve("pointer",chi,True,s) for s in (365,367,368)]))
    # Post-primary pair with dynamical ready states. Compare to the lowest
    # action on the entire complementary class using an explicit trial.
    delta = -1/16
    best_complementary = BARE+1/8
    trial = BARE+2*delta**2+(1+np.sin(delta))**2/8
    gain = best_complementary-trial
    paired = dict(action_improvement=float(gain),
                  minimum_defect_bound=float(1-np.sqrt(1-8*gain)),
                  fits=[solve("ready_pair",1.,True,s) for s in (365,367,368)])
    return dict(segments=8, lam=2, seeds=[365,367,368], direct=direct, pointer=pointer,
                ready_pair_followup=paired)


if __name__ == "__main__":
    print(json.dumps(probe(), indent=2, allow_nan=False))
