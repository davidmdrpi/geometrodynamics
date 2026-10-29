"""PR #344's incident packet on the glued Klein surface of #338, instead of a bare S².

The PR reads the throat drive from u(A, t)/peak on an unholed sphere, with gain .18.
Here the same S² initial data at t = 0 (far from both mouths) is evolved with the
conformal operator -Δ + 1/4 on (i) the bare sphere and (ii) the sphere with mouths
at ±x joined by #338's Klein neck r(s) = sqrt(s² + b²), using the repository's
intrinsic cotangent operators. We record what actually reaches the throat.

    PYTHONPATH=.:scripts python docs/reviews/pr344/glued_drive.py     # from a checkout of #344
"""
import json, sys
import numpy as np
from geometrodynamics.viz.klein_handles import build_surface, HandleSpec
from klein_release_model import KleinRelease

m = KleinRelease()
N_SPHERE = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
T_END, N_OUT = 2.6, 521

def run(handles):
    s = build_surface(handles=handles, n_sphere=N_SPHERE, n_ring=32, ds_target=.01)
    K, M = s.operators()
    sph = s.on_sphere
    mu = s.points[:, 0] / np.linalg.norm(s.points, axis=1)
    a0, p0 = m.incident(0.)
    u = np.where(sph, m.basis(mu) @ a0, 0.); v = np.where(sph, m.basis(mu) @ p0, 0.)
    acc = lambda w: -(K @ w)/M - .25*w
    lam = 0.; x = np.random.default_rng(0).standard_normal(len(u))
    for _ in range(80):
        y = -acc(x); lam = np.linalg.norm(y)/np.linalg.norm(x); x = y/np.linalg.norm(y)
    out_dt = T_END/(N_OUT-1); sub = int(np.ceil(out_dt/(1.8/np.sqrt(1.05*lam)))); dt = out_dt/sub
    neck = s.tube_of >= 0
    tri_neck = neck[s.tris].any(axis=1)
    area, KT = s.local_stiffness()
    nearB = sph & (mu < -np.cos(1.))                  # sphere within 1 rad of mouth B's centre
    tri_B = nearB[s.tris].all(axis=1)
    def energies(u, v):
        pot_t = .5*np.einsum('ti,tij,tj->t', u[s.tris], KT, u[s.tris])
        kin = .5*M*v*v + .5*.25*M*u*u
        return float(kin.sum()+pot_t.sum()), float(kin[neck].sum()+pot_t[tri_neck].sum()), float(kin[nearB].sum()+pot_t[tri_B].sum())
    pole = int(np.argmax(np.where(sph, mu, -2)))       # sphere vertex closest to +x (rim on the glued surface)
    waist = neck & (np.abs(np.nan_to_num(s.s, nan=9.)) <= (np.abs(s.s[neck]).min()+1e-12 if neck.any() else -1))
    rec = dict(t=[], pole=[], waist=[], E=[], E_neck=[], E_B=[])
    # velocity Verlet
    a = acc(u)
    for k in range(N_OUT):
        E, En, EB = energies(u, v)
        rec['t'].append(k*out_dt); rec['pole'].append(float(u[pole])); rec['E'].append(E); rec['E_neck'].append(En); rec['E_B'].append(EB)
        rec['waist'].append(float(np.mean(u[waist])) if waist.any() else np.nan)
        for _ in range(sub):
            v += .5*dt*a; u += dt*v; a = acc(u); v += .5*dt*a
    return s, {k: np.array(w) for k, w in rec.items()}, dict(pole_mu=float(mu[pole]), waist_vertices=int(waist.sum()), dt=dt, n=s.n_vertices)

if __name__ == '__main__':
    bare, rb, ib = run([])
    klein, rk, ik = run([HandleSpec((1, 0, 0), (0, -1, 0), True)])
    t = rb['t']; modal = m.mouth(t)*m.peak
    i = np.argmax(np.abs(rb['pole']))
    res = dict(
        mesh=dict(bare=ib, klein=ik, n_sphere=N_SPHERE),
        bare_pole_peak_vs_modal=float(np.max(rb['pole'])/m.peak), bare_pole_peak_time=float(t[i]),
        bare_energy=[float(rb['E'][0]), float(np.max(np.abs(rb['E']-rb['E'][0])))],
        klein_energy=[float(rk['E'][0]), float(np.max(np.abs(rk['E']-rk['E'][0])))],
        neck_energy_max=float(rk['E_neck'].max()), neck_energy_max_time=float(t[np.argmax(rk['E_neck'])]),
        neck_energy_fraction_max=float(rk['E_neck'].max()/rk['E'][0]),
        waist_peak_over_s2_peak=float(np.nanmax(rk['waist'])/m.peak), waist_peak_time=float(t[np.nanargmax(rk['waist'])]),
        waist_trough_over_s2_peak=float(np.nanmin(rk['waist'])/m.peak),
        energy_near_B_by_t=dict(t=float(t[-1]), bare=float(rb['E_B'][-1]), klein=float(rk['E_B'][-1]),
                                klein_max=float(rk['E_B'].max()), klein_max_time=float(t[np.argmax(rk['E_B'])])),
        pr344_readout_waist_peak_time=float(m.focus+m.length/2),
        contact_gain_needed_at_waist=float((m.neck.par.waist-m.neck.level_radius)),
    )
    res['effective_gain_at_waist'] = .18*res['waist_peak_over_s2_peak']
    print(json.dumps(res, indent=2))
    np.savez('glued_drive.npz', t=t, modal=modal, peak=m.peak, **{f'bare_{k}': v for k, v in rb.items()}, **{f'klein_{k}': v for k, v in rk.items()})
