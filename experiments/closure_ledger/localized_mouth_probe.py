"""Frozen localized-quartet constraints with retained successes and failures."""
import argparse
import gzip
import json
import hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from geometrodynamics.waves import localized_mouth as m
from .evidence_archive import read_bytes

GATES=('action','momentum','hamiltonian','seam','physical','localization','controls','evidence')


def peak(x):
    x=np.asarray(x,dtype=float)
    if not x.size or not np.isfinite(x).all():raise ValueError('empty/nonfinite evidence')
    return float(np.max(abs(x)))


def finite(x):
    if isinstance(x,dict):return all(finite(v) for v in x.values())
    if isinstance(x,(list,tuple)):return all(finite(v) for v in x)
    if isinstance(x,(float,int)):return bool(np.isfinite(x))
    return True


def run():
    out=dict(prereg=m.PREREG,gate_schema=list(GATES),symbolic=list(m.symbolic()),profiles={},
             solutions=[],physical=[],diagnostics=[],seams=[],failures=[])
    for L in m.LENGTHS:
        theta,a=m.profile(L)
        out['profiles'][str(L)]=dict(theta=m.serialize(theta),tensor=m.serialize(a))
        for eta in m.AMPLITUDES:
            for n,tol in m.SCHEDULE:
                try:
                    r=m.solve(L,eta,n,tol);out['solutions'].append(r)
                    print(f'L={L} eta={eta} n={n}: status={r["status"]}, nodes={r["nodes"]}',flush=True)
                    if n==513:
                        d=m.Data(r,out['profiles'][str(L)])
                        out['diagnostics'].append(dict(L=L,eta=eta,data=d.diagnostics()))
                        out['seams'].append(dict(L=L,eta=eta,data=m.seam(d)))
                except Exception as e:
                    out['failures'].append(dict(L=L,eta=eta,n=n,error=str(e)))
    r=next((r for r in out['solutions'] if (r['L'],r['eta'],r['n_initial'])==(5.5,.3,513)),None)
    if r is not None:
        d=m.Data(r,out['profiles']['5.5'])
        for h in (1e-3,5e-4,2.5e-4):
            cases=[m.physical_constraints(d,(frac*5.5,t,.37),h) for frac in (.07,.23,.51,.79,.93) for t in (.43,.91,1.47,2.13)]
            out['physical'].append(dict(h=h,cases=cases))
    out['reversal']=m.solve(4.5,-.3,257,1e-8)
    s=np.linspace(-4.5,4.5,4001);tp=m.theta_prime(abs(s),4.5)
    evenp=.3*m.sech(4.5)**3*np.cos(np.pi*s/(2*4.5))
    out['compatibility']=dict(s=s.tolist(),even_source=(-m.ALPHA*evenp*tp).tolist(),
        odd_source=(-m.ALPHA*m.momentum(s,4.5,.3)*tp).tolist())
    return out


def score(out,derivative_order=2):
    gates=dict.fromkeys(GATES,False);metrics={}
    try:
        integrity=finite(out) and out['prereg']==m.PREREG and out['gate_schema']==list(GATES) and not out['failures']
        gates['action']=out['symbolic']==list(m.symbolic()) and all(v=='0' for v in out['symbolic']) and m.F>0
        records=out['solutions'];expected=[(L,e,n,t) for L in m.LENGTHS for e in m.AMPLITUDES for n,t in m.SCHEDULE]
        keys=[(r['L'],r['eta'],r['n_initial'],r['tolerance']) for r in records]
        integrity &= keys==expected
        lookup={(r['L'],r['eta'],r['n_initial']):r for r in records}
        integrity &= set(out['profiles'])=={str(L) for L in m.LENGTHS}
        # Profile data must encode the registered initial data, not a chosen source.
        profile_error=0.
        for L in m.LENGTHS:
            theta,a=m.profile(L);stored=out['profiles'][str(L)];s=np.linspace(0,L,1001)
            for name,actual in [('theta',theta),('tensor',a)]:
                q=m.restore(stored[name]);profile_error=max(profile_error,peak(q(s)-actual(s)))
            integrity &= profile_error<1e-12
        mres=[];corrections=[];currents=[];Hres=[];boundaries=[];minimum=[];differences=[]
        statuses=[];rows=[]
        for L in m.LENGTHS:
            for eta in m.AMPLITUDES:
                s=np.linspace(0,L,1001);ds=[]
                for n,tol in m.SCHEDULE:
                    r=lookup[L,eta,n];d=m.Data(r,out['profiles'][str(L)]);ds.append(d)
                    statuses.append(r['success'] and r['status']==0)
                    minimum.append(float(np.min(d.sol(s)[0])))
                    boundaries.append(peak(d.sol(np.array([0.,L]))[1]))
                fine=ds[-1];H,M=fine.residuals(s)
                Hres.append(peak(H));mres.append(peak(M))
                difference=peak((ds[1].sol(s)[0]-ds[2].sol(s)[0])/ds[2].sol(s)[0]);differences.append(difference)
                if eta:
                    corrections.append(peak(eta*fine.tensor(s)-eta*m.sech(L)**3/8))
                    currents.append(peak(m.ALPHA*m.momentum(s,L,eta)*m.theta_prime(s,L)))
                rows.append(dict(L=L,eta=eta,H_residual=peak(H),M_residual=peak(M),
                                 minimum_psi=float(np.min(fine.sol(s)[0])),relative_refinement=difference))
        comp=out['compatibility'];s=np.asarray(comp['s']);tp=m.theta_prime(abs(s),4.5)
        even=-m.ALPHA*.3*m.sech(4.5)**3*np.cos(np.pi*s/9)*tp
        odd=-m.ALPHA*m.momentum(s,4.5,.3)*tp
        integrity &= peak(np.asarray(comp['even_source'])-even)<1e-13 and peak(np.asarray(comp['odd_source'])-odd)<1e-13
        even_integral=float(simpson(even,x=s));odd_integral=float(simpson(odd,x=s))
        gates['momentum']=max(mres)<1e-9 and min(corrections)>0 and min(currents)>0 and abs(even_integral)>1e-7 and peak(odd)>1e-7 and abs(odd_integral)<1e-12
        gates['hamiltonian']=all(statuses) and min(minimum)>0 and max(boundaries)<1e-9 and max(Hres)<1e-7 and max(differences)<1e-6
        metrics.update(cases=rows,momentum_max=max(mres),H_offgrid_max=max(Hres),boundary_max=max(boundaries),
                       relative_refinement_max=max(differences),even_source_integral=even_integral,
                       odd_source_integral=odd_integral,omitted_momentum_residual=peak(odd))

        seam_max=0.;missing_sign=0.
        integrity &= [(r['L'],r['eta']) for r in out['seams']]==[(L,e) for L in m.LENGTHS for e in m.AMPLITUDES]
        for row in out['seams']:
            L,e=row['L'],row['eta'];d=m.Data(lookup[L,e,513],out['profiles'][str(L)])
            computed=m.seam(d)
            integrity &= len(row['data'])==len(computed)
            for stored,actual in zip(row['data'],computed):
                for key in actual:integrity &= peak(np.asarray(stored[key])-actual[key])<1e-11
                seam_max=max(seam_max,*(peak(actual[key]) for key in actual if key not in ('t','missing_sign')))
                missing_sign=max(missing_sign,peak(actual['missing_sign']))
        gates['seam']=seam_max<1e-7 and missing_sign>.1
        metrics.update(seam_max=seam_max,missing_scalar_sign=missing_sign)

        # Recompute independent physical contractions from archived coordinate jets.
        d=m.Data(lookup[5.5,.3,513],out['profiles']['5.5'])
        points=[(frac*5.5,t,.37) for frac in (.07,.23,.51,.79,.93) for t in (.43,.91,1.47,2.13)]
        Herr=[];Merr=[];wrong=[];null_values=[]
        steps=[1e-3,5e-4,2.5e-4] if derivative_order==2 else [.008,.004,.002]
        integrity &= derivative_order in (2,4) and [r['h'] for r in out['physical']]==steps
        for group in out['physical']:
            hs=[];ms=[]
            integrity &= [tuple(c['geometry']['point']) for c in group['cases']]==points
            for c in group['cases']:
                gdata=c['geometry'];x=gdata['point'];h=group['h']
                integrity &= gdata['h']==h
                g,K=d.metric_tensor(x);phi,Pi=d.fields(x)
                integrity &= peak(np.asarray(gdata['metric'])-g)<1e-12 and peak(np.asarray(gdata['extrinsic_curvature'])-K)<1e-12
                integrity &= peak(np.asarray(c['fields'])-phi)<1e-12 and peak(np.asarray(c['Pi_J'])-Pi)<1e-12
                if derivative_order==2:
                    gradient=np.column_stack([(d.fields(np.asarray(x)+np.eye(3)[i]*h)[0]-d.fields(np.asarray(x)-np.eye(3)[i]*h)[0])/(2*h) for i in range(3)])
                else:
                    offsets,weights=(-2,-1,1,2),(1/12,-2/3,2/3,-1/12)
                    gradient=np.column_stack([sum(w*d.fields(np.asarray(x)+np.eye(3)[i]*h*o)[0] for o,w in zip(offsets,weights))/h for i in range(3)])
                integrity &= peak(np.asarray(c['gradient'])-gradient)<1e-12
                kinetic=Pi@Pi;spatial=np.einsum('Ai,ij,Aj',gradient,np.linalg.inv(g),gradient);current=Pi@gradient
                mixed=np.linalg.solve(g,K);K2=np.trace(mixed@mixed);tr=np.trace(mixed)
                integrity &= abs(c['kinetic']-kinetic)<1e-12 and abs(c['spatial']-spatial)<1e-12
                integrity &= peak(np.asarray(c['current'])-current)<1e-12
                integrity &= abs(gdata['K2']-K2)<1e-12 and abs(gdata['trace_K']-tr)<1e-12
                geom=gdata['R']+tr*tr-K2
                div=np.asarray(gdata['momentum_terms']).sum(axis=0)
                integrity &= peak(div-np.asarray(gdata['momentum']))<1e-12
                H=m.F*geom-kinetic-spatial-3;M=m.F*div+current
                integrity &= abs(H-c['H'])<1e-12 and peak(M-c['M'])<1e-12
                Hscale=max(1.,m.F*(abs(gdata['R'])+tr*tr+K2)+kinetic+spatial+3)
                Mscale=max(1.,m.F*sum(np.linalg.norm(t) for t in gdata['momentum_terms'])+np.linalg.norm(current))
                hs.append(abs(H)/Hscale);ms.append(np.linalg.norm(M)/Mscale)
                wrong.append(abs(geom-kinetic-spatial-3)/Hscale)
                G=np.eye(4)/m.F+np.outer(phi,phi)/(6*m.F*m.F)
                # All +/- coordinate orthonormal null directions; no assumed NEC flag.
                for i in range(3):
                    for sign in (-1,1):
                        v=Pi/np.sqrt(m.F)+sign*gradient[:,i]/np.sqrt(m.F*g[i,i])
                        null_values.append(float(v@G@v))
            Herr.append(max(hs));Merr.append(max(ms))
        ratios=[float(v[-2]/v[-1]) if v[-1] else None for v in (Herr,Merr)]
        low,high=(2.5,5.5) if derivative_order==2 else (8,24)
        convergence=all(v[-2]<=1e-8 or (r is not None and low<=r<=high) for v,r in zip((Herr,Merr),ratios))
        gates['physical']=Herr[-1]<1e-5 and Merr[-1]<1e-5 and convergence and max(wrong)>10*max(Herr[-1],Merr[-1])
        metrics.update(physical_H=Herr,physical_M=Merr,physical_ratios=ratios,wrong_f_H=max(wrong),null_min=min(null_values))

        localization=True;diagnostics=[];bulk={e:[] for e in m.AMPLITUDES}
        integrity &= [(r['L'],r['eta']) for r in out['diagnostics']]==[(L,e) for L in m.LENGTHS for e in m.AMPLITUDES]
        for row in out['diagnostics']:
            L,e=row['L'],row['eta'];d=m.Data(lookup[L,e,513],out['profiles'][str(L)])
            actual=d.diagnostics()
            for key in actual:integrity &= peak(np.asarray(row['data'][key])-actual[key])<1e-11
            error=peak(actual['bulk_relative']);bulk[e].append(error)
            r=actual['radii'];ratio=r[-1]/r[0]
            if L==5.5:localization &= error<.1 and 0<ratio<.15 and r[-3]>r[-1] and r[-2]>r[-1]
            diagnostics.append(dict(L=L,eta=e,bulk_error=error,radii=r,neck_ratio=ratio,
                                    theta_plus=actual['theta_plus'][-1],theta_minus=actual['theta_minus'][-1]))
        localization &= all(values[0]>values[1]>values[2] for values in bulk.values())
        gates['localization']=localization;metrics['localization']=diagnostics
        rev=out['reversal'];integrity &= (rev['L'],rev['eta'],rev['n_initial'],rev['tolerance'])==(4.5,-.3,257,1e-8)
        r=m.Data(rev,out['profiles']['4.5']);f=m.Data(lookup[4.5,.3,257],out['profiles']['4.5'])
        s=np.linspace(0,4.5,1001);metricdiff=peak(r.sol(s)[0]-f.sol(s)[0])
        reversal=metricdiff<1e-8
        for frac,t,_ in [(v/5.5,t,p) for v,t,p in points]:
            x=(frac*4.5,t,.37)
            reversal &= peak(r.metric_tensor(x)[1]+f.metric_tensor(x)[1])<1e-8
            reversal &= peak(r.fields(x)[1]+f.fields(x)[1])<1e-8
        dr=r.diagnostics();df=f.diagnostics()
        reversal &= peak(np.asarray(dr['theta_plus'])+df['theta_minus'])<1e-8
        reversal &= peak(np.asarray(dr['theta_minus'])+df['theta_plus'])<1e-8
        gates['controls']=reversal and min(null_values)>=-1e-14
        metrics['reversal_scalar_difference']=metricdiff
        gates['evidence']=integrity
    except (KeyError,TypeError,ValueError,IndexError,ArithmeticError) as error:
        gates['evidence']=False;metrics['evidence_error']=str(error)
    gates={k:bool(v) for k,v in gates.items()}
    data=all(gates[k] for k in GATES if k!='localization')
    return dict(gates=gates,passed=sum(gates.values()),total=8,metrics=metrics,
        verdicts={'FOUR_SCALAR_HANDLE_CONSTRAINT_DATA':data,'LOCALIZED_BULK_MOUTH_INITIAL_DATA':data and gates['localization']},
        unestablished=['traversability','crossing_evolution','momentum_transfer_events','sector_selection','discrete_action','quantum_statistics'])


def write_result(output,data,result):
    raw=json.dumps(data,sort_keys=True,indent=2,allow_nan=False).encode()+b'\n'
    (output/'probe.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    (output/'verdict.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    (output/'probe.md').write_text('# Localized four-scalar handle experiment\n\nFreeze: `'+m.PREREG+'`.\n\n'+f"Frozen gates: {result['passed']}/8.\n\n"+'```json\n'+json.dumps(result,indent=2,allow_nan=False)+'\n```\n')
    print(json.dumps({'passed':result['passed'],'gates':result['gates'],'metrics':result['metrics'],'raw_sha256':hashlib.sha256(raw).hexdigest()},indent=2),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--rescore',type=Path)
    args=p.parse_args();output=args.output_dir;output.mkdir(parents=True,exist_ok=True)
    (output/'probe.md').write_text('# Incomplete run\n\nNo affirmative verdict.\n')
    try:
        if args.rescore:
            raw=read_bytes(args.rescore);data=json.loads(gzip.decompress(raw) if args.rescore.suffix=='.gz' else raw)
        else:data=run()
        result=score(data);write_result(output,data,result)
        return 0 if all(result['gates'].values()) else 1
    except Exception as e:
        (output/'probe.md').write_text(f'# Failed run\n\nNo affirmative verdict.\n\n{e}\n')
        (output/'verdict.json').write_text(json.dumps({'error':str(e),'verdicts':{'FOUR_SCALAR_HANDLE_CONSTRAINT_DATA':False,'LOCALIZED_BULK_MOUTH_INITIAL_DATA':False}})+'\n')
        raise


if __name__=='__main__':raise SystemExit(main())
