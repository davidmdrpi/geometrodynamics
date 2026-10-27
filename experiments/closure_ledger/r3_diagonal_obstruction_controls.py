"""Pre-production controls; historical LRS inputs are read from a pinned git commit."""
import argparse
import hashlib
import json
from multiprocessing import Pool
from pathlib import Path
import subprocess
import time
import numpy as np
import gmpy2 as gp
from gmpy2 import mpfr as mp
from geometrodynamics.waves import r3_diagonal_taylor as dt
from geometrodynamics.waves import r3_high_precision_shooting as sh
from geometrodynamics.waves import r3_diagonal_symmetry as ds
from geometrodynamics.waves import r3_breaking as b

ROOT=Path(__file__).resolve().parents[2]
REF='e7afdeb46e274799bd1c0ca749e85e52f70dc50c'
BASE='experiments/closure_ledger/runs/20261010_r3_ladder/'
PINS={'hp_3_7.json':'ea1f54af3181d67c8f6f98879f0faf758791b74e786e3441ef23cdf764b854c1',
      'scan_3_7.json':'a91756fa9a55768013e4b2149e88431e5fb24cf9999090bb1872a84688151b9a'}


def historical(name):
    raw=subprocess.check_output(['git','show',REF+':'+BASE+name],cwd=ROOT)
    if hashlib.sha256(raw).hexdigest()!=PINS[name]:raise ValueError('control archive fingerprint')
    return json.loads(raw)


def lrs(j):
    hp=historical('hp_3_7.json')['rows'][j];point=historical('scan_3_7.json')['points'][j]
    results=[];start=time.time()
    for cfg,tol in ((dt.PRIMARY,'1e-45'),(dt.CHECK,'1e-65')):
        def P(z):
            v=dt.clock_map(list(z)+[0,0],**cfg)['z'];return v[:4]
        result=sh.solve(P,point,cfg['bits'],tol,start=hp if not results else results[-1])
        if not result['ok']:raise ArithmeticError('LRS control did not converge')
        results.append(result)
    with gp.context(gp.get_context(),precision=336):
        err=abs(mp(results[1]['lam'])-mp(hp['lam']));difference=abs(mp(results[0]['lam'])-mp(results[1]['lam']))
    return dict(j=j,archive_error=str(err),formulation_difference=str(difference),primary=results[0],check=results[1],seconds=time.time()-start,ok=bool(err<mp('1e-32') and difference<mp('1e-40')))


def synthetic(epsilon='1e-24'):
    """Kick-drift symplectic map with exact lambda=-epsilon sin(12 phi).

    Coordinates (theta,I); fixed theta slice, g=(0,-1), q=1.
    The solved I=-epsilon sin(12 phi); lambda has 24 simple zeros if epsilon!=0.
    """
    errors=[];lams=[]
    with gp.context(gp.get_context(),precision=224):
        e=mp(epsilon)
        for j in range(192):
            phi=2*gp.const_pi()*j/192+mp('.013');s=gp.sin(12*phi);der=12*e*gp.cos(12*phi)
            def P(z):
                th,I=z;after=I+e*gp.sin(12*th);return [th+after,after]
            point=dict(nodes=[[str(phi),'0']],lam='0',g=[0,-1],t0=[1,0],c0=[str(phi),'0'],
                       J=[[float(der),1,0],[float(der),0,1],[1,0,0]])
            out=sh.solve(P,point,224,'1e-60');errors.append(abs(mp(out['lam'])+e*s));lams.append(float(out['lam']))
            if not out['ok']:raise ArithmeticError('synthetic control failed')
    return dict(epsilon=epsilon,max_lambda_error=str(max(errors)),sign_changes=b.sign_changes(lams),
                dominant_harmonic=int(np.argmax(abs(np.fft.rfft(lams))[1:])+1) if e else None,
                ok=bool(max(errors)<mp('1e-55') and (not e or b.sign_changes(lams)==24)))


def off_family():
    z=['1.0004','-.002','.05','.03','-.04','.1'];out=[]
    for cfg in (dt.PRIMARY,dt.CHECK):out.append(dt.clock_map(z,**cfg))
    with gp.context(gp.get_context(),precision=336):err=max(abs(a-c) for a,c in zip(out[0]['z'],out[1]['z']))
    reference=ds.clock_map(np.array(z,float),half=False,matrix=True)['z']
    matrix_error=float(max(abs(reference-np.array(out[0]['z'],float))))
    return dict(formulation_difference=str(err),matrix_error=matrix_error,
                constraints=[str(r['constraint']) for r in out],ok=bool(err<mp('1e-50') and matrix_error<1e-10))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    with Pool(3) as pool:rows=pool.map(lrs,[0,17,41])
    result=dict(reference_commit=REF,archive_pins=PINS,lrs=rows,synthetic=synthetic(),null=synthetic('0'),off_family=off_family())
    result['ok']=all(r['ok'] for r in rows+[result['synthetic'],result['null'],result['off_family']])
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='lrs'},indent=2))


if __name__=='__main__':main()
