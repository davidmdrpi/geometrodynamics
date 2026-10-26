"""Audit exact focus/contact ordering in pinned PR #345, without editing it.

Requires git object 25ce28cc8e33ce7541b7d4dfac034a5574f6dcff (fetch that PR's
branch first). Run as a separate process: its temporary package path selects
the two pinned upstream modules. Remaining dependencies must match this HEAD.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[2]
PIN='25ce28cc8e33ce7541b7d4dfac034a5574f6dcff'
OUT=ROOT/'experiments/closure_ledger/runs/20261025_quotient_flash'


def main():
    provenance={}
    def blob(name):
        return subprocess.check_output(['git','show',PIN+':geometrodynamics/waves/'+name],cwd=ROOT)
    for name in ['moving_mouth_pairs.py','antipodal_kick.py']:
        data=blob(name)
        if data!=(ROOT/'geometrodynamics/waves'/name).read_bytes():
            raise RuntimeError('Unpinned dependency differs: '+name)
        provenance[name]=hashlib.sha256(data).hexdigest()
    with tempfile.TemporaryDirectory() as temp:
        for name in ['pair_birth.py','antipodal_pair_birth.py']:
            data=blob(name);(Path(temp)/name).write_bytes(data)
            provenance[name]=hashlib.sha256(data).hexdigest()
        import geometrodynamics.waves
        geometrodynamics.waves.__path__.insert(0,temp)
        from geometrodynamics.waves.antipodal_pair_birth import AntipodalPairBirth
        from geometrodynamics.waves.pair_birth import Parameters
        model=AntipodalPairBirth(Parameters(cutoff=12,footprint=.2))
        rows=[]
        for k in [1,2]:
            uk=(np.pi-.6)/((2*k+1)*np.pi)
            for delta in [-1e-9,0.,1e-9]:
                target=uk+delta
                kick=brentq(lambda v:np.linalg.norm(model.velocity(model.initial(v)),axis=1)[0]-target,.02,.95,xtol=1e-14)
                for dt in [.02,.01,.005]:
                    initial=model.initial(kick);speed=float(np.linalg.norm(model.velocity(initial),axis=1)[0])
                    _,info=model.run(initial,(np.pi-.6)/speed+.05,cycles=1,dt=dt,output_dt=.5)
                    events=info['events'];contact=next(e for e in events if e['kind']=='contact')
                    captures=[e for e in events if e['kind']=='capture' and e['front'].startswith('flash')]
                    row=dict(k=k,target_speed=target,delta=delta,actual_speed=speed,dt=dt,
                        contact_time=contact['time'],contact_minus_focus=contact['time']-(2*k+1)*np.pi,
                        capture_times=[e['time'] for e in captures],
                        event_order=[e['kind'] for e in events if abs(e['time']-(2*k+1)*np.pi)<.01],
                        energy_drift=info['energy_drift'],remaining_fronts=len(info['final_state'].fronts))
                    rows.append(row);print(json.dumps(row),flush=True)
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/'pinned_endpoint_order.json').write_text(json.dumps(dict(commit=PIN,source_hashes=provenance,rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
