"""Deflection, leading junction source, and near-null stress budget."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import scipy
from geometrodynamics.tangherlini.characteristic_packet import solve_pair
from geometrodynamics.tangherlini.spacelike_seam import deflected_seam,field_derivative_source_control


def sample(d,index=-1):
    keys=['V','x','y','H','common_radius_shift','jump_K_length','jump_K_angle',
          'S_length_leading','S_angle_leading','densitized_S_angle','densitized_null_generator',
          'scalar_normal_jump_leading','scalar_delta_source_leading']
    return {k:float(d[k][index]) for k in keys}


def run(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for N in (64,128,256,512):
        packet=solve_pair(N)
        d=deflected_seam(packet['V'],packet['field'][-1,:,0])
        control=field_derivative_source_control(packet)
        window=packet['V']>=2.5
        exact=d['scalar_normal_jump_leading'][-1];measured=control['scalar_normal_jump'][-1]
        rows.append(dict(N=N,endpoint=sample(d),scalar_jump_2d=float(measured),
            scalar_jump_relative_discrepancy=float(abs(measured/exact-1)),
            radius_match_max=float(np.max(np.abs(d['radius_match_residual']))),
            temporal_match_max=float(np.max(np.abs(d['temporal_match_residual']))),
            min_H_after_peak=float(np.min(d['H'][window])),
            codazzi_leading_max_after_peak=float(np.max(np.abs(d['codazzi_leading_residual'][window]))),
            gauss_leading_max_after_peak=float(np.max(np.abs(d['gauss_leading_residual'][window])))))
    # Scaling uses the existing linear trace, not a claimed nonlinear amplitude experiment.
    amplitudes=[]
    for amplitude in (.005,.01,.02):
        s=deflected_seam(packet['V'],packet['field'][-1,:,0]*(amplitude/.01))
        amplitudes.append(dict(A=amplitude,**sample(s),
            proper_normal_width_budgets=[dict(width=w,
                necessary_peak_angular_stress_lower_bound=float(abs(s['S_angle_leading'][-1])/(2*w)))
                for w in (.01,.1,1.)]))
    jerrors=[r['scalar_jump_relative_discrepancy'] for r in rows]
    gates=dict(induced_metric_match=max(r['radius_match_max'] for r in rows)<1e-12 and max(r['temporal_match_max'] for r in rows)<1e-12,
               spacelike_after_peak=all(r['min_H_after_peak']>0 for r in rows),
               scalar_source_independent_control=jerrors[-1]<1e-2 and all(b<a for a,b in zip(jerrors[:-1],jerrors[1:])),
               leading_conservation_identities=rows[-1]['codazzi_leading_max_after_peak']<1e-12 and rows[-1]['gauss_leading_max_after_peak']<1e-12,
               null_density_limit=all(abs(r['densitized_S_angle']-.5)<1e-12 and abs(r['densitized_null_generator']+24)<1e-12 for r in amplitudes))
    keys=['V','trace','x','y','H','common_radius_shift','jump_K_length','jump_K_angle',
          'S_length_leading','S_angle_leading','scalar_normal_jump_leading']
    # NaN deliberately denotes undefined unit-normal quantities on the null prefix.
    np.savetxt(output/'seam_N512.csv',np.column_stack([d[k] for k in keys]),
               delimiter=',',header=','.join(keys),comments='')
    root=Path(__file__).resolve().parents[2]
    paths=['geometrodynamics/tangherlini/spacelike_seam.py',
           'geometrodynamics/tangherlini/characteristic_packet.py',
           'experiments/throat_dynamics/spacelike_seam_probe.py',
           'tests/test_spacelike_seam.py','docs/spacelike_seam_protocol.md',
           'docs/spacelike_seam_junction.md',
           'experiments/throat_dynamics/runs/20261014_reciprocal_packet/result.json']
    result=dict(label='SPACELIKE_MATCH_REQUIRES_SURFACE_STRESS_AND_SCALAR_SOURCE' if all(gates.values()) else 'NUMERICALLY_UNRESOLVED',
                numerical_gates=gates,refinement=rows,amplitude_diagnostics=amplitudes,
                conclusion='induced metric can match, but source-free scalar matching fails and surface stress remains required',
                normalization='b=c=kappa5=1; spacelike unit normal future-directed; scalar constant on S3',
                orders=dict(embedding='O(A²)',K_length='O(1/A), remainder O(A)',K_angle='O(A), remainder O(A³)',
                            S_length='O(A), remainder O(A³)',S_angle='O(1/A), remainder O(A)',scalar_source='O(A²), remainder O(A⁴)'),
                uncomputed=['surface action','full subleading conservation','finite-thickness metric solution',
                            'nonlinear stability','global continuation through the signature-change locus'],
                sources={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths},
                artifacts={'seam_N512.csv':hashlib.sha256((output/'seam_N512.csv').read_bytes()).hexdigest()},
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    (output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('label','numerical_gates','amplitude_diagnostics')},indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=Path('experiments/throat_dynamics/runs/20261014_spacelike_seam'))
    run(parser.parse_args().output)
