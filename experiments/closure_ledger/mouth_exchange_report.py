#!/usr/bin/env python
"""Evaluate frozen gates and draw the first-transit exchange diagnostics."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('experiments/closure_ledger/runs/20261017_mouth_exchange')

def key(lam, dt=.002, cutoff=32):
    return f'L{cutoff}_dt{dt}_a{lam}'

def evaluate(result):
    report = {}
    for name, case in result['cases'].items():
        runs = case['moving']
        lam = [1., .5, .25, .125]
        r = np.array([runs[key(a)]['normalized_final'] for a in lam])
        c = (4*r[1:]-r[:-1])/3
        reference = np.asarray(case['reference']['32']['pair_mixed_impulse'])
        ref = reference[-1]
        budget = max(.01*np.max(abs(c[-1,:,1])), 1e-8)
        change = max(np.max(abs(np.asarray(runs[key(a, dt, cutoff)]['normalized_final'])[:,1]-r[i,:,1]))
                     for i,a in [(2,.25),(3,.125)] for dt,cutoff in [(.001,32),(.002,40)])
        field = np.array([runs[key(a)]['mixed_field'][-1][1] for a in lam])
        composite = np.array([runs[key(a)]['composite_change'][-1] for a in lam])
        field_noise = max(abs(runs[key(a,dt,cutoff)]['mixed_field'][-1][1]-runs[key(a)]['mixed_field'][-1][1])
                          for a in [.25,.125] for dt,cutoff in [(.001,32),(.002,40)])
        composite_noise = max(abs(runs[key(a,dt,cutoff)]['composite_change'][-1]-runs[key(a)]['composite_change'][-1])
                          for a in [.25,.125] for dt,cutoff in [(.001,32),(.002,40)])
        # Absolute uncertainty also includes every-run normalized Noether residual times its scale.
        ledger_noise = max(run['stats']['max_angular_relative']*max(sum(np.sum(v) for v in run['records'][0]['energy'].values()),
                           np.linalg.norm(np.array(run['records'][0]['field'])+np.array(run['records'][0]['angular']).sum(axis=0)),1e-12)
                           for triple in runs.values() for run in triple['runs'])
        field_noise = max(field_noise, 3*ledger_noise)
        ratios = field[:-1]/field[1:]
        comp_ratios = composite[:-1]/composite[1:]
        maxe = max(run['stats']['max_energy_relative'] for triple in runs.values() for run in triple['runs'])
        maxj = max(run['stats']['max_angular_relative'] for triple in runs.values() for run in triple['runs'])
        ref_change = np.max(abs(np.asarray(case['reference']['40']['pair_mixed_impulse'])[-1,:,1]-ref[:,1]))
        resolved_field = min(abs(field[-2:])) > field_noise
        resolved_composite = min(abs(composite[-2:])) > composite_noise
        gates = dict(
            complete=not any(t['incomplete'] for t in runs.values()),
            extrapolation=bool(np.max(abs(c[-1,:,1]-c[-2,:,1])) <= budget),
            reciprocal=bool(abs(c[-1,:,1].sum()) <= budget),
            independent_reference=bool(np.max(abs(c[-1,:,1]-ref[:,1])) <= budget),
            moving_refinement=bool(change <= budget),
            reference_refinement=bool(ref_change <= budget),
            reference_ledger=bool(max(x['max_mixed_ledger'] for x in case['reference'].values()) < 1e-8),
            field_quartic=bool(np.all(abs(ratios[-2:]/16-1) <= .2)) if resolved_field else None,
            composite_quartic=bool(np.all(abs(comp_ratios[-2:]/16-1) <= .2)) if resolved_composite else None,
            energy=bool(maxe < 1e-4), angular_ledger=bool(maxj < 1e-8))
        report[name] = dict(gates=gates, budget=budget, normalized_final=r.tolist(),
            richardson=c.tolist(), reference=ref.tolist(), field_mixed=field.tolist(),
            field_ratios=ratios.tolist(), field_uncertainty=field_noise,
            composite_ratios=comp_ratios.tolist(), composite_uncertainty=composite_noise,
            max_energy_relative=maxe, max_angular_relative=maxj,
            max_moving_refinement=change, reference_cutoff_change=float(ref_change),
            reference_source_half=(reference[1]-reference[0]).tolist(),
            reference_target_full=(reference[4]-reference[2]).tolist())
    return report


def main():
    result = json.loads((ROOT/'results.json').read_text())
    report = evaluate(result)
    (ROOT/'gates.json').write_text(json.dumps(report, indent=2)+'\n')
    # Compact archive retains all recorded states, energies, momenta and events.
    (ROOT/'results.json').write_text(json.dumps(result, separators=(',', ':'))+'\n')
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for col,(name,case) in enumerate(result['cases'].items()):
        r = report[name]
        ax = axes[0,col]
        x = np.array([1., .5, .25, .125])**2
        normalized = np.array(r['normalized_final'])
        for h,color in [(0,'#166c9e'),(1,'#cb6633')]:
            ax.plot(x, normalized[:,h,1]*1e4, 'o-', color=color, label=f'Pair {h+1}')
            ax.axhline(r['reference'][h][1]*1e4, color=color, ls='--', lw=1)
        ax.set(title=name.title(), xlabel=r'$\lambda^2$', ylabel=r'Mixed $L_y/\lambda^2$ ($10^{-4}$)')
        ax.legend(frameon=False)
        ax = axes[1,col]
        ref = case['reference']['32']
        impulse = np.asarray(ref['pair_mixed_impulse'])
        ax.plot(ref['t'],impulse[:,0,1]*1e4,'o-',color='#166c9e',label='Pair 1')
        ax.plot(ref['t'],impulse[:,1,1]*1e4,'o-',color='#cb6633',label='Pair 2')
        ax.axvspan(np.pi-.6,np.pi+.6,color='gray',alpha=.15,label='Target focus window')
        ax.set(xlabel='Time (unit sphere)',ylabel=r'Fixed-axis mixed $L_y$ ($10^{-4}$)')
        ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Leading classical exchange: direct bilinear torque vs moving-mouth limit\nDashed lines: independent fixed-axis response; event samples joined for guidance')
    dest = Path('docs/images/mouth_exchange'); dest.mkdir(parents=True,exist_ok=True)
    fig.savefig(dest/'exchange.svg')
    fig.savefig('/tmp/mouth_exchange.png',dpi=130)
    files = sorted(ROOT.glob('*.json'))
    (ROOT/'sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.name!='sha256.json'},indent=2)+'\n')
    print(json.dumps(report, indent=2))

if __name__ == '__main__':
    main()
