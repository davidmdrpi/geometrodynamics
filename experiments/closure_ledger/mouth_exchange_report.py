#!/usr/bin/env python
"""Evaluate frozen gates and draw the first-transit exchange diagnostics."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'mouth-exchange'
import matplotlib.pyplot as plt
from mouth_exchange_probe import archive_text

ROOT = Path('experiments/closure_ledger/runs/20261017_mouth_exchange')

def key(lam, dt=.002, cutoff=32):
    return f'L{cutoff}_dt{dt}_a{lam}'

def evaluate(result):
    report = {}
    for name, case in result['cases'].items():
        runs = case['moving']
        if any(t['incomplete'] for t in runs.values()):
            report[name] = dict(gates=dict(complete=False), status='INCOMPLETE_COLLISION_STOP')
            continue
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
        ratios = np.divide(field[:-1], field[1:], out=np.full(3, np.nan), where=field[1:] != 0)
        comp_ratios = np.divide(composite[:-1], composite[1:], out=np.full(3, np.nan), where=composite[1:] != 0)
        composite_noise = max(composite_noise, 32*np.finfo(float).eps)
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
            composite_ratios=[float(x) if np.isfinite(x) else None for x in comp_ratios], composite_uncertainty=composite_noise,
            max_energy_relative=maxe, max_angular_relative=maxj,
            max_moving_refinement=change, reference_cutoff_change=float(ref_change),
            reference_source_half=(reference[1]-reference[0]).tolist(),
            reference_target_full=(reference[4]-reference[2]).tolist())
    return report



def write_numerics(result, report):
    lines = ['# First-transit exchange: numerical results', '',
             'All quantities use the unit-sphere model units. These are mixed coefficients, not full recoil momenta.', '',
             '| Preparation | Direct C1y | Direct C2y | Finest extrapolated C1y | Finest extrapolated C2y |',
             '|---|---:|---:|---:|---:|']
    for name,r in report.items():
        if not r['gates']['complete']:
            lines.append(f'| {name} | incomplete | incomplete | incomplete | incomplete |')
            continue
        ref, c = np.array(r['reference']), np.array(r['richardson'])[-1]
        lines.append(f'| {name} | {ref[0,1]:.10g} | {ref[1,1]:.10g} | {c[0,1]:.10g} | {c[1,1]:.10g} |')
    lines += ['', '![Exchange coefficient and time history](images/mouth_exchange/exchange.svg)', '',
              '## Frozen gates', '', '| Gate | Directional | Symmetric |', '|---|---|---|']
    for gate in sorted(set().union(*(r['gates'] for r in report.values()))):
        values = [{True:'PASS',False:'FAIL',None:'UNRESOLVED'}[report[n]['gates'].get(gate)] for n in ['directional','symmetric']]
        lines.append(f'| {gate} | {values[0]} | {values[1]} |')
    lines += ['', 'The symmetric single-packet composite separation is zero to floating-point precision. Its quartic ratio is undefined, so that gate is UNRESOLVED rather than counted as a pass.', '',
              'The two smallest adjacent amplitude ratios of retained mixed field Jy (expected 16) are:', '']
    for name,r in report.items():
        if not r['gates']['complete']: continue
        ratios = r['field_ratios'][-2:]
        lines.append(f"- {name}: {ratios[0]:.6f}, {ratios[1]:.6f}; estimated absolute numerical uncertainty {r['field_uncertainty']:.3g}.")
    lines += ['', '## Where the exchange accumulates', '',
              'Signed direct mixed impulse of pair 1 in four non-overlapping intervals; pair 2 has opposite y impulse.', '',
              '| Preparation | Source half [0,.6] | Between [.6,pi-.6] | Target [pi-.6,pi+.6] | After [pi+.6,pi+1.2] |',
              '|---|---:|---:|---:|---:|']
    for name,case in result['cases'].items():
        v = np.array(case['reference']['32']['pair_mixed_impulse'])[:,0,1]
        windows = [v[1]-v[0],v[2]-v[1],v[4]-v[2],v[5]-v[4]]
        lines.append('| '+name+' | '+' | '.join(f'{x:.9g}' for x in windows)+' |')
    lines += ['', '**The target-focus window partly reverses an exchange that has already accumulated.**',
              'The passing reciprocity test does not establish a localized antipodal exchange event.',
              'A prospective localization test must measure signed and absolute torque inside and outside the focus window;',
              'the window-resolved result here gives no basis to describe the whole endpoint impulse as an antipodal kick.', '',
              '## Numerical checks', '']
    for name,r in report.items():
        if not r['gates']['complete']: continue
        lines.append(f"- {name}: maximum relative energy error {r['max_energy_relative']:.3g}; normalized angular ledger error {r['max_angular_relative']:.3g}; largest timestep/cutoff change in normalized mixed y {r['max_moving_refinement']:.3g}; comparison budget {r['budget']:.3g}.")
    lines += ['', 'The archive contains 48 moving trajectories (16 three-packet comparisons), plus four direct references.',
              'The last two amplitudes are independently refined in timestep and spectral cutoff.',
              'No admitted run crosses the .65-radian stopping boundary if the complete gate passes.',
              'The source solver also supplies a dense interpolated reference curve for the figure; event accounting uses exact requested times.', '',
              '[Equations, interpretation and reproduction](mouth_exchange_results.md).',
              '[Frozen protocol](mouth_exchange_protocol.md).',
              '[Machine-readable gates and uncertainties](../experiments/closure_ledger/runs/20261017_mouth_exchange/gates.json).', '']
    Path('docs/mouth_exchange_numerics.md').write_text('\n'.join(lines))


def main():
    result = json.loads((ROOT/'results.json').read_text())
    report = evaluate(result)
    (ROOT/'gates.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    write_numerics(result, report)
    # Compact archive retains all recorded states, energies, momenta and events.
    (ROOT/'results.json').write_text(archive_text(result))
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for col,(name,case) in enumerate(result['cases'].items()):
        r = report[name]
        if not r['gates']['complete']:
            for ax in axes[:,col]:
                ax.text(.5,.5, name+': stopped before horizon', ha='center')
            continue
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
        impulse = np.asarray(ref['curve_pair_mixed_impulse'])
        ax.plot(ref['curve_t'],impulse[:,0,1]*1e4,'-',color='#166c9e',label='Pair 1')
        ax.plot(ref['curve_t'],impulse[:,1,1]*1e4,'-',color='#cb6633',label='Pair 2')
        ax.axvspan(np.pi-.6,np.pi+.6,color='gray',alpha=.15,label='Target focus window')
        ax.set(xlabel='Time (unit sphere)',ylabel=r'Fixed-axis mixed $L_y$ ($10^{-4}$)')
        ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Leading classical exchange: direct bilinear torque vs moving-mouth limit\nDashed lines: independent fixed-axis response; lower panels: adaptive reference solution')
    dest = Path('docs/images/mouth_exchange'); dest.mkdir(parents=True,exist_ok=True)
    fig.savefig(dest/'exchange.svg', metadata={'Date': None})
    svg = dest/'exchange.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    files = sorted(ROOT.glob('*.json'))
    (ROOT/'sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.name!='sha256.json'},indent=2)+'\n')
    print(json.dumps(report, indent=2, allow_nan=False))

if __name__ == '__main__':
    main()
