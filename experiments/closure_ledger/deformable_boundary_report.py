"""Replay saved states, plot measured fields, and bind the PR #349 evidence."""
import hashlib
import json
from pathlib import Path

import numpy as np

from geometrodynamics.waves.alice_boundary import time_deck
from geometrodynamics.waves.deformable_mode_conversion import DeformableDefect,Parameters
from .deformable_boundary_probe import ROOT,OUT


def read(name):
    return json.loads((OUT/(name+'.json')).read_text())


def separation(row):
    roots=row['centers']
    return roots[-1]-roots[0] if len(roots)==2 else np.nan


def main():
    checks=[]
    for case in read('schedule'):
        name=case['name'];data=read(name);state=OUT/(name+'.npz')
        assert hashlib.sha256(state.read_bytes()).hexdigest()==data['state_sha256'],name
        with np.load(state) as a:
            model=DeformableDefect(Parameters(**data['spec']['parameters']))
            initial=model.energy(a['u_initial'],a['v_initial'])
            final=model.energy(a['u'],a['v'])
            residual=abs(final-data['rows'][-1]['energy'])
            assert residual<2e-12,name
            assert abs(initial-data['background_energy']-data['seed_energy'])<2e-12,name
            assert np.max(abs(model.momenta(a['u'],a['v'])-data['rows'][-1]['momenta']))<2e-12,name
            roots=model.centers(a['u'])
            np.testing.assert_allclose(roots,data['rows'][-1]['centers'],atol=1e-12)
            checks.append(dict(case=name,final_energy_replay_error=residual,
                               relative_recorded_energy_error=data['maximum_energy_error']/initial))
    q=read('nonlinear_time_quotient')
    state=OUT/'nonlinear_time_quotient.npz'
    assert hashlib.sha256(state.read_bytes()).hexdigest()==q['state_sha256']
    with np.load(state) as a:
        np.testing.assert_allclose(a['u_past'],time_deck(a['u_future']),atol=1e-12)
        np.testing.assert_allclose(a['v_past'],-time_deck(a['v_future']),atol=1e-12)
    comparison={}
    for name in ['pair_no_wave','pair_heavy','pair_heavy_fine','pair_heavy_time_fine','pair_light','pair_offset']:
        data=read(name);rows=data['rows'];at120=min(rows,key=lambda r:abs(r['time']-120))
        first=next((r['time'] for r in rows if len(r['centers'])!=2 or separation(r)<2),None)
        comparison[name]=dict(separation_at_120=separation(at120),first_separation_below_2=first,
                              seed_over_background=data['seed_energy']/data['background_energy'])
    single=read('single_heavy_finer')
    late=[r for r in single['rows'] if r['time']>=65]
    velocity=float(np.polyfit([r['time'] for r in late],[r['centers'][0] for r in late],1)[0])
    with np.load(OUT/'single_heavy_finer.npz') as a:
        dx=float(a['x'][1]-a['x'][0]);center=single['rows'][-1]['centers'][0]
        gradient=(np.roll(a['u'][0],-1)-np.roll(a['u'][0],1))/(2*dx)
        mask=abs(a['x']-center)<5
        local_momentum=float(-dx*np.sum(a['v'][0,mask]*gradient[mask]))
    translation=dict(fitted_center_velocity=velocity,localized_phi_momentum=local_momentum,
                     kink_mass_times_gamma_velocity=4/3*velocity/np.sqrt(1-velocity*velocity))
    replay=dict(state_replays=checks,quotient_verified=True,comparison=comparison,translation_check=translation)
    (OUT/'replay.json').write_text(json.dumps(replay,indent=2,allow_nan=False)+'\n')
    plot()
    sources=[Path(__file__),ROOT/'experiments/closure_ledger/deformable_boundary_probe.py',
             ROOT/'geometrodynamics/waves/alice_boundary.py',ROOT/'geometrodynamics/waves/deformable_mode_conversion.py',
             ROOT/'geometrodynamics/waves/pair_birth.py',ROOT/'tests/test_alice_boundary.py',
             ROOT/'tests/test_deformable_mode_conversion.py',ROOT/'docs/deformable_alice_boundary.md',
             ROOT/'docs/images/history_throat/deformable_boundary.png']
    files=sorted([p for p in OUT.iterdir() if p.name!='manifest.json'])+sources
    manifest=dict(files=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
                              sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files])
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(states_verified=len(checks)+1,manifest_files=len(files),comparisons=comparison),indent=2))


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,2,figsize=(12,12),constrained_layout=True)
    ax=axes[0,0]
    for name,label in [('single_heavy_finer','Higher-mass incident mode'),('single_light','Opposite mode'),('single_disconnected','Coupling off')]:
        data=read(name)
        ax.plot([r['time'] for r in data['rows']],[r['centers'][0] for r in data['rows']],label=label)
    ax.axhline(0,color='.7',lw=.7);ax.set(title='(a) Defect motion follows the evolved field',xlabel='Time',ylabel='Position')
    ax.legend(fontsize=8)
    ax=axes[0,1]
    for name,label in [('pair_no_wave','No wave'),('pair_heavy','Inward branch'),('pair_heavy_fine','Refined inward branch'),('pair_light','Opposite mode')]:
        rows=read(name)['rows'];ax.plot([r['time'] for r in rows],[separation(r) for r in rows],label=label)
    ax.set(title='(b) Initially separating pair',xlabel='Time',ylabel='Zero-crossing separation',ylim=(0,50))
    ax.legend(fontsize=8)
    ax=axes[1,0];rows=read('single_closed_return')['rows']
    ax.plot([r['time'] for r in rows],[r['momenta'][0] for r in rows],color='#9550a0')
    ax.axhline(0,color='.7',lw=.7);ax.set(title='(c) Closed return nearly cancels first recoil',xlabel='Time',ylabel='Defect-field momentum')
    ax=axes[1,1]
    gaps=read('directional_gap_screen')
    for section,label in [('circle_m1','Circle m=1'),('sphere_l1','Sphere l=1')]:
        rows=sorted([r for r in gaps if r['angle']==.3 and r['doppler_factor']==1 and r['cross_section']==section],
                    key=lambda r:r['neck_over_mouth'])
        ax.semilogy([r['neck_over_mouth'] for r in rows],[r['flash_energy_above_gap'] for r in rows],'o-',label=label)
    ax.set(title='(d) Product-waveguide screen, a=0.3',xlabel='Neck / mouth radius',ylabel='Birth energy above directional gap')
    ax.legend(fontsize=8)
    ax=axes[2,0]
    with np.load(OUT/'pair_heavy_fine.npz') as a:
        select=abs(a['snapshot_x'])<25
        mesh=ax.pcolormesh(a['snapshot_x'][select],a['snapshot_t'],a['snapshots'][:,0,select],shading='auto',cmap='coolwarm',vmin=-1,vmax=1)
    ax.set(title='(e) Measured deformation field, refined pair',xlabel='Position',ylabel='Time')
    fig.colorbar(mesh,ax=ax,label='φ',shrink=.8)
    ax=axes[2,1];rows=read('nonlinear_time_quotient')['rows']
    ax.plot([-r['time'] for r in rows[::-1]],[r['past_cell_centers'][0] for r in rows[::-1]],label='Independently evolved past')
    ax.plot([r['time'] for r in rows],[r['future_cell_centers'][0] for r in rows],label='Future')
    ax.set(title='(f) One consistent time-reversing quotient history',xlabel='Cover time',ylabel='Cell defect position')
    ax.legend(fontsize=8)
    fig.suptitle('Deformable wave-conversion benchmark and Alice boundary checks\nReference action and quotient; BAM metric and focal matching remain open',fontsize=14)
    path=ROOT/'docs/images/history_throat/deformable_boundary.png'
    fig.savefig(path,dpi=150);plt.close(fig)


if __name__=='__main__':main()
