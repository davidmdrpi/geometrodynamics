"""Render release evidence from the checked-in benchmarks; no trajectory art."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometrodynamics.waves.s3_maxwell_encounter import DipolePulse, radial_basis

DATA = ROOT/'experiments/closure_ledger/runs'


def read(path):
    return json.loads((DATA/path).read_text())


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4), constrained_layout=True)
    pulse = DipolePulse(64, .2)
    chi = np.linspace(.0001, np.pi-.0001, 650)
    h, deriv = radial_basis(pulse.n, chi)
    f = h*np.sin(chi)**2
    fp = deriv*np.sin(chi)
    for time, color, label in [(0., '#226b98', 'Prepared pulse, t=0'),
                               (np.pi, '#c15127', 'Antipodal focus, t=π')]:
        q, p = pulse.state(time)
        u, us, ut = q@f, q@fp, p@f
        density = .5*(ut**2+us**2+2*u*u/np.sin(chi)**2)
        axes[0, 0].plot(chi/np.pi, density, color=color, lw=2.2, label=label)
    axes[0, 0].set(title='A. Exact S³ Maxwell refocusing', xlabel='Polar distance / π',
                   ylabel='Angular-integrated energy density')
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[0, 0].text(.5, .48, 'Full-field refocus defect\n8.66 × 10⁻¹⁵', ha='center',
                    transform=axes[0, 0].transAxes)
    rows = read('20261028_maxwell_throat/scattering.json')
    for depth, color in [(0., '#226b98'), (.4, '#55834b'), (.8, '#c15127')]:
        sub = [r for r in rows if r['angle']==.3 and r['ell']==1
               and r['charge_fraction']==0 and r['lapse_depth']==depth]
        axes[0, 1].semilogy([r['omega'] for r in sub], [r['transmission'] for r in sub],
                           'o-', color=color, label=f'Neck lapse {1-depth:.1f}')
    axes[0, 1].set(title='B. Geometry determines throat filtering', xlabel='Frequency (R=c=1)',
                   ylabel='Finite-port transmission fraction')
    axes[0, 1].legend(frameon=False, fontsize=9)
    for name, label, color in [('pair_heavy_fine', 'Prepared waves', '#c15127'),
                                ('pair_no_wave', 'No-wave control', '#226b98')]:
        data = read('20261023_deformable_boundary/'+name+'.json')['rows']
        sub = [r for r in data if r['time']<=120 and len(r['centers'])==2]
        axes[1, 0].plot([r['time'] for r in sub],
                        [r['centers'][1]-r['centers'][0] for r in sub], color=color, lw=2, label=label)
    axes[1, 0].set(title='C. Action-derived attraction in a field analogue',
                   xlabel='Analogue time', ylabel='Defect separation')
    axes[1, 0].legend(frameon=False, fontsize=9)
    axes[1, 0].text(.04, .05, 'Prepared radiation ≈ 42% of rest energy\nAt BAM scale (0.6%), the pair keeps separating',
                    transform=axes[1, 0].transAxes, fontsize=9)
    rows = read('20261028_einstein_throat/histories.json')
    for fraction, color in [(0., '#226b98'), (.9, '#c15127')]:
        row = next(r for r in rows if r['angle']==.3 and r['charge_fraction']==fraction)
        history = row['history']
        axes[1, 1].plot(np.array(history['tau'])/row['endpoint_proper_time'],
                        np.array(history['radius'])/row['initial_neck'], color=color, lw=2,
                        label='Neutral' if fraction==0 else 'Charged, q/qmax=0.9')
        axes[1, 1].scatter([1.], [row['inner_radius']/row['initial_neck']], color=color, s=20)
    axes[1, 1].set(title='D. Einstein evolution supplies radius dynamics',
                   xlabel='Central proper time / (πM)', ylabel='Neck radius / initial radius')
    axes[1, 1].legend(frameon=False, fontsize=9)
    axes[1, 1].text(.04, .05, 'Neutral endpoint: singularity\nCharged endpoint: ideal inner horizon',
                    transform=axes[1, 1].transAxes, fontsize=9)
    fig.suptitle('BAM foundations · verified building blocks for the coupled model', fontsize=17)
    target = ROOT/'docs/images/releases/bam_foundations.png'
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
