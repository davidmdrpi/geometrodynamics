"""Build the antipodal-only routing controls and connected-mouth replay.

These are prescribed worldtube tests. They do not use the superseded scalar
port forces, integrate a recoil law, or rerun the archived fourteen cases.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

from geometrodynamics.waves.history_throat import AntipodalTransport, inertial_mouth
from geometrodynamics.viz.klein_handles import build_surface, default_handles

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO/'experiments/closure_ledger/runs/20261021_antipodal_transport'
TEMPLATE = REPO/'experiments/closure_ledger/antipodal_throat_viewer.html'
VIEWER = REPO/'docs/images/history_throat/replay.html'


def controls():
    def mouth(q, v=0., **kw):
        return inertial_mouth(q, v, horizon=20., radius=.24, **kw)
    north = np.array([0., 0., 1.])
    specs = [
        ('through', 'Antipodal mouth → throat', [mouth(0.), mouth(np.pi, chart_side=-1)], {},
         'Fixed antipodal mouths. Each surface focus enters B, crosses the same throat, and exits A.'),
        ('empty', 'Empty antipode → reflection', [mouth(0.), mouth(np.pi-.8, chart_side=-1)], {},
         'B is displaced from the antipodal focus. The front reflects at the empty antipode and later returns to A.'),
        ('moving_hit', 'Moving mouth at arrival', [mouth(0., .04), mouth(np.pi+.08*np.pi, -.08, chart_side=-1)], {},
         'Prescribed inertial B arrives at the first focus at t=π. Later captures are decided anew from arrival positions.'),
        ('moving_miss', 'Moving mouth misses arrival', [mouth(0., .04), mouth(np.pi, -.12, chart_side=-1)], {},
         'B starts at the antipode but has left its aperture by arrival. Its old position cannot capture the front.'),
        ('reversed', 'Imposed loop · channel arrows reset', [mouth(0.), mouth(np.pi, chart_side=-1, orientation=-1, anchor=16.)], dict(emission=1.),
         'The supplied clock map and channel-forward resets impose a closed null route. Two caustics give gain −1, so only the zero amplitude satisfies its source-free loop condition.'),
        ('reversed_precaustic', 'Imposed loop · pre-caustic capture', [mouth(0.), mouth(np.pi, chart_side=-1, orientation=-1, anchor=16.)], dict(emission=1.),
         'The same channel-forward loop, but open-mouth capture precedes the caustic. Its front gain is +1; arbitrary homogeneous amplitudes are consistent, without selecting a driven field solution.'),
        ('reversed_oriented', 'Reversed time · transported arrow', [mouth(0.), mouth(np.pi, chart_side=-1, orientation=-1, anchor=16.)], dict(emission=1.),
         'The same clock map with orientation carried through each attachment. The exit surface leg runs backward in exterior time when required; the four-leg loop no longer closes.'),
        ('neck_reflection', 'Closed outlet → throat reflection', [mouth(0.), mouth(np.pi, chart_side=-1, opened=False, end_condition='dirichlet')], dict(emission=1., channel='throat'),
         'An initial front inside A travels to a closed B outlet, reflects inside the throat, and returns to A. This is a separate boundary control.'),
    ]
    result = {}
    for key, title, mouths, settings, description in specs:
        options = dict(emission=0., origin=north, max_legs=16)
        options.update(settings)
        r = AntipodalTransport(mouths, time_policy='channel_forward' if key in ('reversed', 'reversed_precaustic')
                              else 'transport_orientation', capture_phase='before_caustic' if key == 'reversed_precaustic'
                              else 'after_caustic').trace(**options)
        r.update(title=title, description=description)
        result[key] = r
    return result


def mesh():
    # Reuse the existing Klein-type glued drawing, including shared rim indices.
    # Later display warps move these same attached vertices. Their drawn arc
    # length is never used in characteristic travel times.
    m = build_surface(a=.24, handles=default_handles(1), n_ring=24, n_sphere=900,
                      ds_target=.012)
    neck = np.isfinite(m.s)
    fraction = np.full(m.n_vertices, -1.)
    fraction[neck] = (m.neck.S-m.s[neck])/(2*m.neck.S)
    fraction[m.rims[0][0]] = 0.
    fraction[m.rims[0][1]] = 1.
    return dict(points=np.round(m.points, 7).tolist(), triangles=m.tris.tolist(),
                neck_fraction=np.round(fraction, 7).tolist(),
                tube_triangles=(m.tri_tube >= 0).astype(int).tolist(),
                rims=[[r.tolist() for r in pair] for pair in m.rims],
                tube_rings=m.tube_ids[0].tolist(), radius=m.a)


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    cases = controls()
    for name, route in cases.items():
        (ROOT/f'{name}.json').write_text(json.dumps(route, separators=(',', ':'), allow_nan=False)+'\n')
    counts = {name: {kind: sum(e['event'] == kind for e in r['edges'])
                    for kind in ['mouth_entry', 'mouth_exit', 'antipodal_reflection', 'throat_reflection']}
              for name, r in cases.items()}
    summary = dict(cases=counts, status={k:r['status'] for k,r in cases.items()},
                   loop_diagnostics={k:r['continuation'] for k,r in cases.items() if r['continuation']},
                   scope='Prescribed-motion antipodal-only characteristic transport; no recoil law or field-energy ledger solved.')
    (ROOT/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    payload = dict(cases=cases, mesh=mesh())
    VIEWER.write_text(TEMPLATE.read_text().replace('/*TRANSPORT_DATA*/{}', json.dumps(payload, separators=(',', ':'), allow_nan=False)))
    paths = list(ROOT.glob('*.json')) + [VIEWER, TEMPLATE, Path(__file__),
            Path('experiments/closure_ledger/check_antipodal_throat_viewer.cjs'),
            Path('geometrodynamics/waves/antipodal_throat_transport.py'),
            Path('geometrodynamics/waves/history_throat.py'),
            Path('tests/test_antipodal_throat_transport.py')]
    paths += [p for p in [REPO/'docs/antipodal_throat_transport.md',
                          VIEWER.with_name('antipodal_transport.png'), ROOT/'README.md'] if p.exists()]
    paths = [p if p.is_absolute() else REPO/p for p in paths]
    hashes = {p.resolve().relative_to(REPO).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(paths) if p.name != 'sha256.json'}
    (ROOT/'sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
