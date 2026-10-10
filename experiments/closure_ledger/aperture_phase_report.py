"""Render the frozen family's saved measurements; never launch propagation."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from experiments.closure_ledger import aperture_phase_probe as p


def coherence_review(result):
    """Post-hoc comparison of already archived quantities; no propagation."""
    ds=result['diagnostics'];cases=[]
    for n,c in p.schedule():
        if not n.endswith('_fine') or c.squash==1:continue
        phase=ds[n]['phase'];coherence=phase['weighted_free_coherence']
        retention=result['retention'][n]
        cases.append(dict(name=n,squash=c.squash,carrier=c.carrier,footprint=c.footprint,
                          proxy=phase['proxy'],retention=retention,free_coherence=coherence,
                          capture_to_coherence=retention/coherence,
                          coherence_to_simple_fresnel=coherence/(np.pi/(4*phase['proxy']))))
    fits={}
    for f in (4.8,7.2):
        subset=[c for c in cases if c['squash']==.8 and c['footprint']==f]
        slope,intercept=np.polyfit(np.log([c['proxy'] for c in subset]),
                                  np.log([c['free_coherence'] for c in subset]),1)
        fits[str(f)]=dict(free_coherence_beta=float(-slope),
                         capture_beta=result['fits'][str(f)]['fine']['beta'])
    ratios=[c['capture_to_coherence'] for c in cases]
    high=[c for c in cases if c['proxy']>=10]
    return dict(interpretation='Post-hoc; no new production or registered gate.',
                manifest_sha=p.digest(p.RUN/'manifest.json'),cases=cases,fits=fits,
                ratio_min=min(ratios),ratio_max=max(ratios),ratio_median=float(np.median(ratios)),
                high_phase_count=len(high),
                high_phase_free_coherence_below_point_one=sum(c['free_coherence']<.1 for c in high))


def render():
    r=json.loads((p.RUN/'result.json').read_text());ds=r['diagnostics']
    review=coherence_review(r)
    (p.ROOT/'docs/aperture_phase_coherence.json').write_text(json.dumps(review,indent=2,allow_nan=False)+'\n')
    coherence_rows=[f'| {c["footprint"]:g} | {c["squash"]:g} | {c["carrier"]} | {c["retention"]:.6f} | {c["free_coherence"]:.6f} | {c["capture_to_coherence"]:.6f} |'
                    for c in review['cases']]
    coherence_fits=[f'| {f} | {v["capture_beta"]:.6f} | {v["free_coherence_beta"]:.6f} |'
                    for f,v in review['fits'].items()]
    colors=plt.cm.viridis(np.linspace(.05,.95,6))
    fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
    for color,b in zip(colors,(.8,.9,.98,1.02,1.1,1.2)):
        names=[p.name(b,w,4.8,'fine') for w in (12,24,48)]
        axes[0].loglog([ds[n]['phase']['proxy'] for n in names],
                       [r['retention'][n] for n in names],'-o',color=color,label=f'b={b:g}')
    names=[p.name(.8,w,7.2,'fine') for w in (12,24,48)]
    axes[0].loglog([ds[n]['phase']['proxy'] for n in names],
                   [r['retention'][n] for n in names],'--s',color='tab:red',label='b=.8, f=7.2')
    axes[0].axhline(.1,color='black',ls=':',label='Retention threshold')
    axes[0].axvspan(10,45,color='grey',alpha=.08)
    axes[0].set(xlabel='Carrier phase proxy Phi (rad)',ylabel='Capture / matching round capture',
                title='Retention decreases across the controlled family')
    axes[0].legend(fontsize=8,ncol=2);axes[0].grid(alpha=.2,which='both')
    for f,color in [(4.8,'tab:blue'),(7.2,'tab:red')]:
        for b,style in [(1.,'-o'),(.8,'--s')]:
            values=[100*ds[p.name(b,w,f,'fine')]['capture'] for w in (12,24,48)]
            axes[1].plot([12,24,48],values,style,color=color,label=f'f={f:g}, b={b:g}')
    axes[1].set(xlabel='Carrier w',ylabel='Captured incident energy (%)',
                title='Absolute capture also depends on the transducer',xticks=[12,24,48])
    axes[1].legend(fontsize=8);axes[1].grid(alpha=.2)
    figure=p.ROOT/'docs/figures/aperture_phase.png';fig.savefig(figure,dpi=160);plt.close(fig)
    main=[];secondary=[]
    for n,c in p.schedule():
        if not n.endswith('_fine'):continue
        d=ds[n]
        row=f'| {c.squash:g} | {c.carrier} | {d["phase"]["proxy"]:.3f} | {d["phase"]["weighted_std"]:.3f} | {100*d["capture"]:.4f} | {r["retention"][n]:.4f} | {100*d["remaining"]:.3f} |'
        (main if c.footprint==4.8 else secondary).append(row)
    fits=[]
    for f,v in r['fits'].items():
        fits.append(f'| {f} | {v["fine"]["beta"]:.6f} | {v["coarse"]["beta"]:.6f} | {v["fine"]["max_log_residual"]:.6f} | {v["fine"]["pass_law"]} |')
    diagnostic_limits={'energy_error':1e-9,'final_energy_error':1e-9,'port_error':1e-8,
                       'reconstructed_energy_error':1e-8,'state_error':1e-8,'early_flux':1e-5,'cap_tail':1e-4}
    diagnostics=[f'| {key} | {max(d[key] for d in ds.values()):.6g} | {limit:g} |'
                 for key,limit in diagnostic_limits.items()]
    provenance=json.loads((p.RUN/'provenance.json').read_text())
    high=[(n,v) for n,v in r['retention'].items() if ds[n]['phase']['proxy']>=10]
    failed=[(n,v) for n,v in high if v<.1]
    base=p.name(.8,48,4.8,'fine');late=ds['extended']['capture']-ds[base]['capture']
    body=f'''# Controlled aperture phase family: measured result

The 57 trajectories were run once after publication of freeze
`{provenance['freeze']}` in [PR #324](https://github.com/davidmdrpi/geometrodynamics/pull/324).
The [preregistered](aperture_phase_prereg.md) numerical and phase-coverage gates
are **{r['numerical_validity']}** and **{r['phase_coverage']}**, respectively.

| Registered claim | Frozen result |
|---|---|
| Every designated high-phase case retains >=10% of round capture | **{r['retention_label']}** |
| Both specified fits obey the approximate inverse-phase law | **{r['inverse_phase_label']}** |

{len(failed)} of {len(high)} designated high-phase cases fall below 10% retention.
The smallest is {min(high,key=lambda x:x[1])[1]:.5f}. This rejects uniform 10%
retention in the declared family. It does not show zero transport, nor does it
exclude other port models or geometries. The inverse-phase result is a
three-carrier, finite-range test at b=.8, not an asymptotic scaling theorem.
Neither result establishes quantum mechanics or self-consistent feedback.
The post-hoc review comparison below shows that this inverse-phase trend
largely tracks free-field spectral dephasing already encoded in the frozen
phase diagnostic. The additional empirical result is how closely integrated
lead capture follows that free coherence under the specified port dynamics.

![Capture and retention](figures/aperture_phase.png)

## What was controlled

Primary a*w=4.8; the second coherent-port family has a*w=7.2. The carrier
values are 12, 24 and 48; gamma/w=2/3 and the compact source contains three
cycles at each carrier. This fixes coordinate aperture/wavelength, relative
pulse bandwidth and dimensionless damping strength. All captures use the
complete incident energy and the same stop time 1.75 in round-transit units.
The horizon is not enlarged after observing a poor capture. Changes in bulk
spectral density and the geometry remain part of the question, rather than
being normalized away. The numerical implementation removes only exact dark
and +/-m degeneracies of the #323 port representation.

These distributed L2-normalized coherent transducers remain assumptions.
They are not physical holes with a derived boundary-matching law. Gamma is
scaled as a controlled parameter, not determined by GR. The corrective-kick
budget is zero. Static field-plus-lead energy is accounted for throughout.
There is no metric pump, dynamical R3 background, moving mouth or reinjection.

## Primary family, a*w=4.8

Capture and remaining energy are percentages of complete incident energy.
Retention is relative to the matching round case. Exact phase spread uses the
registered round free-source modal weights, not measured receiver amplitudes.

| b | w | Proxy Phi | Exact weighted phase std | Capture (%) | Retention | Remaining bulk (%) |
|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(main)}

## Second aperture family, a*w=7.2

| b | w | Proxy Phi | Exact weighted phase std | Capture (%) | Retention | Remaining bulk (%) |
|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(secondary)}

## Preregistered inverse-phase fits

Fit log(retention)=intercept-beta*log(Phi), only b=.8, at both footprints.
The proxy spans 10.60–42.41 rad. The registered acceptance region is
beta in [.5,1.5] and maximum absolute log residual <=.25 for **both** fits.
No alternative b slice, peak-flux exponent or time window replaces this test.

| a*w | Fine beta | Coarse beta | Fine max log residual | Within declared law bounds |
|---:|---:|---:|---:|---|
{chr(10).join(fits)}

These deterministic fits have only three points each. A passing law is only
supported in this finite interval. The exact source-weighted phase standard
deviation passes the separately stated coverage gate, but a broad asymptotic
regime and extrapolation to throat-scale resolution remain unestablished.
The near-round deformation cases report behavior outside the fit range;
they are not added as convenient points to improve the fitted law.

## Review follow-up: free coherence versus integrated capture

The [independent review](https://github.com/davidmdrpi/geometrodynamics/pull/324#issuecomment-6099053228)
correctly identifies the leading dephasing mechanism. This comparison is
**post-hoc interpretation** of quantities already archived by the frozen run,
not a newly registered test or a new simulation. Both original labels remain.
The [machine-readable comparison](aperture_phase_coherence.json) contains all
21 nonround fine cases and binds the original manifest.

Let C_free=|sum rho_lm exp(i delta_phi_lm)|^2, using the registered normalized
round free-source energy weights. This depends only on the spectrum and
source/aperture weights, without evolving the coupled ports. Let R denote
the measured nonround/round integrated lead capture. Across all 21 cases,
R/C_free lies in **[{review['ratio_min']:.6f}, {review['ratio_max']:.6f}]**, with
median **{review['ratio_median']:.6f}**.

| a*w | b | w | Measured retention R | Free coherence C_free | R/C_free |
|---:|---:|---:|---:|---:|---:|
{chr(10).join(coherence_rows)}

On the same preselected b=.8 slices, the free-coherence fits already closely
predict the integrated-capture exponents:

| a*w | Registered capture beta | Post-hoc free-coherence beta |
|---:|---:|---:|
{chr(10).join(coherence_fits)}

All {review['high_phase_count']} designated high-phase cases have C_free<.1.
The one measured retention above .1 is b=.8,w=12,a*w=7.2: free coherence
about .0944 is multiplied by R/C_free about 1.125, lifting retention to .1062.
Thus the free diagnostic predicted suppression before port evolution in
principle. It was not used as a prospective capture predictor or ratio gate
in this protocol. The principal additional information from the trajectories
is the measured port/propagation correction R/C_free and its energy ledger,
not discovery of an otherwise unknown inverse-phase mechanism.

### What the Fresnel asymptotic does and does not establish

For a single high-l block, uniform m weights and the quadratic phase
approximation give the continuum factor

    A(Phi) = integral_0^1 exp(i Phi x^2) dx
           = Phi^(-1/2) integral_0^sqrt(Phi) exp(i u^2) du,
    |A(Phi)|^2 ~ pi/(4 Phi) as Phi -> infinity.

The sign-reversed phase has the same squared magnitude. The standard Fresnel
limits and asymptotic expansions are given in
[DLMF 7.5](https://dlmf.nist.gov/7.5) and
[DLMF 7.12(ii)](https://dlmf.nist.gov/7.12#ii).
There is therefore an analytic inverse-phase asymptotic for this idealized
free factor. The earlier caution about no established asymptotic exponent
applies to **integrated capture in this full packet/port family**, not to
the Fresnel integral itself.

The archived C_free uses the exact square-root frequencies and a weighted
sum over l as well as discrete m; Phi uses the carrier alone. Consequently
pi/(4 Phi) is not its exact finite-range normalization. At b=.8 the ratios
C_free/[pi/(4 Phi)] are 1.222, 1.148, 1.148 for a*w=4.8 and
1.275, 1.446, 1.447 for a*w=7.2 as w increases. In particular the reported
band on R/C_free cannot be copied unchanged onto R/[pi/(4 Phi)].
The latter reaches about 1.82 in the b=.8,w=24,a*w=7.2 case.

Whether R/C_free stays bounded above and away from zero at much larger phase
is not settled by these samples. A subsequent protocol could put a prospective
band on that ratio, explicitly excluding or handling near-zero C_free, and
check the packet-weighted analytic approximation separately. The current
observed [0.987,1.367] band is not retroactively made a success criterion.

### Portable re-scoring test

The review also found that the evidence test required exact dictionary
equality after recomputing np.polyfit. Least-squares/LAPACK roundoff made
Python 3.10 CI fail although its fit values differed only in the last digits;
the Python 3.12 job passed. The revised test admits absolute 1e-12, with zero
relative tolerance, **only for beta and maximum log residual**. Fit decisions,
all labels and other re-scored fields remain exact, as do file/source hashes.
Regression tests accept last-bit fit changes but reject material fit changes,
changed labels, changed fit decisions and changes to other recorded values.
This edits the post-freeze evidence test, not the frozen scorer, protocol,
sources, thresholds or simulation data. CPU-dispatch testing alone does not
cover differences across NumPy/LAPACK versions.

Local validation uses Python 3.12.14 with NumPy 2.5.3 and separately 2.2.6.
The latter reproduces failure of the former exact-equality assertion, with
maximum fit drift 2.67e-15; all 64 focused tests pass with the revised test
in both environments. The Python 3.10 CI job remains a separate check.

## Verification and energy accounting

| Diagnostic | Largest observed | Frozen limit |
|---|---:|---:|
{chr(10).join(diagnostics)}

Largest paired/refinement comparison is {max(r['comparisons'].values()):.6g}
in source-normalized two-port output L2, against .03. Separate time, mode and
extension comparisons are {r['comparisons']['time_refine']:.6g},
{r['comparisons']['mode_refine']:.6g} and {r['comparisons']['extended']:.6g}.
The extension adds {late:.6g} of incident energy to B capture after 1.75 in
the one designated case. Its later flux is reported, not folded into the
primary measurement or either verdict. Reflected A energy plus captured B
energy plus remaining bulk energy closes the ledger in every trajectory.

Independent reconstruction uses archived forces incoming-outgoing, rebuilding
all field states and energies without re-solving the coupled port evolution.
Source formulas allow absolute roundoff <=1e-14 while grid, compact support
and inactive port remain exact. Stored incoming samples are used in physics
ledgers after validation. The frozen thresholds and source files are intact.

The original 58 focused tests passed before publication. A complete authenticated replay with AVX2, FMA3
and AVX512F NumPy features disabled reproduces both labels, both exponents
and every archived capture measurement. All 57 NPZ archives were inspected
with pickle disabled: 342 finite float64 arrays, with the declared configuration
metadata, and every manifest hash verified. CI runs the focused tests and full
archive reconstruction before the repository-wide suite.

## Evidence and chronology

- Run start: `{provenance['started_utc']}`.
- Run finish: `{provenance['finished_utc']}`.
- Runtime: Python {provenance['python']}, NumPy {provenance['numpy']}, SciPy {provenance['scipy']}.
- [All 57 raw trajectories, result, manifest and timestamped provenance](../experiments/closure_ledger/runs/20261010_aperture_phase/).
- Manifest SHA-256: `{p.digest(p.RUN/'manifest.json')}`.

Per-case start and finish UTC are recorded directly; no NPZ ZIP timestamp is
used as a run clock. The manifest authenticates the original arrays before
reconstruction. Replay checks all source hashes, measurements and both
verdicts. It never launches new production trajectories.

```sh
OPENBLAS_NUM_THREADS=1 python -m experiments.closure_ledger.aperture_phase_probe --manifest-sha {p.digest(p.RUN/'manifest.json')}
OPENBLAS_NUM_THREADS=1 python -m pytest -q tests/test_aperture_phase.py tests/test_aperture_transfer.py tests/test_mty_packet.py
python -m experiments.closure_ledger.aperture_phase_report
```

## Consequence for the next physical test

This family measures how much the static result depends on phase accumulation
under declared source and port scaling. It supplies no reason to extrapolate
static suppression unchanged to an evolving triaxial R3 geometry: anisotropy
can change during a transit, and mode coupling and metric work must then be
included. A future R3 study needs a new freeze, independently sampled geometry
phases, an explicit metric/port-work ledger and a pre-stated capture failure
criterion. The present study does not claim to have carried out that test.
'''
    (p.ROOT/'docs/aperture_phase.md').write_text(body)


if __name__=='__main__':render()
