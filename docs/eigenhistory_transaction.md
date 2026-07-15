# The eigenhistory transaction (PR #218)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. #217 solved the *driven* loop and
> found G_eff diverging exactly at the completed transaction. That
> divergence is a message: a completed transaction is not a response to
> an external source. This PR formulates it as what it is — a
> **homogeneous, globally constrained eigenhistory** whose amplitude is
> fixed by energy and state closure — and establishes the theorem: **a
> finite, source-inclusive, energy-conserving, self-consistent wormhole
> transaction exists.** The companion probe machine-checks every claim
> (~2 s).

## 0. From driven loop to eigenhistory

#217: F = F₀ + ΛF, G_eff = F₀/(1−Λ), pole at Λ → 1. At the pole the
driven description fails — because the completed transaction has no F₀.
It is a solution of the *homogeneous* equation

```
F = Λ_tot(ω, I) · F ,        Λ_tot(ω*, I*) = 1 ,
```

a self-sustaining history threading source → antipodal mouth → throat →
past mouth → source, forever — an **eigenhistory**. Existence is a
null-space statement: det(I − M_tot) = 0 for the full transfer system.
And the amplitude I* = |F|², which linearity cannot fix, is fixed by
the two closure constraints.

## 1. The source joins the loop

The source is no longer external: it sits *in* the loop as a scatterer
with internal state — an energy-conserving **reactive** element
(|s| = 1 exactly, zero net absorbed power) whose phase is pulled by the
field it sits in:

```
s(I) = e^{iφ_s(I)} ,      φ_s(I) = β·I / (1 + I/I_sat)
```

— the classical analog of the internal-state (level) shift of the
source during the transaction, with saturation as its finite capacity.
The full loop eigenvalue is Λ_tot(ω, I) = t_net(ω)·e^{iφ_s(I)} (time
closure absorbing the geometric phases per #217).

## 2. Energy closure: |Λ_tot| = 1

Every element must be lossless at the working point:

- exterior legs and mouth offset: unit modulus (geometry);
- the source: reactive by construction (|s| = 1 to 10⁻¹⁶);
- **the throat: lossless *exactly on* its interior resonance** — an
  algebraic unitarity identity for identical mouths: at resonance the
  loop factor r²e^{2iωτ} is real positive |r|², so
  |t_net| = |t|²/(1−|r|²) = T/T = **1 exactly**.

Machine-checked at two tiers: the exact tier (unitarity-projected
ports) gives |t_net| − 1 = 10⁻¹⁶; the physical Tangherlini ports give
deficit 2.0×10⁻⁴ = 10× the solver's own port flux error
(finesse-amplified) — the physical throat realizes the identity to
solver precision. Off resonance the loop is strictly passive
(|Λ| < 0.9), so the eigenhistory can live *only* on the resonance comb
— #217's fixed-point result, now forced by energy closure alone.

## 3. State closure: arg Λ_tot = 0 fixes the amplitude

At the resonance the throat still carries its scattering phase
χ = arg t_net(ω_res) = 2·arg t (measured 2.8631). The source must
supply the cancellation:

```
φ_s(I*) = −χ   (mod 2π)
```

**Existence by the intermediate value theorem**: φ_s sweeps
[0, β·I_sat) = [0, 12) > 2π continuously from zero, so the target is
always reachable — a root I* exists, finite (below saturation) and
nonzero (χ ≠ 0 generically). Solved: **I* = 1.5944**, phase residual
8×10⁻¹⁴, giving **Λ_tot − 1 = 3×10⁻¹⁶**.

**The null space opens exactly there**: the smallest singular value of
(I − M_tot) is 10⁻¹⁶ at (ω*, I*) versus 0.14 detuned in intensity and
0.9 detuned in frequency. The homogeneous solution exists at exactly
one point per branch — and

**the branches are discrete**: φ_s(I) = −χ + 2πm gives an amplitude
spectrum I*₀ = 1.594, I*₁ = 16.90, … — discrete transaction amplitudes
from a classical closure condition (branch structure, *not*
quantization — scope §6).

## 4. The theorem, dynamically verified

- **energy conserved pass by pass, exactly**: the loop map
  z ← Λ_tot(|z|²)·z preserves |z| to 10⁻¹⁶ per pass (|Λ| = 1);
- **persistence**: from z* = √I*, 10⁴ passes leave the amplitude
  drifted by 10⁻¹⁶ and the phase by 10⁻¹¹ — the eigenhistory is a
  genuine periodic self-consistent history;
- **source-inclusive**: the null eigenvector carries the source's
  internal-state shift φ_s(I*) = 3.42 and the mouths' resonantly
  enhanced interior amplitude (finite) — source state, mouth state, and
  field are parts of *one* closed history;
- **perturbed amplitudes dephase** at exactly dφ_s/dI (measured =
  predicted to 10⁻¹²): state closure *selects* I*;
- **the physical tier** decays at exactly its solver deficit
  (|t_net|^N matched to 10⁻⁶ over 200 passes) — the residual damping is
  numerical truncation, not physics.

## 5. The driven ↔ homogeneous correspondence

With the source frozen at I*, #217's driven response |1/(1−Λ_tot)|
exceeds 10¹⁵ at the eigenhistory and is O(10) detuned: **the driven
pole is the eigenhistory** — the marginal Novikov point (Λ = 1, flagged
passive-but-unpopulated in #217) is now populated by an explicit
history. Weak driving shadows it: the nonlinear driven steady state
converges to z* as F₀ → 0.

The Wheeler–Feynman/Cramer arc closes conceptually: #213 derived the
propagator from complete histories; #216 gave the advanced half a
causal mechanism; #217 resummed the loop; #218 exhibits the completed
transaction itself — not a limit, not a pole of a response function,
but a finite self-consistent object with its own fixed amplitude.

## 6. Honest scope

- **ℏ is not derived.** The amplitude scale comes from the source's
  nonlinear parameters (β, I_sat — the classical saturation analog);
  the discreteness of the branch spectrum is closure structure.
  Connecting the eigenhistory scale to the quantum of action is the
  successor question.
- **Monochromatic skeleton**: carrier closure only; a localized packet
  eigenhistory needs the group condition too (#217's Wigner-corrected
  closure) — named successor.
- The losslessness identity is exact in the unitary model class;
  Tangherlini realizes it to solver precision.
- Single-channel reactive source: no radiation into other modes,
  no radiative corrections to the source state.
- **Marginal stability**: the eigenhistory neither grows nor decays;
  selection against neighbors is by dephasing, not attraction. A weakly
  dissipative registration mechanism (the #209 opens) would make it
  attracting.
- Frozen geometry, classical zonal scalar, MTY network history posited.

## 7. What would falsify this

- A phase target outside the source's range — no root, no existence.
  (Checked: range 12 > 2π covers every χ.)
- A nonzero smallest singular value at the solved point — no
  homogeneous solution. (Checked: 10⁻¹⁶ vs 0.14/0.9 detuned.)
- Energy drift of the eigenhistory. (Checked: 10⁻¹⁶ per pass, 10⁻¹⁶
  over 10⁴ passes.)
- A continuous amplitude family — closure would fix nothing. (Checked:
  discrete branches, neighbors dephase at the predicted rate.)

## 8. Companion probe

`experiments/closure_ledger/eigenhistory_transaction_probe.py` (T1–T8,
~2 s): the losslessness identity at both tiers; the IVT existence
solve; the null-space scan; the 10⁴-pass persistence; the branch
spectrum and dephasing rates; the driven-pole correspondence.

**Verdict:**
`A_FINITE_SOURCE_INCLUSIVE_ENERGY_CONSERVING_SELF_CONSISTENT_WORMHOLE_TRANSACTION_EXISTS_THE_EIGENHISTORY_AMPLITUDE_FIXED_BY_ENERGY_AND_STATE_CLOSURE`
