# The geon-statistics adjudication: is the BAM exchange −1 a theorem? (PR #196)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This document is **mathematics against a
> specific literature** (Friedman–Sorkin; Sorkin–Surya; Dowker–Sorkin;
> Hendriks; Giulini), not a simulation. A companion probe machine-checks the
> finite/algebraic steps; the argument lives here.

## 0. The verdict, up front

**The spin-statistics *correlation* is a theorem for BAM throats** — the
throat prime passes all three hypotheses of the Dowker–Sorkin theorem
(prime ✓, non-chiral ✓, abelian ✓), which is the strongest topological
proof point available to the program.

**The *sign* is framing-dependent, and the honest label of the arc changes:**

- In **bare metric GR** (quantizing Diff of the 3-manifold, no spinor
  structure), the BAM throat prime RP³ is **non-spinorial** (Hendriks: the
  2π rotation is *isotopic to the identity*), so the same Dowker–Sorkin
  theorem that establishes the correlation **selects Bose (+1)** for
  pair-created throats. In bare canonical GR without topology change,
  statistics is a θ-sector *choice* (Bose, Fermi, and a continuum of
  indefinite-statistics sectors all exist). Either way, the −1 is **not** a
  bare-GR theorem.
- In **Pin⁻-framed GR** — the state space BAM actually has, because the
  antipodal quotient forces the non-orientable RP² brane slice (#169), RP²
  admits Pin⁻ *only* (#170), and the throat's particle modes are charged
  pin spinors (#195) — the 2π rotation, though isotopic to the identity as
  a diffeomorphism, **lifts to −1 on the pin frame bundle** (the rotation
  loop is the generator of π₁(SO(3)) ≅ Z₂, and RP³ ≅ SO(3)). The throat is
  *spinorial with framing*, and the same SSC theorem then **selects Fermi
  (−1)** — consistent with the holonomy measured in #188.

So: **"Pauli from GR" must honestly be labeled "Pauli from GR + the Pin⁻
framing that BAM's own antipodal quotient forces."** The framing is not an
import bolted on to get the answer — it is forced by the construction
(§1, §5) — but it *is* a step beyond vacuum Diff quantization, and the
thesis text (#170/#171) currently overstates the bare-topology case (§3).

---

## 1. Lemma 1 — the throat prime is RP³ (fixed by the repo's own construction)

The #169 quotient acts on the 5D Tangherlini geometry by
`J: (χ,θ,φ) ↦ (π−χ, π−θ, φ+π)`, restricting on the brane (χ = π/2) to the
S² antipodal map. Near the throat the brane's spatial slice is the
Einstein–Rosen neighborhood S² × ℝ, and J acts on it by
`(x, y) ↦ (−x, −y)` (antipodal on the sphere, reflection through the
throat). This involution is *free* and *orientation-preserving*:
`det(dJ) = det(dA|_{TS²}) · (−1) = (−1)(−1) = +1` (the S² antipodal map
reverses the 2-sphere's orientation; the reflection of the ℝ factor
reverses once more). The quotient is therefore an **orientable** 3-manifold:
the twisted ℝ-line bundle over RP², which is `RP³ ∖ {point}` — the
**RP³ geon** of Friedman–Sorkin and Louko–Marolf, whose minimal surface is
the one-sided **RP² cross-cap**: exactly the repo's "non-orientable throat
slice inside an orientable bulk" (#167–#169).

Consequently the spatial manifold of the brane with n throats is

```
Σₙ  =  S³ # RP³ # ⋯ # RP³  =  #ⁿ RP³ ,
```

and the BAM throat prime is `P = RP³` — *the* canonical example of the
entire geon-statistics literature. This is a genuinely favorable accident:
no other prime has been analyzed as thoroughly.

## 2. Lemma 2 — the three Dowker–Sorkin hypotheses hold

**(a) Prime.** RP³ = S³/Z₂ is an elliptic space form, hence irreducible
(every embedded 2-sphere bounds a ball), hence prime.

**(b) Non-chiral (amphichiral).** A lens space L(p,q) admits an
orientation-reversing self-diffeomorphism iff `q² ≡ −1 (mod p)`. For
RP³ = L(2,1): `1 ≡ −1 (mod 2)` ✓. Explicitly: any reflection R ∈ O(4),
det R = −1, commutes with the deck map −I₄, hence descends to an
orientation-reversing diffeomorphism of RP³. **The throat is its own
mirror image** — so the mirror partner created in pair production (§6) is
diffeomorphic to the throat itself, and "identical-particle statistics" is
well-posed for the created pair.

**(c) Abelian.** π₁(RP³) = Z₂, abelian.

These are precisely the hypotheses under which Dowker–Sorkin prove their
spin-statistics theorem for topological geons ("non-chiral abelian geons
… satisfy a spin-statistics correlation if … described by a functional
integral over metrics on a four-manifold describing a topology-changing
process which creates a pair of geons"). **The BAM throat is in the
covered class.** The answers to the adjudication's three questions are
yes, yes, yes.

## 3. Lemma 3 — RP³ is NOT spinorial in bare Diff (a correction to #170/#171)

Hendriks' classification (with Friedman–Sorkin's physical reading, and
Giulini's tabulation): a closed oriented 3-manifold is **non-spinorial**
iff its prime decomposition consists entirely of **lens spaces and handles
(S¹×S²)**; in particular *every prime with cyclic fundamental group is
non-spinorial*. RP³ = L(2,1) is a lens space with π₁ = Z₂:

> **The 2π rotation of a BAM throat, parallel to its connecting sphere, is
> isotopic to the identity in Diff.** In bare metric GR the throat carries
> only integer angular momentum; the Friedman–Sorkin route to spin-½ from
> pure spatial topology is **not available** for RP³.

This corrects the thesis text at #170/#171, which asserts "the single
geon's 2π rotation acts as −I (spinorial; … Friedman–Sorkin's spin-½)."
As a statement about π₀Diff (bare large diffeomorphisms) that is **false**
for RP³ — for lens spaces the isotopy trivializing the rotation is
explicit. What *is* true is the **pin-lifted** statement of §5. The #188
probe in fact measured the lift, not the bare diffeomorphism class — its
result stands, but its interpretation sharpens.

## 4. Lemma 4 — the two-throat mapping class group and its sectors

For two RP³ geons (asymptotically flat setting; the transfer to the closed
S³ background only adds background isometries, which act trivially on the
sector structure — McCullough's generation theorem localizes the MCG to
the summands), Sorkin–Surya give the mapping class group explicitly:

```
G  =  ⟨ s₁, s₂, E  |  s₁² = s₂² = E² = 1,  E s₁ E = s₂ ⟩
   =  (Z₂ ∗ Z₂) ⋊ S₂ ,
```

where s₁, s₂ are the **slides** of each prime through the Z₂ loop of the
other (`s² = e` because π₁(RP³) = Z₂ — the Fouxe-Rabinovitch relations),
E is the exchange, internal diffeos contribute nothing (π₀Diff(RP³) is
trivial for orientation-preserving classes), and the rotation is trivial
by Lemma 3 (consistently, E² = 1 appears directly in the presentation).
Note the slide subgroup is the **free** product Z₂ ∗ Z₂ ≅ D_∞ — slides do
*not* commute — so G is infinite.

**Unitary sector structure** (the θ-sectors of canonical quantization):

- **Four one-dimensional sectors**: s₁ = s₂ = ±1 (forced equal, being
  E-conjugate), E = ±1. These are Bose and Fermi, each in a slide-trivial
  and a slide-twisted version. **The Fermi sector exists** — the #171
  claim of *consistency* survives adjudication.
- **A continuum of two-dimensional sectors** (s₁, s₂ two reflections at
  relative angle θ; E a reflection interchanging them): here exchange does
  not commute with slides — geons with **indefinite statistics**, the
  Sorkin–Surya anomalous sectors, realized concretely on BAM's own prime.

**Conclusion of the canonical analysis:** bare frozen-topology GR does
*not* select the −1. Statistics is a superselection choice among Bose,
Fermi, and indefinite — exactly Dowker–Sorkin's "no obstruction to
anomalous pairings without topology change."

## 5. Theorem — the sign under topology change, and the pin lift

**(i) The SOH selection.** In the sum-over-histories framework with
topology change, the weights come from **abelian (one-dimensional)
representations** of the mapping class group (Sorkin–Surya) — the
indefinite-statistics continuum of §4 is eliminated — and for geons
**pair-created** on a suitable four-manifold, Dowker–Sorkin's theorem
forces the **spin-statistics correlation**: exchange phase = 2π-rotation
phase. All its hypotheses hold for the BAM throat (Lemma 2).

**(ii) Bare GR evaluates the rotation phase to +1.** By Lemma 3 the
rotation is isotopically trivial: tensorial. SSC then gives **Bose**.
Pair-created BAM throats in bare metric GR are bosons. *If BAM were bare
geometrodynamics, the exchange arc would be refuted here.*

**(iii) BAM is not bare geometrodynamics: the Pin⁻ framing is forced.**
Three repo-level facts, none optional:
1. the antipodal quotient makes the brane slice RP² non-orientable
   (#169) — any spinor transport on it requires a *pin* structure;
2. RP² admits **Pin⁻ only** (w₂ + w₁² = 0; #170) — the structure is
   unique, not chosen;
3. the throat's particle modes *are* charged pin spinors — the #195 index
   mechanism (the k zero modes of the charge-k/2 monopole Dirac operator)
   is what makes the electron level natural at all.

So the physical state space carries the Pin⁻ frame bundle, and the
symmetry group acting on states is the group of diffeomorphisms **lifted
to the pin bundle**.

**(iv) The pin lift of the trivial rotation is −1.** RP³ ≅ SO(3), and the
2π rotation about an axis is the loop `t ↦ R(t)`, t ∈ [0, 2π] — the
generator of π₁(SO(3)) ≅ Z₂. Its lift to the double (spin/pin) cover SU(2)
is `t ↦ exp(−i t σ·n̂/2)`, ending at **−1**. Concretely: the isotopy that
trivializes the rotation in Diff(RP³) traces this noncontractible loop, so
the induced bundle map on the Pin⁻ frame bundle is the deck element −1.
On the pin-framed state space the throat **is spinorial**: the 2π rotation
acts as −I. (This is exactly the holonomy #188 measured; the measurement
was correct, and Lemma 3 shows it was measuring the lift, not the bare
mapping class.)

**(v) The sign.** Applying the same SSC theorem on the pin-framed state
space: exchange phase = pin-rotation phase = **−1. Fermi.** The Pauli
exclusion of #185/#187 follows.

**Status of (v):** the correlation is Dowker–Sorkin's theorem; the
evaluation of the rotation phase on the pin-framed states is the SU(2)
lift above (elementary and machine-checked in the companion probe); the
one step taken on trust is that the SOH framework extends to pin-framed
state spaces with the correlation intact — standard in the literature's
own treatment of geons with fermionic matter, but not re-proved here.

## 6. The #58 nucleation cobordism against the Dowker–Sorkin condition

The #58 channel nucleates a **throat–antithroat pair** (Hopf charges ±1,
threshold 2m_ec²). Against the theorem's requirements:

- **Mirror-pair identity.** Pair creation produces a geon and its mirror.
  By Lemma 2(b) the mirror of RP³ is RP³ — the created pair are
  *identical* geons, as the statistics question requires. ✓
- **Existence of the cobordism.** Ω₃^{SO} = 0: every closed orientable
  3-manifold bounds, so 4-manifolds interpolating S³ → S³ # RP³ # RP³
  exist. For the pin-framed version: RP³ is spin (orientable,
  parallelizable — it is a Lie group), Σ = RP³ # RP³ is spin, and
  Ω₃^{Spin} = 0, so the spin structure extends over a bounding 4-manifold.
  The structure needed for §5(iii–v) is available at the bordism level. ✓
- **The specific Dowker–Sorkin 4-manifold.** Their theorem is proved for
  the wave function defined by *their* pair-creation four-manifold (built
  over the ℝ³ background). The BAM background is S³, but connected sums
  are local: the creation region is a ball, and gluing their construction
  into a ball of S³ (rather than of ℝ³) is the standard transplant. What
  is *not* verified here: that the transplanted cobordism satisfies the
  finer causal/Morse-structure conditions (elementary cobordism, causal
  continuity) that the SOH weighting arguments use, on the closed
  background. **This is the honest open item** — a construction, not an
  obstruction: no topological invariant blocks it (all the bordism-level
  checks pass), but the explicit BAM 4-manifold has not been written down.

  **[CLOSED in PR #200.]** The explicit manifold is exhibited in
  `docs/pair_creation_cobordism_capstone.md`: two 2-handles attached to
  S³×I along an unlink with framings ±2 — boundary S³ → RP³ # RP³, spin
  (even framings), signature 0, only index-2 Morse points (causally
  continuous, avoiding the proven-discontinuous classes) — built
  directly over the S³ background, no ℝ³ transplant needed.

## 7. The adjudication table

| Question | Answer | Source |
|---|---|---|
| Is the throat prime? | **Yes** — RP³, irreducible elliptic | Lemma 2a |
| Is it non-chiral? | **Yes** — L(2,1): q² ≡ −1 (mod 2); explicit reflection descends | Lemma 2b |
| Is π₁ abelian? | **Yes** — Z₂ | Lemma 2c |
| Is the throat spinorial in bare Diff? | **No** — lens spaces are non-spinorial (Hendriks) | Lemma 3 |
| Does bare canonical GR fix the statistics? | **No** — Bose, Fermi, and indefinite sectors all exist ((Z₂∗Z₂)⋊S₂) | Lemma 4 |
| Does the SSC theorem apply to BAM throats? | **Yes** — all three hypotheses hold | Lemma 2, §5(i) |
| Sign in bare metric GR (SOH, pair-created) | **+1 (Bose)** | §5(ii) |
| Sign in Pin⁻-framed GR (BAM's state space) | **−1 (Fermi)** — the SU(2)/pin lift of the trivialized rotation | §5(iii–v) |
| Does the #58 channel satisfy the 4-manifold condition? | **Structurally yes** (mirror-identity ✓, bordism existence ✓, spin extension ✓); the explicit transplanted 4-manifold remains to be constructed | §6 |

**Outcome:** neither pure Outcome A nor pure Outcome B — something more
precise than either. BAM **owns the spin-statistics correlation as a
theorem** on its topology (the strongest available form of "Pauli from
GR"), and owns the **−1 conditionally**: conditional on the Pin⁻ framing,
which its own quotient forces, and on the SOH transplant of §6. The
anomalous sector is **closed by the SOH abelianness**, not open — but the
*bare-GR* reading of the arc is not survivable: bare geometrodynamics
would make the throats bosons. The thesis label changes from "Pauli from
GR" to **"Pauli from GR + the forced Pin⁻ framing"**, and the #170/#171
spinoriality sentence is corrected per Lemma 3.

## References

- J. L. Friedman, R. D. Sorkin, *Spin ½ from gravity*, PRL 44 (1980) 1100.
- R. D. Sorkin, S. Surya, *An analysis of the representations of the
  mapping class group of a multi-geon three-manifold*, IJMPA 13 (1998)
  3749, gr-qc/9605050. [MCG presentation; slides s²=e for RP³; slide
  subgroup Z₂∗Z₂ for two geons; abelian weights in the SOH; anomalous
  sectors.]
- H. F. Dowker, R. D. Sorkin, *A spin-statistics theorem for certain
  topological geons*, CQG 15 (1998) 1153, gr-qc/9609064. [SSC for
  non-chiral abelian geons pair-created from the background; no
  obstruction without topology change.]
- H. Hendriks, Bull. Soc. Math. France Mém. 53 (1977). [Non-spinorial
  classification: lens spaces and handles.]
- D. Giulini, *Mapping-class groups of 3-manifolds in canonical quantum
  gravity*, math-ph/0606066; *Matter from space*, 0910.2574. [Tabulation:
  cyclic-π₁ primes non-spinorial; internal MCG of RP³ trivial.]
- J. Louko, D. Marolf, *Inextendible Schwarzschild black hole with a
  single exterior*, PRD 58 (1998) 024007. [The RP³ geon spacetime.]
- C. Aneziris, A. P. Balachandran, M. Bourdeau, S. Jo, T. R. Ramadas,
  R. D. Sorkin, *Statistics and general relativity*, IJMPA 4 (1989) 5459.
- D. McCullough, *Mappings of reducible 3-manifolds*, in Geometric and
  Algebraic Topology, Banach Center Publ. 18 (1986). [Generation theorem:
  internals, exchanges, slides; Fouxe-Rabinovitch relations.]
- Amphichirality of lens spaces: L(p,q) reverses orientation iff
  q² ≡ −1 (mod p) (standard; see e.g. Hodgson–Rubinstein).

## Reproduce (the machine-checkable steps)

```bash
python -m experiments.closure_ledger.geon_statistics_adjudication_probe
# Verdict: SSC_THEOREM_HYPOTHESES_HOLD_SIGN_IS_FRAMING_DEPENDENT
#          _BARE_GR_SELECTS_BOSE_PIN_MINUS_FRAMING_SELECTS_FERMI
```
