# Entanglement swapping is bridge surgery: linking never-co-nucleated throats (PR #207)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR opens the dynamical half of
> #206's register consequence: how throats that never shared a
> nucleation become entangled. The companion probe measures both sides
> (algebra and lattice dynamics) in ~8 s.

## 0. The mechanism

#206 derived the entangled sector for throats that *share* a
nucleation: one bridge, one shared fiber, the singlet from the Pin⁻
transport, CHSH = 2√(1+C²) from bridge content. The open half was
N-body reach: by #206's own theorem (the LHV cap, enumerated), throats
that never shared a nucleation cannot become Bell-correlated through
any 3-space-local dynamics. The mechanism that closes the gap:
**throats interact**.

> A proximity pair can pinch off in pair annihilation — the #58/#200
> topology-change channel run in reverse, returning the pair's mass to
> wave fronts — and the pinch **relinks the bridges** of two previously
> unrelated pairs in the bulk, leaving the two distant survivors
> connected exactly as if they shared a nucleation.

That is **entanglement swapping**, realized as bridge surgery: pairs
(1,2) and (3,4); mouths 2 and 3 meet and annihilate (a *local* event);
the composite bridge 1–4 remains. This PR makes the correspondence
quantitative.

## 1. The composition law is the quantum swapping law

**QM side (16-dim, machine-checked).** For pair states
Ψ_a(1,2) ⊗ Ψ_b(3,4) (with Ψ_φ = (|+−⟩ + e^{iφ}|−+⟩)/√2, the
charge-zero Bell circle), projecting mouths (2,3) onto Ψ_c leaves (1,4)
in **Ψ_{a+b+c}** — fidelity 1 − 10⁻¹², each of the four outcomes with
probability exactly ¼.

**Bulk side (measured, §2).** Bridge surgery composes transports; the
surgered chain's (1,4) pair phase obeys

```
φ₁₄ = φ(s₁₂) + φ(s₃₄) + π s_j / 2  =  φ_a + φ_b + φ_c ,
```

the *same* a+b+c law — with the Bell-outcome phase φ_c supplied by the
**junction holonomy s_j of the pinch**: quantum swapping's measurement
outcome is a topological datum of the local annihilation event.

**Superselection.** Annihilation conserves winding — the (2,3) pair
must carry zero net charge — so the surgery orbit is confined to the
charge-zero (Ψ) sector: the charged Bell states are not annihilation
outcomes, exactly as a charge-conserving Bell measurement is
sector-restricted in QM. And the **unconditioned mixture** over the
four holonomies is **separable** (negativity 0, machine-checked):
without the outcome record there is no usable (1,4) correlation —
swapping's classical-communication requirement and no-signaling
(#204), in one identity.

## 2. The composite bridge, live (the #206 lattice, extended)

Four mouths on a 192-site ring × 8-site fiber; bridges (1↔2) and (3↔4)
with holonomies s₁₂ = s₃₄ = 2; the junction glues (m₂, χ) ↔
(m₃, (s_j − χ) mod 8) — a **local** operation (the mouths are 8 sites
apart). Measured:

- **End-to-end conjugation:** winding k = +1 into mouth 1 arrives at
  mouth 4 as k = −1 with channel purity **0.999832** — the composite
  bridge conserves k₁ + k₄ = 0 (the superselection of §1, live).
- **The four-outcome orbit:**

| s_j | φ₁₄ measured | predicted | state fidelity |
|---:|---:|---:|---:|
| 0 | 0.0000 | 0.0000 | 1.00000 |
| 1 | 1.5874 | 1.5708 | 0.99991 |
| 2 | 3.1650 | 3.1416 | 0.99983 |
| 3 | 4.7288 | 4.7124 | 0.99992 |

The four Bell outcomes of quantum swapping = the four fiber holonomies
at which the proximity pair can pinch.

## 3. The annihilation event

Two *populated, disjoint* bridges — (1,2) carrying Alice's preparation
phase α₁, (3,4) prepared independently. At t_j the proximity pair (2,3)
is annihilated: their fibers glue and their binding wells vanish (the
pinch-off). Measured:

- **Before:** the distant pair (1,4) shows no response to α₁ —
  differential 3×10⁻⁶ (the ballistic Lieb–Robinson tail).
  **After:** 1.3×10⁻² — a factor **4650**.
- **The response lands in exactly the swapped channel** (α₁ rides m₁'s
  k = −1 mode, whose surgery image is m₄'s k = +1): channel power ratio
  ≈ 290.
- **The linkage carries Alice's phase linearly:** the transited
  amplitude obeys D(α) ∝ (e^{iα} − 1) — magnitude ratio 1.99 vs
  predicted 1.99, phase step 0.100 vs predicted 0.100.
- **The pair returns to wave fronts:** the middle bound population
  falls 0.502 → 0.285 → 0.20 (t = 2 → 5 → 8) while the
  no-annihilation control keeps ≈ 0.54, the freed debris dispersing
  into propagating waves (population outside every mouth neighborhood
  rises 0.000 → 0.045 and climbing as the fronts clear the windows).

A local event; a nonlocal-in-projection correlation;
**never-co-nucleated mouths, linked** — mouths 1 and 4 shared no gluing
at any time.

## 4. The swapped pair is a Bell pair; the network is reachable

- The extracted (1,4) effective states across the holonomy orbit reach
  **CHSH ≥ 2.8280** (Tsirelson 2.8284) — maximal Bell violation between
  never-co-nucleated throats.
- **The repeater:** surgeries compose associatively (the a+b+c law
  applied twice on a three-bridge chain reproduces the total-sum
  prediction exactly; end pair at CHSH = 2√2): a chain of local
  annihilations links **any** two throats in the network — pairwise
  entanglement distribution over the throat network is
  **mechanism-complete**.

## 5. Monogamy is the matching

Bridges pair mouths: the bridge configuration of a many-throat universe
is a **perfect matching**, and #206 gives one Bell pair per bridge.
Machine-checked on a three-bridge (six-mouth) state: a mouth with its
bridge partner — negativity **1/2** (maximal); a mouth with any
non-partner — negativity **0** (exactly unentangled). The monogamy of
maximal entanglement is not an inequality here; it is the statement
that a mouth has exactly one bridge. Surgery is the only rewiring move
(two pairs in → one pair + radiation out). What the matching *cannot*
reach is genuinely multipartite entanglement (GHZ/W) — that would
require junctions joining three or more mouths into one bridge: a
sharp, named open construction.

## 6. What is and is not established (honest scope)

**Established:** the composition law = the QM swapping law (both sides
machine-checked/measured); the junction holonomy as the Bell outcome;
the winding superselection and the separable unconditioned mixture
(no-signaling + the classical-communication requirement); the live
surgery event (local pinch, distant linkage, radiated middle); Tsirelson
saturation for never-co-nucleated mouths; monogamy as the matching;
repeater associativity.

**Not established (the conditions):**
1. The pinch is modeled by its **topological content** — fiber gluing
   plus well removal; the full 5D pinch-off dynamics (the #58/#200
   topology-change channel in reverse, with its Sorkin causal-continuity
   structure) is not solved here.
2. The **outcome distribution** over the four holonomies is not
   derived: which holonomy a given annihilation realizes is set by the
   local beables of the event; QM's ¼-each is recovered under the
   equilibrium hypothesis (#198/#204, dBB grade).
3. Still the **internal** (fiber/winding) sector; spatial entanglement
   and measurement transport remain open (#206 scope).
4. **Multipartite (≥ 3-mouth) junctions** — the GHZ sector — are the
   named open construction.

**The register moves:** #206's dynamical half narrows to (i) the
spatial/measurement sector and (ii) multipartite junctions — both
named, both sharp. Pairwise, the entangled sector of BAM is now
mechanism-complete: nucleation makes Bell pairs, annihilation swaps
them, monogamy is the matching, and the Bell outcome is a holonomy.

## References

- M. Żukowski, A. Zeilinger, M. A. Horne, A. K. Ekert, PRL 71 (1993)
  4287. [Entanglement swapping.]
- H.-J. Briegel, W. Dür, J. I. Cirac, P. Zoller, PRL 81 (1998) 5932.
  [Quantum repeaters.]
- J. Maldacena, L. Susskind, Fortsch. Phys. 61 (2013) 781. [ER=EPR;
  here the swapping face, classical and quantitative.]
- The BAM chain: PR #206 (entanglement is bridge topology), #200 (the
  pair-creation/annihilation cobordism), #58 (the pair threshold),
  #198/#204 (equilibrium; no-signaling).

## Reproduce

```bash
python -m experiments.closure_ledger.bridge_surgery_swapping_probe
# Verdict: ENTANGLEMENT_SWAPPING_IS_BRIDGE_SURGERY_THE_JUNCTION_HOLONOMY_IS
#          _THE_BELL_OUTCOME_LOCAL_ANNIHILATION_LINKS_NEVER_CONUCLEATED
#          _THROATS_MONOGAMY_IS_THE_MATCHING
```
