# Configuration-space emergence: entanglement is bridge topology (PR #206)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR attacks #198's condition 2 at
> its sharpest point — the emergence of the entangled sector — and makes
> "classical ER=EPR" quantitative. The companion probe machine-checks
> every step and runs the lattice adjudication.

## 0. The sharp target

The honest statement of #198's condition 2 ("subsystem effective
wavefunctions; the general nonlinear measurement theory is open")
conceals its strongest form, which should be stated as the target:

> **A single classical field on 3-space with local dynamics and local
> readout is a local-hidden-variable model, and cannot violate CHSH.
> That is Bell's theorem itself.**

The field configuration is λ; detectors respond locally; without
anything else, the dBB guidance of two independent throats factorizes.
The probe enumerates all 16 deterministic local strategies: max CHSH
= 2, exactly. The quantum sector sits at 2√2 (Tsirelson — computed from
the singlet). **BAM's entangled sector cannot emerge from field
correlations alone** — no refinement of the 3-space dynamics can cross
that gap.

So it must come from the one nonlocal element the theory legitimately
owns: **the bridge**. The two mouths of a throat are one object through
the bulk — the #168 embedding (the throat as the regular 5D Killing
horizon of the Tangherlini bulk), the #169 J-quotient, and
`bell.bulk_identity`'s kinematic reading ("one continuous 5D defect
observed from two boundary frames"). Adjacent in the bulk, distant on
the brane: nonlocal *in the 3-space projection only*. That is the
budget this PR spends, and the concrete problem is:

> Derive the effective ψ(x₁, x₂) of a throat pair from the universal
> 3-space field plus the bulk identification, and show the bulk
> connectivity is exactly what carries the non-factorizable
> correlation.

## 1. The identification, live on a lattice

The sandbox: a universal field on (x, χ) — a 128-site 3-space ring
carrying an 8-site fiber — with **local dynamics everywhere except the
two mouth cells, which are glued through the bulk**:
(x_A, χ) ↔ (x_B, (s − χ) mod 8). Local in the full topology; nonlocal
only when projected to the ring. Measured:

1. **Charge conjugation.** Winding k = +1 sent into mouth A arrives at
   mouth B in channel k = −1 with purity 1.000000 — the C-conjugate
   pair (winding = charge, #42–#44; Σc₁ = 0, #58/#200), exactly as the
   orientation-reversing identification requires.
2. **Phase locking.** The inter-channel phase at B tracks the phase at
   A with slope 1.00000000: the two mouths hold **one** phase, not two.
3. **The holonomy is the Bell-state selector.** The gluing shift s — a
   π fiber holonomy of the handle, a *topological datum* — offsets the
   locked phase by exactly π (s = 2): it selects the singlet sign
   versus the triplet-0 sign.
4. **One object.** The B-side/A-side readout ratio of the shared mode
   is identical across channels to 10⁻⁸: both mouths read the *same*
   variable.
5. **The cut control.** Remove the handle: B receives 10⁻³⁰. Everything
   above is carried by the bridge and nothing else.

## 2. The emergence, as three machine-checked lemmas

**Lemma A (the tensor structure).** One shared fiber mode, read at both
mouths through the bridge transport T, defines the map
`W|k⟩ = |k⟩_A ⊗ T|k⟩_B` — verified an **isometry**. A single shared
circle described by two local observers *is* a maximally correlated
two-party state. The tensor-product structure of the pair is not
postulated: it is the two-frame description of one bulk object. This is
where configuration space comes from.

**Lemma B (which state).** The repo's derived non-orientable throat
transport is T = iσ_y with **T² = −1** — the same Pin⁻ minus sign that
forces Fermi statistics (#195/#196). Applied to the symmetric bridge:

```
ψ_eff = (I ⊗ T)|Φ⁺⟩ = the singlet        (fidelity 1.000000000000)
```

— which is exactly the pair state `bell.bulk_identity` *postulated*
from topology. What that module assumed is now **derived** from field
plus identification; two independent paths, one state. The exchange
sign and the entanglement sign have a single origin: the bridge
transport that squares to −1.

**Lemma C (the law).** The derived state's correlation is
E(a, b) = −cos(a − b), matching the `bulk_identity` module pointwise to
2×10⁻¹⁶.

## 3. The quantitative ER=EPR law

Sweeping the bridge preparation (the nucleation stand-in: bridge-mode
superposition cos β |+1⟩ + sin β |−1⟩) and reading the effective
two-observer state off the evolved field:

| β | sin 2β | C extracted | S_ent | CHSH | 2√(1+C²) |
|---|---:|---:|---:|---:|---:|
| 0 | 0 | 0.000 | 0.000 | 2.000 | 2.000 |
| π/12 | 0.500 | 0.500 | 0.246 | 2.236 | 2.236 |
| π/8 | 0.707 | 0.707 | 0.417 | 2.449 | 2.449 |
| π/6 | 0.866 | 0.866 | 0.562 | 2.646 | 2.646 |
| π/4 | 1.000 | 1.000 | 0.693 | 2.828 | 2.828 |

- **Schmidt weights = bridge-mode amplitudes** (extracted concurrence
  tracks sin 2β to 2×10⁻³ across the sweep); entanglement entropy =
  bridge participation entropy.
- **CHSH(ψ_eff) = 2√(1 + C²) exactly** (Horodecki criterion, verified
  by direct optimization over measurement settings), running from the
  Bell bound 2 at zero bridge coherence to Tsirelson 2√2 at the
  symmetric bridge — where the extracted state *is* the singlet
  (fidelity 1.000000).
- **The cut control:** two throats without a bridge give a product
  state (entropy 0) and CHSH = 2.000000000 — on the LHV side of the
  gap, as Bell demands of any 3-space-local remainder.

**Entanglement — and its Bell violation — is a measured function of
bridge topology.**

## 4. The nonlocality budget, audited

- **Marginals:** Bob's reduced state of the derived ψ_eff is invariant
  under every Alice-side unitary (10⁻¹⁶); the module-level marginals
  are setting-independent (10⁻¹⁶). The bridge supplies Bell
  correlations, not a telegraph.
- **Statistics:** the operational probabilities ride on quantum
  equilibrium — Born-distributed beables on the effective configuration
  space (#198, dBB grade), whose equilibrium signal-locality (and
  non-equilibrium Valentini signal) **#204 measured on this very
  structure**.
- **Non-traversability, honestly:** the physical bridge is
  non-traversable — the identification is a property of the geometry,
  not a transport channel. The correlation is *imprinted at
  nucleation* — the two mouths are born as the boundary of **one
  2-handle** (#200's explicit cobordism) — and conserved topologically
  thereafter. The lattice handle is a modeling stand-in for that
  identification; what it demonstrates (conjugation, locking, holonomy,
  one-object readout) are properties of the identification itself.

The triple {bulk-geometric correlation, no superluminal signal, Bell
violation at the effective level} is exactly the triple quantum
mechanics itself exhibits.

## 5. What is and is not established (honest scope)

**Established:**
- the LHV cap for the bridge-free sector (enumerated exactly) and the
  quantum gap;
- the identification's field-level signature (conjugation, locking,
  holonomy, one-object readout; all killed by the cut);
- the tensor-product structure and the singlet, derived — the
  previously postulated `bulk_identity` state recovered from field +
  identification;
- the quantitative law: Schmidt = bridge amplitudes,
  CHSH = 2√(1 + C²) from bridge content.

**Not established (the conditions):**
1. The derivation covers the pair's **internal** (fiber/winding =
   spin/charge) sector — where Bell lives. The spatial part of
   ψ(x₁,x₂) and the full measurement dynamics (Stern–Gerlach-type
   transport on the effective space) are not derived; Born statistics
   enter at dBB grade (#198).
2. **N bodies:** one shared variable per bridge; many-pair states are
   tensor products over bridges. Entangling throats that never shared a
   nucleation (bridge networks / swapping) is an open construction —
   the mechanism is exhibited, not the general case.
3. The lattice handle is traversable (a coupling) — a stand-in for the
   nucleation-imprinted identification of the non-traversable physical
   bridge.
4. N_χ = 8 fiber sites; the k = ±1 doublet is the scalar reduction of
   the #195/#197 spinor structure (T = iσ_y is the repo's derived
   transport).

**The register consequence:** #198's condition 2 **splits**. Its
Bell-sector half — where does the non-factorizable structure come
from? — is *discharged*: from the bridge, quantitatively. Its dynamical
half (spatial sector, N-body networks, measurement transport) remains
the standing item, now sharply bounded.

## References

- J. S. Bell, Physics 1 (1964) 195; B. S. Cirel'son, Lett. Math. Phys.
  4 (1980) 93. [The bound and the quantum maximum.]
- R. Horodecki, P. Horodecki, M. Horodecki, Phys. Lett. A 200 (1995)
  340. [CHSH maximum for two-qubit states.]
- J. Maldacena, L. Susskind, Fortsch. Phys. 61 (2013) 781. [ER=EPR —
  here in a classical-geometric, quantitative form.]
- The BAM chain: PR #167–#169 (the 5D embedding and the J-quotient),
  `bell.bulk_identity` (the kinematic one-object reading; the T = iσ_y
  transport), PR #195/#196 (the Pin⁻ sign), PR #198/#204 (equivariance,
  equilibrium signal-locality), PR #200 (the pair-nucleation 2-handle).

## Reproduce

```bash
python -m experiments.closure_ledger.configuration_space_emergence_probe
# Verdict: CONFIGURATION_SPACE_EMERGES_FROM_THE_BRIDGE_SCHMIDT_WEIGHTS_ARE
#          _BRIDGE_MODE_AMPLITUDES_CHSH_2_WITHOUT_2SQRT2_WITH_CLASSICAL
#          _ER_EPR_QUANTITATIVE
```
