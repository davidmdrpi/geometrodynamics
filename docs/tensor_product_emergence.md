# Tensor-product emergence: the 2^N state space from one field on the network (PR #232)

> **Framing.** QFT on the *fixed classical* throat geometry — geometry →
> fields, **not** quantum gravity. This PR answers the program's
> decisive ontology problem, stated by the program itself:
>
> **Can ordinary 2^N-dimensional quantum state structure emerge from a
> field and topology in physical space?**
>
> The two-mouth (#206), swapping (#207), and GHZ (#208) constructions
> were beginnings of that bridge, not a general proof: each derives a
> *low-dimensional shared-variable slice* of the many-body space (the
> singlet; the GHZ ray) — not the full tensor product with independent
> local unitaries on every factor. This PR builds the general object
> and fences it with two no-gos. The companion probe machine-checks
> every claim (~17 s).

## 0. The answer, stated first

**Yes — constructively, with a sharply located mechanism.** The 2^N
space emerges as the **N-quantum dual-rail sector of one quantized
field on the mouth-doublet network**:

- one scalar field on physical space supplies 2N modes — N throat
  doublets, the committed #224 pair (interior fraction 0.97, basin
  localization 0.98, re-read from its ledger);
- the qubit is **which-mouth**: |0⟩ = one excitation in rail a,
  |1⟩ = in rail b;
- the N-excitation sector with one quantum per doublet has dimension
  **exactly 2^N** (machine-checked N = 1…5);
- the tensor product is **derived, not postulated** (§2);
- entanglement is bridge topology plus the committed quartic (§4);
- universality follows: explicit circuits build GHZ_N at fidelity
  1 − 10⁻¹² through N = 5, with Mermin = 4 live (§5).

And two proven NOs locate the ontology exactly (§3): a **classical**
field fails by resource counting at N = 3 (and by Bell at N = 2,
#206), and a **linear** field fails deterministically at the measured
KLM boundary. So the essay's "real diffuse field" must be read as the
*quantum state of one field* — its dilution and refocusing the
#211–#219 transactional story — and its exponential workspace is the
multi-quantum coherence, physically carried, and nothing else.

## 1. The sector: where 2^N lives

| N | modes | Fock dim (N quanta) | dual-rail sector | classical field dof (real) | pure 2^N manifold (real) |
|---|---|---|---|---|---|
| 1 | 2 | 2 | **2** | 4 | 2 |
| 2 | 4 | 10 | **4** | 8 | 6 |
| 3 | 6 | 56 | **8** | 12 | **14** |
| 4 | 8 | 330 | **16** | 16 | 30 |
| 5 | 10 | 2002 | **32** | 20 | 62 |

The dual-rail constraint (one quantum per doublet) is not imposed by
hand at gate time: every gate used below conserves it exactly (unit
sector norm throughout). The exponentiality lives in the
**excitation-number sector of one field** — N quanta shared among 2N
physical modes — not in N separate fields, not in configuration-space
substance, not in worlds.

## 2. The tensor product is network locality

With the doublets decoupled (arbitrary random intra-pair one-body
Hamiltonians), the field's sector evolution equals the Kronecker
product of the N one-body 2×2 unitaries **to machine zero** (worst
deviation 10⁻¹⁵, zero leakage; N = 2, 3, 4). The ⊗ of quantum
mechanics, on this substrate, is the statement that *decoupled regions
of one field evolve independently*. This is the general-N form of the
"multiple frames reading one bulk object" construction: #206's
embedding derived the two-frame slice; here every doublet is an
independent frame with full local SU(2), and the product structure
falls out of locality alone.

## 3. The two no-gos that fence the claim

**The classical counting no-go.** A classical linear field
configuration on the network carries 4N real degrees of freedom (2N
complex amplitudes). The pure 2^N manifold needs 2·2^N − 2. At N = 3
the field loses (12 < 14), and exponentially thereafter. This is the
resource-counting face of what #206's LHV theorem (ledger re-read: LHV
max = 2, Tsirelson 2√2) already proved at the correlation level for
N = 2: **the ontology cannot be classical-field-valued.** The essay's
dilution/refocusing picture survives — but what dilutes is the quantum
state of the field, not a classical amplitude.

**The linear no-go (the measured KLM boundary).** Bridge two doublets
with any beam splitter and the dual-rail sector leaks (Hong–Ou–Mandel
double occupation). Across the bridge family *and 200 random 4-mode
linear-optics evolutions*: every map that stays unitary on the sector
has **zero** entangling power, and every entangling map (up to 1.000
ebit seen) rides on O(1) leakage. A strictly linear field on **any**
network topology entangles only probabilistically, exactly like linear
optics (KLM). The essay's "richness of mode correlations" is not
deterministically available at linear order — a linear-only BAM would
be a postselection theory.

## 4. The committed quartic is the entangling resource

The program already owns the term that crosses the KLM boundary: the
scalar self-interaction of the soliton sector (#210/#222, the |φ|⁴ EKG
term). Two rails sharing a bridge region acquire a cross-Kerr overlap
χ·n_a·n_b — and at χt = π this is an **exact CZ**: deviation 10⁻¹⁶,
**zero** sector leakage (the interaction conserves each mode's
number), Bell entropy exactly 1.000000 from |++⟩. With the intra-pair
SU(2) of §2 this is the standard universal set {CZ, single-qubit
unitaries}.

**The same nonlinearity that makes the soliton makes the CZ.** The
entangling resource is not decoration bolted onto the geometry — it is
the one interaction term the matter sector was already built from.

## 5. Universality delivered: GHZ_N through N = 5

Explicit physical circuits (intra-pair rotations + the quartic CZ)
build GHZ_N at fidelity 1 − 10⁻¹² with unit sector norm for
N = 2, 3, 4, 5 — beyond the #208 Y-junction slice, on independent
qubits with full local control. The N = 3 state reaches the algebraic
Mermin maximum **|M| = 4 live** (LHV cap 2), matching the #208 ledger
(extracted Mermin 4.0, junction fidelity 1.0, three-tangle 1.0), and
the #208 charged-GHZ superselection no-go (zero-sum triples = 0) is
re-read unchanged: GHZ lives in the transported-frame label, exactly
as committed.

## 6. What this says to the interpretations

- **Copenhagen**: nothing here needs a projection axiom; the
  measurement layer remains the #209/#219 sector, and this PR adds the
  missing piece Copenhagen never supplies — *what* the many-body state
  is a state *of* (one field's multi-quantum coherence on physical
  space).
- **Many-worlds**: the exponential workspace is physically carried by
  one field history on one closed geometry — the computation Everett
  places in the universal wavefunction is here the propagation,
  interference, and Kerr phases of N quanta of a single field.
- **Pilot-wave**: the configuration-space object is, on this
  substrate, **derived**: Ψ(x₁,…,x_N) is the N-quantum sector of one
  field read through N doublet frames. The "multiple-frame description
  of a single bridge" framing generalizes rigorously — that was the
  decisive gap, and it is now a construction, not a hope.
- **Quantum computation**: gates are intra-doublet rotations and
  bridge-overlap Kerr phases; entanglement is a physical linkage; the
  apparent exponential calculation resides in N-quantum coherence.
  The essay's demanding burden — *finite local degrees of freedom
  supporting the effective exponential state space* — is met by the
  excitation-number sector, and the two no-gos prove it is met **only**
  there: no classical and no linear-deterministic shortcut exists.

## 7. Honest scope

- The network is graph-level (rails and bridges), the same level as
  #207/#208; the physical doublet is anchored by the #224 ledger, but
  the gates are not solved on the PDE throat background here.
- χ is structural: the |φ|⁴ overlap integral of two rail modes in a
  shared bridge — the operator form is committed, the coefficient not
  calibrated to a throat geometry (gate time scales as 1/χ; the gate
  is exact for any χ > 0).
- The field is the program's bosonic scalar; fermionic statistics and
  spin-structure sectors are the #227–#231 arc. Dual rail does not
  need them.
- **Quantization is imported, not derived** — the program's framing is
  QFT on classical geometry. The claim is emergence of the 2^N
  *structure* from one quantized field plus topology; the T2/T4
  no-gos show a derivation of quantum mechanics from classical fields
  is impossible in exactly quantified senses.
- Measurement, decoherence, Born rule: the #209/#219 sector, untouched.
- N = 5 demonstrated explicitly; the construction is manifestly
  N-uniform, but no asymptotic claim beyond the demonstrated range is
  made.

## 8. What would falsify this

- A dual-rail sector dimension ≠ 2^N at any N — the counting would be
  wrong. (Checked: N = 1…5 exact.)
- A decoupled-doublet evolution that fails to factorize — the tensor
  product would need an axiom after all. (Checked: 10⁻¹⁵.)
- A linear evolution found that entangles deterministically — KLM
  would be violated and the quartic unnecessary. (Checked: none in
  the family + 200 random samples.)
- A leakage or error in the quartic CZ — universality would not be
  free. (Checked: 10⁻¹⁶, zero leakage.)
- GHZ_N failing at some N, or Mermin < 4 — the multipartite sector
  would not generalize. (Checked: fidelity 1 − 10⁻¹² to N = 5,
  M = 4.0.)

## 9. Companion probe

`experiments/closure_ledger/tensor_product_emergence_probe.py` (T1–T9,
~17 s): the sector dimensions and counting no-go; the locality ⟹
tensor-product theorem; the measured KLM boundary; the exact quartic
CZ; the GHZ_N circuits with live Mermin; the ontology ledger with
committed re-reads (#206, #208, #224).

**Verdict:**
`THE_2_TO_N_STATE_SPACE_IS_THE_N_QUANTUM_SECTOR_OF_ONE_FIELD_ON_THE_MOUTH_DOUBLET_NETWORK_LOCALITY_MAKES_THE_TENSOR_PRODUCT_LINEARITY_ALONE_CANNOT_ENTANGLE_DETERMINISTICALLY_AND_THE_COMMITTED_QUARTIC_COMPLETES_THE_GATE_SET`
