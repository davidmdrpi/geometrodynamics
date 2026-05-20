# Tree-level QED from BAM geometry — status / release note

**Thread:** PR #35 → PR #46 (the "QFT-event reinterpretation"
amplitude programme).

**Status:** the scalar-intensity structure of every tree-level
`2 → 2` QED process — Compton, Breit–Wheeler, pair annihilation,
Bhabha, Møller — is now reproduced from a small, fixed set of
BAM-geometric primitives, with no remaining QED-specific ansatz at
tree level. Each ingredient (vertex factor, propagator, spinor
trace, Fermi-statistics sign, photon polarisation) is derived from
S³ closure, the Hopf bundle, the non-orientable throat transport,
or the S³ Green function. Every claim is backed by a probe with a
machine-precision numerical test.

This note collects the result. It is a synthesis, not a new probe;
the per-PR research plans under `docs/` carry the detailed
derivations and the run archives under
`experiments/closure_ledger/runs/` carry the numerical evidence.

---

## 1. The headline

> All tree-level `2 → 2` QED photon-mediated scalar intensities are
> reproducible from BAM geometry. The closed-form Compton vertex
> factor, its crossing partners, the 4-fermion Dirac-trace
> numerators, the Fermi-statistics interference signs, and the
> photon propagator (scalar and full Lorentz-tensor) all derive from
> the same geometric substrate.

The substrate, in full:

| BAM-geometric primitive | QED ingredient it produces | PR |
|---|---|---|
| Antipodal `S³` Green function (pole) | Compton propagator pole | #35 |
| Transverse Hopf-fibre polarisation | Thomson angular factor `(1+c²)/2` | #35, #38 |
| Single throat action (closure 2π + antipodal symmetry + stationary action) | closed-form vertex factor `F²(x, c)` | #41 |
| Mandelstam crossing of `F²` | Breit–Wheeler + pair annihilation | #36, #37 |
| Harmonic-mean throat-rate (equal-action) | kinematic Padé `K(x) = 2x/(1+x)` | #39 |
| Hopf-fibre helicity spinor (equal spin-action) | polarisation factor `Q(x, c)` | #40 |
| SU(2) Hopf-bundle Pauli/Weyl traces | 4-fermion numerators `(s²+u²)`, `(u²+t²)`, `(s²+t²)` | #43 |
| Non-orientable throat transport `T = iσ_y = ε` | Fermi-statistics interference signs | #44 |
| `S³` scalar Green function (flat limit) | photon propagator magnitude `1/q²` | #45 |
| Hopf-bundle U(1) connection | photon propagator Lorentz tensor `−η^{μν}/q²` | #46 |

---

## 2. The two sub-threads

### 2a. Compton and its crossing triangle (PRs #35–#41)

The Compton amplitude is reproduced by a BAM construction with a
single closed-form **vertex factor**:

```
F²(x, c) = 4·x³·(x² + 1 − x·sin²θ) / [(1 + c²)·(1 + x)²]
```

with `x = ω'/ω = 1/(1 + ε(1 − cos θ))`, `c = cos θ`. The BAM
amplitude `f_BAM_baseline · F²` reproduces Klein–Nishina **exactly
at all orders in ε up to ε ~ 2** (machine precision, PR #35).

The factor decomposes geometrically (PR #38):

```
F²(x, c) = K(x)² · Q(x, c)
K(x)     = 2x/(1+x)                          (caustic / throat-rate Padé)
Q(x, c)  = x² + x·(1−x)²/(1+c²)              (Hopf-fibre helicity channels)
```

via the exact identity `x² + 1 − x·sin²θ ≡ (1−x)² + x·(1+c²)`.

Both factors are *derived*, not fitted:

  - **`K(x)`** (PR #39) is the harmonic mean of the in/out photon
    frequencies, the unique flux-continuous throat-rate; alternative
    means (arithmetic, geometric, linear) leave `Q` non-polynomial.
  - **`Q(x, c)`** (PR #40) is the squared norm of a two-component
    Hopf-fibre helicity spinor with helicity-preserving amplitude
    `A_pres = x` and helicity-flipping amplitude
    `A_flip = √x·(1−x)/√(1+c²)`. The `(1+c²)` normalisation is the
    Wigner-d¹ helicity sum `cos⁴(θ/2) + sin⁴(θ/2)`.

Both equal-action splittings (energy → `K`, spin → `Q`) follow from
a single **BAM throat action** (PR #41) under three principles:
closure quantum `S = 2π` (`action_base`), `S³` antipodal symmetry
`σ(p) = −p`, and stationary action. Alternative principles are
rejected.

Under standard Mandelstam crossing the same `F²` reproduces
**Breit–Wheeler** `γγ → e⁺e⁻` (PR #36) and **pair annihilation**
`e⁺e⁻ → γγ` (PR #37); the Compton/BW/annihilation crossing triangle
closes at machine precision, the loop being the identity.

### 2b. The 4-fermion processes (PRs #42–#46)

Bhabha (`e⁺e⁻ → e⁺e⁻`) and Møller (`e⁻e⁻ → e⁻e⁻`) sit *outside* the
Compton/BW/annihilation crossing orbit — they have four external
fermions, not two fermions and two photons. PR #42 documented this
explicitly: the scalar Compton kernel alone cannot produce them.
Four further ingredients close the gap.

  - **Diagonal numerators** (PR #43). The parity-symmetric SU(2)
    Hopf-bundle Pauli/Weyl trace product

    ```
    T_BAM(p_a,p_b,p_c,p_d) = 32·[(p_a·p_c)(p_b·p_d) + (p_a·p_d)(p_b·p_c)]
    ```

    equals the QED Dirac trace product. Channel-specific kinematic
    pairings select the QED numerators: `(s²+u²)` for t-channel,
    `(u²+t²)` for s-channel, `(s²+t²)` for u-channel.

  - **Interference signs** (PR #44). The non-orientable throat
    transport `T = iσ_y = [[0,1],[−1,0]]` *is* the SU(2)
    antisymmetric ε tensor. One fermion-leg transposition acts as
    `−1` on the antisymmetric singlet — reproducing the Bhabha
    s↔t Wick sign and the Møller t↔u Pauli sign automatically.

  - **Propagator magnitude** (PR #45). The `S³` scalar Green
    function `G(ψ) = ((π−ψ)cot ψ − ½)/(4π²R)` reduces in the
    flat-space limit to the Coulomb potential `1/(4π d)`, whose
    Fourier transform is `1/q²` — the photon propagator, recovered
    as a geometric exchange kernel rather than a virtual particle.

  - **Propagator tensor structure** (PR #46). The Hopf-bundle U(1)
    connection gives the full Lorentz structure
    `D_F^{μν}(q) = −η^{μν}/q²`; in Feynman gauge this factors as
    `−η^{μν}` times the PR #45 scalar kernel. Gauge-mode differences
    decouple from physical amplitudes via the Ward identity
    `q_μ Tr[γ^μ p̸_1 γ^ν p̸_2] = 0`.

With these, end-to-end Bhabha and Møller `|M̄|²/(8e⁴)` match the QED
textbook formulas to machine precision.

---

## 3. Probe ledger

Each probe writes JSON + markdown archives under
`experiments/closure_ledger/runs/<timestamp>_<probe>/`.

| PR | probe | verdict |
|---|---|---|
| #35 | `compton_vertex_resummation_probe` | EXACT_RESUMMATION — closed-form `F²` reproduces KN to ε~2 |
| #36 | `breit_wheeler_cross_process_probe` | PROCESS_GENERAL_UNDER_CROSSING |
| #37 | `pair_annihilation_crossing_probe` | TRIANGLE_CLOSES |
| #38 | `throat_nucleation_caustic_derivation_probe` | DERIVATION_COMPLETE — `F² = K²·Q` |
| #39 | `two_mouth_flux_action_probe` | PADÉ_DERIVED — `K = 2x/(1+x)` |
| #40 | `hopf_helicity_transport_probe` | Q_DERIVED — helicity spinor |
| #41 | `throat_action_derivation_probe` | ACTIONS_DERIVED — both splittings from one action |
| #42 | `bhabha_moller_interference_probe` | DIAGONAL_FAILS_TOO — gap identified |
| #43 | `dirac_trace_geometry_probe` | DIAGONALS_DERIVED_INTERFERENCE_MAGNITUDE_OK — SU(2) Pauli traces |
| #44 | `mobius_exchange_sign_probe` | SIGNS_DERIVED — `T = iσ_y = ε` |
| #45 | `bam_exchange_kernel_probe` | PROPAGATOR_FROM_GEOMETRY — `1/q²` |
| #46 | `hopf_vector_exchange_kernel_probe` | VECTOR_KERNEL_FROM_HOPF_BUNDLE |

All verdicts are PASS-class (PR #42 is a deliberate, well-characterised
negative result that scopes the 4-fermion gap subsequently closed by
PRs #43–#46).

---

## 4. What this is and is not

**Is:**

  - A demonstration that the *kinematic and algebraic structure* of
    tree-level `2 → 2` QED scalar intensities is reproducible from a
    fixed BAM-geometric substrate, process by process, to machine
    precision.
  - A set of falsifiable probes: each derivation has a numerical test
    that would fail if the geometric identification were wrong.
    Several plausible alternatives are explicitly rejected (non-harmonic
    throat rates in #39; alternative flip amplitudes in #40; alternative
    action principles in #41; non-Möbius exchange conventions in #44).

**Is not:**

  - A first-principles BAM *Lagrangian* for QED. The probes show that
    QED tree structure is consistent with — and reconstructible from —
    BAM geometry, but a single action that generates all of QED by
    variation is not claimed.
  - A loop-level result. Everything here is tree level.
  - An absolute-normalisation claim. The probes work at the level of
    ratios and structure (`|M̄|²/(8e⁴)`); the overall coupling `e`
    (and `ℏ`) enters through the separate closure-ledger anchoring
    (see `docs/hbar_origin_status.md`).

---

## 5. Open problems (next-thread targets)

  - **Loop corrections.** Vertex correction, self-energy, vacuum
    polarisation. In BAM these would couple to closed throat-fibre
    loops on `S³` and the bulk radial channel, and would require a
    Faddeev-Popov / ghost treatment of the gauge-covariant Laplacian
    (flagged in PR #46).

  - **Non-Abelian extension (QCD tree level).** The 4-fermion
    machinery used SU(2)/U(1) Hopf-bundle structure. A non-Abelian
    Hopf-bundle generalisation (SU(2)/SU(3) connections on `S³`)
    would target `qq → qq`, `q q̄ → q q̄`, and gluon vertices. BAM
    already has a separate quark-sector programme; connecting the two
    is open.

  - **First-principles BAM action for the full vertex.** PR #41
    derived the equal-action splittings from a throat action, but a
    single covariant action generating the entire QED vertex +
    propagator + spinor structure by variation is not yet written.

  - **High-`q²` / Planck-scale curvature regime.** The propagator
    derivations (PRs #45, #46) take the flat-space limit `R → ∞`. At
    `q² ∼ 1/R²` the `S³` curvature corrections (`2/R²` Ricci mass for
    the vector kernel) become non-negligible — a potential
    BAM-specific deviation from flat QED, and a possible falsification
    handle.

---

## 6. Cross-references

Research plans (one per PR):

  - `docs/compton_vertex_resummation_research_plan.md` (#35)
  - `docs/breit_wheeler_cross_process_research_plan.md` (#36)
  - `docs/pair_annihilation_crossing_research_plan.md` (#37)
  - `docs/throat_nucleation_caustic_derivation_research_plan.md` (#38)
  - `docs/two_mouth_flux_action_research_plan.md` (#39)
  - `docs/hopf_helicity_transport_research_plan.md` (#40)
  - `docs/throat_action_derivation_research_plan.md` (#41)
  - `docs/bhabha_moller_interference_research_plan.md` (#42)
  - `docs/dirac_trace_geometry_research_plan.md` (#43)
  - `docs/mobius_exchange_sign_research_plan.md` (#44)
  - `docs/bam_exchange_kernel_research_plan.md` (#45)
  - `docs/hopf_vector_exchange_kernel_research_plan.md` (#46)

Foundational thesis context: `docs/THESIS.md` (QFT-event
reinterpretation section). Dimensional / `ℏ` anchoring:
`docs/hbar_origin_status.md`.
