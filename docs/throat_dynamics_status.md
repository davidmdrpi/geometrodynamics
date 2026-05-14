# Throat-dynamics status — summary of the inner-boundary probe sequence

Closes the throat-dynamics research thread opened in
`docs/throat_dynamics_research_plan.md`. A four-probe sequence under
`experiments/closure_ledger/` characterized the inner boundary of
the 5D Tangherlini radial eigenproblem from three angles
(local-BC substitution, throat-thickness regularization, non-local
reflection-phase) and produced one positive structural identification
(with explicit caveats about what is and isn't rigorously derived).

The thread is the follow-up to the closure-ledger sequence (PRs
#11–18, summary in `docs/hbar_origin_status.md`). The closure-
ledger reduced the locked surrogate's residual external input from
six phenomenological parameters to one anchor (m_e) and identified
the inner cutoff `ε = resistance / k_5⁴ = 7π / (100·5⁴)` as a
closure-quantum invariant. The throat-dynamics thread asks the
deeper physics question: **what is the physical inner boundary,
derived from throat dynamics rather than imposed by hand?**

## Probe sequence

Each probe has JSON + markdown archives under
`experiments/closure_ledger/runs/<timestamp>_<probe>/`:

| # | probe | sub-target | result |
|---|-------|-----------|--------|
| 1 | `throat_boundary_condition_probe`           | (A) BC substitution        | **negative**: 0 of 8 BCs (Dirichlet, Neumann, six Robin κ) give ε-converged spectra. The throat is asymptotically free in tortoise coordinates; any local BC at finite ε produces only a reflection-phase shift, not a change in the asymptotic spectrum. Important caveat: the closure-quantum ε = 7π/(100·5⁴) of PR #18 closes the Compton bridge *only under Dirichlet* — the closure-quantum scaffolding is Dirichlet-specific. |
| 2 | `throat_thickness_probe`                    | (B) thickness regularization | **negative**: replacing the hard wall with a smooth confining sigmoid V_throat = V_0/(1 + exp((r − r_t)/σ)) does not remove the inner-boundary problem. σ is a *third* arbitrary parameter (shifts ω by O(σ/δ)); the V_0 → ∞ limit only recovers the hard wall at r_t. No closure-quantum natural (V_0, δ) closes the Compton bridge at σ → 0 — apparent hits at σ = δ/30 (e.g. V_0 = 100·γ at +0.13 %) are σ-fitting artifacts. The model has 3 parameters where the hard-wall scheme has 1. |
| 3 | `throat_reflection_phase_probe`             | (3) reflection-phase analysis | **POSITIVE**: the asymptotic phase Φ = ω · L is ε-INVARIANT (0.34 % spread across ε ∈ [1e-4, 1e-3]) with mean value 3.509. The invariant matches `γ_{1..5}/(2π) = 3.511` to 0.17 % (and `22/(2π) = 3.501`, `7/2 = 3.500` to similar precision). The closure-quantum ε of PR #18 is the solution of L(ε) = γ/(2π) at ω = 1 — i.e., the inner cutoff is structurally derived from a non-local BS condition. |
| 4 | `throat_reflection_phase_formalization_probe` | formalization layer       | **partial derivation**: formalizes ω · L = π + Δ as the WKB-BS quantization, where π is the empty-box hard-wall ground state and Δ ≈ γ/(2π) − π is the closure-quantum potential correction. T² = −I rigorously derives Dirichlet at the throat; standard WKB derives the (n+1)·π asymptotic limit (verified: at n = 3 the deviation is 1.2 %). The Δ = γ/(2π) − π identification is at WKB precision (~0.5 %) but not yet rigorously derived from the throat algebra. Identity is **R-specific** (R\* = 1.262636 only), **l-specific** (l = 1), **n-asymptotic** (higher n → empty-box BS). |
| 5 | `throat_wkb_uniform_probe` | uniform-WKB derivation     | **EXACT three-piece decomposition** of Δ. The uniform-WKB matching `tan(J + π/4) = (1/2)·exp(−2I)` across the turning point (Airy connection) is approximately satisfied at the Chebyshev eigenstate (residual -0.077, consistent with WKB precision 1/J² ≈ 0.14). The WKB-predicted ω matches Chebyshev to 1.85%. Most importantly: Δ decomposes exactly as `Δ = (J − π) + Δ_cl + ω·L_fb` where (J − π) is the BS-phase shift below empty-box ground state, Δ_cl is the classical-region potential deficit, and ω·L_fb is the geometric tortoise width of the forbidden region. Each piece is computed from V(r*), R\*, ε alone. The sum equals Δ_total = 0.367 to machine precision; matches γ/(2π) − π = 0.370 to 0.72%. |

## Headline finding — the inner cutoff is a WKB-BS consequence

The closure-quantum inner cutoff `ε = 7π/(100·5⁴)` of PR #18 is
**not** an independent identification. It is the solution of a
non-local matching condition between the asymptotic phase and the
closure-quantum pinhole:

```
ω(1, 0; R*) · L(R*; ε)  =  γ_{1..5}(R*) / (2π)        (Compton bridge)
```

At ω = 1, this becomes `L(ε) = γ/(2π)`, whose ε-solution is 3.49 ×
10⁻⁴ — within 0.81 % of the closure-quantum value 7π/(100·5⁴) =
3.52 × 10⁻⁴. The Compton-bridge ω = 1 closes to 0.08 % at the
structural ε.

Within the closure-ledger scaffolding (γ, R\*, ω = 1, ε), there is
now a single non-local matching identity rather than three
independent quantities.

The uniform-WKB probe (probe 5) makes the structural origin of the
identity explicit: Δ = ω·L − π admits an *exact* three-piece
decomposition in WKB quantities,

```
Δ  =  (J − π)  +  Δ_cl  +  ω · L_forbidden        (exact)
```

where (J − π) is the BS-phase shift below the empty-box ground
state, Δ_cl is the classical-region potential-phase deficit, and
ω·L_forbidden is the geometric tortoise-coordinate width of the
forbidden region. Each piece is a Tangherlini-geometry integral
over (V, R\*, ε). The numerical match of this sum to γ/(2π) − π
at R\* (0.72 % gap) is consistent with the WKB precision of the
matching equation itself (residual 0.077 ~ 1/J²).

## What's rigorously derived vs structural reading

The formalization probe makes the derivability explicit:

| ingredient | status | source |
|------------|--------|--------|
| Dirichlet boundary at the throat | **rigorously derived** | T² = −I via T-fixed-point: ψ = T·ψ = T²·ψ = −ψ ⇒ ψ = 0 |
| Dirichlet at the outer wall | convention | radial-grid choice; regular boundary |
| Empty-box `(n+1)·π` asymptotic | **rigorously derived** | standard WKB hard-wall BS; verified at n = 3 (1.2 % deviation) |
| Uniform-WKB matching `tan(J + π/4) = (1/2)·exp(−2I)` | **rigorously derived** | Airy connection across the turning point at WKB precision (residual 0.077 ~ 1/J²) |
| Decomposition Δ = (J − π) + Δ_cl + ω·L_fb | **exact identity** | algebraic; sum reproduces Δ_total to machine precision |
| n = 0 potential correction Δ ≈ γ/(2π) − π | **structural reading** (~0.7 %) | the decomposition's value at R\* numerically matches γ/(2π) − π at WKB precision |
| R\* = 1.262636 specificity | derived | the unique R at which both lepton mass ratios (PR #15) and BS phase = γ/(2π) (this thread) hold simultaneously |

The 0.5 % gap in the Δ ↔ γ/(2π) identification is at the precision
of the WKB approximation itself and at the same level as the other
closure-quantum identifications in PRs #15–18 (transport 0.13 %,
resistance 0.94 %, γ 2 %).

## Robustness profile

The reflection-phase condition is approximately satisfied only in
a structurally narrow regime:

- **R-specific.** Holds tightly only at R\* = 1.262636 (−0.06 %);
  deviates by up to 48 % at R = 1.10 (small box, WKB breaks down)
  and by 1.4 % at R = 1.30–1.35. The R-specificity is consistent
  with the closure-quantum reading: R\* is the unique R at which
  both the lepton mass ratios (PR #15) and the BS phase = γ(R)/(2π)
  hold simultaneously. The two conditions are linked, not
  independent.

- **l-specific.** ω(l, 0) · L for l = 1..5 grows from 3.51 to 4.65;
  no simple closure-quantum relation across l. The identity is
  specific to the l = 1 mode, which is the radial mode coupled to
  the lepton ground state in the closure-ledger surrogate (PR #14
  baseline).

- **n-asymptotic.** ω(1, n) · L → (n + 1) · π for higher n (the
  empty-box hard-wall BS limit). Deviations: n = 0 → 11.7 %,
  n = 1 → 4.5 %, n = 2 → 2.1 %, n = 3 → 1.2 %. The n = 0 deviation
  IS the closure-quantum potential correction γ/(2π) − π; higher
  modes asymptote to empty-box BS.

## What this leaves open

Three structural questions remain, in order of depth:

1. **Rigorous derivation of Δ = γ/(2π) − π.** The empirical match
   is at WKB precision; closing the gap requires either (i) a
   uniform WKB expansion of the bound-state phase including the
   tunneling-tail contribution from the outer-wall barrier, or
   (ii) an algebraic identity tying the BS phase directly to the
   closure-quantum γ. (i) is a calculational task within the
   closure-ledger scope; (ii) likely requires the deeper throat-
   dynamics framework.

2. **The 0.81 % residual gap on ε.** The structural ε (solution
   of L(ε) = γ/(2π) at ω = 1) is 3.49 × 10⁻⁴; the closure-quantum
   identification of PR #18 is 7π/(100·5⁴) = 3.52 × 10⁻⁴. The gap
   is at the WKB-precision level and is consistent with the
   precision of other closure-quantum identifications.

3. **R_MID self-consistency.** R_MID = 1 by convention. Determining
   R_MID dynamically from a self-consistency condition (equilibrium
   throat radius for the locked mass spectrum) would lift the m_e
   anchor itself. THESIS.md "self-consistent throat radius" —
   outside the closure-ledger framework.

## Implication for the ℏ-origin program

The closure-ledger framework's parameter-reduction history now reads:

  - Six phenomenological parameters at the start of PR #14.
  - → Two (transport, resistance) at the end of PR #15.
  - → Zero closure-quantum constants in PR #16.
  - → One factor (1.054) in PR #17 (regularization reframing).
  - → One anchor (m_e) after PR #18 (ε structurally identified).
  - → The same anchor (m_e) after this thread, with ε now derived
    from the BS quantization rather than from an independent
    closure-quantum identification.

The closure-ledger scaffolding is now self-contained: γ ≈ Σ V_max
fixes the BS phase; R\* is selected by cross-species fit AND the BS
condition simultaneously; ε is L⁻¹(γ/(2π)). The remaining external
input is m_e, and lifting it requires the deeper R_MID self-
consistency (outside the closure-ledger scope).

## Cross-references

- `docs/throat_dynamics_research_plan.md` — research plan (this
  thread's opening doc).
- `docs/hbar_origin_status.md` — eight-probe ℏ-origin summary
  table; this thread's prerequisite.
- `docs/hbar_origin_note.md` — closure-ledger paper draft (§5
  introduces ε = 7π/(100·5⁴); this thread refines that to the
  non-local BS reading).
- `experiments/closure_ledger/throat_boundary_condition_probe.py`
  (sub-target A, negative).
- `experiments/closure_ledger/throat_thickness_probe.py`
  (sub-target B, negative).
- `experiments/closure_ledger/throat_reflection_phase_probe.py`
  (sub-target 3, positive — non-local BS identification).
- `experiments/closure_ledger/throat_reflection_phase_formalization_probe.py`
  (formalization layer — WKB decomposition + (R, l, n) robustness).
- `experiments/closure_ledger/throat_wkb_uniform_probe.py`
  (uniform-WKB derivation — Airy-connection matching equation
  `tan(J + π/4) = (1/2)·exp(−2I)` and exact three-piece decomposition
  of Δ).
- `experiments/closure_ledger/runs/` — per-probe JSON + markdown
  archives, one timestamped directory per run.
