# Pre-registration: can any existing BAM field keep the throat open?

**Status: frozen before the audit module exists.** Committed alone, ahead of
the code, for the same reason as `docs/finite_mouth_prereg.md`.

---

## The binary question

> **Does any classical field already present in BAM have a stress tensor
> capable of supplying `T_kk < 0` at the finite mouth's neck, without changing
> the fundamental theory?**

The number every candidate must confront is not a fit:

```
8πG₅ T_ab k^a k^b |₀  =  −3/b²  =  −3/(R² sin⁴a)
```

At `R = 1, a = 0.3` that is `−393.343997824`.

---

## The requirement is a flare-out theorem, not a feature of `N = 1`

Verified symbolically before freezing, in the ultrastatic finite-mouth metric
with the radial null vector `k = (1,1,0,0,0)`:

| step | result |
|--|--|
| `k` is null | exact |
| `R_ab k^a k^b` | `−3b²/(s²+b²)²`, giving `−3/b²` at the neck |
| `θ = 3f′/f` at the neck | `0` |
| `dθ/dλ` at the neck | `+3/b²` |
| Raychaudhuri `dθ/dλ = −θ²/3 − σ² − R_kk` | residual **exactly 0**, with `σ = 0` |

The screen is three-dimensional in `D = 5`, so the `θ²` coefficient is `1/3`.
The `Λ` term drops under null contraction, so Einstein gives
`8πG₅T_kk = R_kk = −3/b²`. This recovers P4 without touching the component
`p_s` equation, and shows it is not an artefact of `N = 1`, of reflection
symmetry, or of the stress decomposition:

```
smooth radial flare-out  +  Einstein gravity   ⟹   T_kk < 0
```

**Attribution.** This is *not* a new result. It is the 5D specialisation of the
standard traversable-wormhole flare-out requirement (Morris–Thorne 1988) in the
local form proved by Hochberg–Visser for a throat defined as a marginally
anti-trapped surface, without symmetry or asymptotic-flatness assumptions.
Reproducing a known theorem is a **validation of the construction**, and it is
recorded as such rather than as a discovery. The audit found the repository has
essentially one external anchor for 324 claim rows; this is a second.

---

## Frozen per-candidate predictions

Each is a prediction about what the module will find when the stress tensor is
built from the **actual action**, not from prose.

| # | candidate | predicted `T_kk` | reason |
|--|--|--|--|
| C1 | minimally coupled scalar `ψ` | `= (k·∇ψ)² ≥ 0` | the `g_ab` terms, `V` included, drop under null contraction |
| C2 | complex GL throat-order `q` | `= κ\|k·∇q\|² ≥ 0` | same, with the repository's positive gradient term |
| C3 | GL potential `V(q)`, quartic and symmetry-breaking | `= 0` | pure `g_ab`; irrelevant to the NEC even where `V < 0` |
| C4 | Maxwell / KK gauge field | `= \|F_{ab}k^b\|² ≥ 0` | `V_a = F_{ab}k^b` satisfies `V·k = 0`, so `V` is spacelike or null |
| C5 | cosmological constant | `= 0` | `T_ab ∝ g_ab` |
| C6 | perfect fluid | `= (ρ+p)(u·k)² ≥ 0` under the NEC | vanishes only for `ρ + p = 0` |
| C7 | classical 5D gravitational waves | `R_kk = 0` exactly in vacuum | Isaacson effective energy is positive |
| C8 | non-orientable identification, wrap sign | no contribution | global boundary data, not local `R_kk` |
| C9 | projected bulk Weyl stress (#167/#168) | no contribution **to the 5D equation** | `T^{(5)}_ab = 0`; the projection is an effective *4D brane* source |
| C10 | conformally improved scalar | `∝ (1−2ξ)(dq/dλ)²` at a node | the only candidate whose sign is not fixed a priori |

**C10 is the only real loophole, and it is predicted to close.** At a node
`q = 0` the improvement term gives `d²q²/dλ² = 2(dq/dλ)²` because `q q″`
vanishes, and the prefactor `1 − 8πG₅ξq²` is exactly `1`. So

```
R_kk ∝ (1 − 2ξ)(dq/dλ)²
```

and a sign flip needs `ξ > 1/2`. The conformal value in `D` dimensions is
`ξ_c = (D−2)/(4(D−1))`, giving

```
1 − 2ξ_c = D / (2(D−1))  >  0     for every D
```

— `2/3` in 4D, `5/8` in 5D, `3/5` in 6D. Conformal coupling **never** flips the
sign at a node, in any dimension. And BAM places the defect core exactly at
`q = 0`, so this is the relevant point.

**Frozen prediction: every one of C1–C10 gives `T_kk ≥ 0`, hence the verdict is
NO.** If that is what the module finds, the honest headline is a negative
result: *the current classical BAM field content cannot support a traversable
particle throat.*

---

## The dynamic escape, and why it is predicted to fail too

The finite-mouth construction used a time-symmetric slice. One might hope
nonzero `K_ij` — the genuine gravitational-wave degrees of freedom — supplies
the support without exotic matter. Raychaudhuri forbids it: if `T_kk ≥ 0` then

```
dθ/dλ = −θ²/3 − σ² − 8πG₅T_kk  ≤  0
```

so `θ` is non-increasing along the congruence and can never turn from negative
to positive. A ray entering the throat converges and cannot flare back out.

**Frozen prediction:** integrating Raychaudhuri numerically with any `R_kk ≥ 0`
profile and any initial `θ < 0` will produce **no** turning point, and `θ` will
run to `−∞` in finite affine parameter (the focusing theorem). This makes the
statement stronger than the static one:

> Not just static BAM, but **any** smooth two-way traversable BAM throat in
> classical Einstein gravity requires null-convergence violation somewhere.

---

## Falsifiers

1. Any candidate whose correctly derived `T_kk` is negative somewhere on a
   physically admissible configuration falsifies its row.
2. A Raychaudhuri integration with `R_kk ≥ 0` that produces `θ: − → +`
   falsifies the dynamic no-go — and would mean a sign error in this document.
3. A nonminimal coupling with `ξ > 1/2` arising from an action already in the
   repository, rather than being inserted to order, falsifies C10's closure.
4. Finding a matter action in `geometrodynamics/` not covered by C1–C10 makes
   the audit incomplete, and it must then be extended rather than reported.

---

## What is deliberately *not* claimed

The theorem's premises are **classical Einstein gravity in 5D** with the listed
matter. Each premise is a real escape hatch, and none is refuted here:

- **higher-curvature gravity** — in `D = 5` the natural term is Gauss–Bonnet,
  which is exactly where "Einstein gravity" fails as a premise; Lovelock
  wormholes can satisfy the matter NEC;
- **ghost / wrong-sign fields** — buy `T_kk < 0` at the cost of stability;
- **quantum stress** (Casimir-type) — then the particle geometry is no longer
  derived from purely classical GR;
- **classical Dirac field** — the one classical matter field that can violate
  the NEC. BAM does **not** currently contain one: `hopf/spinor.py` is SU(2)
  holonomy transport with no stress tensor, so `U_spin` is a transport object,
  not a matter source. Recorded because its absence is informative;
- **reinterpretation** — particle exchange need not require a traversable
  throat at all.

If the verdict is NO, those five are the remaining branches, and choosing among
them is a physics decision this audit does not make.
