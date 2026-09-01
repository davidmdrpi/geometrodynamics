# Can any classical field already in BAM keep the throat open?

**5/5 checks pass — and the result is negative.**

Every prediction was frozen in `docs/source_audit_prereg.md` before this module existed.

> ## The verdict

> **NO -- the current classical BAM field content cannot support a traversable particle throat**

> The neck demands `8πG₅T_kk = -393.343998` at `R = 1, a = 0.3`. Nine candidates give exactly zero or a manifest square; the tenth — a nonminimally coupled scalar, the only one whose sign is not fixed a priori — closes at the defect core in **every** dimension.

---

## S0 — the requirement is a flare-out theorem

The radial null congruence has a three-dimensional screen in `D = 5`, so `θ = 3f′/f` and

```
dθ/dλ = −θ²/3 − σ² − R_ab k^a k^b
```

At the neck `θ = 0.0e+00` and `σ = 0` by spherical symmetry, while `dθ/dλ = 393.343998` — the flare-out — against the expected `3/b² = 393.343998`. The Raychaudhuri identity holds to `1.1e-13`, so

```
8πG₅T_kk = R_kk = -393.343998 = −3/b²
```

> Derived from theta and Raychaudhuri alone, with no reference to p_s, to N = 1, or to reflection symmetry. Smooth radial flare-out plus Einstein gravity forces T_kk < 0, whatever the lapse and whatever the matter.

> **NOT a new theorem: this is the 5D specialisation of the standard traversable-wormhole flare-out requirement (Morris-Thorne 1988) in the local form of Hochberg-Visser, who proved NEC violation at a throat defined as a marginally anti-trapped surface. Recovering it here VALIDATES the finite-mouth construction; it is not a discovery, and the repository has too few external anchors to waste one by miscrediting it.**

## S1 — every candidate, from its actual stress tensor

| | candidate | min `T_kk` | sign | reason |
|--|--|--|--|--|
| `C1` | minimally coupled scalar psi | `+1.35e-03` | **NON-NEGATIVE** | T_kk = (k.grad psi)^2; all g_ab terms drop under the null contraction |
| `C2` | complex GL throat-order field q | `+3.42e-02` | **NON-NEGATIVE** | T_kk = kappa |k.Dq|^2 with the repository's positive gradient term; this is the field introduced to REPRESENT the throat, and it cannot pay for it |
| `C3` | GL potential V(q), incl. symmetry-breaking | `+4.13e-14` | **EXACTLY ZERO** | contracting the full T_ab for V = -9, -1, 0, +4 leaves (k.grad q)^2 to 4.1e-14: the potential drops identically, so V < 0 regions are irrelevant to the NEC |
| `C4` | Maxwell / KK gauge field | `+5.29e+01` | **NON-NEGATIVE** | T_kk = |F.k|^2 with V.k = 0 to 0.0e+00 over 200 random field strengths, so V is spacelike or null |
| `C5` | cosmological constant | `+0.00e+00` | **EXACTLY ZERO** | T_ab ∝ g_ab, and g_ab k^a k^b = 0 |
| `C6` | perfect fluid obeying the NEC | `+8.88e-04` | **NON-NEGATIVE** | T_kk = (rho+p)(u.k)^2, zero only for rho + p = 0 |
| `C7` | classical 5D gravitational waves | `+0.00e+00` | **EXACTLY ZERO** | R_ab = 0 in vacuum, so R_kk = 0 identically; the Isaacson effective energy is positive |
| `C8` | non-orientable identification / wrap sign | `+0.00e+00` | **NO CONTRIBUTION** | changes global boundary data, not the local R_kk. hopf/spinor.py is SU(2) holonomy transport with no stress tensor, so U_spin is a transport object, not a source |
| `C9` | projected bulk Weyl stress (#167/#168) | `+0.00e+00` | **NO CONTRIBUTION TO THE 5D EQUATION** | there T^(5)_ab = 0 exactly; the projection is an effective 4D BRANE source and cannot appear on the right-hand side of the 5D Einstein equation, which is what the neck needs |

Candidates with negative null stress: **none**.

> **C2 is the field the repository introduced AS the throat-order degree of freedom. With the ordinary positive gradient term its null stress is a modulus squared, so it cannot support the object it was introduced to represent. That is stronger than saying its constants were never derived.**

## S2 — the nonminimal loophole

At a node `q = 0` the prefactor `1 − 8πG₅ξq²` is exactly `1` and `d²(q²)/dλ² = 2(dq/dλ)²`, so the sign is `sign(1 − 2ξ)`.

| coupling | `ξ` | `1 − 2ξ` | null stress |
|--|--|--|--|
| minimal xi = 0 | `0.0000` | `+1.0000` | positive |
| conformal 4D xi = 1/6 | `0.1667` | `+0.6667` | positive |
| conformal 5D xi = 3/16 | `0.1875` | `+0.6250` | positive |
| conformal 6D xi = 1/5 | `0.2000` | `+0.6000` | positive |
| threshold xi = 1/2 | `0.5000` | `+0.0000` | zero |
| beyond threshold xi = 0.8 | `0.8000` | `-0.6000` | negative |

And the closed form `1 - 2 xi_c = D/(2(D-1)) > 0 for every D >= 2`:

| `D` | `ξ_c` | `1 − 2ξ_c` |
|--|--|--|
| 3 | `0.125000` | `0.750000` |
| 4 | `0.166667` | `0.666667` |
| 5 | `0.187500` | `0.625000` |
| 6 | `0.200000` | `0.600000` |
| 10 | `0.222222` | `0.555556` |
| 26 | `0.240000` | `0.520000` |

> At q = 0 the prefactor 1 - 8 pi G xi q^2 is exactly 1 and d^2(q^2)/dlam^2 = 2 (dq/dlam)^2, because the q q'' term carries a factor q. BAM places the defect core AT the node, so the sign there is sign(1 - 2 xi).

> **A flip needs xi > 1/2. Conformal coupling gives D/(2(D-1)), which is 5/8 in five dimensions and never reaches zero in ANY dimension. So even a conformally improved throat-order field remains NEC-positive at a smooth q = 0 core.**

## S3 — the dynamic escape

If `T_kk ≥ 0` every term on the right of Raychaudhuri is non-positive, so `θ` is non-increasing and a ray that enters converging cannot flare back out. Each ray is integrated past its own analytic caustic bound, `lam_caustic <= (D-2)/|theta_0| = 3/|theta_0|`:

| profile | `θ₀` | bound | measured caustic | turned `−→+`? |
|--|--|--|--|--|
| vacuum R_kk = 0 | `-0.05` | `60.00` | `60.19` | **no** |
| vacuum R_kk = 0 | `-0.50` | `6.00` | `6.02` | **no** |
| vacuum R_kk = 0 | `-2.00` | `1.50` | `1.50` | **no** |
| constant R_kk = 0.5 | `-0.05` | `60.00` | `3.44` | **no** |
| constant R_kk = 0.5 | `-0.50` | `6.00` | `2.66` | **no** |
| constant R_kk = 0.5 | `-2.00` | `1.50` | `1.32` | **no** |
| gaussian bump | `-0.05` | `60.00` | `2.83` | **no** |
| gaussian bump | `-0.50` | `6.00` | `2.29` | **no** |
| gaussian bump | `-2.00` | `1.50` | `1.34` | **no** |
| oscillating, non-negative | `-0.05` | `60.00` | `2.12` | **no** |
| oscillating, non-negative | `-0.50` | `6.00` | `1.78` | **no** |
| oscillating, non-negative | `-2.00` | `1.50` | `1.12` | **no** |
| decaying tail | `-0.05` | `60.00` | `1.80` | **no** |
| decaying tail | `-0.50` | `6.00` | `1.50` | **no** |
| decaying tail | `-2.00` | `1.50` | `1.00` | **no** |

Turning points found: **none**.

> An earlier draft integrated to a fixed lam <= 6 and reported that not every ray focused -- which was an artefact of the window, not a failure of the theorem: at theta0 = -0.05 the vacuum caustic sits at lam = 3/|theta0| = 60. Each ray is now integrated past its own analytic bound. The no-turning-point result never depended on the window; the focusing sub-claim did.

> **Not just static BAM: ANY smooth two-way traversable BAM throat in classical Einstein gravity requires null-convergence violation somewhere. A numerical-relativity round looking for a vacuum dynamic rescue is looking for something the null equations forbid.**

## S4 — the ledger

| claim | verdict | evidence |
|--|--|--|
| the neck NEC price is an artefact of N = 1 or of symmetry | **NO -- IT IS A FLARE-OUT THEOREM** | derived from theta and Raychaudhuri alone, residual 1.1e-13, with no reference to p_s |
| that theorem is new | **NO -- IT IS MORRIS-THORNE / HOCHBERG-VISSER IN 5D** | recovering a known theorem validates the finite-mouth construction; it is a second external anchor, not a discovery |
| the GL throat-order field q can support the throat | **NO** | T_kk = kappa |k.Dq|^2 >= 0 with the repository's own positive gradient term; the field introduced to represent the throat cannot pay for it |
| the GL potential's symmetry-breaking region helps | **NO -- IT DROPS IDENTICALLY** | every g_ab term vanishes under null contraction, so V < 0 is irrelevant to the NEC |
| a conformally improved scalar rescues it at the core | **NO -- 1 - 2 xi_c = D/(2(D-1)) > 0 IN EVERY DIMENSION** | 5/8 in five dimensions; a flip needs xi > 1/2, and BAM puts the defect core exactly at the q = 0 node |
| the #167/#168 Weyl mechanism can keep this throat open | **NO -- WRONG EQUATION** | there T^(5)_ab = 0; a projected Weyl tensor is an effective 4D brane source and cannot appear in the 5D Einstein equation the neck lives in |
| nonzero K_ij (gravitational waves) can rescue it | **NO** | 15 Raychaudhuri integrations with non-negative R_kk and shear produce zero turning points; theta is non-increasing and focuses |

**What this forces.** If every existing BAM action gives T_kk >= 0, the honest headline is a negative result about the CURRENT field content, not about wormholes in general.

**The remaining branches.**

- **1 accept the horizon** — take the Tangherlini branch N(0) = 0 as the particle, and abandon MTY traversability
- **2 add a ghost** — a wrong-sign field buys T_kk < 0 and brings an obvious stability problem
- **3 leave Einstein gravity** — higher curvature -- in D = 5 the natural term is Gauss-Bonnet, which is exactly where the theorem's premise fails; Lovelock wormholes can satisfy the matter NEC
- **4 quantum stress** — Casimir-type support, meaning the particle geometry is no longer derived from classical GR
- **5 reinterpret** — particle exchange does not require a traversable throat at all

**What the audit does not do.** It does not choose among those five, and it refutes none of them. Each is a premise of the theorem rather than a consequence.
