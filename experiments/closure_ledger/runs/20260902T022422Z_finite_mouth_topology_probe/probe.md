# Finite-mouth topology probe — 4/5

Pre-registration: `docs/finite_mouth_topology_prereg.md @ d9d85bc`. Verdict: **`FINITE_MOUTH_ADMITS_BUT_DOES_NOT_SELECT_THE_BAM_LIFT`**

## H1 — the gluing is free; antipodal gluing is orientable

**FAILS** (*topological theorem*)

> Darmois matching is O(4)-invariant, so the seam maps are free and only the loop monodromy m matters. det(-I_4) = +1: the antipodally glued two-mouth handle is ORIENTABLE in the bulk; what it twists is the brane S^2-handle and the brane's normal line bundle. The (-1)^l is the computed action of -I on degree-l harmonics; a reflection does not act as a scalar at all.

## H2 — the unique free involution

**HOLDS** (*topological theorem, conditional on P_B = -P_A (chosen)*)

> -I_5 is the only fixed-point-free involution in O(5), and through the tube it continues uniquely as (s, Omega) -> (-s, -Omega). That map is free, an isometry for both lapses, and reverses the bulk orientation while preserving the brane slice's. The non-orientable object is the QUOTIENT M/iota = RP^4 # RP^4, not the two-mouth handle; the RP^2 the repository names is its brane neck. And the Pin types do not match: RP^4 # RP^4 is Pin+ only, RP^2 is Pin- only.

## H3 — the scalar sector and the oracle

**HOLDS** (*analytic identity, conditional on H2's choices and on eta*)

> On M/iota a scalar obeys psi_l(-s) = eta (-1)^l psi_l(s), so the neck is Neumann or Dirichlet by the parity of l. PR #129's even/odd structure is the eta = +1 sector, obtained at the finite ultrastatic neck with no horizon. The half-tube admittance is the (1, ±1) sector of PR #277's two-mouth oracle, reproduced by an independent second-order solve. Replacing the antipode by the identity fixes the whole neck and erases the grading.

## H4 — no horizon limit; the time-reversal datum

**HOLDS** (*analytic identity*)

> f and both lapses are even in s, so the static operator commutes with s -> -s for the ultrastatic and Tangherlini branches alike: the sector labels are lapse-independent and there is no limit to take. PR #129's (U,V) -> (-U,-V) is the PT-type extension of iota (5D det +1, time orientation reversed); the P-type extension (U,V) -> (V,U) exists too (5D det -1, time orientation kept). Both restrict to iota on the t = 0 slice. Which one is physical is a further choice, recorded as the datum T.

## H5 — what J is

**HOLDS** (*analytic identity (S2, S3); representation choice (the K)*)

> The lift of the neck reflection is ±e_s in Pin^±(4): a real matrix that exchanges chiralities. The lift of -I_4 is ± the volume element: the chirality sign. J = i sigma_y K is left multiplication by -j, a ROTATION with det +1; it is antilinear only with respect to the Hopf complex structure L_i and is the linear SU(2) matrix i sigma_y with respect to R_i. No Pin lift is antilinear, so J is not the lift of any gluing map. Its K is the reversal of the Hopf U(1) - charge conjugation, a bundle datum - and its i sigma_y is a spin rotation by pi, neither of which the mouth gluing supplies.

### Gluing classes

| gluing | det m | det m_3 | eps | w_1 on loop |
|--|--|--|--|--|
| identity | `+1` | `+1` | `+1` | `+1` |
| antipodal | `+1` | `-1` | `-1` | `+1` |
| reflection (normal) | `-1` | `+1` | `-1` | `-1` |
| reflection (brane) | `-1` | `-1` | `+1` | `-1` |

### Half-tube admittance against the pre-registered targets

| ℓ | condition | target | solve | rel. err |
|--|--|--|--|--|
| 0 | Neumann | `0.000000000` | `0.000000000` | `0.0e+00` |
| 0 | Dirichlet | `0.157587622` | `0.157587622` | `4.1e-10` |
| 1 | Neumann | `1.797266559` | `1.797266798` | `1.3e-07` |
| 1 | Dirichlet | `1.804461992` | `1.804462240` | `1.4e-07` |
| 2 | Neumann | `3.524607516` | `3.524608499` | `2.8e-07` |
| 2 | Dirichlet | `3.524854051` | `3.524855035` | `2.8e-07` |
| 3 | Neumann | `5.248595411` | `5.248597926` | `4.8e-07` |
| 3 | Dirichlet | `5.248602920` | `5.248605435` | `4.8e-07` |

Convergence ratios: Neumann ℓ=2: `4.00`, `4.00`, `4.00`; Dirichlet ℓ=1: `4.00`, `4.00`, `4.00`

## Audit table

| Question | Result | Evidence class |
|--|--|--|
| Is the physical mouth non-orientable? | The two-mouth handle with antipodal gluing: NO (bulk w_1 = 0; brane S^2-handle and normal bundle twisted). The quotient M/iota: YES (RP^4 # RP^4). | topological theorem |
| What object carries w_1 != 0? | On the handle: the brane S^2 x~ S^1 and the brane's normal line bundle (eps = -1). On the quotient: the bulk itself and its brane neck RP^2. Never the bulk mouth S^3. | topological theorem |
| Is J = i sigma_y K a valid lift? | NO, of any gluing map: every Pin lift is complex-linear; J is antilinear for the Hopf structure. J = (spin rotation by pi) o (U(1) reversal). | analytic identity |
| Is it unique? | Not applicable as a lift. Among fibre-reversing Hopf isometries it is one point of the component (-j e^{i alpha}, g_R); the lifts of iota are four (±e_s in Pin^±) and of -I_4 two (± volume). | analytic identity |
| Is PR #129's antipodal BC derived? | CONDITIONALLY: it is the eta = +1 sector of the unique free involution, given P_B = -P_A (chosen) and the quotient rather than the double cover (chosen). Neither choice is forced by the geometry. | conditional on an unproved physical identification |
| Does its horizon limit reproduce even/odd BCs? | There is no limit: the sector labels are lapse-independent and hold at the finite ultrastatic neck. Reproduced numerically against the PR #277 oracle to 4.8e-07, second order. | numerically converged |
| Which inputs remain postulates? | P_B = -P_A; quotient vs cover; eta; the extension in time (P or PT); the Hopf complex structure on the state space; charge conjugation C; any Pin type and sign. The geometry fixes none of them. | definition / chosen |

## Dependency ledger

* `antipodal_BC` = BC( P_B=-P_A [chosen], quotient-not-cover [chosen], eta [chosen], iota [derived given P_B=-P_A], (-1)^l [identity], a, R [geometry], P-or-PT [chosen] )
* `J` = J( Hopf complex structure on C^2 [chosen], C = U(1) reversal [chosen], -j in SU(2)_L [gauge/convention] )
* `spinor transport of iota` = ±e_s ( Pin type [chosen], sign [chosen], neck frame [gauge] )
* `handle topology` = mapping torus of m ( m in O(4) [chosen]; (det m_3, eps) [derived from m] )

## Refinement of the trichotomy

The trichotomy is well posed but its middle option needs one sharpening: the geometry does not merely fail to select the BAM lift, it contradicts two of the words attached to it. (i) The antipodally glued two-mouth handle is orientable in the bulk; non-orientability lives on the quotient, and there the Pin type (Pin+ on RP^4 # RP^4) is not the Pin- the repository assigns to the RP^2 mouth. (ii) J is a rotation, not a reflection, and its antilinearity is a choice of complex structure. Neither the non-orientability nor the K can be read off the finite mouth.
