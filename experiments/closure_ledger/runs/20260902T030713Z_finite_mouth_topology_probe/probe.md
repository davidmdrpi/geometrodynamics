# Finite-mouth topology probe — 6/6

Pre-registration: `docs/finite_mouth_topology_prereg.md @ d9d85bc; docs/mouth_pin_holonomy_prereg.md @ 7f46fff`. Verdict: **`FINITE_MOUTH_ADMITS_BUT_DOES_NOT_SELECT_THE_BAM_LIFT`**

## H1 — the gluing is free; antipodal gluing is orientable

**HOLDS** (*topological theorem*)

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

## H6 — the mouth Pin holonomy (after review)

**HOLDS** (*topological theorem (P1-P3), analytic identity (P4, P6, P7), numerically converged (P5), chosen (P8)*)

> The deck generator of the neck RP^2 reverses both normals (d_s and the brane normal) and reflects the tangent plane. The ambient Pin+ structure restricts to the neck (w_2(TM|_N) = 0) and induces the intrinsic Pin- one; the twisted tangent generators e_t e_s e_n generate Cl^-(2) = H. The holonomy lift ±e_s e_n e_t2 squares to -1 in Pin+ -- the square root of the spin holonomy -1 of a 2 pi tangent rotation around the neck's great circle, computed by unwrapping the angle to 2 pi -- and to +1 in Pin-, which is excluded. In the intrinsic module the holonomy is left multiplication by a unit quaternion of the (i, j) plane: the fibre-reversing component of Pin(2) in SU(2) that contains sigma = L_{-j}, up to a Spin(2) conjugation and a sign. J^2 = -1 is geometric; -j and the sign are not selected.

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
| Is J = i sigma_y K a valid lift? | YES: the Spin(4) lift (-j, 1) of the rotation sigma = L_{-j}; and the induced Pin- holonomy of the RP^2 mouth lies in the same fibre-reversing component of Pin(2) in SU(2), with the same square -1. Not a lift of iota (±e_s) or of -I_4 (± volume). [R2 corrects the first draft] | analytic identity |
| Is it unique? | NO: determined up to Spin(2) = U(1) conjugation (the tangent direction of the generator path) and a sign (the Pin structure). Outcome B. | analytic identity |
| Where does J^2 = -1 come from? | From the ambient Pin+ structure and the two twisted normal lines: (e_s e_n)^2 = -1, and the holonomy lift squares to the spin holonomy -1 of the 2 pi tangent rotation around the neck's great circle. Forced by the quotient, not inserted. [R1: the Pin types are not mismatched] | analytic identity; numerically converged |
| Is PR #129's antipodal BC derived? | CONDITIONALLY: it is the eta = +1 sector of the unique free involution, given P_B = -P_A (chosen) and the quotient rather than the double cover (chosen). Neither choice is forced by the geometry. | conditional on an unproved physical identification |
| Does its horizon limit reproduce even/odd BCs? | There is no limit: the sector labels are lapse-independent and hold at the finite ultrastatic neck. Reproduced numerically against the PR #277 oracle to 4.8e-07, second order. | numerically converged |
| Which inputs remain postulates? | P_B = -P_A; quotient vs cover; eta; the extension in time (P or PT); the identification of the Hopf fibre phase with the mouth's Spin(2) rotation; the U(1) direction and the sign of the mouth holonomy. The Pin type is no longer a postulate (Pin+ ambient, Pin- intrinsic, both forced). | definition / chosen |

## Dependency ledger

* `antipodal_BC` = BC( P_B=-P_A [chosen], quotient-not-cover [chosen], eta [chosen], iota [derived given P_B=-P_A], (-1)^l [identity], a, R [geometry], P-or-PT [chosen] )
* `J` = J( Hopf complex structure on C^2 [chosen], C = U(1) reversal [chosen], -j in SU(2)_L [gauge/convention] )
* `mouth holonomy` = H~( iota [derived given P_B=-P_A], nu = lambda+lambda [derived], Pin+ [derived: w_1^2 != 0 and P5], t2 [gauge: Spin(2)], sign [chosen: Pin structure] )
* `H~ = J` = ( H~ [above], Cl^-(2) = H [canonical], H = Hopf C^2 [chosen complex structure], t2 = -j [chosen within U(1)], sign [chosen] )
* `spinor transport of iota` = ±e_s ( Pin type [chosen], sign [chosen], neck frame [gauge] )
* `handle topology` = mapping torus of m ( m in O(4) [chosen]; (det m_3, eps) [derived from m] )

## Refinement of the trichotomy

The trichotomy is well posed; its middle option is sharpened in both directions after review. Not selected: the seam gluing of the two-mouth handle (antipodal gluing is orientable in the bulk), the quotient over the cover, eta, the U(1) direction and the sign of the mouth holonomy. Supplied, given the quotient: the neck RP^2 has normal bundle lambda + lambda, the ambient Pin+ induces the intrinsic Pin- (not a mismatch: R1), and the deck-generator holonomy squares to -1 -- the spin holonomy of a 2 pi rotation on the round neck -- inside the same component of Pin(2) in SU(2) that contains sigma = L_{-j} (which IS the spin lift of a rotation: R2). J^2 = -1 is geometric; -j and the sign are not.
