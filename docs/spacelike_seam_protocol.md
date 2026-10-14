# Deflected spacelike seam: validation protocol

Written after analytic exploration of the deflection and stress scalings, before
implementation of the producer. This is a validation protocol, not a blind
preregistration. Preserve the packet protocol, solver and archives of #335.

Use the same canonical spherical scalar and b=c=kappa5=1. Work with the
O(A²) embeddings U_in=x(V), U_out=-1/V, V_out=y(V), obtained by matching
both components of the first fundamental form. Compute the leading nonzero
orders of the second fundamental form, with future-directed timelike unit
normals, and the spacelike Israel source. The normal has norm -1; no timelike
shell energy-density interpretation is assigned to its spatial eigenvalues.

1. Verify radius and temporal-metric matching analytically, including equal
   initial embeddings in the unperturbed past. Check that the seam is spacelike
   wherever the positive packet-source integral is nonzero.
2. Derive both extrinsic-curvature eigenvalues. Validate the longitudinal
   curvature using the embedding acceleration and normal, not only a difference
   formula. Validate angular curvature with n(r)/r including the leading bulk
   radius response. Label uncomputed subleading terms explicitly.
3. Obtain scalar transverse derivatives in two ways: by integrating the horizon
   characteristic equation and by one-sided differentiation of the #335 2D
   fields. When calculating the source, enforce common scalar trace through
   O(A³); do not misinterpret uncorrected trace displacement as a jump.
4. Test leading Codazzi and Gauss cancellations symbolically. Full next-order
   conservation needs the next metric/scalar response and is not claimed.
5. Run N=64,128,256,512 packet controls and amplitudes A=.005,.01,.02. Separate
   exact perturbative scaling from grid convergence. At V=4 require relative
   differences between the two scalar-source estimates to decrease with
   refinement and be <1e-2 at the finest grid. Check analytic metric-matching
   residuals <1e-12 and positive induced h_VV after the packet's peak.
6. Report both unit-normal-frame surface stresses and coordinate-densitized
   stresses. A 1/A coefficient is not by itself a divergent null distribution.
   For a specified proper-normal smoothing width L_perp, record the necessary
   peak-stress budget |S_angle|/(2 L_perp), and the scaling needed to keep that
   budget bounded. This is a bound/diagnostic, not construction of a smooth GR
   solution or proof that every finite-thickness model fails.

Admission requires a common scalar trace and normal derivative if no surface
scalar coupling is present, and a supplied material law for the required
surface stress. A resolved nonzero normal-derivative jump rejects source-free
scalar continuation for this deflection. Do not call the induced-metric match
alone a dynamical completion. The onset of nonzero spacelike character has a
degenerate induced metric; do not extend unit-normal formulas into the null
segment by dividing by zero or threshold-clipping the geometry.
