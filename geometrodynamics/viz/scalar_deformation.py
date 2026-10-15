"""Display a solved scalar as a deformation of the handled sphere.

The intrinsic wave operator is untouched. On the sphere this is the v41
radial graph r = 1 + excursion*tanh(gain*u/u_ref). In a handle it is a
continuous ambient-vector displacement, NOT a globally chosen surface normal
(which a Klein surface cannot have). The vector blends rim radial directions
into transverse breathing of each drawn tube ring. It can have tangential
components. This is a field visualization, not an evolved metric.
"""

from dataclasses import dataclass

import numpy as np

from .klein_handles import CutMesh, HandledSphere


@dataclass
class ScalarDeformation:
    surface: HandledSphere
    cut: CutMesh
    excursion: float = 0.22

    def __post_init__(self):
        if not np.isfinite(self.excursion) or not 0 <= self.excursion < 0.5:
            raise ValueError("excursion must be finite and in [0, 0.5)")
        surf = self.surface
        vectors = self.excursion * surf.points.copy()
        for ids in surf.tube_ids:
            rings = surf.points[ids]
            transverse = rings - rings.mean(axis=1, keepdims=True)
            transverse /= np.linalg.norm(transverse, axis=2, keepdims=True)
            s = surf.s[ids[:, 0]]
            q = np.abs(s) / surf.neck.S
            w = (q*q*(3 - 2*q))[:, None, None]
            # ids already contains the actual far-rim gluing permutation.
            rim = np.where((s >= 0)[:, None, None], rings[0], rings[-1])
            breathing = self.excursion * surf.neck.r(s)[:, None, None] * transverse
            vectors[ids[1:-1]] = (w*self.excursion*rim + (1-w)*breathing)[1:-1]
        self.reference = surf.points[self.cut.orig].copy()
        self.vectors = vectors[self.cut.orig]

    def coordinates(self, field, amplitude: float, gain: float = 6.0):
        """Use external-time cut-mesh samples and a fixed launch normalization.

        Positive and negative values produce opposite displacements. No
        per-frame normalization or synthetic arrival/ring geometry is used.
        Duplicated waist vertices remain distinct when their clocks differ.
        """
        field = np.asarray(field, dtype=float)
        if field.shape != (len(self.reference),) or not np.isfinite(field).all():
            raise ValueError("field must be finite and sampled on the cut mesh")
        if not np.isfinite(amplitude) or amplitude <= 0:
            raise ValueError("launch amplitude must be finite and positive")
        if not np.isfinite(gain) or gain < 0:
            raise ValueError("display gain must be finite and nonnegative")
        return self.reference + np.tanh(gain*field/amplitude)[:, None]*self.vectors
