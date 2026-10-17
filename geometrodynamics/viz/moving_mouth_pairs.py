"""Move #338's handle drawings with solved pair axes, keeping shared rims.

This coordinate map is NOT a fluid solution or a moving-metric wave operator.
"""
import numpy as np


def axis_rotation(start, end):
    """Shortest proper rotation between unit axes; no arbitrary neck spin."""
    start, end = np.asarray(start, float), np.asarray(end, float)
    if not (np.isclose(start@start, 1., atol=1e-10, rtol=0)
            and np.isclose(end@end, 1., atol=1e-10, rtol=0)):
        raise ValueError('Expected unit axes')
    dot = np.clip(start@end, -1., 1.)
    if dot < -1+1e-8:
        raise ValueError('Antipodal drawing rotation is ambiguous')
    v = np.cross(start, end)
    skew = np.array([[0., -v[2], v[1]], [v[2], 0., -v[0]], [-v[1], v[0], 0.]])
    return np.eye(3)+skew+skew@skew/(1+dot)


class MovingPairGeometry:
    def __init__(self, surface, collar=.68):
        self.surface = surface
        self.axes = np.array([d.c for d in surface.drawn])
        self.sphere = np.flatnonzero(surface.on_sphere)
        self.reference = surface.points.copy()
        if not surface.a < collar < np.pi/4:
            raise ValueError('Collars must be outside rims and disjoint in the reference layout')
        distance = np.arccos(np.clip(abs(self.reference[self.sphere]@self.axes.T), 0, 1))
        t = np.clip((distance-surface.a)/(collar-surface.a), 0, 1)
        self.weights = 1-t**3*(10-15*t+6*t*t)
        if np.any(np.sum(self.weights > 0, axis=1) > 1):
            raise ValueError('Reference collars overlap')

    def coordinates(self, axes):
        surf = self.surface
        axes = np.asarray(axes)
        if axes.shape != self.axes.shape:
            raise ValueError('One axis per handle required')
        if len(axes) == 2 and np.arccos(np.clip(abs(axes[0]@axes[1]), 0, 1)) <= 2*surf.a+.05:
            raise ValueError('Moving mouths overlap')
        points = self.reference.copy()
        original = self.reference[self.sphere]
        for k, (start, end) in enumerate(zip(self.axes, axes)):
            rotation = axis_rotation(start, end)
            points[self.sphere] += self.weights[:, k, None]*(original@rotation.T-original)
            ids = surf.tube_ids[k][1:-1].ravel()
            points[ids] = self.reference[ids]@rotation.T
        points[self.sphere] /= np.linalg.norm(points[self.sphere], axis=1)[:, None]
        return points

    def deform(self, points, sphere_field, amplitude, excursion=.15):
        """Radial field display on sphere; taper shared rim offsets into neutral necks."""
        if amplitude <= 0 or not np.isfinite(amplitude):
            raise ValueError('Positive finite fixed field scale required')
        sphere_field = np.asarray(sphere_field)
        if sphere_field.shape != (len(self.sphere),) or not np.isfinite(sphere_field).all():
            raise ValueError('One finite field sample per moved sphere vertex required')
        out = points.copy()
        displacement = np.zeros_like(points)
        displacement[self.sphere] = points[self.sphere]*(excursion*np.tanh(sphere_field/amplitude))[:, None]
        for ids in self.surface.tube_ids:
            # Actual far-rim permutation is already encoded by ids[-1].
            s = self.surface.s[ids[:, 0]]/self.surface.neck.S
            w = abs(s)**2*(3-2*abs(s))
            shifts = np.where((s >= 0)[:, None, None], displacement[ids[0]], displacement[ids[-1]])
            displacement[ids[1:-1]] = (w[:, None, None]*shifts)[1:-1]
        return out+displacement

    def minimum_sphere_orientation_ratio(self, points):
        """Signed tetrahedral volume ratio; positive means no spherical triangle flip."""
        triangles = self.surface.tris[self.surface.sphere_tris]
        def volume(p):
            a, b, c = np.moveaxis(p[triangles], 1, 0)
            return np.einsum('ij,ij->i', a, np.cross(b, c))
        return float(np.min(volume(points)/volume(self.reference)))
