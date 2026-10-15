#!/usr/bin/env python
"""Klein-handle wormhole explorer (PyVista), with mouth clock offsets.

A unit sphere with one or two antipodal finite-mouth handles, each glued either
orientation-reversingly (Klein-type, the default) or orientation-preservingly
(torus-type), and each carrying a Morris–Thorne–Yurtsever clock offset Δ. A
scalar wave launched on the sphere runs on the *intrinsic* glued surface (true
neck length 2S) and deforms the displayed surface with its signed amplitude,
in external (Killing) time. The intrinsic metric remains fixed. Geometry, topology, clocks and waves come
from ``geometrodynamics.viz.klein_handles``; this file is only the viewer.

Run
    python scripts/klein_wormhole_explorer.py
    python scripts/klein_wormhole_explorer.py --offsets=-0.8,0.6
    python scripts/klein_wormhole_explorer.py --handles klein,torus --a 0.3
    python scripts/klein_wormhole_explorer.py --screenshot klein.png --time 1.1
    python scripts/klein_wormhole_explorer.py --gif klein_waves.gif

Controls
    Space         play / pause            Slider (bottom left)   time t
    1, 2          handle 1 / 2: Klein <-> torus     Slider (bottom right)  mouth angle a
    h             one or two handles      Sliders (top right)    clock offsets Δ1, Δ2
    d             deformation / fixed geometry     f   first-arrival contour
    c             wave colours / neutral shading   Gain slider: display exaggeration
    m             colour: wave / first arrival / slicing tilt / curvature
    z             set both offsets to zero
    g             see-through sphere      o   orientation loop through handle 1
    Right-click   put the wave source at that point on the sphere and relaunch

Reading the picture
    * The reference sphere and mouth rims are exact. The drawn neck is an
      immersion with a stretched middle and a reshaped Klein exit.
    * Shape is driven by the computed signed scalar, with fixed launch
      normalization: sphere radius = 1 + 0.22*tanh(gain*u/u_ref). In handles
      a continuous ambient vector gives transverse breathing and matches
      sphere displacement at the rims. It is not a global normal field.
      This is a display of a fixed-metric wave, not gravitational backreaction.
    * Default: opaque deformed surface, no arrival contour. Press d to compare
      fixed geometry; c removes wave colour to inspect shape alone. No feature
      is placed at a prescribed radius or arrival time.
    * A Klein-type handle passes through the sphere wall and joins its far mouth
      from inside. That crossing is a drawing artifact.
    * Clock offset Δ = o_B − o_A (the convention of transaction.network): a wave
      crossing a neck from A to B emerges Δ later in external time, one crossing
      B to A emerges Δ earlier. The jump sits at the gold ring on the neck
      waist, so the colour is discontinuous there. Δ < −2S makes a wave leave
      mouth B before it has reached mouth A.
    * Chronology: a closed causal curve exists when some closed loop gains at
      least its own length in clock offset. The status line gives the largest
      gain/length ratio and the loop that attains it. Two offsets that are each
      safe alone can close a loop through both handles (the Roman ring). The
      "tilt" colour shows the slope of the slicing used for the evolution, and is
      brightest along the loop that is closest to closing.
    * The wave colour is log-compressed: a narrow neck reflects most of an
      incoming wave, and what crosses is a few percent of the launch amplitude.
    * The optional white curve (f) is a first-threshold-crossing contour: where the field first exceeds
      0.5% of its launch amplitude by the current external time. Before the
      launch (t < 0) a front can appear only beyond a mouth with an advanced
      clock.

Requires numpy, scipy and pyvista (pip install pyvista); --gif also needs imageio.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

try:
    import pyvista as pv
except ImportError:  # pragma: no cover
    sys.exit("This explorer needs PyVista: pip install pyvista")

from geometrodynamics.viz.klein_handles import (
    ChronologyError,
    HandleWaveSim,
    build_surface,
    chronology_ratio,
    default_handles,
    euler_characteristic,
    is_orientable,
    loop_lengths,
    orientation_loop,
    source_clock,
    symlog_curvature,
)

from geometrodynamics.viz.scalar_deformation import ScalarDeformation

MODES = ("wave", "arrival", "tilt", "curvature")
HELP = ("Space play/pause   1/2 handle type   h handles   m colour   z zero offsets   "
        "d deform   f front   c colour\ng see-through   o reference loop   right-click source")
OFFSET_RANGE = 3.3
GOLD = "#ffd166"


def _faces(tris: np.ndarray) -> np.ndarray:
    return np.hstack([np.full((len(tris), 1), 3, dtype=np.int64), tris.astype(np.int64)]).ravel()


class Explorer:
    def __init__(self, a=0.3, kinds=(True, True), source=(-0.55, 0.62, 0.56), offsets=(0.0, 0.0),
                 off_screen=False, window_size=(1400, 900), deform=True, gain=6.0, front=False,
                 neutral=False):
        self.a = float(a)
        self.kinds = list(kinds)
        self.offsets = [float(v) for v in (list(offsets) + [0.0, 0.0])[:2]]
        self.source = np.asarray(source, float) / np.linalg.norm(source)
        self.mode = "wave"
        self.glass = False
        if not np.isfinite(gain) or gain < 0:
            raise ValueError("gain must be finite and nonnegative")
        self.deform, self.gain = bool(deform), float(gain)
        self.show_front, self.neutral = bool(front), bool(neutral)
        self.playing = False
        self.show_loop = False
        self.frame = 0
        self.off_screen = off_screen
        self.pl = pv.Plotter(off_screen=off_screen, window_size=list(window_size),
                             title="Klein-handle wormholes")
        self.pl.set_background("#0e1320", top="#1b2436")
        self.build()
        self.pl.camera_position = [(3.6, -4.2, 2.7), (0.25, -0.05, 0.0), (0.0, 0.0, 1.0)]
        self.pl.reset_camera()  # keep the view direction, fit both handles in frame
        self.pl.camera.zoom(1.4)
        self._add_widgets()

    # ── model ────────────────────────────────────────────────────────────────
    @property
    def active_offsets(self):
        return self.offsets[:len(self.kinds)]

    def build(self):
        """Geometry, then clocks and waves."""
        handles = default_handles(len(self.kinds), self.kinds) if self.kinds else ()
        self.surf = surf = build_surface(self.a, handles)
        self.cut = surf.cut()
        self.deformation = ScalarDeformation(surf, self.cut)
        self.loops = loop_lengths(surf)
        pts = surf.points[self.cut.orig]
        ns = len(surf.sphere_tris)
        sphere_tris = self.cut.tris[:ns].copy()
        tri_points = pts[sphere_tris]
        normal = np.cross(tri_points[:, 1]-tri_points[:, 0], tri_points[:, 2]-tri_points[:, 0])
        reverse = np.einsum("ij,ij->i", normal, tri_points.mean(axis=1)) < 0
        sphere_tris[reverse] = sphere_tris[reverse][:, [0, 2, 1]]
        self.sphere_mesh = pv.PolyData(pts, _faces(sphere_tris))
        self.tube_mesh = pv.PolyData(pts, _faces(self.cut.tris[ns:]))
        self.full_mesh = pv.PolyData(pts, _faces(self.cut.tris))
        self.curv = symlog_curvature(surf)[self.cut.orig]
        self.chi = euler_characteristic(surf.tris)
        self.orientable = is_orientable(surf.tris)
        self.simulate(draw=False)
        self._draw_static()
        self.update()

    def simulate(self, draw=True):
        """Clock field and wave evolution for the current offsets and source."""
        surf = self.surf
        self.ratio, self.binding = chronology_ratio(self.loops, self.active_offsets)
        self.clock = source_clock(surf, self.source, self.active_offsets, cut=self.cut)
        tilt = np.zeros(len(self.cut.orig))
        np.maximum.at(tilt, self.cut.tris.ravel(), np.repeat(self.clock.slope, 3))
        self.tilt = tilt
        self.sim, self.sim_error = None, ""
        try:
            self.sim = HandleWaveSim(surf, source=self.source, clock=self.clock).run()
        except ChronologyError as err:
            self.sim_error = str(err)
        if self.sim is not None:
            self.times = self.sim.times
            arrival = self.sim.arrival.copy()  # unobserved vertices remain undefined
        else:
            self.times = np.array([0.0])
            arrival = np.full(len(self.cut.orig), np.inf)
        self.full_mesh.point_data["arrival"] = arrival
        self.frame = self._frame_at(0.0) if self.sim is not None else 0
        slider = getattr(self, "time_slider", None)
        if slider is not None:
            rep = slider.GetRepresentation()
            rep.SetMinimumValue(float(self.times[0]))
            rep.SetMaximumValue(float(max(self.times[-1], self.times[0] + 1e-3)))
            rep.SetValue(float(self.times[self.frame]))
        if draw:
            self._draw_static()
            self.update()

    def _frame_at(self, t):
        return self.sim.frame_at(t) if self.sim is not None else 0

    def _scalars(self):
        mode = self.mode
        if self.sim is None and mode in ("wave", "arrival"):
            mode = "tilt"
        if mode == "wave":
            values = np.zeros(len(self.cut.orig)) if self.neutral else self.sim.compressed(self.frame)
            return values, "RdBu_r", (-1.0, 1.0), "u (log-compressed)"
        if mode == "arrival":
            hit = self.sim.arrival[np.isfinite(self.sim.arrival)]
            lo, hi = (float(hit.min()), float(hit.max())) if hit.size else (0.0, self.sim.t_max)
            return self.full_mesh.point_data["arrival"], "viridis", (min(lo, 0.0), hi), "first arrival t"
        if mode == "tilt":
            return self.tilt, "magma", (0.0, 1.0), "slicing tilt"
        return self.curv, "RdBu", (-1.0, 1.0), "curvature (signed log)"

    def _draw_static(self):
        pl, surf = self.pl, self.surf
        values, cmap, clim, title = self._scalars()
        for old in list(pl.scalar_bars.keys()):  # replacing an actor leaves its bar behind
            pl.remove_scalar_bar(old)
        for name, mesh, opacity in (("sphere", self.sphere_mesh, 0.38 if self.glass else 1.0),
                                    ("tubes", self.tube_mesh, 1.0)):
            mesh.point_data["c"] = values
            if mesh.n_cells == 0:
                pl.remove_actor(name)
                continue
            actor = pl.add_mesh(mesh, scalars="c", cmap=cmap, clim=clim, name=name, opacity=opacity,
                        smooth_shading=False, specular=0.15,
                        show_scalar_bar=(name == ("tubes" if self.tube_mesh.n_cells else "sphere")
                                         and not (self.neutral and self.mode == "wave")),
                        scalar_bar_args=dict(title=title, vertical=True, position_x=0.86, position_y=0.25,
                                             height=0.45, width=0.04, color="white"))
            if name == "sphere":
                actor.GetProperty().SetInterpolationToPhong()
        wlab, wpos = [], []
        for k, delta in enumerate(self.active_offsets):
            if delta:
                p = surf.points[surf.tube_ids[k][surf.waist]]
                wlab.append(f"Δ{k + 1} = {delta:+.2f}")
                wpos.append(p.mean(axis=0) + 0.12 * surf.drawn[k].e1)
        labels, lpos = [], []
        for k, d in enumerate(surf.drawn):
            for tag, sgn in (("A", 1), ("B", -1)):
                labels.append(f"{tag}{k + 1}")
                lpos.append(1.18 * sgn * d.c)
        for name, pos, text, color in (("mouth_labels", lpos, labels, "white"), ("waist_labels", wpos, wlab, GOLD)):
            if text:
                pl.add_point_labels(np.array(pos), text, name=name, font_size=12, text_color=color,
                                    shape=None, show_points=False, always_visible=True)
            else:
                pl.remove_actor(name)
        self._draw_loop()

    def _draw_loop(self):
        if not (self.show_loop and self.surf.drawn):
            for n in ("loop", "loop_arrows", "loop_start"):
                self.pl.remove_actor(n)
            return
        pts, nrm, flipped = orientation_loop(self.surf, 0)
        self.pl.add_mesh(pv.lines_from_points(pts), color=GOLD, line_width=3, name="loop")
        idx = np.linspace(0, len(pts) - 1, 18).astype(int)
        cloud = pv.PolyData(pts[idx])
        cloud["n"] = nrm[idx]
        self.pl.add_mesh(cloud.glyph(orient="n", scale=False, factor=0.16, geom=pv.Arrow()),
                         color=GOLD, name="loop_arrows")
        start = pv.PolyData(pts[:1])
        start["n"] = nrm[:1]
        self.pl.add_mesh(start.glyph(orient="n", scale=False, factor=0.2, geom=pv.Arrow()),
                         color="#7ae582", name="loop_start")
        self.loop_flipped = flipped

    def update(self):
        values, _, _, _ = self._scalars()
        points = self.deformation.reference
        if self.deform and self.sim is not None and self.mode == "wave":
            points = self.deformation.coordinates(self.sim.frames[self.frame], self.sim.amplitude, self.gain)
        # Update the actual mesh vertices, including the contour's input mesh.
        for mesh in (self.sphere_mesh, self.tube_mesh, self.full_mesh):
            mesh.points = points.copy()
        # Smooth lighting only: vertices and the computed field are unchanged.
        # The sphere band is orientable; no global normals are imposed on a Klein handle.
        self.sphere_mesh.compute_normals(cell_normals=False, point_normals=True,
                                         consistent_normals=False, inplace=True)
        self._moving_markers(points)
        for mesh in (self.sphere_mesh, self.tube_mesh):
            mesh.point_data["c"][:] = values
            mesh.Modified()
        t = float(self.times[self.frame])
        lead = self.show_front and self.sim is not None and self.mode in ("wave", "arrival")
        front = self.full_mesh.contour([t], scalars="arrival") if lead and t > self.times[0] else None
        if front is not None and front.n_points:
            self.pl.add_mesh(front, color="white", line_width=3, name="front")
        else:
            self.pl.remove_actor("front")
        self._status(t)

    def _moving_markers(self, points):
        for name, rings, color, width in (
            ("rims", [r for pair in self.surf.rims for r in pair], "#8fa3c2", 1.0),
            ("waists", [ring for k, ids in enumerate(self.surf.tube_ids)
                        if self.active_offsets[k] != 0 for ring in
                        (ids[self.surf.waist], np.flatnonzero(self.cut.handle == k))], GOLD, 2.0),
        ):
            lines = [pv.lines_from_points(np.vstack([points[r], points[r[:1]]])) for r in rings]
            if lines:
                self.pl.add_mesh(pv.merge(lines), color=color, line_width=width, name=name)
            else:
                self.pl.remove_actor(name)
        sphere_ids = np.flatnonzero(self.surf.on_sphere)
        j = sphere_ids[np.argmax(self.surf.points[sphere_ids] @ self.source)]
        self.pl.add_mesh(pv.Sphere(radius=0.025, center=points[j]), color=GOLD, name="source")

    def _status(self, t):
        surf, nk = self.surf, self.surf.neck
        kinds = ", ".join("Klein" if k else "torus" for k in self.kinds) or "none"
        lines = [f"handles: {kinds}    chi = {self.chi}    "
                 f"{'orientable' if self.orientable else 'non-orientable'}",
                 f"a = {self.a:.3f}   mouth radius {nk.sa:.3f}   waist b = {nk.b:.3f}   "
                 f"neck length 2S = {2 * nk.S:.3f}   (chord {2 * nk.ca:.3f})"]
        if self.kinds:
            offs = "   ".join(f"Δ{k + 1} = {d:+.2f}" for k, d in enumerate(self.active_offsets))
            lines.append(f"clock offsets (A->B emerges Δ later): {offs}")
            if self.binding is None:
                lines.append("chronology: no offsets, static")
            else:
                lp = self.binding
                bound = " (lower bound: path passes other mouths)" if lp.obstructed else ""
                verdict = "no closed causal curve" if self.ratio < 1 else "CLOSED CAUSAL CURVES"
                lines.append(f"chronology: max gain/length = {self.ratio:.2f} on {lp.label} "
                             f"(L = {lp.length:.3f}{bound})  ->  {verdict}")
            refused = self.sim_error.split(":")[0].split(";")[0]
            lines.append(f"slicing tilt (largest slope of theta) = {self.clock.tilt:.2f}" +
                         (f"    no evolution: {refused}" if self.sim is None else ""))
        shape = "deformed" if self.deform and self.sim is not None and self.mode == "wave" else "reference"
        lines.append(f"t = {t:.2f}    colour: {self.mode}    shape: {shape}    gain: {self.gain:.1f}")
        if self.sim is not None:
            raw = self.sim.frames[self.frame] / self.sim.amplitude
            lines.append(f"u/u_ref range: [{raw.min():+.3f}, {raw.max():+.3f}]    "
                         f"tau-energy drift: {self.sim.energy_drift:.1e}")
        lines.append("Shape = display of solved scalar; fixed intrinsic metric")
        if self.show_loop and surf.drawn:
            lines.append("reference-geometry loop through handle 1: normal returns "
                         + ("REVERSED (one-sided)" if self.loop_flipped else "unchanged (two-sided)"))
        self.pl.add_text("\n".join(lines), position="upper_left", font_size=10, color="white", name="status")
        self.pl.add_text(HELP, position=(12, 8), font_size=8, color="#aab4c8", name="help")

    # ── interaction ──────────────────────────────────────────────────────────
    def _slider(self, callback, rng, value, title, a, b, fmt):
        kw = dict(rng=rng, value=value, title=title, fmt=fmt, pointa=a, pointb=b, style="modern",
                  color="white", title_height=0.018)
        try:
            widget = self.pl.add_slider_widget(callback, interaction_event="end", **kw)
        except TypeError:  # older PyVista without interaction_event
            widget = self.pl.add_slider_widget(callback, **kw)
        self._style_slider(widget)
        return widget

    @staticmethod
    def _style_slider(widget):
        rep = widget.GetRepresentation()
        rep.SetSliderWidth(0.018)
        rep.SetSliderLength(0.018)
        rep.SetTubeWidth(0.006)
        rep.SetEndCapLength(0.008)
        rep.SetEndCapWidth(0.018)
        rep.SetTitleHeight(0.016)
        rep.SetLabelHeight(0.018)

    def _add_widgets(self):
        pl = self.pl
        self.time_slider = pl.add_slider_widget(
            self._on_time, rng=[float(self.times[0]), float(max(self.times[-1], self.times[0] + 1e-3))],
            value=float(self.times[self.frame]), title="time t", fmt="%.2f",
            pointa=(0.05, 0.1), pointb=(0.42, 0.1), style="modern", color="white")
        self._style_slider(self.time_slider)
        self._slider(self._on_a, [0.15, 0.45], self.a, "mouth angle a", (0.55, 0.1), (0.85, 0.1), "%.3f")
        self.offset_sliders = [
            self._slider(lambda v, k=k: self._on_offset(k, v), [-OFFSET_RANGE, OFFSET_RANGE], self.offsets[k],
                         f"clock offset Δ{k + 1}", (0.62, 0.93 - 0.1 * k), (0.95, 0.93 - 0.1 * k), "%+.2f")
            for k in range(2)]
        pl.add_checkbox_button_widget(self._on_play, value=False, position=(12, 140), size=32,
                                      color_on="#7ae582", color_off="#3a4458")
        pl.add_text("play", position=(52, 146), font_size=9, color="white", name="play_label")
        pl.add_key_event("space", lambda: self._on_play(not self.playing))
        pl.add_key_event("1", lambda: self._toggle_kind(0))
        pl.add_key_event("2", lambda: self._toggle_kind(1))
        pl.add_key_event("h", self._toggle_count)
        pl.add_key_event("m", self._cycle_mode)
        pl.add_key_event("z", self._zero_offsets)
        pl.add_key_event("d", self._toggle_deformation)
        pl.add_key_event("f", self._toggle_front)
        pl.add_key_event("c", self._toggle_neutral)
        self._slider(self._on_gain, [0.0, max(20.0, self.gain)], self.gain, "shape gain (display only)",
                     (0.55, 0.23), (0.85, 0.23), "%.1f")
        pl.add_key_event("g", self._toggle_glass)
        pl.add_key_event("o", self._toggle_loop)
        try:
            pl.enable_surface_point_picking(callback=self._on_pick, show_point=False, show_message=False,
                                            left_clicking=False)
        except Exception:  # pragma: no cover - picking API differs across versions
            pass
        if not self.off_screen:
            pl.add_timer_event(max_steps=10 ** 9, duration=40, callback=self._tick)

    def _set_frame(self, k):
        self.frame = int(np.clip(k, 0, len(self.times) - 1))
        slider = getattr(self, "time_slider", None)  # absent while the widget is being created
        if slider is not None:
            slider.GetRepresentation().SetValue(float(self.times[self.frame]))
        self.update()

    def _on_time(self, value):
        self._set_frame(self._frame_at(value))

    def _tick(self, _step=None):
        if not self.playing or self.sim is None:
            return
        self._set_frame((self.frame + 1) % len(self.times))
        self.pl.render()

    def _on_play(self, state):
        self.playing = bool(state)

    def _busy(self):
        self.pl.add_text("computing ...", position="upper_right", font_size=11, color=GOLD, name="busy")
        self.pl.render()

    def _done(self):
        self.pl.remove_actor("busy")
        self.pl.render()

    def _rebuild(self):
        self._busy()
        cam = self.pl.camera_position
        self.build()
        self.pl.camera_position = cam
        self._done()

    def _resimulate(self):
        self._busy()
        self.simulate()
        self._done()

    def _on_a(self, value):
        if abs(value - self.a) > 1e-6:
            self.a = float(value)
            self._rebuild()

    def _on_offset(self, k, value):
        if abs(value - self.offsets[k]) > 1e-6:
            self.offsets[k] = float(value)
            if k < len(self.kinds):
                self._resimulate()

    def _zero_offsets(self):
        self.offsets = [0.0, 0.0]
        for s in self.offset_sliders:
            s.GetRepresentation().SetValue(0.0)
        self._resimulate()

    def _toggle_kind(self, i):
        if i < len(self.kinds):
            self.kinds[i] = not self.kinds[i]
            self._rebuild()

    def _toggle_count(self):
        self.kinds = self.kinds[:1] if len(self.kinds) == 2 else (self.kinds + [True])[:2]
        self._rebuild()

    def _cycle_mode(self):
        self.mode = MODES[(MODES.index(self.mode) + 1) % len(MODES)]
        self._draw_static()
        self.update()
        self.pl.render()

    def _toggle_deformation(self):
        self.deform = not self.deform
        self.update()
        self.pl.render()

    def _toggle_front(self):
        self.show_front = not self.show_front
        self.update()
        self.pl.render()

    def _toggle_neutral(self):
        self.neutral = not self.neutral
        self._draw_static()
        self.update()
        self.pl.render()

    def _on_gain(self, value):
        self.gain = float(value)
        self.update()
        self.pl.render()

    def _toggle_glass(self):
        self.glass = not self.glass
        self._draw_static()
        self.update()
        self.pl.render()

    def _toggle_loop(self):
        self.show_loop = not self.show_loop
        self._draw_loop()
        self.update()
        self.pl.render()

    def _on_pick(self, point, *args):
        x = np.asarray(point, float)
        # Picking is in displaced display coordinates; recover the intrinsic
        # sphere vertex rather than rejecting points whose radius is not one.
        j = int(np.argmin(np.linalg.norm(self.full_mesh.points - x, axis=1)))
        orig = self.cut.orig[j]
        if not self.surf.on_sphere[orig]:
            return
        self.source = self.surf.points[orig].copy()
        try:
            self._resimulate()
        except Exception:
            return

    # ── output ───────────────────────────────────────────────────────────────
    def show(self):
        self.pl.show()

    def screenshot(self, path, t):
        self._set_frame(self._frame_at(t))
        self.pl.show(screenshot=path, auto_close=True)

    def gif(self, path, stride=2):
        self.pl.open_gif(path)
        for k in range(0, len(self.times), stride):
            self._set_frame(k)
            self.pl.write_frame()
        self.pl.close()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--a", type=float, default=0.3, help="mouth angular radius (0.15 to 0.45)")
    ap.add_argument("--handles", default="klein,klein",
                    help="comma list of 'klein'/'torus', one or two entries, or 'none'")
    ap.add_argument("--offsets", default="0,0", help="clock offsets Δ1,Δ2 (Δ = o_B − o_A)")
    ap.add_argument("--source", default="-0.55,0.62,0.56", help="source direction x,y,z on the sphere")
    ap.add_argument("--mode", choices=MODES, default="wave")
    ap.add_argument("--flat", action="store_true", help="compare fixed display geometry")
    ap.add_argument("--gain", type=float, default=6.0, help="shape exaggeration; fixed launch normalization")
    ap.add_argument("--front", action="store_true", help="show first threshold crossing contour")
    ap.add_argument("--neutral", action="store_true", help="neutral shading to inspect shape without wave colour")
    ap.add_argument("--loop", action="store_true", help="show the orientation loop through handle 1")
    ap.add_argument("--screenshot", help="render one frame to this PNG and exit")
    ap.add_argument("--time", type=float, default=1.1, help="time for --screenshot")
    ap.add_argument("--gif", help="write an animation to this GIF and exit")
    args = ap.parse_args(argv)
    kinds = [] if args.handles.strip() == "none" else [k.strip().lower() == "klein" for k in args.handles.split(",")]
    headless = bool(args.screenshot or args.gif)
    ex = Explorer(a=args.a, kinds=kinds, source=[float(v) for v in args.source.split(",")],
                  offsets=[float(v) for v in args.offsets.split(",")], off_screen=headless,
                  deform=not args.flat, gain=args.gain, front=args.front, neutral=args.neutral)
    if args.mode != "wave":
        ex.mode = args.mode
        ex._draw_static()
    if args.loop:
        ex.show_loop = True
        ex._draw_loop()
    ex.update()
    if args.screenshot:
        ex.screenshot(args.screenshot, args.time)
    elif args.gif:
        ex.gif(args.gif)
    else:
        ex.show()


if __name__ == "__main__":
    main()
