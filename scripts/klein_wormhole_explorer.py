#!/usr/bin/env python
"""Klein-handle wormhole explorer (PyVista).

A unit sphere with one or two antipodal finite-mouth handles, each glued either
orientation-reversingly (Klein-type, the default) or orientation-preservingly
(torus-type). A scalar wave launched on the sphere runs on the *intrinsic* glued
surface (true neck length 2S), and is shown on the drawing as colour plus its
leading wavefront. Geometry, topology and waves come from
``geometrodynamics.viz.klein_handles``; this file is only the viewer.

Run
    python scripts/klein_wormhole_explorer.py
    python scripts/klein_wormhole_explorer.py --handles klein,torus --a 0.3
    python scripts/klein_wormhole_explorer.py --screenshot klein.png --time 1.1
    python scripts/klein_wormhole_explorer.py --gif klein_waves.gif

Controls
    Space         play / pause            Slider (bottom left)  time
    1, 2          handle 1 / 2: Klein <-> torus     Slider (bottom right) mouth angle a
    h             one or two handles      m   colour: wave / first arrival / curvature
    g             see-through sphere      o   orientation loop through handle 1
    Right-click   put the wave source at that point on the sphere and relaunch

Reading the picture
    * Sphere band, mouth rims, neck radius along its length and the flares at
      each mouth are exact. The thin middle of each neck is stretched (no 3D
      surface can show 2S < 2R cos a at true length).
    * A Klein-type handle passes through the sphere wall and joins its far mouth
      from inside. That crossing is a drawing artifact.
    * The wave colour is log-compressed: a narrow neck reflects most of an
      incoming wave, and what crosses is a few percent of the launch amplitude.
    * The white curve is the leading wavefront: where the field first exceeds
      0.5% of its launch amplitude by the current time.
    * Mouth clock offsets (Morris–Thorne–Yurtsever) are not modelled; the
      surface is static.

Requires numpy, scipy and pyvista (pip install pyvista); --gif also needs imageio.
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

try:
    import pyvista as pv
except ImportError:  # pragma: no cover
    sys.exit("This explorer needs PyVista: pip install pyvista")

from geometrodynamics.viz.klein_handles import (
    HandleWaveSim,
    build_surface,
    default_handles,
    euler_characteristic,
    is_orientable,
    orientation_loop,
    symlog_curvature,
)

MODES = ("wave", "arrival", "curvature")
HELP = ("Space play/pause   1/2 handle type   h handles   m colour   g see-through   "
        "o orientation loop   right-click source")


def _faces(tris: np.ndarray) -> np.ndarray:
    return np.hstack([np.full((len(tris), 1), 3, dtype=np.int64), tris.astype(np.int64)]).ravel()


class Explorer:
    def __init__(self, a=0.3, kinds=(True, True), source=(-0.55, 0.62, 0.56), off_screen=False,
                 window_size=(1400, 900)):
        self.a = float(a)
        self.kinds = list(kinds)
        self.source = np.asarray(source, float) / np.linalg.norm(source)
        self.mode = "wave"
        self.glass = True
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
    def build(self):
        handles = default_handles(len(self.kinds), self.kinds) if self.kinds else ()
        self.surf = build_surface(self.a, handles)
        self.sim = HandleWaveSim(self.surf, source=self.source).run()
        self.frame = min(self.frame, self.sim.n_frames - 1)
        surf = self.surf
        ns = len(surf.sphere_tris)
        self.sphere_mesh = pv.PolyData(surf.points, _faces(surf.tris[:ns]))
        self.tube_mesh = pv.PolyData(surf.points, _faces(surf.tris[ns:]))
        arrival = np.where(np.isnan(self.sim.arrival), self.sim.t_max, self.sim.arrival)
        self.full_mesh = pv.PolyData(surf.points, _faces(surf.tris))
        self.full_mesh.point_data["arrival"] = arrival
        self.curv = symlog_curvature(surf)
        self.chi = euler_characteristic(surf.tris)
        self.orientable = is_orientable(surf.tris)
        self._draw_static()
        self.update()

    def _scalars(self):
        if self.mode == "wave":
            return self.sim.compressed(self.frame), "RdBu_r", (-1.0, 1.0), "u (log-compressed)"
        if self.mode == "arrival":
            hit = self.sim.arrival[~np.isnan(self.sim.arrival)]
            top = float(hit.max()) if hit.size else self.sim.t_max
            return self.full_mesh.point_data["arrival"], "viridis", (0.0, top), "first arrival t"
        return self.curv, "RdBu", (-1.0, 1.0), "curvature (signed log)"

    def _draw_static(self):
        pl, surf = self.pl, self.surf
        values, cmap, clim, title = self._scalars()
        for old in list(pl.scalar_bars.keys()):  # replacing an actor leaves its bar behind
            pl.remove_scalar_bar(old)
        for name, mesh, opacity in (("sphere", self.sphere_mesh, 0.38 if self.glass else 1.0),
                                    ("tubes", self.tube_mesh, 1.0)):
            mesh.point_data["c"] = values
            pl.add_mesh(mesh, scalars="c", cmap=cmap, clim=clim, name=name, opacity=opacity,
                        smooth_shading=True, specular=0.25, show_scalar_bar=(name == "tubes"),
                        scalar_bar_args=dict(title=title, vertical=True, position_x=0.86, position_y=0.25,
                                             height=0.5, width=0.04, color="white"))
        rims = []
        for ra, rb in surf.rims:
            for ring in (ra, rb):
                p = surf.points[ring]
                rims.append(pv.lines_from_points(np.vstack([p, p[:1]])))
        if rims:
            pl.add_mesh(pv.merge(rims), color="#8fa3c2", line_width=2.0, name="rims")
        else:
            pl.remove_actor("rims")
        pl.add_mesh(pv.Sphere(radius=0.035, center=self.source), color="#ffd166", name="source")
        self._draw_loop()

    def _draw_loop(self):
        if not (self.show_loop and self.surf.drawn):
            for n in ("loop", "loop_arrows", "loop_start"):
                self.pl.remove_actor(n)
            return
        pts, nrm, flipped = orientation_loop(self.surf, 0)
        self.pl.add_mesh(pv.lines_from_points(pts), color="#ffd166", line_width=3, name="loop")
        idx = np.linspace(0, len(pts) - 1, 18).astype(int)
        cloud = pv.PolyData(pts[idx])
        cloud["n"] = nrm[idx]
        self.pl.add_mesh(cloud.glyph(orient="n", scale=False, factor=0.16, geom=pv.Arrow()),
                         color="#ffd166", name="loop_arrows")
        start = pv.PolyData(pts[:1])
        start["n"] = nrm[:1]
        self.pl.add_mesh(start.glyph(orient="n", scale=False, factor=0.2, geom=pv.Arrow()),
                         color="#7ae582", name="loop_start")
        self.loop_flipped = flipped

    def update(self):
        values, _, _, _ = self._scalars()
        for mesh in (self.sphere_mesh, self.tube_mesh):
            mesh.point_data["c"][:] = values
            mesh.Modified()
        t = float(self.sim.times[self.frame])
        if self.mode in ("wave", "arrival") and t > 0.0:
            front = self.full_mesh.contour([t], scalars="arrival")
            if front.n_points:
                self.pl.add_mesh(front, color="white", line_width=3, name="front")
            else:
                self.pl.remove_actor("front")
        else:
            self.pl.remove_actor("front")
        self._status(t)

    def _status(self, t):
        surf, nk = self.surf, self.surf.neck
        kinds = ", ".join("Klein" if k else "torus" for k in self.kinds) or "none"
        top = (f"handles: {kinds}    chi = {self.chi}    "
               f"{'orientable' if self.orientable else 'non-orientable'}\n"
               f"a = {self.a:.3f}   mouth radius {nk.sa:.3f}   waist b = {nk.b:.3f}   "
               f"neck length 2S = {2 * nk.S:.3f}   (chord {2 * nk.ca:.3f})\n"
               f"t = {t:.2f}    colour: {self.mode}")
        if self.show_loop and surf.drawn:
            top += ("\nloop through handle 1: normal returns "
                    + ("REVERSED (one-sided)" if self.loop_flipped else "unchanged (two-sided)"))
        self.pl.add_text(top, position="upper_left", font_size=10, color="white", name="status")
        self.pl.add_text(HELP, position="lower_edge", font_size=8, color="#aab4c8", name="help")

    # ── interaction ──────────────────────────────────────────────────────────
    def _add_widgets(self):
        pl = self.pl
        self.time_slider = pl.add_slider_widget(
            self._on_time, rng=[0.0, self.sim.t_max], value=0.0, title="time t", fmt="%.2f",
            pointa=(0.05, 0.1), pointb=(0.42, 0.1), style="modern", color="white")
        try:
            pl.add_slider_widget(self._on_a, rng=[0.15, 0.45], value=self.a, title="mouth angle a", fmt="%.3f",
                                 pointa=(0.55, 0.1), pointb=(0.85, 0.1), style="modern", color="white",
                                 interaction_event="end")
        except TypeError:  # older PyVista without interaction_event
            pl.add_slider_widget(self._on_a, rng=[0.15, 0.45], value=self.a, title="mouth angle a", fmt="%.3f",
                                 pointa=(0.55, 0.1), pointb=(0.85, 0.1), style="modern", color="white")
        pl.add_checkbox_button_widget(self._on_play, value=False, position=(12, 140), size=32,
                                      color_on="#7ae582", color_off="#3a4458")
        pl.add_text("play", position=(52, 146), font_size=9, color="white", name="play_label")
        pl.add_key_event("space", lambda: self._on_play(not self.playing))
        pl.add_key_event("1", lambda: self._toggle_kind(0))
        pl.add_key_event("2", lambda: self._toggle_kind(1))
        pl.add_key_event("h", self._toggle_count)
        pl.add_key_event("m", self._cycle_mode)
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
        self.frame = int(np.clip(k, 0, self.sim.n_frames - 1))
        slider = getattr(self, "time_slider", None)  # absent while the widget is being created
        if slider is not None:
            slider.GetRepresentation().SetValue(float(self.sim.times[self.frame]))
        self.update()

    def _on_time(self, value):
        self._set_frame(round(value / self.sim.t_max * (self.sim.n_frames - 1)))

    def _tick(self, _step=None):
        if not self.playing:
            return
        k = (self.frame + 1) % self.sim.n_frames
        self._set_frame(k)
        self.pl.render()

    def _on_play(self, state):
        self.playing = bool(state)

    def _rebuild(self):
        cam = self.pl.camera_position
        self.build()
        self.pl.camera_position = cam
        self.pl.render()

    def _on_a(self, value):
        if abs(value - self.a) > 1e-6:
            self.a = float(value)
            self._rebuild()

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
        if abs(np.linalg.norm(x) - 1.0) > 0.02:
            return  # only points on the sphere band can host the source
        self.source = x / np.linalg.norm(x)
        self.frame = 0
        try:
            self.sim = HandleWaveSim(self.surf, source=self.source).run()
        except Exception:
            return
        arrival = np.where(np.isnan(self.sim.arrival), self.sim.t_max, self.sim.arrival)
        self.full_mesh.point_data["arrival"] = arrival
        self._draw_static()
        self.update()
        self.pl.render()

    # ── output ───────────────────────────────────────────────────────────────
    def show(self):
        self.pl.show()

    def screenshot(self, path, t):
        self._set_frame(round(t / self.sim.t_max * (self.sim.n_frames - 1)))
        self.pl.show(screenshot=path, auto_close=True)

    def gif(self, path, stride=2):
        self.pl.open_gif(path)
        for k in range(0, self.sim.n_frames, stride):
            self._set_frame(k)
            self.pl.write_frame()
        self.pl.close()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--a", type=float, default=0.3, help="mouth angular radius (0.15 to 0.45)")
    ap.add_argument("--handles", default="klein,klein",
                    help="comma list of 'klein'/'torus', one or two entries, or 'none'")
    ap.add_argument("--source", default="-0.55,0.62,0.56", help="source direction x,y,z on the sphere")
    ap.add_argument("--mode", choices=MODES, default="wave")
    ap.add_argument("--loop", action="store_true", help="show the orientation loop through handle 1")
    ap.add_argument("--screenshot", help="render one frame to this PNG and exit")
    ap.add_argument("--time", type=float, default=1.1, help="time for --screenshot")
    ap.add_argument("--gif", help="write an animation to this GIF and exit")
    args = ap.parse_args(argv)
    kinds = [] if args.handles.strip() == "none" else [k.strip().lower() == "klein" for k in args.handles.split(",")]
    headless = bool(args.screenshot or args.gif)
    ex = Explorer(a=args.a, kinds=kinds, source=[float(v) for v in args.source.split(",")], off_screen=headless)
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
