#!/usr/bin/env python
"""Render the same scalar solution with fixed and moving display geometry.

Run from the repository root with PYTHONPATH=. (plus dependencies as needed).
Needs PyVista, Pillow, imageio and a working VTK rendering backend.
"""
import argparse
from pathlib import Path

import numpy as np
from PIL import Image

from scripts.klein_wormhole_explorer import Explorer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='docs/images/scalar_deformation')
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    ex = Explorer(source=(.5, -.7, .5), off_screen=True)
    ex.pl.show(auto_close=False, interactive=False)
    ex._set_frame(ex._frame_at(1.1))
    for deformed, name in [(False, 'fixed.png'), (True, 'deformed.png')]:
        ex.deform = deformed
        ex.update()
        ex.pl.screenshot(str(out/name))
    # Fixed camera and fixed launch normalization throughout the movie.
    frames = []
    for k in np.unique([ex._frame_at(t) for t in np.linspace(0, 3.4, 48)]):
        ex._set_frame(k)
        pic = ex.pl.screenshot(return_img=True)
        frames.append(Image.fromarray(pic).resize((700, 450)))
    # One palette sampled across the whole movie preserves later red/blue
    # crests while keeping the static background byte-identical between frames.
    swatches = Image.new('RGB', (8*175, 6*112))
    for i, frame in enumerate(frames):
        swatches.paste(frame.resize((175, 112)), ((i % 8)*175, (i // 8)*112))
    palette = swatches.quantize(colors=256)
    frames = [frame.quantize(palette=palette, dither=Image.Dither.NONE) for frame in frames]
    frames[0].save(out/'wave.gif', save_all=True, append_images=frames[1:], duration=110, loop=0, optimize=True)
    ex.pl.close()
    print(f'Rendered fixed.png, deformed.png, wave.gif to {out}')


if __name__ == '__main__':
    main()
