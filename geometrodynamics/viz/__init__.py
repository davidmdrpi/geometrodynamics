"""Visualisation helpers (animation, plotting)."""

from geometrodynamics.viz.antipodal_crossing import (
    AntipodalCrossingSim,
    Absorption,
    Crossing,
    run_animation,
)
from geometrodynamics.viz.antipodal_focusing import (
    PlaneWaveSim,
    SphereWaveSim,
    RefocusResult,
    focusing_factor,
    measure_refocus,
    draw_contrast,
    draw_focus_object,
    draw_focusing_factor,
    plot_contrast_panel,
    plot_focusing_panel,
    plot_object_panel,
    run_animation as run_focusing_animation,
)
from geometrodynamics.viz.geometry_panels import (
    draw_green,
    draw_handshake,
    draw_hopf,
    draw_throat,
    plot_green_panel,
    plot_handshake_panel,
    plot_hopf_panel,
    plot_throat_panel,
    run_dashboard,
)

__all__ = [
    "AntipodalCrossingSim",
    "Absorption",
    "Crossing",
    "run_animation",
    # antipodal focusing
    "PlaneWaveSim",
    "SphereWaveSim",
    "RefocusResult",
    "focusing_factor",
    "measure_refocus",
    "draw_contrast",
    "draw_focus_object",
    "draw_focusing_factor",
    "plot_contrast_panel",
    "plot_focusing_panel",
    "plot_object_panel",
    "run_focusing_animation",
    # geometry panels
    "draw_green",
    "draw_handshake",
    "draw_hopf",
    "draw_throat",
    "plot_green_panel",
    "plot_handshake_panel",
    "plot_hopf_panel",
    "plot_throat_panel",
    "run_dashboard",
]
