"""Photoemission: positive charging, photoelectron cloud and the shadow on the panel."""
import numpy as np
from _common import Run, simulate, save, media, plt, style
from picsp_post import plots


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(f"{d}/output")
    h = r.history()
    V = h["V_panel+mast"]
    print(f"floating potential {V[-50:].mean():.2f} V")
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.4))
    a = axs[0]
    a.plot(h["time"] * 1e6, V, color=style.series_color(0))
    a.set_xlabel("time [µs]")
    a.set_ylabel("potential [V]")
    a.set_title("positive charging by photoemission")
    plots.plot_slice(r, "n_photoelectron", "y", 0.0, -1, ax=axs[1], log=True, title="photoelectron cloud (side view)")
    fig.tight_layout()
    save(fig, d, "11_photoemission.png")
    import pyvista as pv
    from picsp_post.render3d import render_scene
    p = render_scene(r, -1, None, surface_scalars="J_emit", slice_field=None, window=(1500, 1000), zoom=1.9,
                     title="photoelectron current density: the mast shadows the panel")
    p.camera_position = [(1.6, -2.4, 2.4), (0, 0, -0.1), (0, 0, 1)]
    p.screenshot(media("11_photoemission_shadow.png"))
    p.close()
