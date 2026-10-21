"""3U CubeSat from STL: surface current map, plasma slices and the ion-wake movie."""
import numpy as np
from _common import Run, simulate, save, media, plt, style, gif_writer, fig_to_rgb, out, render3d
from picsp_post import plots


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(out(d))
    from picsp_post.render3d import render_scene
    p = render_scene(r, -1, None, surface_scalars="J_net", slice_field="n_O+", slice_normal="z",
                     slice_origin=(0.04, 0, -0.03), window=(1800, 1100), zoom=2.0, show_box=not style.PAPER,
                     title="3U CubeSat: net surface current density and O+ density (z = -3 cm)")
    p.camera.azimuth = -35
    p.camera.elevation = 10
    render3d(p, d, "12_cubesat_render.png", width_mm=183, zoom=0.55)
    p.close()
    fig, axs = plt.subplots(2, 2, figsize=(14, 8.4))
    plots.plot_slice(r, "n_O+", "z", 0.03, -1, ax=axs[0, 0], title="O+ density, z = 3 cm (through the bus)")
    plots.plot_slice(r, "n_electron", "z", 0.03, -1, ax=axs[0, 1], title="electron density, z = 3 cm")
    plots.plot_slice(r, "phi", "y", 0.0, -1, ax=axs[1, 0], streamlines=True, title="potential, y = 0 (side view)")
    plots.plot_slice(r, "n_O+", "y", 0.0, -1, ax=axs[1, 1], title="O+ density, y = 0 (side view, flow along +x)")
    fig.tight_layout()
    save(fig, d, "12_cubesat_slices.png")
    h = r.history()
    print(f"spacecraft {h['V_spacecraft'][-50:].mean():.3f} V, probe {h['V_probe'][-50:].mean():.3f} V")
    w = gif_writer(media("12_cubesat_wake.gif"), fps=6)
    for i in range(len(r.series("fields"))):
        f, ax = plt.subplots(figsize=(7.2, 5.6), dpi=80)
        plots.plot_slice(r, "n_O+", "z", -0.03, i, ax=ax, vmin=0, vmax=1.0e10,
                         title=f"O+ density below the panels   t = {r.time('fields', i) * 1e6:5.1f} µs")
        f.tight_layout()
        w.append_data(fig_to_rgb(f))
        plt.close(f)
    w.close()
    print("wrote 12_cubesat_wake.gif")
