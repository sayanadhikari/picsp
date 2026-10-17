"""Ion wake behind a sphere in mesothermal flow: animation of the ion density and the
final density / potential maps."""
import numpy as np
from _common import Run, simulate, save, media, plt, style, gif_writer, fig_to_rgb
from picsp_post import plots


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(f"{d}/output")
    n0 = 1e12
    w = gif_writer(media("08_wake.gif"), fps=8)
    for i in range(len(r.series("fields"))):
        f, ax = plt.subplots(figsize=(7.6, 3.9), dpi=90)
        plots.plot_slice(r, "n_ion", "z", 0.0, i, ax=ax, vmin=0, vmax=2.0 * n0,
                         title=f"ion density (flow →)  t = {r.time('fields', i) * 1e6:.2f} µs")
        f.tight_layout()
        w.append_data(fig_to_rgb(f))
        plt.close(f)
    w.close()
    print("wrote 08_wake.gif")
    fig, axs = plt.subplots(1, 2, figsize=(14, 3.9))
    plots.plot_slice(r, "n_ion", "z", 0.0, -1, ax=axs[0], vmin=0, vmax=2.0 * n0, title="ion density: wake and focusing")
    plots.plot_slice(r, "phi", "z", 0.0, -1, ax=axs[1], streamlines=True, title="potential")
    fig.tight_layout()
    save(fig, d, "08_wake.png")
