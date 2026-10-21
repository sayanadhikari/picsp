"""CAD import: raw export vs placed object on the grid, and the charging of the imported body."""
import os
import subprocess
import sys

import numpy as np
from _common import Run, simulate, save, media, plt, style, out, screenshot
from matplotlib.collections import PolyCollection
from picsp_post.mesh import Mesh


def simulate_example(d):
    if not os.path.exists(os.path.join(d, "cad", "cubesat_mm.stl")):
        subprocess.run([sys.executable, os.path.join(d, "cad", "make_cad_export.py")], check=True)
    # the tool's view of the import (report and three-view preview)
    subprocess.run([sys.executable, "-m", "picsp_post", "mesh", "info", "cad/cubesat_mm.stl", "--units", "mm"], cwd=d,
                   check=True)
    simulate(d)


def _draw(ax, tris, i, j, color):
    k = 3 - i - j
    order = np.argsort(tris[:, :, k].mean(axis=1))
    ax.add_collection(PolyCollection(tris[order][:, :, [i, j]], facecolors=color, edgecolors=(0, 0, 0, 0.15),
                                     linewidths=0.2))
    ax.autoscale_view()
    ax.set_aspect("equal")


def analyze(d):
    raw = Mesh.load(os.path.join(d, "cad", "cubesat_mm.stl"))
    placed = raw.placed(units="mm", anchor="center", rotate=(0, 0, -30), center=(0, 0, 0))
    r = Run(out(d))
    g = r.fields(-1)
    o, sp, dims = np.asarray(g.origin), np.asarray(g.spacing), np.asarray(g.dimensions) - 1
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6))
    a = axs[0]
    _draw(a, raw.tris, 0, 1, style.series_color(1))
    a.plot(0, 0, "+", color=style.INK, ms=12)
    a.set_title("CAD export: millimetres, offset origin, tilted")
    a.set_xlabel("x [mm]")
    a.set_ylabel("y [mm]")
    a.set_xlim(min(-30, a.get_xlim()[0]), a.get_xlim()[1])
    a.set_ylim(min(-30, a.get_ylim()[0]), max(30, a.get_ylim()[1]))
    a = axs[1]
    for q in range(dims[0] + 1):
        a.plot([o[0] + q * sp[0]] * 2, [o[1], o[1] + dims[1] * sp[1]], color="0.88", lw=0.4, zorder=0)
    for q in range(dims[1] + 1):
        a.plot([o[0], o[0] + dims[0] * sp[0]], [o[1] + q * sp[1]] * 2, color="0.88", lw=0.4, zorder=0)
    _draw(a, placed.tris, 0, 1, style.series_color(0))
    sdf = r.field_array(-1, "sdf")
    kz = int(round(-o[2] / sp[2]))
    X, Y = np.meshgrid(o[0] + sp[0] * np.arange(dims[0] + 1), o[1] + sp[1] * np.arange(dims[1] + 1), indexing="ij")
    ins = sdf[:, :, kz] <= 0
    a.plot(X[ins], Y[ins], "o", ms=2.2, color=style.series_color(1), label="grid nodes inside (solver)")
    a.set_xlim(o[0], o[0] + dims[0] * sp[0])
    a.set_ylim(o[1], o[1] + dims[1] * sp[1])
    a.set_title("after import: metres, centred, axis-aligned")
    a.set_xlabel("x [m]")
    a.set_ylabel("y [m]")
    a.legend(loc="lower right", fontsize=8)
    h = r.history()
    a = axs[2]
    a.plot(h["time"] * 1e6, h["V_cubesat"], color=style.series_color(0))
    a.set_xlabel("time [µs]")
    a.set_ylabel("potential [V]")
    a.set_title("the imported spacecraft charges to its floating potential")
    fig.tight_layout()
    save(fig, d, "16_mesh_import.png")
    print(f"floating potential {h['V_cubesat'][-20:].mean():.2f} V")
    from picsp_post.render3d import render_scene
    p = render_scene(r, -1, None, surface_scalars="potential", slice_field="phi", slice_normal="z",
                     window=(1500, 950), zoom=1.6, title="imported CAD mesh: surface potential and plasma potential")
    p.camera.azimuth = -40
    p.camera.elevation = 18
    screenshot(p, media("16_mesh_import_render.png"))
    p.close()
