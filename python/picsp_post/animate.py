"""Animations (mp4 / gif) of field slices or full 3D scenes over all output steps."""
from __future__ import annotations

import io as _io
import os

import numpy as np
import matplotlib.pyplot as plt

from .io import Run
from . import plots


def _writer(out, fps):
    import imageio.v2 as imageio
    ext = os.path.splitext(out)[1].lower()
    if ext == ".gif":
        return imageio.get_writer(out, mode="I", duration=1.0 / fps, loop=0)
    return imageio.get_writer(out, fps=fps, codec="libx264", quality=8, macro_block_size=8)


def _fig_to_rgb(fig):
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=fig.dpi, bbox_inches=None)
    buf.seek(0)
    import imageio.v2 as imageio
    return imageio.imread(buf)[..., :3]


def animate_slice(run: Run, name: str = "phi", axis: str = "z", pos=None, out: str = "slice.mp4", fps: int = 8,
                  comp=None, log=False, streamlines=False, fixed_limits=True, stride: int = 1):
    """Animate a field slice over time. Colour limits are fixed over the movie (robust range)."""
    series = run.series("fields")[::stride]
    n = len(series)
    if n == 0:
        raise FileNotFoundError("no field output")
    vmin = vmax = None
    if fixed_limits:
        lo, hi = [], []
        for i in range(0, len(run.series("fields")), stride):
            _, _, v, _, _ = plots.slice_data(run, i, name, axis, pos, comp)
            lo.append(np.nanpercentile(v, 0.5))
            hi.append(np.nanpercentile(v, 99.5))
        vmin, vmax = float(np.min(lo)), float(np.max(hi))
    w = _writer(out, fps)
    for k in range(n):
        idx = k * stride
        fig, ax = plt.subplots(figsize=(7.2, 5.8), dpi=120)
        plots.plot_slice(run, name, axis, pos, idx, comp, ax=ax, vmin=vmin, vmax=vmax, log=log, streamlines=streamlines)
        fig.tight_layout()
        w.append_data(_fig_to_rgb(fig))
        plt.close(fig)
    w.close()
    return out


def animate_scene(run: Run, out: str = "scene.mp4", fps: int = 8, window=(1280, 800), orbit: bool = False, **kw):
    """Animate the 3D scene (surfaces, slices, particles) over all field dumps.
    With orbit=True the camera rotates 360 degrees around the scene during the movie."""
    import pyvista as pv
    from .render3d import render_scene
    n = len(run.series("fields"))
    if n == 0:
        raise FileNotFoundError("no field output")
    w = _writer(out, fps)
    for i in range(n):
        p = pv.Plotter(off_screen=True, window_size=window)
        render_scene(run, i, plotter=p, **kw)
        if orbit:
            p.camera.azimuth = 360.0 * i / max(n, 1)
        img = p.screenshot(return_img=True)
        p.close()
        w.append_data(img[..., :3])
    w.close()
    return out


def orbit_movie(run: Run, index: int = -1, out: str = "orbit.mp4", frames: int = 90, fps: int = 24,
                window=(1280, 800), **kw):
    """360-degree camera orbit around a single time step."""
    import pyvista as pv
    from .render3d import render_scene
    p = pv.Plotter(off_screen=True, window_size=window)
    render_scene(run, index, plotter=p, **kw)
    w = _writer(out, fps)
    for f in range(frames):
        p.camera.azimuth = 360.0 / frames
        img = p.screenshot(return_img=True)
        w.append_data(img[..., :3])
    p.close()
    w.close()
    return out
