"""3D rendering with pyvista (off-screen capable): object surfaces coloured by
potential or current density, field slice planes, iso-surfaces and particles."""
from __future__ import annotations

import numpy as np

from . import style
from .io import Run

_SCALAR_TITLES = {
    "potential": "object potential [V]",
    "J_net": "net current density [A/m²]",
    "J_emit": "emitted current density [A/m²]",
    "sun_cos": "cos(sun incidence)",
}


def _pv():
    import pyvista as pv
    pv.global_theme.font.family = "arial"
    pv.global_theme.font.color = "#f0efec"
    return pv


def render_scene(run: Run, index: int = -1, out: str | None = None, surface_scalars: str = "potential",
                 slice_field: str | None = "phi", slice_normal: str = "y", slice_origin=None,
                 iso_field: str | None = None, iso_values=None, particles: bool = False, particle_species=None,
                 window=(1600, 1000), camera: str = "iso", zoom: float = 1.0, clim=None, slice_clim=None,
                 show_box: bool = True, title: str | None = None, plotter=None, dark: bool = True,
                 interactive: bool = False):
    """Render one output step. Returns the pyvista Plotter (screenshot saved if `out`).
    Rendering is off-screen unless `interactive=True` (then call ``.show()`` on the result)."""
    pv = _pv()
    own = plotter is None
    p = plotter or pv.Plotter(off_screen=not interactive, window_size=window)
    if dark:
        p.set_background("#101418", top="#2a3340")
        text_color = "#f0efec"
    else:
        p.set_background("#fcfcfb")
        text_color = "#0b0b0b"
    sbar = dict(color=text_color, title_font_size=18, label_font_size=15, n_labels=5, fmt="%.3g", vertical=True,
                width=0.05, height=0.32, shadow=False)
    grid = None
    try:
        grid = run.fields(index)
    except FileNotFoundError:
        pass

    # ---- object surfaces
    surf = None
    try:
        surf = run.surface(index if len(run.series("surface")) > 1 else 0)
    except FileNotFoundError:
        if grid is not None and "object_id" in grid.point_data:
            surf = grid.contour([-0.5], scalars="object_id").extract_surface()
            surface_scalars = None
    if surf is not None and surf.n_points:
        if surface_scalars and surface_scalars in surf.cell_data:
            data = np.asarray(surf.cell_data[surface_scalars])
            sb = dict(sbar, title=_SCALAR_TITLES.get(surface_scalars, surface_scalars), position_x=0.88, position_y=0.58)
            if surface_scalars.startswith("J_"):
                lim = clim or (np.percentile(np.abs(data), 98) or 1e-30)
                lim = lim if np.isscalar(lim) else max(abs(lim[0]), abs(lim[1]))
                p.add_mesh(surf, scalars=surface_scalars, cmap=style.CMAP_DIV, clim=(-lim, lim), smooth_shading=True,
                           scalar_bar_args=sb, specular=0.3)
            else:
                cl = clim or (float(data.min()), float(data.max()) if data.max() > data.min() else float(data.min()) + 1e-9)
                p.add_mesh(surf, scalars=surface_scalars, cmap="plasma", clim=cl, smooth_shading=True,
                           scalar_bar_args=sb, specular=0.4, specular_power=20)
        else:
            p.add_mesh(surf, color="#c9ccd1", smooth_shading=True, specular=0.5, specular_power=25)

    # ---- field slice
    if grid is not None and slice_field and slice_field in grid.point_data:
        sl = grid.slice(normal=slice_normal, origin=slice_origin or grid.center)
        name = slice_field
        vals = np.asarray(sl.point_data[name])
        if vals.ndim == 2:
            sl.point_data[name + "_mag"] = np.linalg.norm(vals, axis=1)
            name = name + "_mag"
            vals = np.asarray(sl.point_data[name])
        div = slice_field in ("phi", "rho")
        if slice_clim is None:
            if div:
                lim = float(np.max(np.abs(vals))) or 1e-30
                slice_clim = (-lim, lim)
            else:
                slice_clim = (float(np.percentile(vals, 1)), float(np.percentile(vals, 99.5)))
        title_map = {"phi": "potential φ [V]", "rho": "charge density [C/m³]"}
        p.add_mesh(sl, scalars=name, cmap=style.CMAP_DIV if div else style.CMAP_SEQ, clim=slice_clim, opacity=0.92,
                   scalar_bar_args=dict(sbar, title=title_map.get(slice_field, slice_field), position_x=0.88,
                                        position_y=0.12))
    # ---- iso-surfaces
    if grid is not None and iso_field and iso_field in grid.point_data:
        vals = np.asarray(grid.point_data[iso_field])
        if iso_values is None:
            lo, hi = np.percentile(vals, 5), np.percentile(vals, 95)
            iso_values = list(np.linspace(lo, hi, 4)[1:-1])
        iso = grid.contour(iso_values, scalars=iso_field)
        if iso.n_points:
            p.add_mesh(iso, scalars=iso_field, cmap=style.CMAP_SEQ, opacity=0.35, show_scalar_bar=False)
    # ---- particles
    if particles:
        try:
            pd = run.particles(index)
            sid = np.asarray(pd.point_data["species"]).astype(int)
            for s in np.unique(sid):
                if particle_species is not None and s not in np.atleast_1d(particle_species):
                    continue
                sub = pd.extract_points(sid == s, adjacent_cells=False)
                p.add_mesh(sub, color=style.series_color(int(s)), point_size=2.5, render_points_as_spheres=True,
                           opacity=0.7)
        except FileNotFoundError:
            pass
    if show_box and grid is not None:
        p.add_mesh(grid.outline(), color="#8b93a1", line_width=1.2)
    t = run.time("fields", index) if grid is not None else None
    p.add_text(title or (f"t = {t:.3e} s" if t is not None else ""), position="upper_left", font_size=12,
               color=text_color)
    if camera == "iso":
        p.camera_position = "iso"
    elif camera in ("xy", "xz", "yz"):
        p.camera_position = camera
    p.camera.zoom(zoom)
    if out:
        p.screenshot(out)
        if own:
            p.close()
    return p
