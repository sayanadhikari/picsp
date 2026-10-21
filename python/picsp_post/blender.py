"""Publication renders of pyvista scenes with Blender (Cycles path tracing).

A scene built for pyvista (``render_scene`` or any ``pyvista.Plotter``) is exported mesh by mesh
with the colours VTK would display, then lit and rendered by Blender in the background:

    from picsp_post import blender
    if blender.available():
        fig = blender.render_figure(plotter, "scene.png")   # render + colour bars as a figure

Data colours are preserved: coloured field slices are rendered unlit (emission), so they match
their colour bars exactly; objects get soft studio lighting and shadows. The image is rendered
on a transparent background and placed on white by ``render_figure``, which also draws the
colour bars with matplotlib in the figure style.

Blender is found from $PICSP_BLENDER or ``blender`` on PATH (tested with 4.2 LTS).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile

import numpy as np

SCENE_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_blender_scene.py")


def executable():
    return os.environ.get("PICSP_BLENDER") or shutil.which("blender")


def available():
    exe = executable()
    return bool(exe) and os.path.exists(exe)


def _write_ply(path, points, faces, rgba):
    """Binary PLY with per-vertex RGBA (what Blender's importer reads as a colour attribute)."""
    n, m = len(points), len(faces)
    head = ("ply\nformat binary_little_endian 1.0\n"
            f"element vertex {n}\nproperty float x\nproperty float y\nproperty float z\n"
            "property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\n"
            f"element face {m}\nproperty list uchar int vertex_indices\nend_header\n")
    v = np.zeros(n, dtype=[("p", "<f4", 3), ("c", "u1", 4)])
    v["p"], v["c"] = points, rgba
    f = np.zeros(m, dtype=[("k", "u1"), ("i", "<i4", 3)])
    f["k"], f["i"] = 3, faces
    with open(path, "wb") as fh:
        fh.write(head.encode())
        fh.write(v.tobytes())
        fh.write(f.tobytes())


NEUTRAL = (0.74, 0.76, 0.79)   # objects coloured by a constant (e.g. one floating potential)


def _degenerate(lo, hi):
    return hi - lo <= 1e-6 * max(abs(lo), abs(hi), 1e-30)


def _colors(mapper, ds, prop):
    """RGBA (uint8) per point or per cell as VTK maps them, or None for a uniform colour."""
    from vtkmodules.util.numpy_support import vtk_to_numpy
    if mapper.GetScalarVisibility() and _degenerate(*mapper.GetScalarRange()):
        prop.SetColor(*NEUTRAL)
        return None, None
    if mapper.GetScalarVisibility():
        # the array the mapper colours by, through its lookup table
        mode, name = mapper.GetScalarModeAsString(), mapper.GetArrayName()
        pd, cd = ds.GetPointData(), ds.GetCellData()
        if "Cell" in mode:
            arr = cd.GetArray(name) if name else cd.GetScalars()
        elif "Point" in mode:
            arr = pd.GetArray(name) if name else pd.GetScalars()
        else:
            arr = pd.GetScalars() or cd.GetScalars()
        if arr is None:
            return None, None
        lut = mapper.GetLookupTable()
        lut.SetRange(*mapper.GetScalarRange())
        mapped = lut.MapScalars(arr, 1, -1)    # VTK_COLOR_MODE_MAP_SCALARS, magnitude of vectors
        if mapped is not None and mapped.GetNumberOfTuples() in (ds.n_points, ds.n_cells):
            c = vtk_to_numpy(mapped).astype(np.uint8)
            if c.shape[1] == 3:
                c = np.c_[c, np.full(len(c), 255, np.uint8)]
            return c, ("point" if len(c) == ds.n_points else "cell")
    return None, None


def export_plotter(p, workdir, line_radius=None):
    """Write every visible 3D actor of a pyvista Plotter as a PLY and return the scene description."""
    import pyvista as pv
    import vtk
    meshes = []
    bounds = np.array(p.bounds, dtype=float).reshape(3, 2)
    diag = float(np.linalg.norm(bounds[:, 1] - bounds[:, 0])) or 1.0
    line_radius = line_radius or 6e-4 * diag
    for k, actor in enumerate(p.renderer.actors.values()):
        if not isinstance(actor, vtk.vtkActor) or not actor.GetVisibility():
            continue
        mapper = actor.GetMapper()
        if mapper is None or mapper.GetInputAsDataSet() is None:
            continue
        ds = pv.wrap(mapper.GetInputAsDataSet())
        if ds.n_points == 0:
            continue
        prop = actor.GetProperty()
        rgba, assoc = _colors(mapper, ds, prop)
        ds = ds.copy()
        if rgba is not None:
            (ds.point_data if assoc == "point" else ds.cell_data)["_rgba"] = rgba
        surf = ds if isinstance(ds, pv.PolyData) else ds.extract_surface()
        is_lines = surf.n_lines > 0 and surf.n_faces_strict == 0
        if is_lines:
            if assoc == "cell":
                surf = surf.cell_data_to_point_data()
            surf = surf.tube(radius=line_radius * max(prop.GetLineWidth(), 0.5), n_sides=10, capping=True)
            assoc = "point" if rgba is not None else None
        surf = surf.triangulate().clean() if not is_lines else surf.triangulate()
        if surf.n_cells == 0:
            continue
        faces = surf.faces.reshape(-1, 4)[:, 1:]
        pts = np.asarray(surf.points)
        if rgba is not None and assoc == "cell":
            # flat colour per cell: give every triangle its own vertices
            col = np.asarray(surf.cell_data["_rgba"]).astype(np.uint8)
            pts = pts[faces].reshape(-1, 3)
            col = np.repeat(col, 3, axis=0)
            faces = np.arange(len(pts)).reshape(-1, 3)
        elif rgba is not None:
            col = np.asarray(surf.point_data["_rgba"]).astype(np.uint8)
        else:
            c = (np.array(prop.GetColor()) * 255).round().astype(np.uint8)
            col = np.tile(np.r_[c, 255], (len(pts), 1)).astype(np.uint8)
        ext = np.ptp(pts, axis=0)
        planar = ext.min() < 1e-6 * diag and rgba is not None
        path = os.path.join(workdir, f"mesh{k:03d}.ply")
        _write_ply(path, pts, faces, col)
        if planar:
            # a thin frame around a data slice, so a pale slice still reads against white
            edge = surf.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                              non_manifold_edges=False)
            if edge.n_points:
                tube = edge.tube(radius=line_radius, n_sides=8).triangulate()
                epath = os.path.join(workdir, f"mesh{k:03d}_edge.ply")
                _write_ply(epath, tube.points, tube.faces.reshape(-1, 4)[:, 1:],
                           np.tile([130, 134, 142, 255], (tube.n_points, 1)).astype(np.uint8))
                meshes.append(dict(path=epath, kind="surface", opacity=1.0, smooth=True, lines=True))
        meshes.append(dict(path=path, kind="flat" if planar else "surface", data=rgba is not None,
                           opacity=1.0 if planar else float(prop.GetOpacity()),
                           smooth=bool(prop.GetInterpolation() != 0) and assoc != "cell" and not planar,
                           lines=bool(is_lines)))
    cam = p.camera
    w, h = p.window_size
    return dict(meshes=meshes, diag=diag, center=list(bounds.mean(axis=1)),
                camera=dict(position=list(cam.position), focal_point=list(cam.focal_point), up=list(cam.up),
                            view_angle=float(cam.view_angle), parallel=bool(cam.parallel_projection),
                            parallel_scale=float(cam.parallel_scale)),
                aspect=w / h)


def scalar_bars(p):
    """Colour bars of a pyvista Plotter: [(title, (vmin, vmax), [(r, g, b, a), ...], log)]."""
    out = []
    for title, bar in getattr(p, "scalar_bars", {}).items():
        lut = bar.GetLookupTable()
        n = lut.GetNumberOfTableValues() if hasattr(lut, "GetNumberOfTableValues") else 0
        if n == 0:
            continue
        cols = [lut.GetTableValue(i) for i in range(n)]
        log = bool(getattr(lut, "GetScale", lambda: 0)() == 1)
        if _degenerate(*lut.GetRange()):
            continue
        out.append((bar.GetTitle() or title, tuple(lut.GetRange()), cols, log))
    return out


def render(p, path, width=3000, samples=256, zoom=1.0, keep=False):
    """Render a pyvista Plotter with Blender to an RGBA PNG (transparent background). `zoom`
    multiplies the plotter's camera zoom (below 1 shows more of the scene)."""
    work = tempfile.mkdtemp(prefix="picsp_blender_")
    if zoom != 1.0:
        p.camera.zoom(zoom)
    try:
        scene = export_plotter(p, work)
        scene.update(output=os.path.abspath(path), width=int(width), height=int(round(width / scene["aspect"])),
                     samples=int(samples))
        js = os.path.join(work, "scene.json")
        with open(js, "w") as f:
            json.dump(scene, f)
        cmd = [executable(), "-b", "--factory-startup", "-noaudio", "-P", SCENE_SCRIPT, "--", js]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(path):
            raise RuntimeError("Blender render failed:\n" + r.stdout[-3000:] + r.stderr[-3000:])
    finally:
        if not keep:
            shutil.rmtree(work, ignore_errors=True)
    return scene


def project(scene, xyz):
    """Pixel position (column, row) of 3D points in the rendered image."""
    C = scene["camera"]
    pos, foc, up = (np.asarray(C[k], float) for k in ("position", "focal_point", "up"))
    fwd = (foc - pos) / np.linalg.norm(foc - pos)
    right = np.cross(fwd, up)
    right /= np.linalg.norm(right)
    upv = np.cross(right, fwd)
    d = np.atleast_2d(xyz) - pos
    xc, yc, zc = d @ right, d @ upv, d @ fwd
    W, H = scene["width"], scene["height"]
    if C["parallel"]:
        f, zc = H / (2 * C["parallel_scale"]), np.ones_like(zc)
    else:
        f = (H / 2) / np.tan(np.radians(C["view_angle"]) / 2)
    return np.c_[W / 2 + f * xc / zc, H / 2 - f * yc / zc]


def _crop(img, margin=0.02):
    """Crop to the non-transparent part; returns the image and the (row, column) offset."""
    alpha = img[..., 3] if img.shape[-1] == 4 else np.ones(img.shape[:2])
    ys, xs = np.nonzero(alpha > 0.01)
    if len(xs) == 0:
        return img, (0, 0)
    m = int(margin * max(img.shape[:2]))
    y0, y1 = max(ys.min() - m, 0), min(ys.max() + m + 1, img.shape[0])
    x0, x1 = max(xs.min() - m, 0), min(xs.max() + m + 1, img.shape[1])
    return img[y0:y1, x0:x1], (y0, x0)


def show(ax, p, png, labels=(), samples=256, zoom=1.0, width=3600):
    """Render a pyvista Plotter with Blender into matplotlib axes. `labels` are (xyz, text) pairs
    written at those 3D points in the figure font."""
    import matplotlib.pyplot as plt
    scene = render(p, png, width=width, samples=samples, zoom=zoom)
    img, (r0, c0) = _crop(plt.imread(png))
    ax.imshow(img, interpolation="lanczos")
    ax.set_axis_off()
    for xyz, text in labels:
        u, v = project(scene, xyz)[0]
        ax.text(u - c0, v - r0, text, ha="center", va="center")
    return img


def _bar_label(title):
    """'net current density (A/m²)' -> 'net current density [A/m²]', as in the 2D figures."""
    t = title.strip()
    if t.endswith(")") and "(" in t:
        i = t.rindex("(")
        t = f"{t[:i].strip()} [{t[i + 1:-1].strip()}]"
    if t.startswith("n_") and "[" not in t:
        t += " [m⁻³]"
    return t


def render_figure(p, png, width=None, samples=256, bars=None, height_mm=None, zoom=1.0, labels=()):
    """Render with Blender and lay the image out with its colour bars as a matplotlib figure
    (returned, for saving like any other figure). `bars` overrides the plotter's colour bars;
    pass [] for none."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, LogNorm, Normalize
    from . import style
    width = width or style.DOUBLE
    bars = scalar_bars(p) if bars is None else bars
    scene = render(p, png, width=3600, samples=samples, zoom=zoom)
    img, (r0, c0) = _crop(plt.imread(png))
    h_img, w_img = img.shape[:2]
    bar_w = 0.16 * width if bars else 0
    fig_h = (height_mm * style.MM) if height_mm else (width - bar_w) * h_img / w_img
    fig = plt.figure(figsize=(width, fig_h))
    fig._journal_width = width
    if bars:
        gs = fig.add_gridspec(len(bars), 2, width_ratios=[(width - bar_w), bar_w], wspace=0.02, hspace=0.35)
        ax = fig.add_subplot(gs[:, 0])
    else:
        ax = fig.add_subplot(1, 1, 1)
    ax.imshow(img, interpolation="lanczos")
    ax.set_axis_off()
    for xyz, text in labels:
        u, v = project(scene, xyz)[0]
        ax.text(u - c0, v - r0, text, ha="center", va="center")
    for i, (title, (lo, hi), cols, log) in enumerate(bars or []):
        frac = 0.6 if len(bars) == 1 else 0.8
        sub = gs[i, 1].subgridspec(3, 2, width_ratios=[0.1, 0.9],
                                   height_ratios=[(1 - frac) / 2, frac, (1 - frac) / 2])
        cax = fig.add_subplot(sub[1, 0])
        norm = LogNorm(lo, hi) if log else Normalize(lo, hi)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=ListedColormap(cols))
        cb = fig.colorbar(sm, cax=cax)
        cb.outline.set_visible(False)
        cb.set_label(_bar_label(title))
    return fig
