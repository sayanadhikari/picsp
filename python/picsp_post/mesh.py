"""Import CAD geometry (STL) into PICSP3D inputs.

The solver reads STL files directly (``shape = "stl"`` in ``[[objects]]``). This module helps
to get a CAD export into a simulation:

* ``inspect`` a mesh: triangles, size, area, volume, watertightness, orientation, a guess of
  the CAD units;
* ``place`` it: units, rescaling, anchor point (centre, base, centroid, ...), rotation and
  position, with the same conventions as the solver;
* ``suggest_grid``: a domain and cell size around the placed object;
* write a ready-to-paste TOML block (``[grid]`` and ``[[objects]]``), optionally a cleaned
  binary STL in metres, and a three-view preview of the object on the grid.

Command line (see ``python -m picsp_post mesh --help``)::

    python -m picsp_post mesh info   probe.stl
    python -m picsp_post mesh import probe.stl --units mm --anchor center --name probe \\
           --cells-across 24 --margin 1.0 -o probe_object.toml --preview probe_preview.png

Python::

    from picsp_post.mesh import Mesh
    m = Mesh.load("probe.stl").scaled(1e-3)
    print(m.report())
"""
from __future__ import annotations

import math
import os
import struct
from dataclasses import dataclass

import numpy as np

UNITS = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "in": 0.0254, "inch": 0.0254, "ft": 0.3048}
ANCHORS = ("origin", "center", "centroid", "min", "bottom")


# ---------------------------------------------------------------------------------------
# reading / writing
# ---------------------------------------------------------------------------------------
def read_stl(path: str) -> np.ndarray:
    """Read an ASCII or binary STL file; returns triangles as an (N, 3, 3) float64 array."""
    with open(path, "rb") as f:
        data = f.read()
    if len(data) >= 84:
        n = struct.unpack_from("<I", data, 80)[0]
        if 84 + 50 * n == len(data):
            return _read_binary(data, n)
    head = data[:512].lstrip()
    if head[:5].lower() == b"solid" and b"facet" in data[:4096]:
        return _read_ascii(data.decode("latin-1"))
    if len(data) >= 84:
        n = struct.unpack_from("<I", data, 80)[0]
        if 84 + 50 * n <= len(data):
            return _read_binary(data, n)
    raise ValueError(f"{path}: not a valid STL file")


def _read_binary(data: bytes, n: int) -> np.ndarray:
    rec = np.dtype([("normal", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")])
    arr = np.frombuffer(data, dtype=rec, count=n, offset=84)
    return arr["v"].astype(np.float64)


def _read_ascii(text: str) -> np.ndarray:
    verts = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("vertex"):
            p = s.split()
            verts.append((float(p[1]), float(p[2]), float(p[3])))
    if not verts or len(verts) % 3:
        raise ValueError("ASCII STL: vertex count is not a multiple of 3")
    return np.asarray(verts, dtype=np.float64).reshape(-1, 3, 3)


def write_stl(path: str, tris: np.ndarray, header: str = "PICSP3D mesh (metres)") -> None:
    """Write triangles (N, 3, 3) as a binary STL file."""
    tris = np.asarray(tris, dtype=np.float64)
    nrm = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = np.divide(nrm, ln, out=np.zeros_like(nrm), where=ln > 0)
    rec = np.zeros(len(tris), dtype=np.dtype([("normal", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")]))
    rec["normal"] = nrm
    rec["v"] = tris
    with open(path, "wb") as f:
        f.write(header.encode()[:80].ljust(80, b" "))
        f.write(struct.pack("<I", len(tris)))
        f.write(rec.tobytes())


# ---------------------------------------------------------------------------------------
# geometry helpers (same conventions as the solver)
# ---------------------------------------------------------------------------------------
def euler_matrix(deg) -> np.ndarray:
    """Rotation matrix Rz @ Ry @ Rx from Euler angles in degrees (x first, then y, then z)."""
    a, b, c = np.radians(np.asarray(deg, dtype=float))
    rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    rz = np.array([[math.cos(c), -math.sin(c), 0], [math.sin(c), math.cos(c), 0], [0, 0, 1]])
    return rz @ ry @ rx


@dataclass
class MeshReport:
    triangles: int
    degenerate: int
    vertices: int
    open_edges: int
    nonmanifold_edges: int
    bbox_min: np.ndarray
    bbox_max: np.ndarray
    area: float
    volume: float          # signed: negative means inward-facing normals
    centroid: np.ndarray
    min_edge: float
    max_edge: float

    @property
    def size(self) -> np.ndarray:
        return self.bbox_max - self.bbox_min

    @property
    def watertight(self) -> bool:
        return self.open_edges == 0

    def unit_guess(self) -> str | None:
        """Suggest CAD units if the size in metres looks implausible for a plasma object."""
        L = float(self.size.max())
        if L > 20:
            return "mm"
        if L < 1e-4:
            return "m (the file may already be scaled down)"
        return None

    def text(self, name: str = "mesh", units: str = "m") -> str:
        s, lo, hi = self.size, self.bbox_min, self.bbox_max
        lines = [
            f"{name}",
            f"  triangles      {self.triangles}" + (f"  ({self.degenerate} degenerate)" if self.degenerate else ""),
            f"  vertices       {self.vertices} (after welding)",
            f"  size           {s[0]:.6g} x {s[1]:.6g} x {s[2]:.6g} {units}",
            f"  bounding box   [{lo[0]:.6g}, {lo[1]:.6g}, {lo[2]:.6g}] .. [{hi[0]:.6g}, {hi[1]:.6g}, {hi[2]:.6g}] {units}",
            f"  centroid       [{self.centroid[0]:.6g}, {self.centroid[1]:.6g}, {self.centroid[2]:.6g}] {units}",
            f"  surface area   {self.area:.6g} {units}^2",
            f"  volume         {abs(self.volume):.6g} {units}^3",
            f"  edge length    {self.min_edge:.4g} .. {self.max_edge:.4g} {units}",
        ]
        if self.watertight and self.nonmanifold_edges == 0:
            lines.append("  topology       closed and manifold")
        else:
            lines.append(f"  topology       {self.open_edges} open edges, {self.nonmanifold_edges} non-manifold edges"
                         " (fine for PICSP3D: inside tests use generalised winding numbers)")
        if self.volume < 0:
            lines.append("  orientation    normals point inwards (harmless for PICSP3D, may matter for other tools)")
        g = self.unit_guess() if units == "m" else None
        if g:
            lines.append(f"  note           {s.max():.4g} m across: is the file in {g}? use --units")
        return "\n".join(lines)


class Mesh:
    """A triangle mesh held as an (N, 3, 3) array of vertex coordinates."""

    def __init__(self, tris: np.ndarray, path: str | None = None):
        self.tris = np.asarray(tris, dtype=np.float64)
        self.path = path

    @classmethod
    def load(cls, path: str) -> "Mesh":
        return cls(read_stl(path), path)

    def save(self, path: str) -> None:
        write_stl(path, self.tris)

    def copy(self) -> "Mesh":
        return Mesh(self.tris.copy(), self.path)

    # -- transforms ---------------------------------------------------------------------
    def scaled(self, s: float) -> "Mesh":
        return Mesh(self.tris * s, self.path)

    def translated(self, t) -> "Mesh":
        return Mesh(self.tris + np.asarray(t, dtype=float), self.path)

    def rotated(self, deg, about=(0.0, 0.0, 0.0)) -> "Mesh":
        R = euler_matrix(deg)
        c = np.asarray(about, dtype=float)
        return Mesh((self.tris - c) @ R.T + c, self.path)

    def anchor_point(self, anchor: str) -> np.ndarray:
        if anchor not in ANCHORS:
            raise ValueError(f"anchor must be one of {ANCHORS}")
        lo, hi = self.bounds()
        if anchor == "origin":
            return np.zeros(3)
        if anchor == "center":
            return 0.5 * (lo + hi)
        if anchor == "min":
            return lo
        if anchor == "bottom":
            return np.array([0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), lo[2]])
        return self.report().centroid

    def placed(self, units="m", scale=1.0, fit=0.0, anchor="origin", rotate=(0, 0, 0), center=(0, 0, 0)) -> "Mesh":
        """World coordinates of the mesh as the solver places it for the equivalent input keys."""
        m = self.scaled(UNITS[units] * scale)
        if fit and fit > 0:
            m = m.scaled(fit / m.size().max())
        m = m.translated(-m.anchor_point(anchor))
        m = m.rotated(rotate)
        return m.translated(center)

    # -- properties ---------------------------------------------------------------------
    def bounds(self):
        v = self.tris.reshape(-1, 3)
        return v.min(axis=0), v.max(axis=0)

    def size(self) -> np.ndarray:
        lo, hi = self.bounds()
        return hi - lo

    def report(self) -> MeshReport:
        t = self.tris
        a, b, c = t[:, 0], t[:, 1], t[:, 2]
        cr = np.cross(b - a, c - a)
        twice_area = np.linalg.norm(cr, axis=1)
        degenerate = int(np.count_nonzero(twice_area <= 1e-300))
        area = 0.5 * twice_area.sum()
        vol_i = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6.0
        vol = vol_i.sum()
        L = float(self.size().max()) or 1.0
        if abs(vol) > 1e-9 * L**3:
            centroid = (vol_i[:, None] * (a + b + c)).sum(axis=0) / (4.0 * vol)
        else:
            w = 0.5 * twice_area
            centroid = (w[:, None] * (a + b + c) / 3.0).sum(axis=0) / max(w.sum(), 1e-300)
        e = np.concatenate([np.linalg.norm(b - a, axis=1), np.linalg.norm(c - b, axis=1), np.linalg.norm(a - c, axis=1)])
        # weld vertices on a relative tolerance and count edge uses
        q = 1e-9 * L
        keys = np.round(t.reshape(-1, 3) / q).astype(np.int64)
        _, vid = np.unique(keys, axis=0, return_inverse=True)
        vid = vid.reshape(-1, 3)
        edges = np.concatenate([vid[:, [0, 1]], vid[:, [1, 2]], vid[:, [2, 0]]])
        edges.sort(axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        lo, hi = self.bounds()
        return MeshReport(
            triangles=len(t), degenerate=degenerate, vertices=int(vid.max() + 1) if len(t) else 0,
            open_edges=int(np.count_nonzero(counts % 2)), nonmanifold_edges=int(np.count_nonzero(counts > 2)),
            bbox_min=lo, bbox_max=hi, area=float(area), volume=float(vol), centroid=centroid,
            min_edge=float(e.min()) if len(e) else 0.0, max_edge=float(e.max()) if len(e) else 0.0)

    def cleaned(self) -> "Mesh":
        """Drop degenerate triangles (zero area)."""
        t = self.tris
        ok = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) > 1e-300
        return Mesh(t[ok], self.path)


# ---------------------------------------------------------------------------------------
# grid suggestion and TOML output
# ---------------------------------------------------------------------------------------
def _round_cells(n: int) -> int:
    """Round a cell count up to a multiple of 4 (multigrid coarsening, MPI splitting)."""
    return max(8, int(math.ceil(n / 4.0)) * 4)


def suggest_grid(lo, hi, dx=None, cells_across=None, margin=1.0, debye=None):
    """Domain around the box [lo, hi]: `margin` times the largest object extent on each side.

    The cell size is `dx`, or the largest extent divided by `cells_across` (default 24),
    optionally capped at 0.5 Debye lengths when `debye` is given. Returns (origin, cells, dx).
    """
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    L = float((hi - lo).max())
    if dx is None:
        dx = L / float(cells_across or 24)
        if debye:
            dx = min(dx, 0.5 * debye)
    pad = margin * L
    span = (hi - lo) + 2 * pad
    cells = np.array([_round_cells(int(math.ceil(s / dx))) for s in span])
    centre = 0.5 * (lo + hi)
    origin = centre - 0.5 * cells * dx
    return origin, cells, float(dx)


def _fmt(v) -> str:
    return "[" + ", ".join(f"{x:.6g}" for x in v) + "]"


def object_toml(name, file, units="m", scale=1.0, fit=0.0, anchor="center", rotate=(0, 0, 0), center=(0, 0, 0),
                bias="floating", potential=0.0) -> str:
    lines = ["[[objects]]", f'name = "{name}"', 'shape = "stl"', f'file = "{file}"']
    if units != "m":
        lines.append(f'units = "{units}"')
    if scale != 1.0:
        lines.append(f"scale = {scale:.6g}")
    if fit and fit > 0:
        lines.append(f"fit = {fit:.6g}                 # largest extent [m] after import")
    lines.append(f'anchor = "{anchor}"')
    lines.append(f"center = {_fmt(center)}      # where the anchor point is placed [m]")
    if any(abs(r) > 0 for r in rotate):
        lines.append(f"rotate = {_fmt(rotate)}      # degrees about x, then y, then z (about the anchor)")
    lines.append(f'bias = "{bias}"')
    if bias != "floating" or potential:
        lines.append(f"potential = {potential:.6g}")
    return "\n".join(lines)


def grid_toml(origin, cells, dx) -> str:
    return "\n".join(["[grid]", f"cells = [{', '.join(str(int(c)) for c in cells)}]", f"dx = {dx:.6g}",
                      f"origin = {_fmt(origin)}"])


# ---------------------------------------------------------------------------------------
# preview
# ---------------------------------------------------------------------------------------
def preview(mesh: Mesh, path: str, origin=None, cells=None, dx=None, title=None) -> None:
    """Three orthographic views of the placed mesh with the domain and the cell grid."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from . import style
    style.apply()
    views = [("top (x-y)", 0, 1), ("front (x-z)", 0, 2), ("side (y-z)", 1, 2)]
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.2))
    t = mesh.tris
    for ax, (lab, i, j) in zip(axs, views):
        # shade by the normal component towards the viewer for some depth cue
        k = 3 - i - j
        n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
        ln = np.linalg.norm(n, axis=1)
        shade = np.abs(n[:, k]) / np.where(ln > 0, ln, 1)
        order = np.argsort(t[:, :, k].mean(axis=1))
        polys = t[order][:, :, [i, j]]
        cols = plt.cm.Blues(0.35 + 0.55 * shade[order])
        ax.add_collection(PolyCollection(polys, facecolors=cols, edgecolors=(0.1, 0.2, 0.35, 0.25), linewidths=0.2))
        if origin is not None:
            o, c = np.asarray(origin), np.asarray(cells)
            x0, x1 = o[i], o[i] + c[i] * dx
            y0, y1 = o[j], o[j] + c[j] * dx
            step_i = max(1, int(c[i] // 48))
            step_j = max(1, int(c[j] // 48))
            for q in range(0, int(c[i]) + 1, step_i):
                ax.plot([x0 + q * dx] * 2, [y0, y1], color="0.85", lw=0.4, zorder=0)
            for q in range(0, int(c[j]) + 1, step_j):
                ax.plot([x0, x1], [y0 + q * dx] * 2, color="0.85", lw=0.4, zorder=0)
            ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], color="0.25", lw=1.2)
            ax.set_xlim(x0 - 0.03 * (x1 - x0), x1 + 0.03 * (x1 - x0))
            ax.set_ylim(y0 - 0.03 * (y1 - y0), y1 + 0.03 * (y1 - y0))
            grid_note = f"  (lines every {step_i} cells)" if step_i > 1 else "  (cell grid)"
        else:
            ax.autoscale_view()
            grid_note = ""
        ax.set_aspect("equal")
        ax.set_xlabel(f"{'xyz'[i]} [m]")
        ax.set_ylabel(f"{'xyz'[j]} [m]")
        ax.set_title(lab + grid_note, fontsize=10)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------------------
def add_parser(sub) -> None:
    p = sub.add_parser("mesh", help="inspect an STL file or import it as a PICSP3D object")
    msub = p.add_subparsers(dest="mesh_cmd", required=True)
    pi = msub.add_parser("info", help="mesh statistics (size, area, volume, watertightness, unit guess)")
    pi.add_argument("stl")
    pi.add_argument("--units", default="m", choices=sorted(UNITS))
    pm = msub.add_parser("import", help="place an STL file and write the [grid] and [[objects]] input blocks")
    pm.add_argument("stl")
    pm.add_argument("--name", default=None, help="object name (default: file name)")
    pm.add_argument("--units", default="m", choices=sorted(UNITS), help="units of the STL coordinates")
    pm.add_argument("--scale", type=float, default=1.0, help="extra scale factor")
    pm.add_argument("--fit", type=float, default=0.0, help="rescale so that the largest extent is FIT metres")
    pm.add_argument("--anchor", default="center", choices=ANCHORS, help="mesh point placed at --at (default center)")
    pm.add_argument("--rotate", type=float, nargs=3, default=(0.0, 0.0, 0.0), metavar=("AX", "AY", "AZ"),
                    help="rotation in degrees about x, then y, then z (about the anchor)")
    pm.add_argument("--at", type=float, nargs=3, default=(0.0, 0.0, 0.0), metavar=("X", "Y", "Z"),
                    help="position of the anchor point in the domain [m]")
    pm.add_argument("--bias", default="floating", choices=("floating", "fixed"))
    pm.add_argument("--potential", type=float, default=0.0, help="fixed bias or initial potential [V]")
    g = pm.add_mutually_exclusive_group()
    g.add_argument("--dx", type=float, default=None, help="cell size [m]")
    g.add_argument("--cells-across", type=float, default=None, help="cells across the largest object extent (default 24)")
    pm.add_argument("--debye", type=float, default=None, help="Debye length [m]: cap the cell size at 0.5 lambda_D")
    pm.add_argument("--margin", type=float, default=1.0,
                    help="free space on each side, in units of the largest object extent (default 1.0)")
    pm.add_argument("--bake", metavar="OUT.stl", default=None,
                    help="write the placed mesh in metres (units, fit, anchor and rotation applied) and reference it")
    pm.add_argument("-o", "--out", default=None, help="write the TOML blocks to this file (default: print)")
    pm.add_argument("--preview", default=None, metavar="PNG", help="three-view preview of the object on the grid")
    pm.add_argument("--relative-to", default=None, metavar="DIR",
                    help="directory of the input file, for the STL path in the TOML (default: the -o directory)")


def run_cli(a) -> int:
    if a.mesh_cmd == "info":
        m = Mesh.load(a.stl)
        print(m.report().text(os.path.basename(a.stl) if a.units == "m" else f"{os.path.basename(a.stl)} (in {a.units})",
                              a.units))
        if a.units != "m":
            s = m.scaled(UNITS[a.units]).size()
            print(f"  size in metres {s[0]:.6g} x {s[1]:.6g} x {s[2]:.6g} m")
        return 0

    src = Mesh.load(a.stl)
    rep_src = src.report()
    name = a.name or os.path.splitext(os.path.basename(a.stl))[0]
    world = src.placed(a.units, a.scale, a.fit, a.anchor, a.rotate, a.at)
    lo, hi = world.bounds()
    origin, cells, dx = suggest_grid(lo, hi, a.dx, a.cells_across, a.margin, a.debye)
    base = a.relative_to or (os.path.dirname(os.path.abspath(a.out)) if a.out else os.getcwd())
    if a.bake:
        # the baked mesh holds the placed geometry in object-local metres (anchor at origin,
        # rotation applied), so the object only needs `center`
        local = src.placed(a.units, a.scale, a.fit, a.anchor, a.rotate, (0, 0, 0)).cleaned()
        local.save(a.bake)
        file_ref = os.path.relpath(os.path.abspath(a.bake), base)
        obj = object_toml(name, file_ref, anchor="origin", center=a.at, bias=a.bias, potential=a.potential)
    else:
        file_ref = os.path.relpath(os.path.abspath(a.stl), base)
        obj = object_toml(name, file_ref, a.units, a.scale, a.fit, a.anchor, a.rotate, a.at, a.bias, a.potential)
    ext = hi - lo
    head = [
        f"# generated by `picsp-post mesh import` from {os.path.basename(a.stl)}",
        f"# {rep_src.triangles} triangles, {'watertight' if rep_src.watertight else f'{rep_src.open_edges} open edges'};"
        f" placed size {ext[0]:.4g} x {ext[1]:.4g} x {ext[2]:.4g} m;"
        f" {ext.max() / dx:.1f} cells across, margin {a.margin:g} x size",
        f"# domain {cells[0] * dx:.4g} x {cells[1] * dx:.4g} x {cells[2] * dx:.4g} m; check dx against the Debye length",
    ]
    text = "\n".join(head) + "\n\n" + grid_toml(origin, cells, dx) + "\n\n" + obj + "\n"
    if a.out:
        with open(a.out, "w") as f:
            f.write(text)
        print(f"wrote {a.out}")
    else:
        print(text)
    if a.preview:
        preview(world, a.preview, origin, cells, dx, title=f"{name}: {ext.max() / dx:.1f} cells across")
        print(f"wrote {a.preview}")
    warn = rep_src.unit_guess() if a.units == "m" and a.fit <= 0 else None
    if warn:
        print(f"note: the mesh is {rep_src.size.max():.4g} units across; if the file is in {warn}, use --units")
    if ext.max() / dx < 4:
        print("note: the object is less than 4 cells across; consider a smaller --dx")
    return 0
