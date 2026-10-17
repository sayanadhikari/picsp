"""Readers for PICSP3D output directories."""
from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from functools import lru_cache

import numpy as np


class Run:
    """An output directory written by picsp3d.

    Contents: ``history.csv``, ``fields.pvd`` (+ .pvti), ``particles.pvd`` (+ .pvtp),
    ``surface.pvd`` (+ .pvtp).
    """

    def __init__(self, path: str):
        self.path = os.path.abspath(path)
        if not os.path.isdir(self.path):
            raise FileNotFoundError(f"no such output directory: {path}")

    # ------------------------------------------------------------------ history
    def history(self):
        """Time history as a numpy structured array (columns from the CSV header)."""
        f = os.path.join(self.path, "history.csv")
        if not os.path.exists(f):
            raise FileNotFoundError(f"{f} not found (set output.history_every > 0)")
        return np.genfromtxt(f, delimiter=",", names=True, deletechars="")

    # ------------------------------------------------------------------ series
    def series(self, name: str):
        """List of (time, absolute file path) from ``<name>.pvd``."""
        f = os.path.join(self.path, name + ".pvd")
        if not os.path.exists(f):
            return []
        root = ET.parse(f).getroot()
        out = []
        for ds in root.iter("DataSet"):
            out.append((float(ds.get("timestep")), os.path.join(self.path, ds.get("file"))))
        return out

    def times(self, name: str = "fields"):
        return np.array([t for t, _ in self.series(name)])

    def _pick(self, name, index):
        s = self.series(name)
        if not s:
            raise FileNotFoundError(f"no '{name}' output in {self.path}")
        return s[index]

    @staticmethod
    @lru_cache(maxsize=8)
    def _read(path):
        import pyvista as pv
        if path.endswith(".pvtp"):
            from vtkmodules.vtkIOXML import vtkXMLPPolyDataReader
            r = vtkXMLPPolyDataReader()
            r.SetFileName(path)
            r.Update()
            return pv.wrap(r.GetOutput())
        return pv.read(path)

    def fields(self, index: int = -1):
        """pyvista.ImageData of field dump `index` (point data: phi, rho, E, n_*, T_*, u_*, object_id)."""
        return self._read(self._pick("fields", index)[1])

    def particles(self, index: int = -1):
        """pyvista.PolyData of sampled particles (point data: velocity, species, weight, energy_eV)."""
        return self._read(self._pick("particles", index)[1])

    def surface(self, index: int = -1):
        """pyvista.PolyData of object surface panels (cell data: potential, J_*, J_net, J_emit, sun_cos...)."""
        return self._read(self._pick("surface", index)[1])

    def time(self, name="fields", index=-1):
        return self._pick(name, index)[0]

    # ------------------------------------------------------------------ arrays
    def field_array(self, index: int, name: str):
        """Field as numpy array of shape (nx+1, ny+1, nz+1[, 3])."""
        g = self.fields(index)
        dims = g.dimensions
        a = np.asarray(g.point_data[name])
        if a.ndim == 2:
            return a.reshape(dims[2], dims[1], dims[0], a.shape[1]).transpose(2, 1, 0, 3)
        return a.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)

    def axes(self, index: int = -1):
        """1D node coordinate arrays (x, y, z)."""
        g = self.fields(index)
        o, s, d = g.origin, g.spacing, g.dimensions
        return tuple(o[a] + s[a] * np.arange(d[a]) for a in range(3))

    def species_names(self):
        g = self.fields(0)
        return [n[2:] for n in g.point_data.keys() if n.startswith("n_")]

    def info(self) -> str:
        lines = [f"PICSP3D run: {self.path}"]
        for s in ("fields", "particles", "surface"):
            ser = self.series(s)
            if ser:
                lines.append(f"  {s:10s}: {len(ser):4d} dumps, t = {ser[0][0]:.4e} .. {ser[-1][0]:.4e} s")
        if os.path.exists(os.path.join(self.path, "history.csv")):
            h = self.history()
            lines.append(f"  history   : {len(h)} rows, columns: {', '.join(h.dtype.names)}")
        if self.series("fields"):
            g = self.fields(0)
            lines.append(f"  grid      : {tuple(int(d) - 1 for d in g.dimensions)} cells, spacing {tuple(g.spacing)}")
            lines.append(f"  arrays    : {', '.join(g.point_data.keys())}")
        return "\n".join(lines)
