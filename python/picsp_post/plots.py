"""2D plots (matplotlib): time histories, field slices, profiles, phase space, VDFs."""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LogNorm

from . import style
from .io import Run

UNITS = {"phi": "V", "rho": "C/m³", "E": "V/m", "object_id": ""}
AXN = {"x": 0, "y": 1, "z": 2}


def _units(name):
    if name in UNITS:
        return UNITS[name]
    if name.startswith("n_"):
        return "m⁻³"
    if name.startswith("T_"):
        return "eV"
    if name.startswith("u_"):
        return "m/s"
    return ""


def _label(name, comp=None):
    u = _units(name)
    pretty = {"phi": "φ", "rho": "ρ"}.get(name, name)
    if comp is not None:
        pretty += "xyz"[comp] if name in ("E",) or name.startswith("u_") else f"[{comp}]"
    elif name == "E" or name.startswith("u_"):
        pretty = f"|{pretty}|"
    return f"{pretty} [{u}]" if u else pretty


def _diverging(name, comp):
    return name in ("phi", "rho") or (comp is not None)


def _time_scale(t):
    tm = np.max(np.abs(t)) if len(t) else 1.0
    for f, u in ((1e-9, "ns"), (1e-6, "µs"), (1e-3, "ms"), (1.0, "s")):
        if tm < 1000 * f:
            return 1 / f, u
    return 1.0, "s"


# ---------------------------------------------------------------------------- history
def plot_history(run: Run, out: str | None = None):
    """Multi-panel overview of history.csv: energies, object potentials, currents, particle counts."""
    h = run.history()
    names = h.dtype.names
    s, tu = _time_scale(h["time"])
    t = h["time"] * s
    species = [n[2:] for n in names if n.startswith("N_")]
    conds = [n[2:] for n in names if n.startswith("V_")]
    objs = sorted({n.split("_")[1] for n in names if n.startswith("I_") and not n.startswith("Iemit")})
    panels = ["energy", "count"] + (["potential"] if conds else []) + (["current"] if objs else [])
    n = len(panels)
    ncol = 2
    nrow = (n + 1) // 2
    fig, axs = plt.subplots(nrow, ncol, figsize=(11, 3.4 * nrow), squeeze=False)
    axs = axs.ravel()
    for ax in axs[n:]:
        ax.set_visible(False)
    for ax, p in zip(axs, panels):
        if p == "energy":
            for i, sp in enumerate(species):
                ax.plot(t, h["KE_" + sp], color=style.series_color(i), label=f"kinetic, {sp}")
            ax.plot(t, h["E_field"], color=style.INK2, ls="--", label="field")
            ax.plot(t, h["E_total"], color=style.INK, lw=2.2, label="total")
            ax.set_ylabel("energy [J]")
            ax.set_title("Energy")
            ax.set_yscale("log")
        elif p == "count":
            for i, sp in enumerate(species):
                ax.plot(t, h["N_" + sp], color=style.series_color(i), label=sp)
            ax.set_ylabel("macro-particles")
            ax.set_title("Particle count")
        elif p == "potential":
            for i, c in enumerate(conds):
                ax.plot(t, h["V_" + c], color=style.series_color(i), label=c)
            ax.set_ylabel("potential [V]")
            ax.set_title("Object potentials")
        elif p == "current":
            k = 0
            for o in objs:
                for sp in species:
                    col = f"I_{o}_{sp}"
                    if col in names:
                        y = h[col]
                        ax.plot(t, _smooth(y), color=style.series_color(k), label=f"{o}: {sp}")
                        k += 1
                if f"Iemit_{o}" in names and np.any(h[f"Iemit_{o}"] != 0):
                    ax.plot(t, _smooth(h[f"Iemit_{o}"]), color=style.series_color(k), ls="--", label=f"{o}: emitted")
                    k += 1
            ax.axhline(0, color=style.INK2, lw=0.8)
            ax.set_ylabel("current [A]")
            ax.set_title("Collected currents (smoothed)")
        ax.set_xlabel(f"time [{tu}]")
        ax.legend(loc="best")
    fig.tight_layout()
    if out:
        fig.savefig(out)
        plt.close(fig)
    return fig


def _smooth(y, n=None):
    y = np.asarray(y, dtype=float)
    if n is None:
        n = max(1, len(y) // 50)
    if n <= 1:
        return y
    k = np.ones(n) / n
    return np.convolve(y, k, mode="same")


# ---------------------------------------------------------------------------- slices
def slice_data(run: Run, index: int, name: str, axis: str = "z", pos: float | None = None, comp=None):
    """Return (x1, x2, values2d, labels) for a plane perpendicular to `axis` through `pos` (default centre)."""
    a = AXN[axis]
    ax = run.axes(index)
    arr = run.field_array(index, name)
    if arr.ndim == 4:
        arr = arr[..., comp] if comp is not None else np.linalg.norm(arr, axis=-1)
    c = ax[a]
    k = len(c) // 2 if pos is None else int(np.clip(np.round((pos - c[0]) / (c[1] - c[0])), 0, len(c) - 1))
    sl = np.take(arr, k, axis=a)
    rest = [i for i in range(3) if i != a]
    return ax[rest[0]], ax[rest[1]], sl, ("xyz"[rest[0]], "xyz"[rest[1]]), c[k]


def plot_slice(run: Run, name: str = "phi", axis: str = "z", pos=None, index: int = -1, comp=None, ax=None,
               vmin=None, vmax=None, log=False, show_objects=True, streamlines=False, cmap=None, title=None):
    """Colour map of a field on an axis-aligned plane, with object outlines and optional E streamlines."""
    x1, x2, v, (l1, l2), p = slice_data(run, index, name, axis, pos, comp)
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
    else:
        fig = ax.figure
    div = _diverging(name, comp) and not log
    if cmap is None:
        cmap = style.CMAP_DIV if div else (style.CMAP_HOT if name.startswith("T_") else style.CMAP_SEQ)
    if div:
        lim = max(abs(np.nanmin(v)) if vmin is None else abs(vmin), abs(np.nanmax(v)) if vmax is None else abs(vmax), 1e-30)
        lo = -lim if vmin is None else vmin
        hi = lim if vmax is None else vmax
        if lo < 0 < hi:
            norm = TwoSlopeNorm(vcenter=0.0, vmin=lo, vmax=hi)
        else:
            norm = plt.Normalize(lo, hi)
    elif log:
        vv = v[v > 0]
        hi_ = vmax or (vv.max() if vv.size else 10)
        # default floor: 4 decades below the maximum; zeros are clipped to it (no masked holes)
        lo_ = vmin or (max(vv.min(), hi_ * 1e-4) if vv.size else 1)
        norm = LogNorm(vmin=lo_, vmax=hi_)
        v = np.clip(np.nan_to_num(v, nan=lo_), lo_, None)
    else:
        norm = plt.Normalize(vmin if vmin is not None else np.nanmin(v), vmax if vmax is not None else np.nanmax(v))
    im = ax.pcolormesh(x1, x2, v.T, shading="gouraud", cmap=cmap, norm=norm, rasterized=True)
    cb = fig.colorbar(im, ax=ax, shrink=0.88, pad=0.02)
    cb.set_label(_label(name, comp))
    cb.outline.set_visible(False)
    if show_objects:
        _draw_objects(run, index, axis, pos, x1, x2, ax)
    if streamlines:
        Ev = run.field_array(index, "E")
        a = AXN[axis]
        rest = [i for i in range(3) if i != a]
        k = np.argmin(np.abs(run.axes(index)[a] - p))
        e1 = np.take(Ev[..., rest[0]], k, axis=a)
        e2 = np.take(Ev[..., rest[1]], k, axis=a)
        # line width follows |E|, so the noise-dominated weak-field region stays in the background
        # and lines stop where the field is too weak to be distinguished from particle noise
        mag = np.hypot(e1, e2).T
        ref = np.percentile(mag, 98) or 1.0
        lw = 0.15 + 0.85 * np.clip(mag / ref, 0, 1)
        weak = mag < 0.04 * ref
        ax.streamplot(x1, x2, np.ma.masked_where(weak, e1.T), np.ma.masked_where(weak, e2.T), color=(0.15, 0.15, 0.15, 0.6), density=1.0, linewidth=lw, arrowsize=0.5)
    ax.set_aspect("equal")
    ax.set_xlabel(f"{l1} [m]")
    ax.set_ylabel(f"{l2} [m]")
    ax.grid(False)
    t = run.time("fields", index)
    ax.set_title(title or f"{_label(name, comp).split(' [')[0]}   {axis} = {p:.3g} m   t = {t:.3e} s")
    if own:
        fig.tight_layout()
    return fig, ax, im


def _draw_objects(run, index, axis, pos, x1, x2, ax):
    """Fill object cross-sections: smooth outline from the signed distance field when
    available, otherwise from the object-node mask."""
    g = run.fields(index)
    if "sdf" in g.point_data:
        _, _, sd, _, _ = slice_data(run, index, "sdf", axis, pos)
        if np.nanmin(sd) < 0:
            ax.contourf(x1, x2, sd.T, levels=[-1e300, 0.0], colors=["#3b3b3b"], alpha=0.9)
            ax.contour(x1, x2, sd.T, levels=[0.0], colors=["#111111"], linewidths=1.0)
    elif "object_id" in g.point_data:
        _, _, ob, _, _ = slice_data(run, index, "object_id", axis, pos)
        if np.any(ob >= 0):
            m = (ob >= 0).T.astype(float)
            ax.contourf(x1, x2, m, levels=[0.5, 1.5], colors=["#3b3b3b"], alpha=0.85)
            ax.contour(x1, x2, m, levels=[0.5], colors=["#111111"], linewidths=1.0)


def plot_profile(run: Run, name: str = "phi", axis: str = "x", through=None, index: int = -1, ax=None, comp=None, **kw):
    """1D profile along `axis` through point `through` (default: domain centre)."""
    a = AXN[axis]
    axes = run.axes(index)
    arr = run.field_array(index, name)
    if arr.ndim == 4:
        arr = arr[..., comp] if comp is not None else np.linalg.norm(arr, axis=-1)
    idx = []
    for i in range(3):
        if i == a:
            idx.append(slice(None))
        else:
            c = axes[i]
            idx.append(len(c) // 2 if through is None else int(np.argmin(np.abs(c - through[i]))))
    y = arr[tuple(idx)]
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot(axes[a], y, **kw)
    ax.set_xlabel(f"{axis} [m]")
    ax.set_ylabel(_label(name, comp))
    return ax


# ---------------------------------------------------------------------------- particles
def _species_mask(run, pd, species):
    sid = np.asarray(pd.point_data["species"]).astype(int)
    if species is None:
        return np.ones_like(sid, dtype=bool)
    if isinstance(species, str):
        names = run.species_names()
        species = names.index(species)
    return sid == species


def plot_phase(run: Run, species=0, x: str = "x", v: str = "vx", index: int = -1, bins=200, ax=None):
    """Phase-space density (2D histogram) from the sampled particles."""
    pd = run.particles(index)
    m = _species_mask(run, pd, species)
    pts = np.asarray(pd.points)[m]
    vel = np.asarray(pd.point_data["velocity"])[m]
    w = np.asarray(pd.point_data["weight"])[m]
    X = pts[:, AXN[x]]
    V = vel[:, AXN[v[-1]]]
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
    H, xe, ye = np.histogram2d(X, V, bins=bins, weights=w)
    H = np.ma.masked_equal(H, 0)
    im = ax.pcolormesh(xe, ye, H.T, cmap=style.CMAP_HOT, norm=LogNorm(), rasterized=True)
    cb = ax.figure.colorbar(im, ax=ax, pad=0.02)
    cb.set_label("f (weighted counts)")
    cb.outline.set_visible(False)
    ax.set_xlabel(f"{x} [m]")
    ax.set_ylabel(f"{v} [m/s]")
    ax.grid(False)
    ax.set_title(f"Phase space ({_name(run, species)}), t = {run.time('particles', index):.3e} s")
    return ax


def plot_vdf(run: Run, species=0, index: int = -1, bins=120, ax=None, log=True):
    """Velocity distributions f(v_x), f(v_y), f(v_z) and the speed distribution."""
    pd = run.particles(index)
    m = _species_mask(run, pd, species)
    vel = np.asarray(pd.point_data["velocity"])[m]
    w = np.asarray(pd.point_data["weight"])[m]
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
    vmax = np.percentile(np.abs(vel), 99.9)
    edges = np.linspace(-vmax, vmax, bins + 1)
    for c in range(3):
        f, _ = np.histogram(vel[:, c], bins=edges, weights=w, density=True)
        ax.plot(0.5 * (edges[1:] + edges[:-1]), f, color=style.series_color(c), label=f"v{'xyz'[c]}")
    if log:
        ax.set_yscale("log")
    ax.set_xlabel("velocity [m/s]")
    ax.set_ylabel("f(v) [s/m]")
    ax.legend()
    ax.set_title(f"Velocity distribution ({_name(run, species)})")
    return ax


def _name(run, species):
    if isinstance(species, str):
        return species
    try:
        return run.species_names()[species]
    except Exception:
        return f"species {species}"
