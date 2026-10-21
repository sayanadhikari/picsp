"""Gyration and E x B drift: 3D trajectories and the measured drift velocity."""
import os
import numpy as np
from _common import simulate, save, plt, style, out


def simulate_example(d):
    simulate(d)


def _blender_panel(ax, tracks, png):
    """Panel (a) for the paper: the 3D orbits rendered by Blender, with an E/B axis triad."""
    import pyvista as pv
    from picsp_post import blender
    p = pv.Plotter(off_screen=True, window_size=(1600, 1000))
    for pts, col in tracks:
        p.add_mesh(pv.lines_from_points(pts).tube(radius=0.004, n_sides=12), color=col, smooth_shading=True)
    lo = np.min([pts.min(axis=0) for pts, _ in tracks], axis=0)
    hi = np.max([pts.max(axis=0) for pts, _ in tracks], axis=0)
    o = np.array([lo[0], lo[1], hi[2]]) - 0.04 * (hi - lo).max()   # below the orbits, at the +B end
    L = 0.12 * (hi - lo).max()
    labels = []
    for k, text in ((0, "x"), (1, "E"), (2, "B")):
        e = np.eye(3)[k]
        p.add_mesh(pv.Arrow(start=o, direction=e, scale=L, tip_length=0.22, tip_radius=0.06, shaft_radius=0.018),
                   color="#3a3a3a")
        labels.append((o + 1.3 * L * e, text))
    # B (z) runs across the image and the drift (x) upwards; tilted so that E (y) and the
    # gyration loops are seen at an angle
    p.camera_position = [tuple((lo + hi) / 2 + np.array([0.5, -0.75, -0.45]) * (hi - lo).max()),
                         tuple((lo + hi) / 2), (1, 0, 0)]
    p.reset_camera()   # keep the view direction, fit everything (orbits and triad) in the frame
    blender.show(ax, p, png, labels=labels)
    p.close()


def analyze(d):
    t = np.genfromtxt(os.path.join(out(d), "tracks.csv"), delimiter=",", names=True)
    from picsp_post import blender
    paper_3d = style.PAPER and blender.available()
    if paper_3d:
        # the orbits are long along B: the 3D view across the page, the projection below
        fig = plt.figure(figsize=(12, 9))
        gs = fig.add_gridspec(2, 1, height_ratios=[1.5, 1])
        ax = fig.add_subplot(gs[0])
    else:
        fig = plt.figure(figsize=(12, 4.8))
        gs = fig.add_gridspec(1, 2)
        ax = fig.add_subplot(gs[0], projection="3d")
    ax2 = fig.add_subplot(gs[1])
    tracks = []
    names = ["electron", "heavy+ (16 mₑ)"]
    for s in (0, 1):
        for i in range(3):
            m = (t["species"] == s) & (t["index"] == i)
            x, y, z = t["x"][m], t["y"][m], t["z"][m]
            # undo periodic wrapping for display
            for arr in (x, y, z):
                jumps = np.diff(arr)
                arr[1:] -= np.cumsum(np.where(jumps > 0.2, 0.4, np.where(jumps < -0.2, -0.4, 0.0)))
            kw = dict(color=style.series_color(s), lw=0.9, alpha=0.9, label=names[s] if i == 0 else None)
            tracks.append((np.c_[x, y, z], kw["color"]))
            if not paper_3d:
                ax.plot(x, y, z, **kw)
            ax2.plot(x, y, **kw)
            if i == 0:
                vd = np.polyfit(t["time"][m], x, 1)[0]
                print(f"{names[s]}: drift velocity {vd:.4g} m/s (E/B = 1e5 m/s)")
    if paper_3d:
        _blender_panel(ax, tracks, os.path.join(out(d), "03_orbits_blender.png"))
    else:
        _axes3d(ax)
    ax.set_title("helical orbits drifting with E × B")
    ax2.set_xlabel("x [m]")
    ax2.set_ylabel("y [m]")
    ax2.set_aspect("equal")
    ax2.set_title("projection perpendicular to B (E along +y)")
    ax2.set_ylim(top=0.2)
    ax2.legend(loc="upper center", ncol=2)
    fig.tight_layout()
    save(fig, d, "03_gyration_ExB.png")


def _axes3d(ax):
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
