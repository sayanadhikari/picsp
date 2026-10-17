"""Gyration and E x B drift: 3D trajectories and the measured drift velocity."""
import numpy as np
from _common import simulate, save, plt, style


def simulate_example(d):
    simulate(d)


def analyze(d):
    t = np.genfromtxt(f"{d}/output/tracks.csv", delimiter=",", names=True)
    fig = plt.figure(figsize=(12, 4.8))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2)
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
            ax.plot(x, y, z, **kw)
            ax2.plot(x, y, **kw)
            if i == 0:
                vd = np.polyfit(t["time"][m], x, 1)[0]
                print(f"{names[s]}: drift velocity {vd:.4g} m/s (E/B = 1e5 m/s)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title("helical orbits drifting with E × B")
    ax2.set_xlabel("x [m]")
    ax2.set_ylabel("y [m]")
    ax2.set_aspect("equal")
    ax2.set_title("projection perpendicular to B (E along +y)")
    ax2.legend(loc="upper right")
    fig.tight_layout()
    save(fig, d, "03_gyration_ExB.png")
