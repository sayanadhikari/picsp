"""Dipole trapping: 3D trajectories with field lines (pyvista) and the bounce/drift motion."""
import numpy as np
from _common import simulate, save, media, plt, style


def simulate_example(d):
    simulate(d)


def field_line(L, n=200):
    lam = np.linspace(-np.arccos(np.sqrt(1.0 / L)), np.arccos(np.sqrt(1.0 / L)), n)  # down to r = 1
    r = L * np.cos(lam) ** 2
    return r * np.cos(lam), r * np.sin(lam)


def analyze(d):
    t = np.genfromtxt(f"{d}/output/tracks.csv", delimiter=",", names=True)
    import pyvista as pv
    p = pv.Plotter(off_screen=True, window_size=(1600, 1000))
    p.set_background("#0d1117", top="#1f2a38")
    p.add_mesh(pv.Sphere(radius=1.0, theta_resolution=80, phi_resolution=80), color="#3d6fa8", smooth_shading=True,
               specular=0.3)
    for phi in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        for L in (2.0, 3.0, 4.0):
            x, z = field_line(L)
            pts = np.c_[x * np.cos(phi), x * np.sin(phi), z]
            p.add_mesh(pv.lines_from_points(pts), color="#8b93a1", opacity=0.35, line_width=1)
    cols = [style.series_color(k) for k in (1, 2, 0)]
    for i in range(3):
        m = (t["species"] == 0) & (t["index"] == i)
        pts = np.c_[t["x"][m], t["y"][m], t["z"][m]]
        if len(pts) > 2:
            p.add_mesh(pv.lines_from_points(pts).tube(radius=0.009), color=cols[i])
    p.camera_position = [(7.5, -10.5, 7.5), (0, 0, 0), (0, 0, 1)]
    p.add_text("protons trapped in a dipole field: gyration, bounce and drift", position="upper_left", font_size=12,
               color="#f0efec")
    p.screenshot(media("05_dipole_trapping.png"))
    p.close()
    # latitude and longitude vs time for the first particle
    m = (t["species"] == 0) & (t["index"] == 0)
    x, y, z, tt = t["x"][m], t["y"][m], t["z"][m], t["time"][m]
    r = np.sqrt(x * x + y * y + z * z)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 3.8))
    a.plot(tt * 1e6, np.degrees(np.arcsin(z / r)), color=style.series_color(0))
    a.set_xlabel("time [µs]")
    a.set_ylabel("magnetic latitude [deg]")
    a.set_title("bounce motion between mirror points")
    b.plot(tt * 1e6, np.degrees(np.unwrap(np.arctan2(y, x))), color=style.series_color(1))
    b.set_xlabel("time [µs]")
    b.set_ylabel("longitude [deg]")
    b.set_title("gradient-curvature drift around the dipole")
    fig.tight_layout()
    save(fig, d, "05_dipole_motion.png")
