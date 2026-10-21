"""Floating sphere: charging history vs orbit-motion-limited theory, radial potential and
density profiles, and the potential map with field lines."""
import numpy as np
from _common import Run, simulate, save, plt, style, E, ME, out
from picsp_post import plots


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(out(d))
    h = r.history()
    R, n, Te, Ti, mi = 3.717e-3, 1e12, 1.0, 1.0, 100 * ME
    A = 4 * np.pi * R**2
    Ie0 = A * n * E * np.sqrt(Te * E / (2 * np.pi * ME))
    Ii0 = A * n * E * np.sqrt(Ti * E / (2 * np.pi * mi))
    lo, hi = -10.0, 0.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if Ie0 * np.exp(m / Te) < Ii0 * (1 - m / Ti) else (lo, m)
    V_oml = m
    late = slice(len(h) // 2, None)
    print(f"floating potential {h['V_sphere'][late].mean():.3f} V, OML {V_oml:.3f} V")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.4))
    a = axs[0]
    a.plot(h["time"] * 1e6, h["V_sphere"], color=style.series_color(0), label="simulation")
    a.axhline(V_oml, color=style.series_color(1), ls="--", label=f"OML theory {V_oml:.2f} V")
    a.set_xlabel("time [µs]")
    a.set_ylabel("sphere potential [V]")
    a.set_title("charging to the floating potential")
    a.legend()
    x, y, z = r.axes()
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    rr = np.sqrt(X**2 + Y**2 + Z**2)
    bins = np.linspace(R, 0.03, 30)
    idx = np.digitize(rr, bins)
    rc = 0.5 * (bins[1:] + bins[:-1])
    a = axs[1]
    for name, col in (("n_electron", 0), ("n_ion", 1)):
        f = r.field_array(-1, name)
        prof = [f[idx == k].mean() for k in range(1, len(bins))]
        a.plot(rc * 1e3, np.array(prof) / n, color=style.series_color(col), label=name[2:])
    phi = r.field_array(-1, "phi")
    pp = [phi[idx == k].mean() for k in range(1, len(bins))]
    a2 = a.twinx()
    a2.plot(rc * 1e3, pp, color=style.INK2, ls=":", label="φ")
    a2.set_ylabel("φ [V]", color=style.INK2)
    a2.grid(False)
    a.set_xlabel("r [mm]")
    a.set_ylabel("density / n₀")
    a.set_title("radial sheath structure")
    a.legend(loc="lower right")
    plots.plot_slice(r, "phi", "z", ax=axs[2], streamlines=True, title="potential and electric field lines")
    fig.tight_layout(w_pad=3)
    save(fig, d, "07_floating_sphere.png")
