"""Ion mobility with charge exchange: drift velocity and distribution vs theory."""
import numpy as np
from _common import Run, simulate, save, plt, style, E, AMU, out


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(out(d))
    h = r.history()
    M, Ef, n, sig = 39.95 * AMU, 1000.0, 3e21, 5e-19
    lam = 1 / (n * sig)
    vd = np.sqrt(2 * E * Ef * lam / (np.pi * M))
    late = h["time"] > 0.5 * h["time"][-1]
    print(f"drift velocity {h['ux_Ar+'][late].mean():.1f} m/s, theory {vd:.1f} m/s")
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 3.7))
    a.plot(h["time"] * 1e6, h["ux_Ar+"], color=style.series_color(0), label="simulation")
    a.axhline(vd, color=style.INK, ls="--", lw=0.9, label=rf"theory, $\sqrt{{2eE\lambda/\pi M}}$ = {vd:.0f} m/s")
    a.set_xlabel("time [µs]")
    a.set_ylabel("drift velocity [m/s]")
    a.set_title(f"Ar⁺ in Ar, E = 1 kV/m, λ = {lam * 1e3:.2f} mm")
    a.legend()
    pd = r.particles(-1)
    vx = np.asarray(pd.point_data["velocity"])[:, 0]
    hist, e = np.histogram(vx, bins=80, range=(-500, 4000), density=True)
    vc = 0.5 * (e[1:] + e[:-1])
    b.plot(vc, hist, color=style.series_color(0), label="simulation")
    s2 = E * Ef * lam / M
    b.plot(vc, np.where(vc > 0, np.exp(-vc**2 / (2 * s2)) / np.sqrt(np.pi * s2 / 2), 0), "--", color=style.INK,
           label="cold-gas theory (half-Gaussian)")
    b.set_xlabel("v_x [m/s]")
    b.set_ylabel("$f(v_x)$")
    b.set_title("velocity distribution along E")
    b.legend()
    fig.tight_layout()
    save(fig, d, "13_ion_mobility.png")
