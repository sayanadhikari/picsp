"""Sampled velocity distributions vs analytic forms."""
import numpy as np
from scipy.special import gamma as G
from _common import Run, simulate, save, plt, style, E, ME


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(f"{d}/output")
    pd = r.particles(0)
    sid = np.asarray(pd.point_data["species"]).astype(int)
    V = np.asarray(pd.point_data["velocity"])
    vt = np.sqrt(E / ME)  # 1 eV
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.0))
    # (a) Maxwellian vs kappa: 1D f(v_x) on a log scale
    a = axs[0]
    for s, lab, col in ((0, "drifting Maxwellian", 0), (2, "kappa, κ = 2.5", 1)):
        v = V[sid == s, 0]
        h, e = np.histogram(v / 1e6, bins=120, range=(-3.5, 3.5), density=True)
        a.semilogy(0.5 * (e[1:] + e[:-1]), h, color=style.series_color(col), lw=1.6, label=lab)
    x = np.linspace(-3.5, 3.5, 400) * 1e6
    a.semilogy(x / 1e6, 1e6 * np.exp(-(x - 3e5) ** 2 / (2 * vt**2)) / np.sqrt(2 * np.pi * vt**2), "--",
               color=style.INK, lw=1, label="analytic")
    k = 2.5
    th2 = (2 * k - 3) / k * vt**2
    fk = G(k) / (np.sqrt(np.pi * k * th2) * G(k - 0.5)) * (1 + x * x / (k * th2)) ** (-k)  # 1D marginal
    a.semilogy(x / 1e6, 1e6 * fk, "--", color=style.INK, lw=1)
    a.set_ylim(1e-4, 2)
    a.set_xlabel("v_x [10⁶ m/s]")
    a.set_ylabel("f(v_x)")
    a.set_title("Maxwellian vs kappa tails")
    a.legend(fontsize=8)
    # (b) bi-Maxwellian, (c) ring: (v_perp, v_par) histograms
    for ax, s, title in ((axs[1], 1, "bi-Maxwellian  T∥ = 4 eV, T⊥ = 0.5 eV"), (axs[2], 3, "ring distribution")):
        v = V[sid == s]
        vperp = np.hypot(v[:, 0], v[:, 1])
        ax.hist2d(v[:, 2] / 1e6, vperp / 1e6, bins=80, cmap=style.CMAP_HOT, range=[[-3, 3], [0, 2]])
        ax.set_xlabel("v∥ [10⁶ m/s]")
        ax.set_ylabel("v⊥ [10⁶ m/s]")
        ax.set_title(title)
        ax.grid(False)
    # (d) beam: v_z distribution and moments of all species
    ax = axs[3]
    names = ["Maxwellian", "bi-Maxwellian", "kappa", "ring", "beam"]
    rows = []
    for s in range(5):
        v = V[sid == s]
        u = v.mean(0)
        T = ME * ((v - u) ** 2).mean(0) / E
        rows.append(f"{names[s]:<13s} T = ({T[0]:.2f}, {T[1]:.2f}, {T[2]:.2f}) eV")
    h, e = np.histogram(V[sid == 4, 2] / 1e6, bins=80, density=True)
    ax.plot(0.5 * (e[1:] + e[:-1]), h, color=style.series_color(2))
    ax.set_xlabel("v_z [10⁶ m/s]")
    ax.set_ylabel("f(v_z)")
    ax.set_title("cold beam (1.5×10⁶ m/s, 0.05 eV)")
    print("\n".join(rows))
    fig.tight_layout()
    save(fig, d, "06_velocity_distributions.png")
