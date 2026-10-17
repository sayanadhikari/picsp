"""Magnetic mirror: trapped fraction vs loss-cone theory, emptying of the loss cone in
pitch-angle space, and an animation of electrons bouncing in the bottle."""
import numpy as np
from _common import Run, simulate, save, media, plt, style, gif_writer, fig_to_rgb, MU0


def loop_B(z, a=0.25, I=2000.0, zc=0.0):
    return MU0 * I * a * a / (2 * (a * a + (z - zc) ** 2) ** 1.5)


def Bz_axis(z):
    return loop_B(z, zc=-0.4) + loop_B(z, zc=0.4)


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(f"{d}/output")
    h = r.history()
    t, N = h["time"], h["N_electron"]
    Bmin, Bmax = Bz_axis(0.0), Bz_axis(0.4)
    # isotropic distribution loaded around z ~ 0: trapped fraction = cos(alpha_lc) averaged
    zs = np.linspace(-0.06, 0.06, 41)
    trapped = np.mean([np.sqrt(1 - Bz_axis(z) / Bmax) for z in zs])
    print(f"mirror ratio {Bmax / Bmin:.2f}; trapped fraction {N[-1] / N[0]:.3f}, theory {trapped:.3f}")
    ser = r.series("particles")

    def pitch(i):
        pd = r.particles(i)
        P, V = np.asarray(pd.points), np.asarray(pd.point_data["velocity"])
        vpar = V[:, 2]  # B is ~axial near the axis
        return P, V, np.degrees(np.arccos(np.clip(vpar / np.linalg.norm(V, axis=1), -1, 1)))

    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.2))
    zz = np.linspace(-0.5, 0.5, 400)
    a1.plot(zz, Bz_axis(zz) * 1e3, color=style.series_color(0))
    a1.set_xlabel("z [m]")
    a1.set_ylabel("B on axis [mT]")
    a1.set_title(f"mirror field, R = B_max/B_min = {Bmax / Bmin:.2f}")
    a2.plot(t * 1e6, N / N[0], color=style.series_color(0), label="simulation")
    a2.axhline(trapped, color=style.series_color(1), ls="--", label=f"loss-cone theory {trapped:.3f}")
    a2.set_xlabel("time [µs]")
    a2.set_ylabel("fraction of electrons retained")
    a2.set_ylim(0, 1.05)
    a2.legend()
    a2.set_title("escape through the loss cone")
    lc = np.degrees(np.arcsin(np.sqrt(Bmin / Bmax)))
    bins = np.linspace(0, 180, 46)
    for k, lab in ((0, "t = 0"), (len(ser) - 1, f"t = {ser[-1][0] * 1e6:.1f} µs")):
        _, _, pa = pitch(k)
        a3.hist(pa, bins=bins, histtype="step", lw=1.8, color=style.series_color(0 if k == 0 else 1), label=lab)
    for x in (lc, 180 - lc):
        a3.axvline(x, color=style.INK2, ls=":", lw=1)
    a3.set_xlabel("local pitch angle [deg]")
    a3.set_ylabel("macro-particles")
    a3.set_title(f"loss cone α_lc = {lc:.1f}° empties")
    a3.legend()
    fig.tight_layout()
    save(fig, d, "04_mirror.png")
    # animation: electrons in the (z, x) plane, coloured by pitch angle
    w = gif_writer(media("04_mirror_bounce.gif"), fps=10)
    for i in range(len(ser)):
        P, V, pa = pitch(i)
        f, ax = plt.subplots(figsize=(7.2, 3.2), dpi=90)
        for zc in (-0.4, 0.4):
            ax.plot([zc, zc], [0.2, 0.25], color=style.INK, lw=6, solid_capstyle="butt")
            ax.plot([zc, zc], [-0.25, -0.2], color=style.INK, lw=6, solid_capstyle="butt")
        sc = ax.scatter(P[:, 2], P[:, 0], c=np.minimum(pa, 180 - pa), cmap="viridis", vmin=0, vmax=90, s=2, lw=0)
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-0.25, 0.25)
        ax.set_xlabel("z [m]")
        ax.set_ylabel("x [m]")
        ax.set_title(f"magnetic bottle   t = {ser[i][0] * 1e6:4.2f} µs")
        cb = f.colorbar(sc, ax=ax, pad=0.01)
        cb.set_label("pitch angle [deg]")
        cb.outline.set_visible(False)
        ax.grid(False)
        f.tight_layout()
        w.append_data(fig_to_rgb(f))
        plt.close(f)
    w.close()
    print("wrote 04_mirror_bounce.gif")
