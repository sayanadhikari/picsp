"""Two-stream instability: growth rate vs cold-beam theory and a phase-space movie."""
import numpy as np
from _common import Run, simulate, save, media, plt, style, gif_writer, fig_to_rgb, E, ME, EPS0, out


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(out(d))
    h = r.history()
    nb, v0, L = 5e13, 1.0e6, 0.0186
    wb = np.sqrt(nb * E * E / (EPS0 * ME))
    k = 2 * np.pi / L
    a = (k * v0) ** 2
    x = a + wb**2 - wb * np.sqrt(wb**2 + 4 * a)
    gamma = np.sqrt(-x)
    t, W = h["time"], h["E_field"]
    m = (W > 30 * W[:20].mean()) & (W < 0.05 * W.max()) & (t < t[np.argmax(W)])
    p = np.polyfit(t[m], np.log(W[m]), 1)
    print(f"growth rate {0.5 * p[0] / wb:.3f} w_b, theory {gamma / wb:.3f} w_b")
    fig, ax = plt.subplots(figsize=(6.4, 4))
    ax.semilogy(t * wb, W, color=style.series_color(0), label="field energy")
    ax.semilogy(t[m] * wb, np.exp(np.polyval(p, t[m])), "--", color=style.INK, label=f"fit γ = {0.5 * p[0] / wb:.3f} ω_b")
    ax.semilogy(t[m] * wb, np.exp(np.polyval(p, t[m][0])) * np.exp(2 * gamma * (t[m] - t[m][0])), ":",
                color=style.series_color(1), lw=2, label=f"theory γ = {gamma / wb:.3f} ω_b")
    ax.set_xlabel("time [1/ω_b]")
    ax.set_ylabel("field energy [J]")
    ax.set_title("Two-stream instability: linear growth")
    ax.legend()
    fig.tight_layout()
    save(fig, d, "02_two_stream_growth.png")
    # phase-space movie
    w = gif_writer(media("02_two_stream_phase.gif"), fps=10)
    ser = r.series("particles")
    for i in range(len(ser)):
        pd = r.particles(i)
        sid = np.asarray(pd.point_data["species"]).astype(int)
        P = np.asarray(pd.points)
        V = np.asarray(pd.point_data["velocity"])
        f, ax = plt.subplots(figsize=(6.4, 4.0), dpi=90)
        for s in (0, 1):
            mm = sid == s
            ax.scatter(P[mm, 0] * 1e3, V[mm, 0] / 1e6, s=0.4, color=style.series_color(s), alpha=0.6, lw=0)
        ax.set_xlim(0, L * 1e3)
        ax.set_ylim(-2.6, 2.6)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("v_x [10⁶ m/s]")
        ax.set_title(f"two-stream phase space   t = {ser[i][0] * wb:5.1f} / ω_b")
        f.tight_layout()
        w.append_data(fig_to_rgb(f))
        plt.close(f)
    w.close()
    print("wrote 02_two_stream_phase.gif")
    # static snapshots for print: linear growth, saturation (first field-energy maximum), late time
    tb = np.array([s[0] for s in ser]) * wb
    first = next(i for i in range(1, len(W) - 1) if W[i] >= W[i - 1] and W[i] >= W[i + 1] and W[i] > 0.3 * W.max())
    t_sat = t[first] * wb
    picks = sorted({int(np.argmin(np.abs(tb - x))) for x in (0.6 * t_sat, t_sat, tb[-1])})
    fig, axs = plt.subplots(1, len(picks), figsize=(4.2 * len(picks), 3.6), sharey=True)
    for ax, i in zip(np.atleast_1d(axs), picks):
        pd = r.particles(i)
        sid = np.asarray(pd.point_data["species"]).astype(int)
        P = np.asarray(pd.points)
        V = np.asarray(pd.point_data["velocity"])
        for s in (0, 1):
            mm = sid == s
            # rasterised so that the PDF stays small with many points
            ax.scatter(P[mm, 0] * 1e3, V[mm, 0] / 1e6, s=0.3, color=style.series_color(s), alpha=0.5, lw=0,
                       rasterized=True)
        ax.set_xlim(0, L * 1e3)
        ax.set_ylim(-2.6, 2.6)
        ax.set_xlabel("x [mm]")
        ax.set_title(f"t = {tb[i]:.1f} / ω_b")
    np.atleast_1d(axs)[0].set_ylabel("v_x [10⁶ m/s]")
    fig.tight_layout()
    save(fig, d, "02_two_stream_phase.png")
