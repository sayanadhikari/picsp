"""Two-stream instability: growth rate vs cold-beam theory and a phase-space movie."""
import numpy as np
from _common import Run, simulate, save, media, plt, style, gif_writer, fig_to_rgb, E, ME, EPS0


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(f"{d}/output")
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
