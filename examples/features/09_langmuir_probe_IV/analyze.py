"""Langmuir probe I-V characteristic and electron temperature from the retardation region."""
import numpy as np
from _common import Run, simulate, save, plt, style, out


def simulate_example(d):
    simulate(d)


def analyze(d):
    h = Run(out(d)).history()
    V = h["V_probe"]
    Ie = -h["I_probe_electron"]
    Ii = h["I_probe_ion"]
    It = Ii - Ie
    bins = np.linspace(V.min(), V.max(), 34)
    k = np.digitize(V, bins)
    sel = [i for i in range(1, len(bins)) if np.any(k == i)]
    Vb = np.array([V[k == i].mean() for i in sel])
    Ib = np.array([It[k == i].mean() for i in sel])
    Ieb = np.array([Ie[k == i].mean() for i in sel])
    vf = Vb[np.argmin(np.abs(Ib))]
    m = (Vb > vf - 1.0) & (Vb < -0.4) & (Ieb > 0)
    p = np.polyfit(Vb[m], np.log(Ieb[m]), 1)
    print(f"floating potential {vf:.2f} V; T_e from slope {1 / p[0]:.2f} eV (input 2.0 eV)")
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 3.7))
    a.plot(Vb, Ib * 1e6, "o-", ms=4, color=style.series_color(0))
    a.axhline(0, color=style.INK2, lw=0.8)
    a.axvline(0, color=style.INK2, lw=0.8, ls=":")
    a.annotate("ion collection (OML)", (Vb[1], Ib[1] * 1e6), xytext=(Vb[1] + 2.5, Ib.max() * 0.85e6),
               arrowprops=dict(arrowstyle="->", color=style.INK2))
    a.set_xlabel("probe bias [V]")
    a.set_ylabel("probe current [µA]")
    a.set_title("swept Langmuir probe: I-V characteristic")
    b.semilogy(Vb, np.clip(Ieb, 1e-12, None) * 1e6, "o", ms=4, color=style.series_color(1), label="electron current")
    b.semilogy(Vb[m], np.exp(np.polyval(p, Vb[m])) * 1e6, "-", color=style.INK,
               label=f"fit: T_e = {1 / p[0]:.2f} eV (input 2.0)")
    b.set_xlabel("probe bias [V]")
    b.set_ylabel("|I_e| [µA]")
    b.set_title("exponential retardation region")
    b.legend()
    fig.tight_layout()
    save(fig, d, "09_probe_IV.png")
