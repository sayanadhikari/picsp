"""Secondary emission: floating potentials with and without SEE, secondary-electron cloud."""
import numpy as np
from _common import Run, simulate, save, plt, style
from picsp_post import plots


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(f"{d}/output")
    h = r.history()
    t = h["time"] * 1e9
    print(f"final potentials: no SEE {h['V_plate_noSEE'][-50:].mean():.1f} V, with SEE {h['V_plate_SEE'][-50:].mean():.2f} V")
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.3))
    a = axs[0]
    a.plot(t, h["V_plate_noSEE"], color=style.series_color(0), label="plate without SEE")
    a.plot(t, h["V_plate_SEE"], color=style.series_color(1), label="plate with SEE (δ_max = 1.6)")
    a.axhline(0, color=style.INK2, lw=0.8)
    a.set_xlabel("time [ns]")
    a.set_ylabel("floating potential [V]")
    a.set_title("charging in a 100 eV electron beam")
    a.legend()
    a = axs[1]
    sm = lambda y: np.convolve(y, np.ones(15) / 15, mode="valid")
    ts = t[7:-7]
    parts = [("beam electrons", h["I_plate_SEE_beam"]), ("plasma electrons", h["I_plate_SEE_electron"]),
             ("returning secondaries", h["I_plate_SEE_secondary"]), ("ions", h["I_plate_SEE_ion"]),
             ("emitted secondaries", h["Iemit_plate_SEE"])]
    total = sum(v for _, v in parts)
    for k, (lab, v) in enumerate(parts):
        a.plot(ts, sm(v) * 1e3, color=style.series_color(k), label=lab)
    a.plot(ts, sm(total) * 1e3, color=style.INK, lw=2, label="net current")
    a.axhline(0, color=style.INK2, lw=0.8)
    a.set_xlabel("time [ns]")
    a.set_ylabel("current to the plate [mA]")
    a.set_title("SEE plate: current balance")
    a.legend(fontsize=8, ncol=2)
    plots.plot_slice(r, "phi", "z", 0.0, -1, ax=axs[2], streamlines=True, title="potential (beam enters from the left)")
    fig.tight_layout()
    save(fig, d, "10_secondary_emission.png")
