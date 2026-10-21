"""Ionisation avalanche: exponential electron multiplication and the electron energy
distribution relative to the inelastic thresholds."""
import numpy as np
from _common import Run, simulate, save, plt, style, out


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(out(d))
    h = r.history()
    t, Ne, Ni = h["time"], h["N_electron"], h["N_Ar+"]
    m = t > 0.35 * t[-1]
    nu = np.polyfit(t[m], np.log(Ne[m]), 1)[0]
    print(f"ionisation growth rate {nu:.3e} /s; electrons x{Ne[-1] / Ne[0]:.1f}; ions created {Ni[-1]:.0f}")
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(16, 4.2))
    a.semilogy(t * 1e9, Ne, color=style.series_color(0), label="electrons")
    a.semilogy(t * 1e9, np.maximum(Ni, 1), color=style.series_color(1), label="Ar⁺ created")
    a.semilogy(t[m] * 1e9, np.exp(np.polyval(np.polyfit(t[m], np.log(Ne[m]), 1), t[m])), "--", color=style.INK,
               label=f"ν_i = {nu:.2e} s⁻¹")
    a.set_xlabel("time [ns]")
    a.set_ylabel("macro-particles")
    a.set_title("Townsend avalanche, E/N = 150 Td")
    a.legend()
    b.plot(t * 1e9, 1.5 * h["T_electron"], color=style.series_color(0), label="mean thermal energy")
    b.plot(t * 1e9, 0.5 * 9.109e-31 * h["ux_electron"] ** 2 / 1.602e-19, color=style.series_color(2), label="drift energy")
    b.set_xlabel("time [ns]")
    b.set_ylabel("energy [eV]")
    b.set_title("swarm relaxes to a steady state")
    b.legend()
    pd = r.particles(-1)
    sid = np.asarray(pd.point_data["species"]).astype(int)
    en = np.asarray(pd.point_data["energy_eV"])[sid == 0]
    hist, e = np.histogram(en, bins=60, range=(0, 50), density=True)
    ec = 0.5 * (e[1:] + e[:-1])
    c.semilogy(ec, hist / np.sqrt(ec), color=style.series_color(0), label="simulation")
    Tm = 2 / 3 * en.mean()
    c.semilogy(ec, 2 / np.sqrt(np.pi) * Tm**-1.5 * np.exp(-ec / Tm), "--", color=style.INK,
               label=f"Maxwellian, same mean energy")
    for th, lab in ((11.55, "excitation"), (15.76, "ionisation")):
        c.axvline(th, color=style.INK2, ls=":", lw=1)
        c.text(th + 0.4, 0.97, lab, rotation=90, va="top", color=style.INK2, fontsize=8, transform=c.get_xaxis_transform())
    c.set_xlabel("electron energy [eV]")
    c.set_ylabel("EEPF  f(ε)/√ε  [eV⁻³ᐟ²]")
    c.set_title("electron energy probability function")
    c.legend(loc="lower left")
    fig.tight_layout()
    save(fig, d, "14_ionization_avalanche.png")
