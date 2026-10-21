"""Coulomb collisions: e-i temperature relaxation and anisotropy relaxation vs theory."""
import numpy as np
from _common import Run, simulate, save, plt, style, E, ME, EPS0


def simulate_example(d):
    simulate(d)
    import os
    import subprocess
    from _common import exe
    env = dict(os.environ, OMP_NUM_THREADS="2")
    cmd = ["mpirun", "-np", "1", "--bind-to", "none", exe(), "input_aniso.toml", "--output", "output_aniso"]
    if os.geteuid() == 0:
        cmd.insert(1, "--allow-run-as-root")
    subprocess.run(cmd, cwd=d, env=env, check=True)


def spitzer(ma, mb, nb, Ta, Tb, lnL=10.0):
    e4 = E**4 / (4 * np.pi * EPS0) ** 2
    return 8 * np.sqrt(2 * np.pi) / 3 * nb * e4 * lnL * np.sqrt(ma * mb) / ((ma * Tb + mb * Ta) * E) ** 1.5


def nu_T(Tpar, Tperp, n, lnL=10.0):
    """Anisotropy relaxation rate of electrons, NRL Plasma Formulary: arctan for A > 0,
    artanh for A < 0; the bracket vanishes as (4/15) A^2 at isotropy (limit used for |A| -> 0)."""
    A = Tperp / Tpar - 1
    if abs(A) < 1e-6:
        return 2 * np.sqrt(np.pi) * E**4 / (4 * np.pi * EPS0) ** 2 * n * lnL / (np.sqrt(ME) * (Tpar * E) ** 1.5) * 4 / 15
    f = -3 + (A + 3) * (np.arctan(np.sqrt(A)) / np.sqrt(A) if A > 0 else np.arctanh(np.sqrt(-A)) / np.sqrt(-A))
    e4 = E**4 / (4 * np.pi * EPS0) ** 2
    return 2 * np.sqrt(np.pi) * e4 * n * lnL / (np.sqrt(ME) * (Tpar * E) ** 1.5) * f / A**2


def analyze(d):
    h = Run(f"{d}/output").history()
    t = h["time"]
    mi, n = 100 * ME, 1e20
    te, ti, tt = [h["T_electron"][0]], [h["T_ion"][0]], np.linspace(0, t[-1], 4000)
    for k in range(1, len(tt)):
        dt = tt[k] - tt[k - 1]
        a, b = te[-1], ti[-1]
        te.append(a + spitzer(ME, mi, n, a, b) * (b - a) * dt)
        ti.append(b + spitzer(mi, ME, n, b, a) * (a - b) * dt)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2))
    a1.plot(t * 1e9, h["T_electron"], color=style.series_color(0), label="electrons (PIC)")
    a1.plot(t * 1e9, h["T_ion"], color=style.series_color(1), label="ions (PIC)")
    a1.plot(tt * 1e9, te, "--", color=style.INK, lw=1, label="Spitzer rate equations")
    a1.plot(tt * 1e9, ti, "--", color=style.INK, lw=1)
    a1.set_xlabel("time [ns]")
    a1.set_ylabel("temperature [eV]")
    a1.set_title("electron-ion relaxation (m_i = 100 m_e)")
    a1.legend()
    try:
        g = Run(f"{d}/output_aniso").history()
        tp, tq = [g["Tz_electron"][0]], [0.5 * (g["Tx_electron"][0] + g["Ty_electron"][0])]
        ts = np.linspace(0, g["time"][-1], 4000)
        for k in range(1, len(ts)):
            dt = ts[k] - ts[k - 1]
            nu = nu_T(tp[-1], tq[-1], n)
            dperp = -nu * (tq[-1] - tp[-1]) * dt
            tq.append(tq[-1] + dperp)
            tp.append(tp[-1] - 2 * dperp)
        a2.plot(g["time"] * 1e9, g["Tz_electron"], color=style.series_color(0), label="T∥ (PIC)")
        a2.plot(g["time"] * 1e9, 0.5 * (g["Tx_electron"] + g["Ty_electron"]), color=style.series_color(1), label="T⊥ (PIC)")
        a2.plot(ts * 1e9, tp, "--", color=style.INK, lw=1, label="NRL isotropisation rate")
        a2.plot(ts * 1e9, tq, "--", color=style.INK, lw=1)
        a2.set_xlim(0, 8)
        a2.set_xlabel("time [ns]")
        a2.set_ylabel("temperature [eV]")
        a2.set_title("anisotropy relaxation by e-e collisions")
        a2.legend()
    except FileNotFoundError:
        a2.set_visible(False)
    fig.tight_layout()
    save(fig, d, "15_coulomb_relaxation.png")
