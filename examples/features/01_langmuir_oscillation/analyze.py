"""Langmuir oscillation: field-energy oscillation frequency vs Bohm-Gross, and the
space-time evolution of the electron density."""
import numpy as np
from _common import Run, simulate, save, plt, style, E, ME, EPS0, out


def simulate_example(d):
    simulate(d)


def analyze(d):
    r = Run(out(d))
    h = r.history()
    n, T = 1e14, 0.5
    wpe = np.sqrt(n * E * E / (EPS0 * ME))
    lam = np.sqrt(EPS0 * T / (n * E))
    L = 64 * 0.5 * lam
    k = 2 * np.pi * 2 / L
    w_bg = np.sqrt(wpe**2 + 3 * k * k * T * E / ME)
    # exact kinetic (Vlasov) root: 1 + (1 + zeta Z(zeta)) / (k lambda_D)^2 = 0, zeta = w / (sqrt(2) k v_th)
    from scipy.special import wofz
    kl = k * lam
    D = lambda w: 1 + (1 + (w / (np.sqrt(2) * kl)) * 1j * np.sqrt(np.pi) * wofz(w / (np.sqrt(2) * kl))) / kl**2
    w = complex(w_bg / wpe, -0.01)
    for _ in range(50):
        h_ = 1e-7
        w = w - D(w) / ((D(w + h_) - D(w - h_)) / (2 * h_))
    w_kin = w.real * wpe
    t, W = h["time"], h["E_field"]
    # field energy oscillates at 2w: estimate from the mean peak spacing (the wave is Landau damped,
    # so take every peak that stands above 5 % of the initial energy, not only the first ones)
    pk = [i for i in range(1, len(W) - 1) if W[i] > W[i - 1] and W[i] >= W[i + 1] and W[i] > 0.05 * W.max()]
    w_meas = np.pi / np.mean(np.diff(t[pk]))
    print(f"w/w_pe measured {w_meas / wpe:.4f}, kinetic theory {w_kin / wpe:.4f} (damping {w.imag:.4f}), "
          f"Bohm-Gross {w_bg / wpe:.4f}")
    # space-time diagram of the electron density along x, averaged over y and z (the problem is 1D,
    # and a single line of cells carries the full particle noise)
    xs = r.axes(0)[0]
    times = r.times()
    ne = np.array([r.field_array(i, "n_electron").mean(axis=(1, 2)) for i in range(len(times))])
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1.1, 1]})
    a.plot(t * wpe, W / W.max(), color=style.series_color(0))
    a.set_xlabel("time [1/ω_pe]")
    a.set_ylabel("field energy (normalised)")
    a.set_title(f"ω = {w_meas / wpe:.3f} ω_pe  (kinetic theory {w_kin / wpe:.3f}, kλ_D = {kl:.2f})")
    im = b.pcolormesh(xs * 1e3, times * wpe, ne / n - 1, cmap=style.CMAP_DIV, shading="auto",
                      vmin=-0.06, vmax=0.06, rasterized=True)
    b.set_xlabel("x [mm]")
    b.set_ylabel("time [1/ω_pe]")
    b.set_title("electron density perturbation δn/n")
    b.grid(False)
    fig.colorbar(im, ax=b, pad=0.02).outline.set_visible(False)
    fig.tight_layout()
    save(fig, d, "01_langmuir.png")
