#!/usr/bin/env python3
"""Langmuir probe I-V characteristic from a swept-bias run (langmuir_probe_sweep.toml).
Electron temperature from the slope of ln(I_e) in the retardation region, compared with
the input temperature; ion saturation current from the most negative biases."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from picsp_post import Run, style

run = Run(sys.argv[1] if len(sys.argv) > 1 else "output_probe")
obj = sys.argv[2] if len(sys.argv) > 2 else "probe"
h = run.history()
names = h.dtype.names
V = h[f"V_{obj}"]
Itot = sum(h[c] for c in names if c.startswith(f"I_{obj}_")) + h[f"Iemit_{obj}"]
Ie = sum(h[c] for c in names if c.startswith(f"I_{obj}_") and "electron" in c)
# bin by bias to beat particle noise
bins = np.linspace(V.min(), V.max(), 49)
idx = np.digitize(V, bins)
Vb = np.array([V[idx == i].mean() for i in range(1, len(bins)) if np.any(idx == i)])
Ib = np.array([Itot[idx == i].mean() for i in range(1, len(bins)) if np.any(idx == i)])
Ieb = np.array([-Ie[idx == i].mean() for i in range(1, len(bins)) if np.any(idx == i)])
Ii_sat = -np.mean(Ib[: max(3, len(Ib) // 8)])
# retardation region: between floating potential and plasma potential (0 V)
vf = Vb[np.argmin(np.abs(Ib))]
m = (Vb > vf) & (Vb < -0.3) & (Ieb > 0)
Te = 1.0 / np.polyfit(Vb[m], np.log(Ieb[m]), 1)[0] if m.sum() > 3 else float("nan")
print(f"floating potential ~ {vf:.2f} V, ion saturation current ~ {Ii_sat:.3e} A, T_e from slope = {Te:.2f} eV")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
a1.plot(Vb, Ib * 1e6, "o-", color=style.series_color(0), ms=4)
a1.axhline(0, color=style.INK2, lw=0.8)
a1.set_xlabel("probe bias [V]")
a1.set_ylabel("collected current [µA]")
a1.set_title("I-V characteristic")
a2.semilogy(Vb, np.clip(Ieb, 1e-12, None), "o", color=style.series_color(1), ms=4, label="electron current")
if m.sum() > 3:
    p = np.polyfit(Vb[m], np.log(Ieb[m]), 1)
    a2.semilogy(Vb[m], np.exp(np.polyval(p, Vb[m])), color=style.INK, ls="--", label=f"fit: T_e = {Te:.2f} eV")
a2.set_xlabel("probe bias [V]")
a2.set_ylabel("|I_e| [A]")
a2.set_title("Electron retardation region")
a2.legend()
fig.tight_layout()
fig.savefig(sys.argv[3] if len(sys.argv) > 3 else "probe_iv.png")
