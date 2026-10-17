#!/usr/bin/env python3
"""Two-stream instability: measured growth rate of the field energy vs the cold-beam
dispersion relation  1 = wb^2/(w - k v0)^2 + wb^2/(w + k v0)^2  (mode 1 of the box)."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from picsp_post import Run, style

e, me, eps0 = 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12
run = Run(sys.argv[1] if len(sys.argv) > 1 else "output_two_stream")
h = run.history()
nb, v0, L = 5e13, 1.0e6, 0.0186                   # must match two_stream.toml
wb = np.sqrt(nb * e * e / (eps0 * me))
k = 2 * np.pi / L
a = (k * v0) ** 2
x = a + wb**2 - wb * np.sqrt(wb**2 + 4 * a)
gamma = np.sqrt(-x) if x < 0 else 0.0
t, W = h["time"], h["E_field"]
# fit the linear phase: between 1e2 and 1e-2 of saturation
Wsat = W.max()
m = (W > 30 * W[:20].mean()) & (W < 0.05 * Wsat) & (t < t[np.argmax(W)])
p = np.polyfit(t[m], np.log(W[m]), 1)
print(f"theory gamma/wb = {gamma / wb:.4f}   measured = {0.5 * p[0] / wb:.4f}   ratio = {0.5 * p[0] / gamma:.3f}")
fig, ax = plt.subplots(figsize=(6.4, 4))
ax.semilogy(t * wb, W, color=style.series_color(0), label="field energy")
ax.semilogy(t[m] * wb, np.exp(np.polyval(p, t[m])), color=style.INK, ls="--", lw=1.4,
            label=f"fit: γ = {0.5 * p[0] / wb:.3f} ω_b")
tt = t[m] * wb
ax.semilogy(tt, np.exp(p[1] + p[0] * t[m][0]) * np.exp(2 * gamma / wb * (tt - tt[0])), color=style.series_color(1),
            ls=":", lw=2, label=f"theory: γ = {gamma / wb:.3f} ω_b")
ax.set_xlabel("time [1/ω_b]")
ax.set_ylabel("field energy [J]")
ax.set_title("Two-stream instability")
ax.legend()
fig.tight_layout()
fig.savefig(sys.argv[2] if len(sys.argv) > 2 else "two_stream_growth.png")
