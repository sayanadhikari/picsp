#!/usr/bin/env python3
"""Compare the floating potential and collected currents of sphere_oml.toml with
orbit-motion-limited theory for a Maxwellian plasma (Te = Ti = 1 eV, m_i = 100 m_e)."""
import sys
import numpy as np
from picsp_post import Run

e, me = 1.602176634e-19, 9.1093837015e-31
R, n, Te, Ti, mi = 3.717e-3, 1e12, 1.0, 1.0, 100 * me
h = Run(sys.argv[1] if len(sys.argv) > 1 else "output_sphere_oml").history()
late = slice(len(h) // 2, None)
V = h["V_sphere"][late].mean()
A = 4 * np.pi * R**2
Ie0 = A * n * e * np.sqrt(Te * e / (2 * np.pi * me))
Ii0 = A * n * e * np.sqrt(Ti * e / (2 * np.pi * mi))
lo, hi = -10.0, 0.0
for _ in range(200):
    mid = 0.5 * (lo + hi)
    lo, hi = (mid, hi) if Ie0 * np.exp(mid / Te) < Ii0 * (1 - mid / Ti) else (lo, mid)
print(f"floating potential: simulation {V:.3f} V, OML {mid:.3f} V")
print(f"ion current      : simulation {h['I_sphere_ion'][late].mean():.3e} A, OML at V_sim {Ii0 * (1 - V / Ti):.3e} A")
print(f"electron current : simulation {-h['I_sphere_electron'][late].mean():.3e} A, OML at V_sim {Ie0 * np.exp(V / Te):.3e} A")
