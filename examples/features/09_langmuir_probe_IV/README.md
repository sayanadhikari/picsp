# 09 Langmuir probe I-V sweep

**Feature:** fixed-bias conductor with a time-ramped potential (`bias = "fixed"`, `ramp`),
per-object and per-species collected currents in `history.csv`.

A small spherical probe is swept from −8 V to +3 V relative to the grounded boundary. The I-V
curve shows the ion-collection branch (which keeps growing with bias for a probe smaller than
the sheath, as in orbit-motion-limited theory), the floating potential, and the exponential
electron-retardation region whose slope gives the electron temperature.

## Result

Floating potential −2.45 V; the exponential fit gives T_e = 2.21 eV for an input of 2.0 eV
(the fit window extends close to the plasma potential, where the sheath slightly flattens the
characteristic).

![result](../../../docs/media/09_probe_IV.png)
