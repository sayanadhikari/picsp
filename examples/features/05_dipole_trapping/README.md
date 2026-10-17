# 05 Dipole trapping

**Feature:** magnetic dipole source (`type = "dipole"`), adiabatic invariants in a strongly
inhomogeneous field, trajectory tracking.

Energetic protons launched at the equator of a dipole (scaled so that the Larmor radius is a
few percent of the field-line radius) gyrate, bounce between magnetic mirror points at higher
latitudes and drift azimuthally (gradient and curvature drift) around the "planet".

## Result

The proton bounces between ±8° magnetic latitude and drifts azimuthally at a constant rate (about 160° in 560 µs), with no secular change of the mirror latitude.

![result](../../../docs/media/05_dipole_motion.png)

Run with `python3 ../run_example.py 05_dipole_trapping` (from this directory) or `python3 run_example.py 05_dipole_trapping` from `examples/features`.
