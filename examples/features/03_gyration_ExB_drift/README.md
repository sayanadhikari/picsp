# 03 Gyration and E × B drift

**Feature:** non-relativistic Boris integrator in a uniform magnetic field with an external
electric field; test-particle mode (`[fields] self_consistent = false`); trajectory output
(`track_every`, `track_count`).

Electrons and heavier positive particles gyrate in opposite directions with different Larmor
radii, while both drift with v = E × B / B² (charge and mass independent). The analysis plots
3D helices and measures the drift velocity (expected 10⁵ m/s).

## Result

Measured drifts 1.000 × 10⁵ m/s (electrons) and 1.03 × 10⁵ m/s (ions, fewer gyro-periods averaged) against 10⁵ m/s.

![result](../../../docs/media/03_gyration_ExB.png)

Run with `python3 ../run_example.py 03_gyration_ExB_drift` (from this directory) or `python3 run_example.py 03_gyration_ExB_drift` from `examples/features`.
