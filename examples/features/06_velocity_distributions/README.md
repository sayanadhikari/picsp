# 06 Velocity distributions

**Feature:** velocity-space loading: drifting Maxwellian, bi-Maxwellian (temperatures parallel
and perpendicular to B or a given axis), kappa (suprathermal power-law tails, sampled as a
multivariate Student-t), ring (perpendicular ring speed with spread) and cold beams.

The analysis histograms the loaded particles and overlays the analytic Maxwellian and kappa
distributions; directional temperatures of all species are printed.

## Result

Sampled temperatures and tails match the analytic distributions (Maxwellian T = 1.00 eV in every direction; bi-Maxwellian T⊥ = 0.50 eV, T∥ = 4.01 eV against 0.5 and 4; kappa T = 1.00 eV with the power-law tail; ring centred on v⊥ = 1.2 × 10⁶ m/s; beam 0.05 eV).

![result](../../../docs/media/06_velocity_distributions.png)

Run with `python3 ../run_example.py 06_velocity_distributions` (from this directory) or `python3 run_example.py 06_velocity_distributions` from `examples/features`.
