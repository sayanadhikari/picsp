# 04 Magnetic mirror

**Feature:** spatially varying magnetic fields from analytic sources (`[[fields.magnet]]`,
type `mirror` = two coaxial current loops evaluated exactly with elliptic integrals),
per-particle Boris rotation with the local B, absorbing walls.

Electrons loaded isotropically near the bottle centre either stay trapped (μ = m v⊥²/2B
conserved, reflection where B = B_min / sin²α) or escape through the loss cone. The retained
fraction approaches the loss-cone prediction and the pitch-angle distribution develops an empty
loss cone.

## Result

Trapped fraction 0.841 against the loss-cone prediction 0.839.

![result](../../../docs/media/04_mirror.png)

Run with `python3 ../run_example.py 04_magnetic_mirror` (from this directory) or `python3 run_example.py 04_magnetic_mirror` from `examples/features`.
