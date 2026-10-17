# 14 Ionisation avalanche

**Feature:** electron-impact ionisation (Lotz cross section) creating electron-ion pairs,
excitation (inelastic energy loss), elastic scattering; energy-dependent cross sections
(`model = "constant" | "step" | "lotz" | "table"`, the latter from inline tables or LXCat-style
files).

An electron swarm in argon at E/N = 150 Td reaches a steady mean energy within a few
nanoseconds and then multiplies exponentially (growth rate about 1.4e8 s⁻¹, every ionisation
creating one Ar⁺ ion and one secondary electron). The right panel shows the electron
energy-probability function together with the excitation and ionisation thresholds.

## Result

Mean energy settles within about 2 ns, then the swarm grows exponentially with ν_i ≈ 1.4 × 10⁸ s⁻¹.

![result](../../../docs/media/14_ionization_avalanche.png)

Run with `python3 ../run_example.py 14_ionization_avalanche` (from this directory) or `python3 run_example.py 14_ionization_avalanche` from `examples/features`.
