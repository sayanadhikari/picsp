# 01 Langmuir oscillation

**Feature:** self-consistent electrostatics in a fully periodic box (singular Poisson
problem solved with zero-mean projection), quiet-start loading with a density perturbation.

A warm electron plasma with a 5% mode-2 density perturbation oscillates at the Bohm-Gross
frequency ω² = ω_pe² + 3k²v_th². The analysis measures the oscillation frequency from the
field energy (which oscillates at 2ω) and plots the density perturbation in space and time.

```bash
python3 ../run_example.py 01_langmuir_oscillation
```

## Result

Measured frequency 1.257 ω_pe; the exact kinetic (Landau) root for this kλ_D is 1.276 ω_pe (Bohm-Gross gives a slightly different value at this wavenumber).

![result](../../../docs/media/01_langmuir.png)
