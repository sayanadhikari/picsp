# 13 Ion mobility with charge exchange (MCC)

**Feature:** Monte Carlo collisions with a neutral background (null-collision method), resonant
charge exchange, external electric field in test-particle mode.

Ar⁺ ions accelerate in a uniform field and lose their momentum in charge-exchange collisions.
For a constant mean free path and cold gas, the distribution along E is a half-Gaussian and the
drift velocity is √(2eEλ/πM); both are reproduced.

## Result

Drift velocity 1017 m/s against 1012 m/s; the distribution along E matches the half-Gaussian.

![result](../../../docs/media/13_ion_mobility.png)

Run with `python3 ../run_example.py 13_MCC_ion_mobility` (from this directory) or `python3 run_example.py 13_MCC_ion_mobility` from `examples/features`.
