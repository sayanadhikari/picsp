# 07 Floating sphere and OML theory

**Feature:** embedded conductor with floating potential (capacitance-matrix method), cut-cell
boundaries, open particle boundaries with exact flux injection, surface current collection.

A small sphere (R = λ_D/2) in a Maxwellian plasma charges until the electron and ion currents
balance. The floating potential is compared with orbit-motion-limited theory; the radial
sheath profiles and the potential with field lines are plotted.

## Result

Floating potential −1.51 V against −1.42 V from orbit-motion-limited theory (the sheath is not infinitely thick, so OML is only approached).

![result](../../../docs/media/07_floating_sphere.png)

Run with `python3 ../run_example.py 07_floating_sphere_OML` (from this directory) or `python3 run_example.py 07_floating_sphere_OML` from `examples/features`.
