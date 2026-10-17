# 15 Coulomb collisions

**Feature:** binary Coulomb collisions (Takizuka-Abe pairing within cells, Nanbu-Yonemura weight
correction), intra- and inter-species, constant or automatic Coulomb logarithm, collision
sub-cycling (`every`).

Two runs: hot electrons and cold ions relax towards a common temperature (compared with the
Spitzer energy-exchange equations), and an anisotropic bi-Maxwellian electron population
isotropises (compared with the NRL anisotropy-relaxation rate).

## Result

Both relaxations follow the theoretical rates; the binary method is about 10% slower than Spitzer at this time step, as expected.

![result](../../../docs/media/15_coulomb_relaxation.png)

Run with `python3 ../run_example.py 15_coulomb_relaxation` (from this directory) or `python3 run_example.py 15_coulomb_relaxation` from `examples/features`.
