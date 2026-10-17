# PICSP3D examples

For one small, analysed example per feature see [`features/`](features). The inputs below are
larger production-style cases.

Run from this directory, e.g.

```bash
export OMP_NUM_THREADS=2
mpirun -np 4 --bind-to none ../build/picsp3d sphere_oml.toml
python3 -m picsp_post history output_sphere_oml
```

| Input | What it shows | Size / cost | Analysis |
|---|---|---|---|
| `langmuir.toml` | Langmuir oscillation in a periodic box (quiet start, density perturbation) | 131k particles, 1200 steps, seconds | `picsp_post history` |
| `two_stream.toml` | Electrostatic two-stream instability, phase-space vortex formation | 147k particles, 1000 steps, ~10 s | `analysis/two_stream_growth.py` |
| `sphere_oml.toml` | Floating sphere (R = λ_D/2) in a Maxwellian plasma, open boundaries | 1M particles, 3000 steps, ~3 min on 2 cores | `analysis/oml_check.py` |
| `cubesat_wake.toml` | 3U CubeSat (STL bus with rails, camera, antenna + solar panels + boom) and an isolated spherical probe in a Mach-4 O+ flow, B field, motional E, photoemission, Sternglass SEE | 2.3M particles, 3500 steps, ~7 min on 2 cores | `picsp_post render/slice/animate` |
| `langmuir_probe_sweep.toml` | Swept-bias spherical Langmuir probe in a magnetised H plasma | 2.7M particles, 24000 steps (cluster) | `analysis/probe_iv.py` |
| `collisional_electrode.toml` | Negatively biased disk electrode in argon with e-Ar elastic and Ar+-Ar charge-exchange collisions, ion-induced SEE | 1.3M particles, 6000 steps (workstation) | `picsp_post slice -f n_Ar+` |

`geometry/make_cubesat_stl.py` generates the CubeSat STL (overlapping closed shells, no boolean
operations needed). Replace it with any CAD export: set `scale = 1e-3` for millimetre files.
