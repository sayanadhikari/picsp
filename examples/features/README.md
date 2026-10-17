# Feature examples

One small run per capability. Each directory holds `input.toml`, `analyze.py` (comparison with
theory, figures and animations written to `output/` and `docs/media/`) and a README.

```bash
python3 run_example.py                          # list
python3 run_example.py 04_magnetic_mirror       # simulate + analyse
python3 run_example.py 04_magnetic_mirror --analyze-only
python3 run_example.py all
```

Runs use one MPI rank with two OpenMP threads by default (see `_common.simulate`). The
post-processing needs the `picsp_post` package (`pip install -e ../../python`).

| # | Feature | Typical cost (2 cores) |
|---|---|---|
| 01 | Langmuir oscillation: periodic electrostatics, quiet start | seconds |
| 02 | Two-stream instability: beams, phase-space movie | ~10 s |
| 03 | Gyration and E × B drift: Boris push, test particles, tracks | seconds |
| 04 | Magnetic mirror: current-loop B sources, loss cone | ~1 min |
| 05 | Dipole trapping: dipole B, bounce and drift | ~3 min |
| 06 | Velocity distributions: all samplers | seconds |
| 07 | Floating sphere: capacitance-matrix charging vs OML | ~2 min |
| 08 | Sphere wake: mesothermal flow, wake movie | ~2 min |
| 09 | Langmuir probe I-V: swept bias | ~6 min |
| 10 | Secondary emission: Sternglass SEE on floating plates | ~5 min |
| 11 | Photoemission with shadowing | ~6 min |
| 12 | Spacecraft from STL: CubeSat with all surface physics | ~6 min |
| 13 | Ion mobility: MCC charge exchange | seconds |
| 14 | Ionisation avalanche: ionisation, excitation | seconds |
| 15 | Coulomb relaxation: Takizuka-Abe collisions | ~3 min |
