# Scaling benchmark

Strong and weak scaling of PICSP3D on one (multi-socket) node, comparing one MPI rank per core,
one rank per socket with OpenMP threads, and a single threaded rank.

```bash
python3 run_scaling.py --quick                       # about a minute: checks the setup
python3 run_scaling.py                               # full sweep, writes results/<cpu-model>/
python3 run_scaling.py --mpirun /usr/bin/mpirun      # pick the MPI launcher explicitly
python3 plot_scaling.py results/<cpu-model>          # figures into docs/media/
```

`python3 run_scaling.py --help` lists the options (core counts, problem size, repeats). Results
for a 4 × 18-core Xeon Platinum 8354H are in `results/xeon-platinum-8354h-4s`; the discussion is in
[docs/performance.md](../../docs/performance.md). Contributions of results from other machines are
welcome.
