#!/usr/bin/env python3
"""Strong and weak scaling benchmark for PICSP3D on a single (multi-socket) node.

    python3 run_scaling.py                        # strong + weak, all modes, results/<cpu-model>/
    python3 run_scaling.py --kind strong --modes mpi hybrid --max-cores 36
    python3 run_scaling.py --quick                # small grid, short runs: checks the setup in ~1 min
    python3 plot_scaling.py results/<cpu-model>    # figures

The test problem is a periodic, thermal hydrogen plasma (electrons 1 eV, ions 0.1 eV, 16 particles
per cell and species) so the step mixes particle work (push, deposit, sort, migration) and the
multigrid Poisson solve the way a production run does. Three ways of using the cores are compared:

    mpi      one MPI rank per core, one thread each
    hybrid   one MPI rank per socket (NUMA domain), one OpenMP thread per core
    openmp   a single rank, one thread per core (across sockets)

Open MPI launch flags are used to pin ranks and threads to physical cores; with another MPI the
runs still work but are not pinned. Results go to results/<cpu-model>/{strong,weak}.csv together
with a system.json describing the machine, build and commit.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

INPUT = """\
title = "scaling benchmark"
seed = 1

[grid]
cells = [{nx}, {ny}, {nz}]
dx_debye = 0.5

[time]
dt_wpe = 0.1
steps = {steps}

[boundaries]
default = {{ field = "periodic", particle = "periodic" }}
{parallel}
[[species]]
name = "electron"
charge = -1
mass_me = 1
density = 1e14
temperature = 1.0
ppc = {ppc}

[[species]]
name = "ion"
charge = 1
mass_amu = 1
density = 1e14
temperature = 0.1
ppc = {ppc}

[output]
dir = "output"
history_every = 0
fields_every = 0
status_every = 0
"""

FIELDS = ["kind", "mode", "cores", "ranks", "threads", "decomposition", "cells", "particles", "steps",
          "mg_levels", "poisson_iters", "ms_per_step", "particle_steps_per_s", "ns_per_particle_step_core",
          "t_push", "t_deposit", "t_field_solve", "t_inject_migrate", "t_sort", "t_other"]


# ---------------------------------------------------------------------------------- machine
def topology():
    """(sockets, physical cores per socket), from lscpu with a fallback to os.cpu_count()."""
    try:
        out = subprocess.run(["lscpu"], capture_output=True, text=True, check=True).stdout
        info = dict(l.split(":", 1) for l in out.splitlines() if ":" in l)
        sockets = int(info["Socket(s)"])
        per_socket = int(info["Core(s) per socket"])
        return sockets, per_socket, info.get("Model name", "").strip()
    except (OSError, KeyError, ValueError, subprocess.CalledProcessError):
        return 1, os.cpu_count() or 1, platform.processor()


def machine_slug(model, sockets):
    """Results directory name from the CPU model, e.g. xeon-platinum-8354h-4s (no hostnames)."""
    words = re.sub(r"\(R\)|\(TM\)|CPU|@.*|Intel|AMD|Processor|\d+-Core", " ", model).lower().split()
    return "-".join(words + [f"{sockets}s"]) if words else f"cpu-{sockets}s"


def mpi_flavour(mpirun):
    try:
        v = subprocess.run([mpirun, "--version"], capture_output=True, text=True, timeout=20)
        text = v.stdout + v.stderr
    except OSError:
        sys.exit(f"cannot run '{mpirun}'; pass --mpirun /path/to/mpirun")
    if "Open MPI" in text or "OpenRTE" in text:
        return "openmpi", text.strip().splitlines()[0]
    return "other", text.strip().splitlines()[0] if text.strip() else "unknown"


def build_info(exe):
    """Compiler and options from the CMake cache next to the executable, if there is one."""
    info = {}
    try:
        with open(os.path.join(os.path.dirname(os.path.abspath(exe)), "CMakeCache.txt")) as f:
            for line in f:
                k, _, v = line.strip().partition("=")
                name = k.split(":")[0]
                if name in ("CMAKE_CXX_COMPILER", "CMAKE_BUILD_TYPE", "PICSP_NATIVE", "PICSP_FAST_MATH"):
                    info[name] = v
    except OSError:
        pass
    cxx = info.get("CMAKE_CXX_COMPILER")
    if cxx:
        v = subprocess.run([cxx, "--version"], capture_output=True, text=True).stdout
        info["compiler"] = v.splitlines()[0] if v else cxx
    return info


# ---------------------------------------------------------------------------------- layouts
def layout(mode, cores, per_socket):
    """(ranks, threads) for a mode and core count."""
    if mode == "mpi":
        return cores, 1
    if mode == "openmp":
        return 1, cores
    ranks = max(1, math.ceil(cores / per_socket))  # hybrid: one rank per socket in use
    return ranks, cores // ranks


def launch_cmd(mpirun, flavour, ranks, threads, mode):
    cmd = [mpirun, "-np", str(ranks)]
    if flavour != "openmpi":
        return cmd
    if mode == "openmp" and threads > 1:
        # one rank whose threads span several sockets: let OMP_PLACES do the pinning
        return cmd + ["--bind-to", "none"]
    if mode == "hybrid" and ranks > 1:
        return cmd + ["--map-by", f"ppr:1:socket:PE={threads}", "--bind-to", "core"]
    return cmd + ["--map-by", f"slot:PE={threads}", "--bind-to", "core"]


# ---------------------------------------------------------------------------------- one run
TIMER = re.compile(r"^(\S.*?)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)%\s*$")


def parse(text):
    r = {}
    m = re.search(r"finished (\d+) steps in [\d.]+ s\s+\(([\d.]+) ms/step, ([\d.eE+-]+) particle-steps/s, ([\d.]+) ns", text)
    if not m:
        raise RuntimeError("no 'finished' line in the output:\n" + text[-3000:])
    r["steps"], r["ms_per_step"] = int(m[1]), float(m[2])
    r["particle_steps_per_s"], r["ns_per_particle_step_core"] = float(m[3]), float(m[4])
    m = re.search(r"MPI ranks: (\d+) \(decomposition (\d+) x (\d+) x (\d+)\)", text)
    r["decomposition"] = f"{m[2]}x{m[3]}x{m[4]}" if m else ""
    m = re.search(r"(\d+) multigrid level", text)
    r["mg_levels"] = int(m[1]) if m else ""
    m = re.search(r"average Poisson iterations per step: ([\d.]+)", text)
    r["poisson_iters"] = float(m[1]) if m else ""
    timers = {}
    for line in text.splitlines():
        t = TIMER.match(line.strip())
        if t and t[1] != "timer":
            timers[t[1].strip()] = float(t[4])  # max over ranks: what the step waits for
    per = lambda k: 1e3 * timers.get(k, 0.0) / r["steps"]
    r["t_push"], r["t_deposit"], r["t_field_solve"] = per("push"), per("deposit"), per("field solve")
    r["t_inject_migrate"], r["t_sort"] = per("inject+migrate"), per("sort")
    r["t_other"] = max(0.0, r["ms_per_step"] - sum(r[k] for k in ("t_push", "t_deposit", "t_field_solve",
                                                                 "t_inject_migrate", "t_sort")))
    return r


def run_case(args, flavour, workdir, kind, mode, cores, cells, ppc, steps, decomposition=None):
    ranks, threads = layout(mode, cores, args.cores_per_socket)
    os.makedirs(workdir, exist_ok=True)
    par = "\n[parallel]\ndecomposition = [{}, {}, {}]\n".format(*decomposition) if decomposition else ""
    with open(os.path.join(workdir, "input.toml"), "w") as f:
        f.write(INPUT.format(nx=cells[0], ny=cells[1], nz=cells[2], steps=steps, ppc=ppc, parallel=par))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_PLACES="cores",
               OMP_PROC_BIND="spread" if mode == "openmp" else "close")
    env.pop("DISPLAY", None)  # hwloc otherwise probes the X server
    cmd = launch_cmd(args.mpirun, flavour, ranks, threads, mode) + [args.exe, "input.toml"]
    best = None
    for _ in range(args.repeats):
        p = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(f"{' '.join(cmd)} failed in {workdir}:\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
        r = parse(p.stdout)
        if best is None or r["ms_per_step"] < best["ms_per_step"]:
            best = r
    shutil.rmtree(os.path.join(workdir, "output"), ignore_errors=True)
    best.update(kind=kind, mode=mode, cores=cores, ranks=ranks, threads=threads,
                cells=cells[0] * cells[1] * cells[2], particles=2 * ppc * cells[0] * cells[1] * cells[2])
    print(f"  {kind:6s} {mode:7s} {cores:4d} cores = {ranks:3d} x {threads:3d}   {best['ms_per_step']:9.2f} ms/step"
          f"   solve {best['t_field_solve']:8.2f} ms   MG {best['mg_levels']}   dec {best['decomposition']}", flush=True)
    return best


# ---------------------------------------------------------------------------------- sweeps
def core_counts(sockets, per_socket, max_cores):
    """1, 2, 4, ... up to one socket, then whole sockets."""
    counts, c = [], 1
    while c < per_socket:
        counts.append(c)
        c *= 2
    counts += [per_socket * s for s in range(1, sockets + 1) if s in (1, 2) or s == sockets or s % 2 == 0]
    return sorted({c for c in counts if c <= max_cores})


def weak_decomposition(cores):
    """Near-cubic factorisation px >= py >= pz of the core count."""
    best = None
    for px in range(1, cores + 1):
        for py in range(1, cores // px + 1):
            if cores % (px * py):
                continue
            pz = cores // (px * py)
            d = sorted((px, py, pz), reverse=True)
            score = d[0] / d[2]
            if best is None or score < best[0]:
                best = (score, tuple(d))
    return best[1]


def main():
    sockets, per_socket, model = topology()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", default=os.environ.get("PICSP3D", os.path.join(ROOT, "build", "picsp3d")))
    ap.add_argument("--mpirun", default=os.environ.get("MPIEXEC", "mpirun"))
    ap.add_argument("--kind", nargs="+", default=["strong", "weak"], choices=["strong", "weak"])
    ap.add_argument("--modes", nargs="+", default=["mpi", "hybrid", "openmp"], choices=["mpi", "hybrid", "openmp"])
    ap.add_argument("--max-cores", type=int, default=sockets * per_socket)
    ap.add_argument("--cores-per-socket", type=int, default=per_socket)
    ap.add_argument("--strong-cells", type=int, default=144, help="cells per axis of the strong-scaling grid")
    ap.add_argument("--weak-cells", type=int, default=48, help="cells per axis per core in the weak-scaling runs")
    ap.add_argument("--ppc", type=int, default=16, help="particles per cell and species")
    ap.add_argument("--steps", type=int, default=40)
    ap.add_argument("--repeats", type=int, default=2, help="runs per point; the fastest is kept")
    ap.add_argument("--quick", action="store_true", help="tiny problem, to check the setup")
    ap.add_argument("--out", default=os.path.join(HERE, "results", machine_slug(model, sockets)))
    args = ap.parse_args()
    if args.quick:
        args.strong_cells, args.weak_cells, args.ppc, args.steps, args.repeats = 48, 24, 4, 10, 1
    if not os.path.isfile(args.exe):
        sys.exit(f"executable not found: {args.exe} (build first, or pass --exe)")

    flavour, mpi_version = mpi_flavour(args.mpirun)
    if flavour != "openmpi":
        print(f"note: {args.mpirun} is not Open MPI; ranks and threads will not be pinned", file=sys.stderr)
    nsock = math.ceil(args.max_cores / args.cores_per_socket)
    counts = core_counts(nsock, args.cores_per_socket, args.max_cores)
    os.makedirs(args.out, exist_ok=True)
    work = os.path.join(args.out, "work")

    commit = subprocess.run(["git", "-C", ROOT, "describe", "--always", "--dirty"], capture_output=True, text=True).stdout.strip()
    system = dict(cpu=model, sockets=sockets, cores_per_socket=per_socket,
                  kernel=platform.release(), mpi=mpi_version, commit=commit, date=dt.date.today().isoformat(),
                  build=build_info(args.exe),
                  settings=dict(strong_cells=args.strong_cells, weak_cells=args.weak_cells, ppc=args.ppc,
                                steps=args.steps, repeats=args.repeats, core_counts=counts))
    with open(os.path.join(args.out, "system.json"), "w") as f:
        json.dump(system, f, indent=2)
    print(f"{model}: {sockets} sockets x {per_socket} cores; core counts {counts}")

    for kind in args.kind:
        rows, done = [], {}
        for mode in args.modes:
            for c in counts:
                ranks, threads = layout(mode, c, args.cores_per_socket)
                if (ranks, threads) in done:  # e.g. 1 core, or hybrid within one socket = openmp
                    r = dict(done[ranks, threads], mode=mode)
                    print(f"  {kind:6s} {mode:7s} {c:4d} cores = {ranks:3d} x {threads:3d}   same layout as {done[ranks, threads]['mode']}")
                elif kind == "strong":
                    n = args.strong_cells
                    r = run_case(args, flavour, os.path.join(work, f"{kind}_{mode}_{c}"), kind, mode, c,
                                 (n, n, n), args.ppc, args.steps)
                else:
                    d = weak_decomposition(c)
                    cells = tuple(args.weak_cells * x for x in d)
                    r = run_case(args, flavour, os.path.join(work, f"{kind}_{mode}_{c}"), kind, mode, c, cells,
                                 args.ppc, args.steps, decomposition=d if ranks == c else None)
                done.setdefault((ranks, threads), r)
                rows.append(r)
                with open(os.path.join(args.out, f"{kind}.csv"), "w", newline="") as f:  # rewritten as we go
                    w = csv.DictWriter(f, fieldnames=FIELDS)
                    w.writeheader()
                    w.writerows({k: x.get(k, "") for k in FIELDS} for x in rows)
    shutil.rmtree(work, ignore_errors=True)
    print(f"results in {args.out}")


if __name__ == "__main__":
    main()
