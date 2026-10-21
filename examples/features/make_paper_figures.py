#!/usr/bin/env python3
"""Publication versions of the feature-example figures: finer grids, more particles, 300 dpi.

    python3 make_paper_figures.py                      # all cases, figures into docs/paper/figures
    python3 make_paper_figures.py 07_floating_sphere_OML 12_spacecraft_STL
    python3 make_paper_figures.py --list               # show the scaling per case
    python3 make_paper_figures.py --analyze-only ...   # re-plot existing output_hires/
    python3 make_paper_figures.py --dry-run ...        # only write input_hires.toml

For every case an input_hires.toml is written next to input.toml with the grid refined (cells
multiplied, cell size divided, so the domain is unchanged), the time step refined with steps and
output intervals scaled to cover the same physical time, and more particles per cell. The run goes
to output_hires/ and the analysis script of the example produces the figures (PNG and PDF).

The figures are laid out in the journal style of picsp_post.style at their printed width (89 mm
for one column, 183 mm for the page). To change only the layout of existing runs, without Blender
and without the animations:

    PICSP_REUSE_RENDERS=1 PICSP_NO_ANIM=1 python3 make_paper_figures.py --analyze-only

The scaling figures of the paper come from the benchmark results:

    PICSP_STYLE=paper python3 ../../benchmarks/scaling/plot_scaling.py \
        ../../benchmarks/scaling/results/xeon-platinum-8354h-4s --out ../../docs/paper/figures --formats pdf
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

# grid: factor per axis on `cells` (cell size divided by it, unless the domain is given as `length`)
# time: time step divided, steps and *_every multiplied; ppc: particles per cell multiplied
# longer: run for this many times the physical time of the example (default 1)
# ranks, threads: parallel layout (small grids cannot be split over many ranks)
CASES = {
    "01_langmuir_oscillation":    dict(grid=(4, 4, 4), time=2, ppc=1, ranks=4, threads=18),
    "02_two_stream_instability":  dict(grid=(4, 1, 1), time=2, ppc=4, ranks=1, threads=18),
    "03_gyration_ExB_drift":      dict(grid=(1, 1, 1), time=1, ppc=1, ranks=1, threads=4),
    "04_magnetic_mirror":         dict(grid=(1, 1, 1), time=1, ppc=2, ranks=1, threads=18),
    "05_dipole_trapping":         dict(grid=(1, 1, 1), time=1, ppc=1, ranks=1, threads=4),
    "06_velocity_distributions":  dict(grid=(1, 1, 1), time=1, ppc=4, ranks=1, threads=18),
    "07_floating_sphere_OML":     dict(grid=(2, 2, 2), time=2, ppc=4, ranks=72, threads=1),
    "08_sphere_wake":             dict(grid=(2, 2, 2), time=2, ppc=4, ranks=72, threads=1),
    "09_langmuir_probe_IV":       dict(grid=(2, 2, 2), time=2, ppc=2, ranks=72, threads=1),
    "10_secondary_emission":      dict(grid=(2, 2, 2), time=2, ppc=4, ranks=72, threads=1),
    "11_photoemission_shadowing": dict(grid=(2, 2, 2), time=2, ppc=2, ranks=72, threads=1, longer=3),
    "12_spacecraft_STL":          dict(grid=(2, 2, 2), time=2, ppc=2, ranks=72, threads=1),
    "13_MCC_ion_mobility":        dict(grid=(1, 1, 1), time=1, ppc=16, ranks=1, threads=18),
    "14_ionization_avalanche":    dict(grid=(1, 1, 1), time=1, ppc=16, ranks=1, threads=18),
    "15_coulomb_relaxation":      dict(grid=(1, 1, 1), time=1, ppc=4, ranks=1, threads=18),
    "16_mesh_import":             dict(grid=(2, 2, 2), time=2, ppc=4, ranks=72, threads=1),
}

NUM = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"


def _fmt(x):
    return str(int(x)) if float(x).is_integer() and abs(x) < 1e15 else repr(float(x))


def scale_input(text: str, grid, tf: int, pf: int, longer: int = 1) -> str:
    """Refine an example input. Only top-level keys of [grid], [time], [output] and per-species
    ppc are touched; physical quantities (domain, rates, objects) are unchanged."""
    out, section = [], ""
    # `length` given in [grid] fixes the domain; `length` keys of objects (cylinders) do not count
    grid_block = re.search(r"^\s*\[grid\][^\n]*\n(.*?)(?=^\s*\[|\Z)", text, re.M | re.S)
    has_length = bool(grid_block) and re.search(r"^\s*length\s*=", grid_block.group(1), re.M) is not None
    for line in text.splitlines():
        m = re.match(r"^\s*\[\[?([\w.]+)\]?\]", line)
        if m:
            section = m.group(1)
        key = re.match(r"^\s*(\w+)\s*=", line)
        k = key.group(1) if key else None

        def sub_num(f):
            return re.sub(r"=\s*" + NUM, lambda mm: "= " + _fmt(f(float(mm.group(1)))), line, count=1)

        if section == "grid" and k == "cells":
            vals = [int(v) for v in re.findall(r"\d+", line.split("=", 1)[1].split("#")[0])]
            line = re.sub(r"\[[^\]]*\]", "[" + ", ".join(str(v * g) for v, g in zip(vals, grid)) + "]", line, count=1)
        elif section == "grid" and k in ("dx", "dx_debye") and not has_length:
            if len(set(grid)) != 1:
                raise ValueError("a scalar cell size cannot be refined along one axis only; give the domain as length")
            line = sub_num(lambda v: v / grid[0])
        elif section == "time" and k in ("dt", "dt_wpe"):
            line = sub_num(lambda v: v / tf)
        elif section == "time" and k == "steps":
            line = sub_num(lambda v: v * tf * longer)
        elif section == "output" and k and k.endswith("_every"):
            line = sub_num(lambda v: v * tf)
        elif section == "output" and k == "particle_stride":
            line = sub_num(lambda v: v * pf)  # keep particle dumps at the same size
        elif section == "species" and k == "ppc":
            line = sub_num(lambda v: v * pf)
        out.append(line)
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases", nargs="*", help="example directories (default: all)")
    ap.add_argument("--media", default=os.path.join(ROOT, "docs", "paper", "figures"))
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--formats", default="png,pdf")
    ap.add_argument("--render-scale", type=int, default=2)
    ap.add_argument("--mpirun", default=os.environ.get("MPIEXEC", "mpirun"))
    ap.add_argument("--analyze-only", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    names = [c.rstrip("/") for c in args.cases] or list(CASES)
    if args.list:
        for n in names:
            c = CASES[n]
            print(f"{n:30s} grid x{c['grid']}  dt /{c['time']}  ppc x{c['ppc']}  {c['ranks']} ranks x {c['threads']} threads")
        return
    os.environ.update(PICSP_INPUT="input_hires.toml", PICSP_OUTPUT_DIR="output_hires", PICSP_MEDIA_DIR=args.media,
                      PICSP_FIG_DPI=str(args.dpi), PICSP_FIG_FORMATS=args.formats,
                      PICSP_RENDER_SCALE=str(args.render_scale), PICSP_STYLE="paper", MPIEXEC=args.mpirun)
    sys.path.insert(0, HERE)
    import run_example  # after the environment is set: _common reads it at import
    for n in names:
        c = CASES[n]
        d = os.path.join(HERE, n)
        with open(os.path.join(d, "input.toml")) as f:
            text = scale_input(f.read(), c["grid"], c["time"], c["ppc"], c.get("longer", 1))
        with open(os.path.join(d, "input_hires.toml"), "w") as f:
            f.write(text)
        if args.dry_run:
            print(f"wrote {n}/input_hires.toml")
            continue
        os.environ.update(PICSP_RANKS=str(c["ranks"]), PICSP_THREADS=str(c["threads"]))
        t0 = time.time()
        print(f"===== {n}", flush=True)
        run_example.run(n, args.analyze_only)
        print(f"===== {n} done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
