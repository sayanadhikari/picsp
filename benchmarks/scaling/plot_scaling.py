#!/usr/bin/env python3
"""Figures from run_scaling.py results.

    python3 plot_scaling.py results/xeon-platinum-8354h-4s            # -> docs/media/scaling_*.png
    python3 plot_scaling.py results/<machine> --out figs --dpi 200

Writes scaling_strong.png (speed-up and parallel efficiency), scaling_breakdown.png (cost per
particle-step and core by stage: flat bars mean perfect scaling), scaling_weak.png and, when a
before_<commit>/ subdirectory with an earlier strong-scaling run exists, scaling_before_after.png
(whole step, and the stage that changed most).
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullLocator, ScalarFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "python"))
from picsp_post import style  # noqa: E402  (also applies the house style)

MODES = {
    "mpi": ("MPI only (1 rank per core)", style.SERIES[0], "o"),
    "hybrid": ("hybrid (1 rank per socket x OpenMP)", style.SERIES[1], "s"),
    "openmp": ("OpenMP only (1 rank)", style.SERIES[2], "^"),
}
if style.PAPER:   # journal figures: short labels, the layouts are explained in the caption
    MODES = {"mpi": ("MPI only", style.SERIES[0], "o"), "hybrid": ("hybrid", style.SERIES[1], "s"),
             "openmp": ("OpenMP only", style.SERIES[2], "^")}
FORMATS = ["png"]


def finish(fig, out, name, dpi):
    """Save in every requested format; in the journal style (PICSP_STYLE=paper) laid out at the
    printed width with panel letters, without the machine line (it goes into the caption)."""
    if style.PAPER:
        if fig._suptitle is not None:
            fig._suptitle.set_visible(False)
        for ax in fig.axes:
            for line in ax.get_lines():
                line.set_markersize(3)
        style.journal_layout(fig, width=style.DOUBLE)
    else:
        fig.tight_layout()
    for fmt in FORMATS:
        fig.savefig(os.path.join(out, f"{name}.{fmt}"), dpi=dpi)
    plt.close(fig)
STAGES = [("t_push", "push"), ("t_deposit", "deposit"), ("t_field_solve", "field solve"),
          ("t_inject_migrate", "migrate"), ("t_sort", "sort"), ("t_other", "other")]


def load(path):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v) if "." in v or "e" in v else int(v)
            except (TypeError, ValueError):
                pass
    return rows


def by_mode(rows):
    out = {}
    for r in rows:
        out.setdefault(r["mode"], []).append(r)
    return {m: sorted(v, key=lambda r: r["cores"]) for m, v in out.items()}


def core_axis(ax, cores, per_socket):
    ax.set_xscale("log", base=2)
    # drop a tick that would crowd the next one (e.g. 16 next to 18)
    ticks = [c for i, c in enumerate(cores) if i + 1 == len(cores) or cores[i + 1] / c >= 1.3]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xlabel("cores")
    for s in range(1, 1 + max(cores) // per_socket):
        if s * per_socket < max(cores):
            ax.axvline(s * per_socket, color=style.GRID if not style.PAPER else "0.82",
                       lw=1.2 if not style.PAPER else 0.5, ls="-" if not style.PAPER else (0, (3, 2)), zorder=0)


def subtitle(system):
    s = system.get("settings", {})
    return (f"{system.get('cpu', '')}, {system.get('sockets', '?')} x {system.get('cores_per_socket', '?')} cores  |  "
            f"{system.get('mpi', '')}  |  commit {system.get('commit', '')}")


def strong(rows, system, out, dpi):
    modes = by_mode(rows)
    per_socket = system.get("cores_per_socket", 18)
    t1 = min(r["ms_per_step"] for r in rows if r["cores"] == 1)
    cores = sorted({r["cores"] for r in rows})
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 3.9))
    a1.plot(cores, cores, color=style.INK2, ls="--", lw=1, label="ideal")
    for m, (label, col, mk) in MODES.items():
        if m not in modes:
            continue
        c = np.array([r["cores"] for r in modes[m]])
        sp = t1 / np.array([r["ms_per_step"] for r in modes[m]])
        a1.plot(c, sp, marker=mk, color=col, label=label)
        a2.plot(c, 100 * sp / c, marker=mk, color=col, label=label)
    for ax in (a1, a2):
        core_axis(ax, cores, per_socket)
    a1.set_yscale("log", base=2)
    a1.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    a1.set_ylabel("speed-up over 1 core")
    a1.set_title("Strong scaling: speed-up")
    a1.legend(loc="upper left")
    a2.set_ylim(0, 110)
    a2.set_ylabel("parallel efficiency [%]")
    a2.set_title("Strong scaling: efficiency")
    n = rows[0]
    a2.text(0.98, 0.97, f"{round(n['cells'] ** (1 / 3))}³ cells, {n['particles'] / 1e6:.0f} M particles\n"
            f"1 core: {t1 / 1e3:.2f} s/step", transform=a2.transAxes, ha="right", va="top", fontsize=9,
            color=style.INK2)
    fig.suptitle(subtitle(system), fontsize=9, color=style.INK2, y=0.995)
    finish(fig, out, "scaling_strong", dpi)


def breakdown(rows, system, out, dpi):
    modes = by_mode(rows)
    show = [m for m in ("mpi", "hybrid") if m in modes]
    fig, axes = plt.subplots(1, len(show), figsize=(6 * len(show), 3.9), sharey=True, squeeze=False)
    for ax, m in zip(axes[0], show):
        rs = modes[m]
        x = np.arange(len(rs))
        bottom = np.zeros(len(rs))
        for i, (key, label) in enumerate(STAGES):
            if not any(r[key] > 0.005 * r["ms_per_step"] for r in rs):
                continue
            # ns per particle-step and core: wall time x cores / particles; constant = ideal scaling
            v = np.array([r[key] * 1e6 * r["cores"] / r["particles"] for r in rs])
            ax.bar(x, v, bottom=bottom, color=style.series_color(i), width=0.7, label=label)
            bottom += v
        ax.set_xticks(x, [str(r["cores"]) for r in rs])
        ax.set_xlabel("cores")
        ax.set_title(MODES[m][0])
        ax.grid(axis="x", visible=False)
        if style.PAPER:
            ax.text(0.03, 0.97, MODES[m][0][0].upper() + MODES[m][0][1:], transform=ax.transAxes, ha="left", va="top",
                    fontweight="bold")
    axes[0][0].set_ylabel("cost per particle-step and core [ns]")
    axes[0][0].set_ylim(0, 1.3 * max(ax.get_ylim()[1] for ax in axes[0]))
    axes[0][-1].legend(loc="upper left", ncol=2, bbox_to_anchor=(0.0, 0.9) if style.PAPER else None)
    fig.suptitle("Where the time goes (strong scaling; equal bar heights = perfect scaling)\n" + subtitle(system),
                 fontsize=9, color=style.INK2)
    finish(fig, out, "scaling_breakdown", dpi)


def weak(rows, system, out, dpi):
    modes = by_mode(rows)
    per_socket = system.get("cores_per_socket", 18)
    cores = sorted({r["cores"] for r in rows})
    t1 = min(r["ms_per_step"] for r in rows if r["cores"] == 1)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 3.9))
    for m, (label, col, mk) in MODES.items():
        if m not in modes:
            continue
        c = np.array([r["cores"] for r in modes[m]])
        t = np.array([r["ms_per_step"] for r in modes[m]])
        a1.plot(c, 100 * t1 / t, marker=mk, color=col, label=label)
    for ax in (a1, a2):
        core_axis(ax, cores, per_socket)
    a1.axhline(100, color=style.INK2, ls="--", lw=1)
    a1.set_ylim(0, 110)
    a1.set_ylabel("weak-scaling efficiency [%]")
    a1.set_title("Weak scaling: efficiency")
    a1.legend(loc="lower left")
    best = "mpi" if "mpi" in modes else next(iter(modes))
    rs = modes[best]
    c = np.array([r["cores"] for r in rs])
    for i, (key, label) in enumerate(STAGES):
        if not any(r[key] > 0.005 * r["ms_per_step"] for r in rs):
            continue
        a2.plot(c, [r[key] for r in rs], marker="o", ms=4, color=style.series_color(i), label=label)
    a2.set_ylabel("time per step [ms]")
    a2.set_title(f"Weak scaling: stages, {MODES[best][0]}")
    if style.PAPER:
        a2.text(0.97, 0.03, f"{MODES[best][0]}, by stage", transform=a2.transAxes, ha="right", va="bottom",
                color=style.INK2)
    a2.legend(loc="upper left", ncol=2)
    n = rows[0]
    a1.text(0.98, 0.86, f"{round(n['cells'] ** (1 / 3))}³ cells and {n['particles'] / 1e6:.1f} M particles per core",
            transform=a1.transAxes, ha="right", va="top", fontsize=9, color=style.INK2)
    fig.suptitle(subtitle(system), fontsize=9, color=style.INK2, y=0.995)
    finish(fig, out, "scaling_weak", dpi)


def before_after(after, before, tag, system, out, dpi):
    a, b = by_mode(after), by_mode(before)
    per_socket = system.get("cores_per_socket", 18)
    cores = sorted({r["cores"] for r in after})
    # second panel: the stage whose total time changed most between the two runs
    def total(rows, key):
        return sum(r[key] for r in rows if r[key] != "")
    stage = max(STAGES, key=lambda s: abs(np.log((total(after, s[0]) + 1e-9) / (total(before, s[0]) + 1e-9))))
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.9))
    for ax, key, title in ((axes[0], "ms_per_step", "whole step"), (axes[1], stage[0], stage[1])):
        for m, (label, col, mk) in MODES.items():
            if m in a and m in b:
                ax.plot([r["cores"] for r in b[m]], [r[key] for r in b[m]], color=col, ls=":", marker=mk, mfc="none")
                ax.plot([r["cores"] for r in a[m]], [r[key] for r in a[m]], color=col, marker=mk, label=label)
        core_axis(ax, cores, per_socket)
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
        ax.set_ylabel("time per step [ms]" if key == "ms_per_step" else f"{title} time per step [ms]")
        ax.set_title(f"Strong scaling, {title}")
        if style.PAPER:
            ax.text(0.03, 0.03, title[0].upper() + title[1:], transform=ax.transAxes, ha="left", va="bottom",
                    fontweight="bold")
    from matplotlib.lines import Line2D
    extra = [Line2D([], [], color=style.INK, ls=":", lw=1, label=f"before ({tag})"),
             Line2D([], [], color=style.INK, lw=1, label="after")]
    h_, l_ = axes[0].get_legend_handles_labels()
    axes[0].legend(h_ + extra, l_ + [e.get_label() for e in extra], loc="upper right")
    fig.suptitle(f"Before (dotted, {tag}) and after (solid)  |  " + subtitle(system), fontsize=9, color=style.INK2, y=0.995)
    finish(fig, out, "scaling_before_after", dpi)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="directory written by run_scaling.py")
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "media"))
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--formats", default="png", help="comma-separated, e.g. png,pdf")
    args = ap.parse_args()
    FORMATS[:] = [f.strip() for f in args.formats.split(",")]
    with open(os.path.join(args.results, "system.json")) as f:
        system = json.load(f)
    os.makedirs(args.out, exist_ok=True)
    s, w = load(os.path.join(args.results, "strong.csv")), load(os.path.join(args.results, "weak.csv"))
    if s:
        strong(s, system, args.out, args.dpi)
        breakdown(s, system, args.out, args.dpi)
    if w:
        weak(w, system, args.out, args.dpi)
    for d in sorted(glob.glob(os.path.join(args.results, "before_*"))):
        b = load(os.path.join(d, "strong.csv"))
        if s and b:
            before_after(s, b, os.path.basename(d).split("_", 1)[1], system, args.out, args.dpi)
    print(f"figures in {args.out}")


if __name__ == "__main__":
    main()
