"""Command line interface:  python -m picsp_post <command> <output_dir> [options]"""
from __future__ import annotations

import argparse
import os
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog="picsp-post", description="PICSP3D post-processing")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("run", help="output directory")
        p.add_argument("-o", "--out", help="output file (png/pdf/svg/mp4/gif)")
        p.add_argument("-i", "--index", type=int, default=-1, help="dump index (default: last)")

    p = sub.add_parser("info", help="summarise an output directory")
    p.add_argument("run")
    p = sub.add_parser("history", help="energies, potentials, currents vs time")
    common(p)
    p = sub.add_parser("slice", help="2D slice of a field")
    common(p)
    p.add_argument("-f", "--field", default="phi")
    p.add_argument("-a", "--axis", default="z", choices="xyz")
    p.add_argument("-p", "--pos", type=float, default=None)
    p.add_argument("-c", "--comp", type=int, default=None, help="vector component 0/1/2 (default magnitude)")
    p.add_argument("--log", action="store_true")
    p.add_argument("--streamlines", action="store_true", help="overlay E-field lines")
    p = sub.add_parser("profile", help="1D profile through the domain centre")
    common(p)
    p.add_argument("-f", "--field", default="phi")
    p.add_argument("-a", "--axis", default="x", choices="xyz")
    p = sub.add_parser("render", help="3D scene (surfaces, slice, particles)")
    common(p)
    p.add_argument("--surface", default="potential", help="surface cell array (potential, J_net, J_<species>, sun_cos)")
    p.add_argument("--slice", default="phi", help="field for the slice plane ('' for none)")
    p.add_argument("--normal", default="y", choices="xyz")
    p.add_argument("--iso", default=None, help="field for iso-surfaces (e.g. n_electron)")
    p.add_argument("--particles", action="store_true")
    p.add_argument("--light", action="store_true", help="light background")
    p.add_argument("--zoom", type=float, default=1.0)
    p = sub.add_parser("animate", help="movie of a slice (default) or of the 3D scene")
    common(p)
    p.add_argument("--scene", action="store_true", help="animate the 3D scene instead of a slice")
    p.add_argument("--orbit", action="store_true", help="rotate the camera during a 3D animation")
    p.add_argument("-f", "--field", default="phi")
    p.add_argument("-a", "--axis", default="z", choices="xyz")
    p.add_argument("-p", "--pos", type=float, default=None)
    p.add_argument("--fps", type=int, default=8)
    p.add_argument("--streamlines", action="store_true")
    p.add_argument("--particles", action="store_true")
    p.add_argument("--surface", default="potential")
    p = sub.add_parser("phase", help="phase space from sampled particles")
    common(p)
    p.add_argument("-s", "--species", default="0")
    p.add_argument("-x", default="x")
    p.add_argument("-v", default="vx")
    p = sub.add_parser("vdf", help="velocity distribution functions")
    common(p)
    p.add_argument("-s", "--species", default="0")

    a = ap.parse_args(argv)
    from .io import Run
    run = Run(a.run)
    if a.cmd == "info":
        print(run.info())
        return 0
    import matplotlib
    matplotlib.use("Agg")
    from . import plots

    def spec(s):
        return int(s) if s.isdigit() else s

    if a.cmd == "history":
        out = a.out or os.path.join(a.run, "history.png")
        plots.plot_history(run, out)
    elif a.cmd == "slice":
        out = a.out or os.path.join(a.run, f"slice_{a.field}_{a.axis}.png")
        fig, _, _ = plots.plot_slice(run, a.field, a.axis, a.pos, a.index, a.comp, log=a.log, streamlines=a.streamlines)
        fig.savefig(out)
    elif a.cmd == "profile":
        import matplotlib.pyplot as plt
        out = a.out or os.path.join(a.run, f"profile_{a.field}_{a.axis}.png")
        ax = plots.plot_profile(run, a.field, a.axis, index=a.index)
        ax.figure.tight_layout()
        ax.figure.savefig(out)
    elif a.cmd == "render":
        from .render3d import render_scene
        out = a.out or os.path.join(a.run, "scene.png")
        render_scene(run, a.index, out, surface_scalars=a.surface, slice_field=a.slice or None, slice_normal=a.normal,
                     iso_field=a.iso, particles=a.particles, dark=not a.light, zoom=a.zoom)
    elif a.cmd == "animate":
        from . import animate
        if a.scene:
            out = a.out or os.path.join(a.run, "scene.mp4")
            animate.animate_scene(run, out, fps=a.fps, orbit=a.orbit, particles=a.particles, surface_scalars=a.surface,
                                  slice_field=a.field)
        else:
            out = a.out or os.path.join(a.run, f"slice_{a.field}_{a.axis}.mp4")
            animate.animate_slice(run, a.field, a.axis, a.pos, out, fps=a.fps, streamlines=a.streamlines)
    elif a.cmd == "phase":
        out = a.out or os.path.join(a.run, "phase.png")
        ax = plots.plot_phase(run, spec(a.species), a.x, a.v, a.index)
        ax.figure.tight_layout()
        ax.figure.savefig(out)
    elif a.cmd == "vdf":
        out = a.out or os.path.join(a.run, "vdf.png")
        ax = plots.plot_vdf(run, spec(a.species), a.index)
        ax.figure.tight_layout()
        ax.figure.savefig(out)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
