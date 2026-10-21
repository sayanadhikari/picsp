"""Shared helpers for the PICSP3D feature examples.

Every example directory contains an `input.toml`, an `analyze.py` and a README. Run one with

    python3 run_example.py 04_magnetic_mirror          # simulate + analyse
    python3 run_example.py 04_magnetic_mirror --analyze-only

Figures and animations are written to the example's `output/` directory and, for the
documentation, to docs/media/.

Environment overrides (used by make_paper_figures.py; the defaults reproduce the examples):
  PICSP_INPUT, PICSP_OUTPUT_DIR   input file and output directory (input.toml, output)
  PICSP_RANKS, PICSP_THREADS      MPI ranks and OpenMP threads per run (1, 2)
  PICSP_MEDIA_DIR                 where published figures go (docs/media)
  PICSP_FIG_DPI, PICSP_FIG_FORMATS  figure resolution and formats, e.g. 300 and png,pdf
  PICSP_RENDER_SCALE              magnification of 3D renders (1)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MEDIA = os.path.join(ROOT, "docs", "media")
sys.path.insert(0, os.path.join(ROOT, "python"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402,F401
from picsp_post import Run, style  # noqa: E402,F401

E = 1.602176634e-19
ME = 9.1093837015e-31
AMU = 1.66053906660e-27
EPS0 = 8.8541878128e-12
MU0 = 1.25663706212e-6

INPUT = os.environ.get("PICSP_INPUT", "input.toml")
OUTPUT = os.environ.get("PICSP_OUTPUT_DIR", "output")


def out(example_dir):
    """The example's simulation output directory."""
    return os.path.join(example_dir, OUTPUT)


def exe():
    for c in (os.environ.get("PICSP3D"), os.path.join(ROOT, "build", "picsp3d")):
        if c and os.path.exists(c):
            return c
    raise FileNotFoundError("picsp3d executable not found; build it or set PICSP3D=/path/to/picsp3d")


def simulate(example_dir, ranks=1, threads=2, extra=()):
    """Run picsp3d on <example_dir>/<INPUT> (output in <example_dir>/<OUTPUT>)."""
    ranks = int(os.environ.get("PICSP_RANKS", ranks))
    threads = int(os.environ.get("PICSP_THREADS", threads))
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    mpirun = os.environ.get("MPIEXEC", "mpirun")
    cmd = [mpirun, "-np", str(ranks), "--bind-to", "none", exe(), INPUT, "--output", OUTPUT, *extra]
    if os.geteuid() == 0:
        cmd.insert(1, "--allow-run-as-root")
    print("$", " ".join(cmd))
    subprocess.run(cmd, cwd=example_dir, env=env, check=True)


def media(name):
    """Path for a published figure: docs/media, or $PICSP_MEDIA_DIR when set (e.g. to render
    high-resolution versions without touching the documentation media)."""
    d = os.environ.get("PICSP_MEDIA_DIR", MEDIA)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


def save(fig, example_dir, name, publish=True):
    os.makedirs(out(example_dir), exist_ok=True)
    fig.savefig(os.path.join(out(example_dir), name))
    if publish:
        stem, ext = os.path.splitext(name)
        for fmt in os.environ.get("PICSP_FIG_FORMATS", ext.lstrip(".")).split(","):
            fig.savefig(media(f"{stem}.{fmt.strip()}"))
    print("wrote", name)


def screenshot(plotter, path):
    """Save a pyvista render, magnified by PICSP_RENDER_SCALE."""
    plotter.screenshot(path, scale=int(os.environ.get("PICSP_RENDER_SCALE", 1)))


class _FrameWriter:
    """GIF writer that also keeps every frame as a lossless PNG (for video editing)."""

    def __init__(self, path, fps):
        import imageio.v2 as imageio
        self.gif = imageio.get_writer(path, mode="I", duration=1.0 / fps, loop=0)
        self.dir = os.path.splitext(path)[0] + "_frames"
        os.makedirs(self.dir, exist_ok=True)
        self.n = 0

    def append_data(self, rgb):
        import imageio.v2 as imageio
        self.gif.append_data(rgb)
        imageio.imwrite(os.path.join(self.dir, f"{self.n:04d}.png"), rgb)
        self.n += 1

    def close(self):
        self.gif.close()


def gif_writer(path, fps=12):
    """Animated GIF writer. With PICSP_FRAMES=1 every frame is also saved as a PNG in
    <name>_frames/ next to the GIF."""
    if os.environ.get("PICSP_FRAMES"):
        return _FrameWriter(path, fps)
    import imageio.v2 as imageio
    return imageio.get_writer(path, mode="I", duration=1.0 / fps, loop=0)


def fig_to_rgb(fig):
    """Rasterise a figure (animation frame). PICSP_ANIM_DPI overrides the figure resolution."""
    import io
    import imageio.v2 as imageio
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=float(os.environ.get("PICSP_ANIM_DPI", fig.dpi)))
    buf.seek(0)
    return imageio.imread(buf)[..., :3]
