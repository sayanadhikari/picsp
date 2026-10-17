"""Shared helpers for the PICSP3D feature examples.

Every example directory contains an `input.toml`, an `analyze.py` and a README. Run one with

    python3 run_example.py 04_magnetic_mirror          # simulate + analyse
    python3 run_example.py 04_magnetic_mirror --analyze-only

Figures and animations are written to the example's `output/` directory and, for the
documentation, to docs/media/.
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


def exe():
    for c in (os.environ.get("PICSP3D"), os.path.join(ROOT, "build", "picsp3d")):
        if c and os.path.exists(c):
            return c
    raise FileNotFoundError("picsp3d executable not found; build it or set PICSP3D=/path/to/picsp3d")


def simulate(example_dir, ranks=1, threads=2, extra=()):
    """Run picsp3d on <example_dir>/input.toml (output in <example_dir>/output)."""
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    cmd = ["mpirun", "-np", str(ranks), "--bind-to", "none", exe(), "input.toml", "--output", "output", *extra]
    if os.geteuid() == 0:
        cmd.insert(1, "--allow-run-as-root")
    print("$", " ".join(cmd))
    subprocess.run(cmd, cwd=example_dir, env=env, check=True)


def media(name):
    os.makedirs(MEDIA, exist_ok=True)
    return os.path.join(MEDIA, name)


def save(fig, example_dir, name, publish=True):
    os.makedirs(os.path.join(example_dir, "output"), exist_ok=True)
    fig.savefig(os.path.join(example_dir, "output", name))
    if publish:
        fig.savefig(media(name))
    print("wrote", name)


def gif_writer(path, fps=12):
    import imageio.v2 as imageio
    return imageio.get_writer(path, mode="I", duration=1.0 / fps, loop=0)


def fig_to_rgb(fig):
    import io
    import imageio.v2 as imageio
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=fig.dpi)
    buf.seek(0)
    return imageio.imread(buf)[..., :3]
