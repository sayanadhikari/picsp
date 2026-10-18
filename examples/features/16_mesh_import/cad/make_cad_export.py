"""Creates cad/cubesat_mm.stl: the 3U CubeSat of examples/geometry as a typical CAD export
would deliver it (millimetres, an arbitrary origin, rotated 30 degrees about z)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "python"))
from picsp_post.mesh import Mesh  # noqa: E402

src = os.path.join(HERE, "..", "..", "..", "geometry", "cubesat_3u.stl")
Mesh.load(src).scaled(1000.0).rotated((0, 0, 30)).translated((250.0, -120.0, 40.0)).save(
    os.path.join(HERE, "cubesat_mm.stl"))
print("wrote cad/cubesat_mm.stl")
