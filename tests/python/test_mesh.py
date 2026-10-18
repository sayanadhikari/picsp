"""Tests for picsp_post.mesh (STL import tool). Run: python -m pytest tests/python"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))
from picsp_post import mesh as M  # noqa: E402


def cube(lo=(0, 0, 0), hi=(1, 1, 1)):
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    c = np.array([[x, y, z] for z in (0, 1) for y in (0, 1) for x in (0, 1)], float)
    v = lo + c * (hi - lo)
    faces = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]  # outward
    tris = []
    for a, b, cc, d in faces:
        tris += [(v[a], v[b], v[cc]), (v[a], v[cc], v[d])]
    return np.asarray(tris)


def test_binary_and_ascii_roundtrip(tmp_path):
    t = cube()
    p = tmp_path / "c.stl"
    M.write_stl(str(p), t)
    assert np.allclose(M.read_stl(str(p)), t)
    a = tmp_path / "a.stl"
    with open(a, "w") as f:
        f.write("solid c\n")
        for tri in t:
            f.write(" facet normal 0 0 0\n  outer loop\n")
            for v in tri:
                f.write(f"   vertex {v[0]} {v[1]} {v[2]}\n")
            f.write("  endloop\n endfacet\n")
        f.write("endsolid c\n")
    assert np.allclose(M.read_stl(str(a)), t)


def test_report_closed_cube():
    r = M.Mesh(cube((0, 0, 0), (2, 1, 1))).report()
    assert r.triangles == 12 and r.watertight and r.nonmanifold_edges == 0
    assert r.area == pytest.approx(10.0) and r.volume == pytest.approx(2.0)
    assert np.allclose(r.centroid, (1.0, 0.5, 0.5))


def test_open_mesh_and_unit_guess():
    r = M.Mesh(cube()[:-2]).report()
    assert not r.watertight and r.open_edges > 0
    assert M.Mesh(cube() * 300).report().unit_guess() == "mm"


@pytest.mark.parametrize("anchor,expected_lo", [("center", (-50e-3, -25e-3, -10e-3)), ("min", (0, 0, 0)),
                                                ("bottom", (-50e-3, -25e-3, 0)), ("origin", (100e-3, 0, 0))])
def test_placement_matches_solver_conventions(anchor, expected_lo):
    m = M.Mesh(cube((100, 0, 0), (200, 50, 20)))  # mm
    lo, hi = m.placed(units="mm", anchor=anchor).bounds()
    assert np.allclose(lo, expected_lo, atol=1e-12)
    assert np.allclose(hi - lo, (0.1, 0.05, 0.02))


def test_fit_rotate_center():
    m = M.Mesh(cube((0, 0, 0), (4, 1, 1)))
    w = m.placed(fit=2.0, anchor="center", rotate=(0, 0, 90), center=(1, 2, 3))
    lo, hi = w.bounds()
    assert np.allclose(hi - lo, (0.5, 2.0, 0.5))           # long axis turned onto y, extent 2
    assert np.allclose(0.5 * (lo + hi), (1, 2, 3))


def test_suggest_grid():
    origin, cells, dx = M.suggest_grid((-0.1, -0.05, -0.05), (0.1, 0.05, 0.05), cells_across=20, margin=1.0)
    assert dx == pytest.approx(0.01)
    assert all(c % 4 == 0 for c in cells) and cells[0] >= 60
    assert np.allclose(origin + 0.5 * cells * dx, 0.0)    # domain centred on the object


def test_cli_import(tmp_path):
    from picsp_post.cli import main
    p = tmp_path / "part.stl"
    M.write_stl(str(p), cube((10, 20, 30), (110, 70, 50)))
    out = tmp_path / "obj.toml"
    assert main(["mesh", "import", str(p), "--units", "mm", "--anchor", "bottom", "--at", "0", "0", "-0.1",
                 "--cells-across", "16", "-o", str(out)]) == 0
    text = out.read_text()
    for key in ('shape = "stl"', 'units = "mm"', 'anchor = "bottom"', "center = [0, 0, -0.1]", "[grid]"):
        assert key in text
    baked = tmp_path / "baked.stl"
    assert main(["mesh", "import", str(p), "--units", "mm", "--anchor", "bottom", "--bake", str(baked),
                 "-o", str(out)]) == 0
    lo, hi = M.Mesh.load(str(baked)).bounds()
    assert np.allclose(lo, (-0.05, -0.025, 0.0), atol=1e-6) and np.allclose(hi, (0.05, 0.025, 0.02), atol=1e-6)
