#!/usr/bin/env python3
"""Generate a 3U CubeSat bus as a binary STL (metres): a 0.34 x 0.10 x 0.10 m body with
four corner rails, a camera barrel and a patch antenna. The parts are separate,
overlapping closed shells - PICSP3D's winding-number inside test handles that, so
there is no need for boolean unions in CAD."""
import struct
import sys
import numpy as np


def box(lo, hi):
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    p = np.array([[hi[0] if c & 1 else lo[0], hi[1] if c & 2 else lo[1], hi[2] if c & 4 else lo[2]] for c in range(8)])
    quads = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    tris = []
    for a, b, c, d in quads:
        tris += [(p[a], p[b], p[c]), (p[a], p[c], p[d])]
    return tris


def cylinder(c0, c1, r, n=32):
    c0, c1 = np.asarray(c0, float), np.asarray(c1, float)
    ax = c1 - c0
    ax /= np.linalg.norm(ax)
    e1 = np.cross(ax, [1, 0, 0] if abs(ax[0]) < 0.9 else [0, 1, 0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(ax, e1)
    ring = [r * (np.cos(t) * e1 + np.sin(t) * e2) for t in np.linspace(0, 2 * np.pi, n, endpoint=False)]
    tris = []
    for i in range(n):
        a, b = ring[i], ring[(i + 1) % n]
        tris += [(c0 + a, c1 + a, c1 + b), (c0 + a, c1 + b, c0 + b)]
        tris += [(c0, c0 + b, c0 + a), (c1, c1 + a, c1 + b)]
    return tris


def write_stl(path, tris):
    with open(path, "wb") as f:
        f.write(b"PICSP3D CubeSat example".ljust(80, b" "))
        f.write(struct.pack("<I", len(tris)))
        for t in tris:
            n = np.cross(t[1] - t[0], t[2] - t[0])
            n = n / (np.linalg.norm(n) or 1)
            f.write(struct.pack("<12fH", *n, *t[0], *t[1], *t[2], 0))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "cubesat_3u.stl"
    L, W = 0.34, 0.10
    tris = box([-L / 2, -W / 2, -W / 2], [L / 2, W / 2, W / 2])
    rail = 0.0085
    for sy in (-1, 1):
        for sz in (-1, 1):
            y0 = sy * (W / 2 - rail / 2)
            z0 = sz * (W / 2 - rail / 2)
            tris += box([-L / 2 - 0.007, y0 - rail / 2, z0 - rail / 2], [L / 2 + 0.007, y0 + rail / 2, z0 + rail / 2])
    tris += cylinder([L / 2 - 0.01, 0, 0], [L / 2 + 0.035, 0, 0], 0.03)             # camera barrel (wake side)
    tris += box([-0.08, -0.035, W / 2 - 0.002], [0.08, 0.035, W / 2 + 0.006])        # patch antenna (+z face)
    write_stl(out, tris)
    print(f"wrote {out}: {len(tris)} triangles")
