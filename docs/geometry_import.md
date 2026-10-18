# Importing CAD geometry (STL)

PICSP3D builds objects from triangulated surfaces in STL format, the export every CAD package
offers (SolidWorks, Fusion 360, FreeCAD, Onshape, CATIA, Blender, OpenSCAD and many more). This
guide covers getting a CAD part into a simulation: checking the file, choosing units and
placement, generating the input blocks and verifying the result on the grid before a run.

A complete worked example is [`examples/features/16_mesh_import`](../examples/features/16_mesh_import).

## The short version

```bash
# 1. look at the file (size, units, watertightness)
python -m picsp_post mesh info part.stl

# 2. place it and generate [grid] + [[objects]] blocks and a preview
python -m picsp_post mesh import part.stl --units mm --anchor center --name probe \
       --cells-across 24 --margin 1.0 -o probe.toml --preview probe.png

# 3. paste probe.toml into your input, add species etc., and check the geometry
picsp3d input.toml --check --output geometry_check
```

`--check` sets up the geometry and the field solver, prints a mesh report and the extent of every
object in cells, warns about objects leaving the domain or crowding a boundary, writes
`geometry_check/fields.pvd` (`object_id`, signed distance `sdf`, `phi`) and
`geometry_check/surface.pvd` (the object surface as the solver resolves it), and exits. Open both
in ParaView to see exactly what the plasma will see.

## STL keys in `[[objects]]`

An STL file is used like any other shape, either directly on the object or as one of its
`[[objects.parts]]` (so it can be combined with primitives by CSG).

```toml
[[objects]]
name = "probe"
shape = "stl"                 # alias: "mesh"
file = "cad/probe.stl"        # relative to the input file
units = "mm"                  # m (default), cm, mm, um, in, ft
anchor = "center"             # which point of the mesh is placed at `center`
center = [0.0, 0.0, 0.0]      # [m]
rotate = [0.0, 0.0, 90.0]     # degrees about x, then y, then z, applied about the anchor
bias = "floating"
```

| key | default | meaning |
|---|---|---|
| `file` | | ASCII or binary STL, path relative to the input file |
| `units` | `m` | unit of the file coordinates: `m`, `cm`, `mm`, `um`, `in` (`inch`), `ft` |
| `scale` | 1 | additional scale factor (multiplies `units`) |
| `fit` | off | rescale the mesh so that its largest extent is `fit` metres (overrides `units`/`scale` in effect) |
| `anchor` | `origin` | the mesh point placed at `center`, see below |
| `center` | `[0, 0, 0]` | position of the anchor point [m] |
| `rotate` | `[0, 0, 0]` | Euler angles [deg], x then y then z, about the anchor |

Anchors:

| `anchor` | point of the mesh placed at `center` |
|---|---|
| `origin` | the CAD origin (file coordinates are kept, only translated by `center`) |
| `center` | centre of the bounding box |
| `centroid` | volume centroid (area centroid for open meshes) |
| `min` | minimum corner of the bounding box |
| `bottom` | centre of the bounding-box face at minimum z (stand an object on a plane) |

The order of operations is: units and scale, `fit`, move the anchor to the origin, rotate,
translate to `center`. Object-level `translate` and `rotate` (for multi-part objects) are applied
after that, to the whole object.

## What the solver prints

At start-up each mesh is reported once:

```
 mesh cad/cubesat_mm.stl
      200 triangles, watertight
      size 0.3708 x 0.2636 x 0.106 m, area 0.2449 m^2, volume 0.00355 m^3, edges 0.00588 .. 0.354 m
      local bbox [-0.1854 -0.1318 -0.053] .. [0.1854 0.1318 0.053] m (anchor 'center')
 object cubesat    spans [-0.1832 -0.04644 -0.053] .. [0.1988 0.05356 0.053] m = 24.0 x 6.3 x 6.7 cells
```

and warnings appear when an object extends outside the domain, sits less than three cells from a
non-periodic boundary, crosses a periodic boundary, or when a mesh is implausibly large or small
(usually a unit mistake).

## The `mesh` tool

`python -m picsp_post mesh` (also installed as `picsp-post mesh`) needs only NumPy and
Matplotlib.

### `mesh info`

```
python -m picsp_post mesh info cad/cubesat_mm.stl
cubesat_mm.stl
  triangles      200
  vertices       114 (after welding)
  size           370.822 x 263.603 x 106 m
  ...
  topology       closed and manifold
  note           370.8 m across: is the file in mm? use --units
```

Reported: triangle and vertex counts, degenerate triangles, bounding box, size, centroid, area,
volume, edge length range, open and non-manifold edges, inverted orientation, and a guess of the
CAD units when the size in metres is implausible. `--units mm` interprets the coordinates in
that unit and adds the size in metres.

### `mesh import`

| option | meaning |
|---|---|
| `--units`, `--scale`, `--fit`, `--anchor`, `--rotate AX AY AZ` | as the input keys above (`--anchor` defaults to `center` here) |
| `--at X Y Z` | position of the anchor in the domain [m] (default origin) |
| `--name`, `--bias`, `--potential` | object settings |
| `--dx` or `--cells-across N` | cell size, directly or as cells across the largest object extent (default 24) |
| `--debye L` | cap the cell size at half a Debye length |
| `--margin F` | free space on every side, in units of the largest object extent (default 1.0) |
| `-o FILE` | write the TOML blocks to a file (default: print) |
| `--preview PNG` | three orthographic views of the placed object with the domain and the cell grid |
| `--bake OUT.stl` | write the placed mesh in metres (units, fit, anchor, rotation applied) and reference that file instead |
| `--relative-to DIR` | directory of the input file, for the STL path written into the TOML |

Cell counts are rounded up to multiples of four (multigrid coarsening and MPI splitting). The
suggested domain is a starting point: enlarge it in the flow direction for wakes, and check the
cell size against the Debye length of your plasma (`--debye`).

### Python

```python
from picsp_post.mesh import Mesh
m = Mesh.load("part.stl")
print(m.report().text("part", units="mm"))
world = m.placed(units="mm", anchor="bottom", rotate=(90, 0, 0), center=(0, 0, -0.2))
world.save("part_in_metres.stl")
```

## Preparing CAD files

* **Units.** Most mechanical CAD exports use millimetres; set `units = "mm"`. A mesh hundreds
  of metres across is the classic symptom.
* **Watertightness is not required.** Inside tests use generalised winding numbers, so
  assemblies exported as overlapping closed bodies, small gaps and inconsistent triangle
  orientation work without boolean clean-up. Large holes (a missing face) do make the inside
  ambiguous near the hole; close those in CAD.
* **Resolution.** Detail below one cell is not resolved by the field solve, except that thin
  sheets (panels, booms, antennas) are kept as one-node-thick conductors. Particle collisions
  with the surface use the exact triangles. Simplify very fine meshes (fillets, threads, text)
  in CAD or a mesh tool: setup time grows with the triangle count, run time does not.
* **Separate conductors.** Parts that are electrically isolated (a probe on a boom, a biased
  electrode) must be separate `[[objects]]`, each from its own STL file. Parts that belong to
  the same conductor can share one file, be listed as `[[objects.parts]]`, or be connected with
  `group`.
* **Mixed geometry.** STL parts combine with primitives by CSG, e.g. an STL bus with box solar
  panels (see example 12), or a primitive with a CAD feature subtracted.
* **Symmetry and periodic axes.** Objects must stay inside the domain on periodic axes.
