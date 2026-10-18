# PICSP3D input reference

Inputs are [TOML](https://toml.io) files. All quantities are SI unless stated otherwise;
temperatures and energies are in electronvolts, charges in units of the elementary charge.
Relative file paths (STL) are resolved relative to the input file.

## Top level

| key | default | meaning |
|---|---|---|
| `title` | `"picsp3d"` | run title |
| `seed` | `12345` | random seed (streams are split per rank and thread) |
| `reference_species` | first negative species with n, T > 0 | species defining ω_p and λ_D for `dt_wpe`, `dx_debye` and the printed diagnostics |

## `[grid]`

| key | meaning |
|---|---|
| `cells = [nx, ny, nz]` | number of cells (nodes are `n+1` per axis) |
| `dx` | cell size in metres (scalar or `[dx, dy, dz]`) |
| `dx_debye` | cell size in reference Debye lengths (alternative to `dx`) |
| `length = [Lx, Ly, Lz]` | domain size (alternative to `dx`) |
| `origin = [x0, y0, z0]` | lower corner (default `[0,0,0]`) |
| `centered = true` | place the domain centre at the origin |

For multigrid efficiency choose cell counts per rank that are divisible by a power of two.

## `[time]`

| key | meaning |
|---|---|
| `dt` / `dt_wpe` | time step in seconds or in units of 1/ω_p of the reference species |
| `steps` / `duration` / `duration_wpe` | number of steps, or total time (s or 1/ω_p) |

## `[boundaries]`

Each face takes a table `{ field = ..., value = ..., particle = ... }`:

* `field`: `periodic`, `dirichlet` (potential `value` in V) or `neumann` (zero normal field);
* `particle`: `periodic`, `absorb`, `reflect` (specular) or `open` (absorb, and inject the
  ambient distribution of every species with `inject = true` with the correct one-sided flux).

Keys: `default`, `x`, `y`, `z` (both faces of an axis), `x_lo`, `x_hi`, `y_lo`, `y_hi`, `z_lo`, `z_hi`.
Periodic must be set on both faces of an axis, for fields and particles.

## `[fields]`

| key | meaning |
|---|---|
| `B0 = [Bx, By, Bz]` | uniform background magnetic field (T) |
| `E0 = [Ex, Ey, Ez]` | uniform external electric field (V/m) |
| `frame_velocity = [vx, vy, vz]` | velocity of the simulation frame (e.g. a spacecraft); adds the motional field E = -v x B0 |
| `self_consistent` | `true`; `false` switches off the Poisson solve (test-particle mode: only E0, motional E and B act) |

### `[[fields.magnet]]`: spatially varying magnetic field

Any number of analytic sources can be added; their fields are superposed on `B0` and evaluated
once on the grid nodes (ghost layers included). Particles gather B with the same shape function as E.

| `type` | keys | field |
|---|---|---|
| `dipole` | `center`, `moment = [mx, my, mz]` (A m²) | point dipole μ0/4π (3(m·r̂)r̂ − m)/r³ |
| `loop` (alias `coil`) | `center`, `axis`, `radius`, `current` | exact circular loop (complete elliptic integrals) |
| `solenoid` | `center`, `axis`, `radius`, `current`, `length`, `loops` | finite solenoid as `loops` coaxial loops |
| `mirror` | `center`, `axis`, `radius`, `current`, `separation` | two coaxial loops `separation` apart (magnetic bottle) |
| `wire` (alias `line`) | `center` (a point on the wire), `axis`, `current` | infinite straight wire μ0 I/(2π ρ) |

Points closer than half a cell to a singular source (dipole centre, wire, loop conductor) are
evaluated at that distance.

```toml
[[fields.magnet]]
type = "mirror"
center = [0.5, 0.5, 0.5]
axis = [0, 0, 1]
radius = 0.12
current = 2.0e4
separation = 0.8
```

## `[sun]`

`direction = [x, y, z]`: unit vector pointing from the objects towards the sun (photoemission).

## `[[species]]` (one table per species)

| key | default | meaning |
|---|---|---|
| `name` | | species name (used in output names) |
| `charge` | | charge number (e.g. `-1`, `1`, `2`) |
| `mass_amu` / `mass_me` / `mass_kg` / `mass` (amu) | | particle mass |
| `density` | 0 | ambient / initial density (m^-3); 0 for emitted-only species |
| `temperature` | 0 | temperature (eV) |
| `drift = [vx, vy, vz]` | 0 | drift velocity (m/s) |
| `distribution` | `maxwellian` | `maxwellian`, `bimaxwellian`, `kappa`, `ring`, `cold` |
| `t_par`, `t_perp` | `temperature` | bi-Maxwellian / ring temperatures (eV) along / across `aniso_axis` |
| `aniso_axis` | direction of B0 (else z) | anisotropy axis |
| `kappa` | 4 | kappa index (> 1.5); `temperature` is the kinetic temperature |
| `ring_speed` | 0 | perpendicular ring speed (m/s) |
| `ppc` | 8 | macro-particles per cell at `density` |
| `weight` | n dV / ppc | explicit macro-particle weight |
| `load` | `random` | `random` or `quiet` (Halton sequence) |
| `load_initial` | `density > 0` | fill the domain at t = 0 |
| `inject` | `density > 0` | inject through `open` boundaries |
| `subcycle` | 1 | push every N steps with N dt (heavy ions) |
| `region = { min = [...], max = [...] }` | whole domain | initial loading box |
| `perturbation = { axis, amplitude, mode, velocity_amplitude }` | | sinusoidal density (and velocity) perturbation for wave studies |

Emitted species (photo- and secondary electrons) use the weight of the first electron species unless
`weight` is given.

## `[[objects]]`

| key | default | meaning |
|---|---|---|
| `name` | | object name |
| `bias` | `floating` | `floating`, `fixed` (alias `grounded`) |
| `potential` | 0 | fixed bias (V), or initial potential of a floating object |
| `ramp` | 0 | linear bias ramp (V/s) for fixed objects, e.g. probe sweeps |
| `group` | none | objects with the same integer group are electrically connected |
| `ext_capacitance` | 0 | extra capacitance to ground (F) |
| `ext_current` | 0 | constant external current into the object (A) |
| `translate`, `rotate` | 0 | rigid transform of the whole object (rotation: Euler angles in degrees, x then y then z) |
| `emit_species` | first electron species | species receiving photo/secondary electrons |
| `see = { model, yield, delta_max, e_max, temperature, ion_yield }` | none | secondary emission: `constant` (yield) or `sternglass` (δ_max at E_max eV); emitted temperature (eV); ion-induced yield |
| `photoemission = { current_density, temperature }` | off | saturation photocurrent density at normal incidence (A/m²) and photoelectron temperature (eV) |

### Shapes

An object is either a single shape (shape keys directly in the object table) or a list of
`[[objects.parts]]` combined with `op = "union" | "subtract" | "intersect"` (the first part must be a union).

| shape | keys |
|---|---|
| `sphere` | `center`, `radius` |
| `box`, `plate` | `center`, `size = [sx, sy, sz]` (full extents), `rounding` |
| `cylinder` | `center`, `axis`, `radius`, `length` |
| `cone` | `center`, `axis`, `radius` (base), `radius2` (top), `length` |
| `torus` | `center`, `axis`, `radius` (major), `radius2` (minor) |
| `ellipsoid` | `center`, `size` (full axes) |
| `halfspace` | `center` (point on plane), `axis` (outward normal) |
| `stl` (alias `mesh`) | `file`, `units` (`m`, `cm`, `mm`, `um`, `in`, `ft`), `scale`, `fit`, `anchor` (`origin`, `center`, `centroid`, `min`, `bottom`), `center` (position of the anchor) |

Every part can be rotated with `rotate = [ax, ay, az]` (degrees, about x then y then z) about its
`center`; on a single-shape object `rotate` turns the shape once about its `center`.
STL files may be ASCII or binary; overlapping shells and inconsistent orientation are handled by the
generalised winding number. `anchor` selects the mesh point placed at `center` (`origin` keeps the
CAD coordinates), `units` converts the file units and `fit` rescales the mesh to a given largest
extent in metres. See [`geometry_import.md`](geometry_import.md) for the import workflow and the
`picsp-post mesh` tool that generates these blocks from a CAD file. Features thinner than a cell (panels, booms) are resolved as one-node-thick
conductor sheets in the field solve and exactly for particle collisions.

## `[collisions]`

| key | meaning |
|---|---|
| `neutral_density` | neutral density (m^-3) |
| `neutral_temperature` | neutral temperature (eV) |
| `neutral_mass_amu` | neutral mass |
| `neutral_drift` | neutral drift velocity (m/s) |
| `[[collisions.process]]` | one table per process, see below |

Each process acts on one species against the neutral background (null-collision Monte Carlo):

| key | meaning |
|---|---|
| `species` | colliding species |
| `type` | `elastic`, `charge_exchange` (`cx`), `excitation`, `ionization` (`ionisation`) |
| `model` | `constant` (default), `step` (σ above `threshold`), `lotz`, `table` |
| `cross_section` | σ (m²) for `constant`/`step`; the Lotz coefficient a (m² eV²) for `lotz`: σ = a ln(E/I)/(E I) |
| `threshold` | energy loss / threshold I (eV) for excitation and ionisation |
| `table = [[E, sigma], ...]` or `file = "..."` | tabulated σ(E) (eV, m²), linear interpolation; files hold two columns, `#` comments allowed |
| `product` | ionisation: ion species created (required) |
| `secondary` | ionisation: species of the secondary electron (default: the colliding species) |

Energies are centre-of-mass energies. Ionisation subtracts the threshold, shares the remaining
energy randomly between the primary and the secondary electron and creates an ion with the
neutral's velocity; new particles carry the primary's weight.

## `[coulomb]`

Binary Coulomb collisions between charged macro-particles in the same cell (Takizuka-Abe).

| key | default | meaning |
|---|---|---|
| `enabled` | true (if the table exists) | |
| `coulomb_log` | 10 | number, or `"auto"`: ln(12π n λ_D³) of the reference plasma (first negative species with a density) |
| `every` | 1 | collide every N steps with an N-times longer time step |
| `pairs = [["a", "b"], ...]` | all pairs | species pairs to collide (`["e", "e"]` for intra-species) |

Unequal weights are handled with the Nanbu-Yonemura acceptance rule (momentum and energy are
conserved on average; exactly for equal weights). The scattering variance requires ν dt ≪ 1;
larger values are clamped to isotropic scattering.

## `[solver]`

| key | default | meaning |
|---|---|---|
| `type` | `mg-pcg` | `mg-pcg` (CG with multigrid preconditioner, BiCGStab fallback), `mg-bicgstab`, `mg` (V-cycles), `jacobi-pcg` |
| `tolerance` | 1e-8 | relative residual |
| `max_iter` | 500 | |
| `cut_cell` | true | Shortley-Weller boundaries at object surfaces |
| `mg_levels`, `mg_pre`, `mg_post`, `mg_coarse_sweeps` | 12, 2, 2, 40 | multigrid controls |

## `[numerics]`

| key | default | meaning |
|---|---|---|
| `shape_order` | 1 | 1 = cloud-in-cell, 2 = triangular-shaped cloud |
| `sort_every` | 20 | particle cell sort interval (steps) |
| `surface_refine` | 2 | sub-sampling per cell for the object surface triangulation |

## `[parallel]`

`decomposition = [px, py, pz]` (default: automatic, minimising halo surface). Each rank needs at
least 3 cells per axis.

## `[output]`

| key | default | meaning |
|---|---|---|
| `dir` | `output` | output directory |
| `fields_every` | 0 | field dumps (VTK) |
| `moments` | true | per-species density, velocity and temperature in field dumps |
| `float32` | true | write field arrays in single precision |
| `particles_every`, `particle_stride` | 0, 100 | particle dumps, every N-th particle |
| `surface_every` | 0 | surface dumps; current densities are averaged over the interval |
| `history_every` | 10 | rows in `history.csv` (currents averaged over the interval) |
| `checkpoint_every` | 0 | checkpoints (always one at the end if > 0) |
| `status_every` | 50 | console status lines |
| `track_every`, `track_count` | 0, 8 | write `tracks.csv` with the trajectories of the first `track_count` particles of each species (rank 0) |
