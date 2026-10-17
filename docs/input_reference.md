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
| `stl` | `file`, `scale` (e.g. `1e-3` for mm), `center` (translation) |

Every part can be rotated with `rotate = [ax, ay, az]` (degrees) about its centre.
STL files may be ASCII or binary; overlapping shells and inconsistent orientation are handled by the
generalised winding number. Features thinner than a cell (panels, booms) are resolved as one-node-thick
conductor sheets in the field solve and exactly for particle collisions.

## `[collisions]`

| key | meaning |
|---|---|
| `neutral_density` | neutral density (m^-3) |
| `neutral_temperature` | neutral temperature (eV) |
| `neutral_mass_amu` | neutral mass |
| `neutral_drift` | neutral drift velocity (m/s) |
| `[[collisions.process]]` | `species`, `type = "elastic" \| "charge_exchange"`, `cross_section` (m²) |

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
