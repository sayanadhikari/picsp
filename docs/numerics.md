# Numerical methods

## Time stepping

Standard leapfrog: positions at integer and velocities at half-integer time levels. Each step:

1. gather **E** at the particles and advance velocities with the Boris scheme (uniform **B**),
   then positions; particles whose start or end cell is flagged as near a surface are tested for
   segment / surface intersection by sphere tracing the signed distance function;
2. absorbed particles charge their conductor and are attributed to the nearest surface panel;
   secondary electrons are emitted (Sternglass or constant yield, ion-induced emission);
3. photoemission from sunlit, unshadowed panels;
4. injection through open boundaries, then outer boundary conditions and MPI migration
   (dimension by dimension, so diagonal moves need no extra messages);
5. Monte Carlo collisions (null-collision method);
6. charge deposition (thread-private buffers, ghost-layer summation across ranks);
7. Poisson solve, conductor potentials, **E** = -∇φ + **E**₀.

Heavy species can be sub-cycled (pushed every N steps with N Δt). Particles are counting-sorted by
cell periodically to keep gathers and deposits cache friendly.

## Poisson equation with embedded conductors

-∇²φ = ρ/ε₀ is discretised with the 7-point stencil on the node grid. Nodes inside objects are
Dirichlet nodes at the conductor potential. With `cut_cell = true` free nodes next to a surface use
the Shortley-Weller stencil, where the neighbour in direction d is replaced by the boundary point
at fractional distance θ (found by bisection on the exact signed distance); nodes with θ < 0.1 are
snapped onto the surface. This makes the boundary treatment second order (third-order potential
convergence observed for the concentric-sphere test), compared with first order for a staircase.

Neumann rows are scaled by the control-volume fraction (1/2 per Neumann face) so that the
staircase operator is symmetric. Fully periodic or all-Neumann problems are singular: the
right-hand side is projected onto the range and the solution has zero mean.

### Solver

Conjugate gradients preconditioned by one geometric multigrid V(2,2) cycle:

* red-black Gauss-Seidel smoothing (symmetric ordering), cut-cell coefficients on the finest level;
* full-weighting restriction R = Pᵀ/8 and trilinear prolongation P;
* rediscretised coarse operators with injected Dirichlet masks;
* coarsening as long as every rank keeps at least two cells per axis; the coarsest level is
  solved with a number of sweeps proportional to its size.

The Shortley-Weller operator is mildly non-symmetric; CG converges well in practice and the solver
switches automatically to BiCGStab if it stalls. Typical cost: 7 to 8 iterations for a 1e-8
relative residual, independent of the grid size; with warm starts from the previous step.

### Floating conductors (capacitance matrix)

With conductors k = 1..K the potential is written as φ = φ₀ + Σₖ Vₖ ψₖ, where φ₀ solves Poisson
with all conductors grounded and ψₖ are vacuum solutions with conductor k at 1 V (computed once).
The induced charge on conductor i is measured with a conservative discrete Gauss law (regular
stencil flux out of the object nodes plus the cut-cell excess of adjacent free nodes), giving
the capacitance matrix C_ik = Q_i[ψ_k]. Charge conservation

  Q_i = Q_ind,i[φ₀] + Σₖ C_ik Vₖ + C_ext,i V_i

is solved for the floating potentials each step (fixed conductors enter the right-hand side), so
only one Poisson solve per step is needed. Q_i accumulates absorbed and emitted charge and
external currents. Objects with the same `group` form one conductor.

### Electric field near surfaces

**E** is obtained by central differences; at free nodes adjacent to a surface the derivative of
the quadratic through the boundary point (consistent with the Shortley-Weller stencil) is used, and
**E** at object nodes bordering free space is extrapolated from the free neighbours so that
particles in cut cells see the correct near-surface field.

## Geometry

Objects are signed distance functions: analytic primitives (Quilez-type exact SDFs), CSG via
min/max, rigid transforms, and triangle meshes. Meshes use a bounding-volume hierarchy for
closest-point queries; inside/outside is the generalised winding number (sum of signed solid angles
over 4π, Van Oosterom-Strackee) evaluated with a dipole far-field approximation on the BVH
(Barill et al. 2018). This is robust to overlapping shells, small gaps and flipped triangles.

Surfaces are triangulated by marching tetrahedra on a refined sampling of the SDF. The iso-surface
is extracted slightly outside the body and projected back onto the zero level, so sheets thinner than
the sampling are still captured. Surface current densities are reported per panel, averaged over
the panels of each grid cell.

Shadowing for photoemission is computed by sphere tracing from every panel towards the sun.

## Particles

* Shapes: cloud-in-cell (linear) or triangular-shaped cloud (quadratic) B-splines, same for
  deposition and gather (momentum conserving).
* Velocity sampling: Maxwellian; bi-Maxwellian and ring distributions about B (or a given axis);
  kappa distributions as multivariate Student-t (ν = 2κ - 1) scaled to the kinetic temperature.
* Open boundaries inject the one-sided flux n⟨v_n H(v_n)⟩. For drifting Maxwellians the
  flux-weighted normal velocity is sampled exactly by a tabulated inverse CDF; other distributions use
  rejection sampling. Injection times are randomised within the step.
* Emitted particles (photo, secondary) are drawn from a flux-weighted half-Maxwellian about the local
  surface normal.

## Spatially varying magnetic field

Analytic sources (`[[fields.magnet]]`) are summed with the uniform B0 and stored once per node
(three arrays for output plus an interleaved copy for the pusher). The circular loop uses the
exact off-axis expressions with complete elliptic integrals K(k) and E(k),

    B_z = μ0 I / (2π √((a+ρ)² + z²)) [K + (a² − ρ² − z²) / ((a−ρ)² + z²) E]
    B_ρ = μ0 I z / (2π ρ √((a+ρ)² + z²)) [−K + (a² + ρ² + z²) / ((a−ρ)² + z²) E]

(checked against direct Biot-Savart quadrature to 1e-13). Solenoids and mirrors are sums of
loops, the dipole and the straight wire use their closed forms. The Boris rotation vectors
t = qB dt/2m and s = 2t/(1+t²) are then formed per particle from the gathered B, which preserves
the magnetic moment in slowly varying fields (the mirror test conserves μ to 0.3% over a bounce).

With `self_consistent = false` the Poisson solve is skipped and particles move in the external
fields only (test-particle mode), which is useful for orbit studies, collision benchmarks and
swarm simulations.

## Monte Carlo collisions

Null-collision method: ν_max = n_g max_g [g Σ σ_k(½ μ g²)] is found by sampling the total
cross section on 4000 relative speeds up to the largest particle speed (plus a margin for the
neutral thermal speed), multiplied by 1.05. A fraction 1 − exp(−ν_max dt) of the particles are
candidates; for each, a neutral partner is drawn from the (drifting) Maxwellian background and
a process is chosen with probability n_g σ_k(E) g / ν_max (the remainder is the null process).

* elastic: isotropic scattering in the centre-of-mass frame
* charge exchange: the ion takes the neutral's velocity
* excitation: the centre-of-mass energy is reduced by the threshold, then isotropic scattering
* ionisation: the energy above the threshold is split randomly between the scattered primary and
  a new secondary electron, an ion is created with the neutral's velocity; both inherit the
  primary's weight. New particles are buffered per thread and appended after the loop.

Cross sections are constant, a step at a threshold, the Lotz form a ln(E/I)/(E I), or a table
with linear interpolation (inline or from a file).

## Coulomb collisions

Binary collisions in the Takizuka-Abe (1977) formulation. Each species is indexed by cell
(counting sort). Within a cell, intra-species collisions pair shuffled particles; an odd count
forms one triplet whose three pairs collide with half the time step. For a species pair the
more numerous population is traversed once and partners are taken cyclically from the other.
For a pair with relative velocity u the scattering angle follows from

    tan(θ/2) = δ,   δ ~ N(0, ⟨δ²⟩),   ⟨δ²⟩ = q_a² q_b² n lnΛ dt / (8π ε0² μ² u³)

with n the smaller of the two densities in the cell, and the azimuth uniformly random. The
relative velocity is rotated and both velocities are updated about the centre of mass, which
conserves momentum and energy exactly for equal weights. For unequal weights each particle of
the pair is updated with probability w_other / max(w_a, w_b) (Nanbu-Yonemura), conserving them
on average. If ⟨δ²⟩ exceeds one the pair is scattered isotropically.

Accuracy: for the e-e and e-i relaxation benchmarks the measured rates lie 5 to 10% below the
Spitzer rates at practical time steps (ν dt ≈ 0.01 to 0.05), a known property of the binary
method at finite time step; the energy drift is at round-off level.

## Parallelisation

3D Cartesian MPI decomposition (automatic factorisation minimising the halo surface). Fields have
two ghost layers; deposition sums ghost contributions into owners, followed by ghost refresh.
Particles migrate dimension by dimension with neighbour exchange. OpenMP threads share each rank's
particles (static schedules, thread-private deposition buffers and random streams). With one
thread per rank, runs are reproducible for a fixed rank count and checkpoints restart bit-exactly;
with several threads, OpenMP sum reductions may combine partial sums in a run-dependent order,
which perturbs results at the round-off level. The cell sort produces the same particle order for
any thread count.

## Output

Parallel VTK XML: one ImageData piece per rank plus a `.pvti` master for fields, PolyData pieces
plus `.pvtp` masters for particles and surfaces, raw-appended binary with 64-bit headers, and
`.pvd` collections for time series.
