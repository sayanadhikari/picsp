// Particle pusher verification with test particles (negligible self-field):
//  1. gyration in a uniform B: Larmor radius and gyro-frequency (Boris accuracy)
//  2. E x B drift velocity
//  3. migration across MPI ranks / periodic wrap conserves particles
//  4. reflecting and absorbing walls
#include <cmath>
#include <cstdio>

#include "Simulation.hpp"
#include "util/Constants.hpp"
#include "util/Log.hpp"

using namespace picsp;
using namespace picsp::constants;

static int failures = 0;
#define CHECK(cond, ...)                     \
  do {                                       \
    if (!(cond)) {                           \
      ++failures;                            \
      if (world_rank() == 0) {               \
        std::printf("  FAIL: " __VA_ARGS__); \
        std::printf("\n");                   \
      }                                      \
    }                                        \
  } while (0)

static Config base(double L, int N) {
  Config c;
  c.cells = {N, N, N};
  c.dx = {L / N, L / N, L / N};
  c.origin = {-L / 2, -L / 2, -L / 2};
  for (auto& f : c.faces) {
    f.field = FieldBC::Periodic;
    f.particle = ParticleBC::Periodic;
  }
  SpeciesConfig e;
  e.name = "electron";
  e.charge = -qe;
  e.mass = me;
  e.density = 0;
  e.weight = 1e-6;  // tiny weight: no self-field
  e.load_initial = false;
  e.inject = false;
  c.species = {e};
  c.output.dir = "test_output";
  c.output.history_every = 0;
  c.output.status_every = 0;
  c.solver.tol = 1e-12;
  return c;
}

// find the particle (it may live on any rank) and broadcast its state
static bool get_particle(Simulation& sim, Vec3& x, Vec3& v) {
  double buf[7] = {0}, out[7];
  auto& p = sim.particles().species[0].p;
  double have = p.size() > 0 ? 1.0 : 0.0;
  if (have > 0) {
    x = p.pos(0);
    v = p.vel(0);
    buf[0] = x.x; buf[1] = x.y; buf[2] = x.z; buf[3] = v.x; buf[4] = v.y; buf[5] = v.z;
  }
  buf[6] = have;
  MPI_Allreduce(buf, out, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  x = {out[0], out[1], out[2]};
  v = {out[3], out[4], out[5]};
  return out[6] == 1.0;
}

static void add_on_owner(Simulation& sim, const Vec3& x, const Vec3& v) {
  const auto& d = sim.decomposition();
  Vec3 lo = d.local_lo(), hi = d.local_hi();
  if (x.x >= lo.x && x.x < hi.x && x.y >= lo.y && x.y < hi.y && x.z >= lo.z && x.z < hi.z)
    sim.particles().add_particle(0, x, v, 1e-6);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  log0("test_particles on %d rank(s)\n", size);

  // ---------------- 1. gyration ----------------
  {
    double B = 1e-3, v0 = 2e5;
    double wc = qe * B / me, rL = v0 / wc;
    Config c = base(8 * rL, 16);
    c.B0 = {0, 0, B};
    c.dt = 0.02 / wc;
    Simulation sim(c, MPI_COMM_WORLD);
    sim.initialise();
    // start at (rL, 0, 0) with v = (0, -v0... ) electron gyrates counter-clockwise about +z:
    // for q<0, omega = +|q|B/m about +z  =>  x = rL cos(wt), y = rL sin(wt) with v = (0, v0, 0)
    // NB: leapfrog velocity is at t - dt/2; we place it directly (small phase offset).
    add_on_owner(sim, {rL, 0, 0}, {0, v0, 0});
    long long nper = (long long)std::llround(constants::twopi / (wc * c.dt));
    double rmin = 1e9, rmax = 0;
    Vec3 x, v;
    double angle = 0, prev = 0;
    for (long long s = 0; s < 4 * nper; ++s) {
      sim.advance();
      get_particle(sim, x, v);
      double r = std::sqrt(x.x * x.x + x.y * x.y);
      rmin = std::min(rmin, r);
      rmax = std::max(rmax, r);
      double a = std::atan2(x.y, x.x);
      double da = a - prev;
      if (da < -pi) da += twopi;
      if (da > pi) da -= twopi;
      angle += da;
      prev = a;
    }
    double w_meas = angle / (4 * nper * c.dt);
    double speed = norm(v);
    log0("[gyration] r/rL in [%.5f, %.5f], |v|/v0 = %.8f, w/wc = %.6f\n", rmin / rL, rmax / rL, speed / v0, w_meas / wc);
    CHECK(std::fabs(speed / v0 - 1) < 1e-10, "Boris conserves |v| in pure B");
    CHECK(std::fabs(w_meas / wc - 1) < 2e-3, "gyro-frequency");
    CHECK(rmax / rL < 1.02 && rmin / rL > 0.98, "Larmor radius");
  }
  // ---------------- 2. E x B drift ----------------
  {
    double B = 1e-3, E = 100.0;
    double wc = qe * B / me, vd = E / B;
    Config c = base(2.0, 16);
    c.B0 = {0, 0, B};
    c.E0 = {0, E, 0};
    c.dt = 0.05 / wc;
    Simulation sim(c, MPI_COMM_WORLD);
    sim.initialise();
    add_on_owner(sim, {0, 0, 0}, {0, 0, 0});
    long long nper = (long long)std::llround(constants::twopi / (wc * c.dt));
    Vec3 x0, v, x1;
    get_particle(sim, x0, v);
    long long n = 10 * nper;
    double unwrapped = 0, last = x0.x;
    for (long long s = 0; s < n; ++s) {
      sim.advance();
      get_particle(sim, x1, v);
      double dxs = x1.x - last;
      if (dxs < -1.0) dxs += 2.0;
      if (dxs > 1.0) dxs -= 2.0;
      unwrapped += dxs;
      last = x1.x;
    }
    double vmeas = unwrapped / (n * c.dt);
    log0("[E x B   ] v_drift/(E/B) = %.5f (particle crossed periodic boundary & ranks)\n", vmeas / vd);
    CHECK(std::fabs(vmeas / vd - 1) < 5e-3, "E x B drift");
    bool found = get_particle(sim, x1, v);
    CHECK(found, "particle conserved across migration");
  }
  // ---------------- 3. walls ----------------
  {
    Config c = base(1.0, 16);
    c.faces[0].field = c.faces[1].field = FieldBC::Dirichlet;
    c.faces[0].particle = ParticleBC::Reflect;
    c.faces[1].particle = ParticleBC::Absorb;
    c.dt = 1e-9;
    Simulation sim(c, MPI_COMM_WORLD);
    sim.initialise();
    add_on_owner(sim, {0.3, 0.1, 0.1}, {-2e7, 0, 0});  // reflects at x=-0.5 then absorbed at x=+0.5
    Vec3 x, v;
    bool reflected = false;
    int steps = 0;
    while (get_particle(sim, x, v) && steps < 1000) {
      if (v.x > 0) reflected = true;
      sim.advance();
      ++steps;
    }
    bool alive = get_particle(sim, x, v);
    log0("[walls   ] reflected=%d, absorbed after %d steps (expected ~%d)\n", int(reflected), steps,
         int(std::ceil((0.8 + 1.0) / (2e7 * 1e-9))));
    CHECK(reflected && !alive, "reflect then absorb");
    CHECK(std::abs(steps - 90) <= 2, "time of flight");
  }
  log0(failures ? "test_particles: %d FAILURE(S)\n" : "test_particles: all checks passed\n", failures);
  MPI_Finalize();
  return failures ? 1 : 0;
}
