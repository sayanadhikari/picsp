// Physics/infrastructure verification:
//  1. velocity samplers: Maxwellian, bi-Maxwellian, kappa, ring moments
//  2. Monte Carlo collisions: momentum relaxation rates (elastic, charge exchange)
//  3. checkpoint / restart reproduces an uninterrupted run
#include <cmath>
#include <cstdio>
#include <filesystem>

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

static Config box(int N, double L) {
  Config c;
  c.cells = {N, N, N};
  c.dx = {L / N, L / N, L / N};
  for (auto& f : c.faces) {
    f.field = FieldBC::Periodic;
    f.particle = ParticleBC::Periodic;
  }
  c.output.dir = "test_output_physics";
  c.output.history_every = 0;
  c.output.status_every = 0;
  return c;
}

static void mean_velocity(Simulation& sim, int s, Vec3& mean, double& w) {
  auto& p = sim.particles().species[s].p;
  double b[4] = {0, 0, 0, 0}, o[4];
  for (std::size_t i = 0; i < p.size(); ++i) {
    b[0] += p.w[i] * p.vx[i];
    b[1] += p.w[i] * p.vy[i];
    b[2] += p.w[i] * p.vz[i];
    b[3] += p.w[i];
  }
  MPI_Allreduce(b, o, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  mean = Vec3{o[0], o[1], o[2]} / o[3];
  w = o[3];
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  log0("test_physics\n");

  // ---------------- 1. velocity samplers ----------------
  {
    Rng rng(99, world_rank());
    const int N = 400000;
    auto moments = [&](const VelocitySampler& vs, Vec3& var, double& kurt) {
      double s2[3] = {0, 0, 0}, s4 = 0;
      for (int n = 0; n < N; ++n) {
        Vec3 v = vs.sample_thermal(rng);
        for (int a = 0; a < 3; ++a) s2[a] += v[a] * v[a];
        s4 += v.x * v.x * v.x * v.x;
      }
      var = Vec3{s2[0], s2[1], s2[2]} / double(N);
      kurt = (s4 / N) / (var.x * var.x);
    };
    SpeciesConfig s;
    s.mass = me;
    s.temperature = 2.0;
    double vt2 = 2.0 * qe / me;
    VelocitySampler m;
    m.init(s, {0, 0, 0});
    Vec3 var;
    double k;
    moments(m, var, k);
    log0("[sampler  ] Maxwellian   <vx^2>/vt^2 = %.4f %.4f %.4f, kurtosis %.3f (3)\n", var.x / vt2, var.y / vt2, var.z / vt2, k);
    CHECK(std::fabs(var.x / vt2 - 1) < 0.01 && std::fabs(k - 3) < 0.05, "Maxwellian moments");
    s.dist = VelocityDist::Kappa;
    s.kappa = 4;
    VelocitySampler kap;
    kap.init(s, {0, 0, 0});
    moments(kap, var, k);
    double kex = 3.0 * (2 * 4.0 - 3) / (2 * 4.0 - 5);  // Student-t kurtosis 3 + 6/(nu-4), nu = 2 kappa - 1
    log0("[sampler  ] kappa=4      <vx^2>/vt^2 = %.4f, kurtosis %.3f (%.3f)\n", var.x / vt2, k, kex);
    CHECK(std::fabs(var.x / vt2 - 1) < 0.03, "kappa temperature");
    CHECK(std::fabs(k / kex - 1) < 0.15, "kappa suprathermal tail (kurtosis)");
    s.dist = VelocityDist::BiMaxwellian;
    s.t_par = 4.0;
    s.t_perp = 1.0;
    VelocitySampler bi;
    bi.init(s, {0, 0, 1e-5});
    moments(bi, var, k);
    log0("[sampler  ] bi-Maxwellian T_perp = %.3f, %.3f eV, T_par = %.3f eV\n", var.x * me / qe, var.y * me / qe, var.z * me / qe);
    CHECK(std::fabs(var.z * me / qe / 4.0 - 1) < 0.02 && std::fabs(var.x * me / qe - 1) < 0.02, "bi-Maxwellian");
  }

  // ---------------- 2. Monte Carlo collisions ----------------
  for (int type = 0; type < 2; ++type) {
    Config c = box(8, 0.08);
    SpeciesConfig sp;
    sp.name = type == 0 ? "electron" : "ion";
    sp.charge = type == 0 ? -qe : qe;
    sp.mass = type == 0 ? me : 40 * amu;
    sp.density = 1e5;  // negligible self-field
    sp.temperature = 0;
    sp.dist = VelocityDist::Cold;
    double v0 = type == 0 ? 1e6 : 1e4;
    sp.drift = {v0, 0, 0};
    sp.ppc = 64;
    sp.inject = false;
    c.species = {sp};
    c.neutrals.enabled = true;
    c.neutrals.density = 1e20;
    c.neutrals.mass = 40 * amu;
    c.neutrals.temperature = type == 0 ? 0.0 : 1e-6;
    double sigma = 1e-19;
    CollisionProcess cp;
    cp.species = sp.name;
    cp.type = type == 0 ? CollisionProcess::Elastic : CollisionProcess::ChargeExchange;
    cp.cross_section = sigma;
    c.neutrals.processes.push_back(cp);
    double nu = c.neutrals.density * sigma * v0;
    c.dt = 0.01 / nu;
    c.steps = 100;
    Simulation sim(c, MPI_COMM_WORLD);
    sim.initialise();
    Vec3 m0, m1;
    double w;
    mean_velocity(sim, 0, m0, w);
    for (int s = 0; s < 100; ++s) sim.advance();
    mean_velocity(sim, 0, m1, w);
    // electrons (heavy neutrals, isotropic): <vx> = v0 exp(-nu t); CX: same law with nu = n sigma v0
    double expected = std::exp(-nu * 100 * c.dt);
    log0("[MCC      ] %-17s <vx>/v0 = %.4f, expected exp(-nu t) = %.4f\n", type == 0 ? "elastic (e-Ar)" : "charge exch. (Ar+)",
         m1.x / m0.x, expected);
    CHECK(std::fabs(m1.x / m0.x - expected) < 0.02, "MCC relaxation rate");
  }

  // ---------------- 3. checkpoint / restart ----------------
  {
    Config c = box(16, 0.02);
    c.faces[0].field = c.faces[1].field = FieldBC::Dirichlet;
    c.faces[0].particle = c.faces[1].particle = ParticleBC::Open;
    SpeciesConfig e;
    e.name = "electron";
    e.charge = -qe;
    e.mass = me;
    e.density = 1e14;
    e.temperature = 1;
    e.ppc = 8;
    SpeciesConfig i = e;
    i.name = "ion";
    i.charge = qe;
    i.mass = amu;
    i.temperature = 0.1;
    c.species = {e, i};
    ObjectConfig o;
    o.name = "ball";
    PartConfig s;
    s.shape = ShapeKind::Sphere;
    s.center = {0.01, 0.01, 0.01};
    s.radius = 0.003;
    o.parts = {s};
    c.objects = {o};
    c.dt = 1e-11;
    c.steps = 40;
    c.output.dir = "test_output_ckpt";
    c.output.checkpoint_every = 20;
    double ke_full, v_full;
    long long n_full;
    {
      Simulation sim(c, MPI_COMM_WORLD);
      sim.run();
      ke_full = sim.particles().kinetic_energy(0);
      n_full = sim.particles().count_global(0);
      v_full = sim.fields().conductors()[0].V;
    }
    // restart from the step-20 checkpoint: overwrite it by a run of 20 steps first
    c.steps = 20;
    c.output.dir = "test_output_ckpt2";
    {
      Simulation sim(c, MPI_COMM_WORLD);
      sim.run();
    }
    c.steps = 40;
    c.restart_from = "test_output_ckpt2/checkpoint";
    c.output.checkpoint_every = 0;
    Simulation sim(c, MPI_COMM_WORLD);
    sim.run();
    double ke = sim.particles().kinetic_energy(0);
    long long n = sim.particles().count_global(0);
    double v = sim.fields().conductors()[0].V;
    log0("[restart  ] N %lld vs %lld, KE rel diff %.2e, V %.6f vs %.6f\n", n, n_full, std::fabs(ke / ke_full - 1), v, v_full);
    CHECK(n == n_full && std::fabs(ke / ke_full - 1) < 1e-9 && std::fabs(v - v_full) < 1e-9, "restart reproduces the run");
  }
  log0(failures ? "test_physics: %d FAILURE(S)\n" : "test_physics: all checks passed\n", failures);
  MPI_Finalize();
  return failures ? 1 : 0;
}
