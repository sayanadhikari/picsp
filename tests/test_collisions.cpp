// Coulomb collisions and ionisation against analytic rates (test-particle mode, no
// self-consistent field):
//  1. temperature equilibration of two electron populations (Spitzer / NRL rate; the binary
//     Takizuka-Abe method reproduces it to ~10% at practical time steps)
//  2. electron-ion temperature relaxation (initial rate)
//  3. electron-impact ionisation rate coefficient for a Maxwellian, step cross section
//  4. tabulated / Lotz cross-section evaluation
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <cstdlib>

#include "Simulation.hpp"
#include "collisions/MCC.hpp"
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

static Config base(double dx) {
  Config c;
  c.cells = {6, 6, 4};
  c.dx = {dx, dx, dx};
  for (auto& f : c.faces) {
    f.field = FieldBC::Periodic;
    f.particle = ParticleBC::Periodic;
  }
  c.electrostatic = false;
  c.output.dir = "test_output_coll";
  c.output.history_every = 0;
  c.output.status_every = 0;
  return c;
}

static SpeciesConfig species(const char* name, double q, double m, double n, double T, double ppc) {
  SpeciesConfig s;
  s.name = name;
  s.charge = q;
  s.mass = m;
  s.density = n;
  s.temperature = T;
  s.ppc = ppc;
  s.inject = false;
  return s;
}

// temperature [eV] and momentum of a species (global)
static double temperature(Simulation& sim, int s, Vec3* mom = nullptr) {
  auto& sp = sim.particles().species[s];
  auto& p = sp.p;
  double b[8] = {0};
  for (std::size_t i = 0; i < p.size(); ++i) {
    b[0] += p.w[i];
    b[1] += p.w[i] * p.vx[i];
    b[2] += p.w[i] * p.vy[i];
    b[3] += p.w[i] * p.vz[i];
    b[4] += p.w[i] * (p.vx[i] * p.vx[i] + p.vy[i] * p.vy[i] + p.vz[i] * p.vz[i]);
  }
  double o[8];
  MPI_Allreduce(b, o, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Vec3 u{o[1] / o[0], o[2] / o[0], o[3] / o[0]};
  if (mom) *mom = Vec3{o[1], o[2], o[3]} * sp.m;
  return sp.m * (o[4] / o[0] - dot(u, u)) / (3 * qe);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  log0("test_collisions\n");
  const double lnL = 10;
  // Spitzer / NRL energy-exchange rate nu^{a|b} [1/s]:  dT_a/dt = nu^{a|b} (T_b - T_a)
  //   nu = 8 sqrt(2 pi)/3 n_b q_a^2 q_b^2 lnL sqrt(m_a m_b) / ((4 pi eps0)^2 (m_a T_b + m_b T_a)^{3/2})
  auto nrl = [&](double ma, double mb, double nb, double Ta, double Tb) {
    double e4 = std::pow(qe, 4) / std::pow(4 * pi * eps0, 2);
    return 8 * std::sqrt(2 * pi) / 3 * nb * e4 * lnL * std::sqrt(ma * mb) / std::pow((ma * Tb + mb * Ta) * qe, 1.5);
  };
  const char* only = std::getenv("ONLY");
  auto run = [&](const char* k) { return !only || std::string(only) == k; };
  // slope of ln(f) over equally spaced samples
  auto log_slope = [](const std::vector<double>& t, const std::vector<double>& f) {
    double tm = 0, lm = 0;
    for (std::size_t k = 0; k < t.size(); ++k) {
      tm += t[k] / t.size();
      lm += std::log(f[k]) / t.size();
    }
    double sxy = 0, sxx = 0;
    for (std::size_t k = 0; k < t.size(); ++k) {
      sxy += (t[k] - tm) * (std::log(f[k]) - lm);
      sxx += (t[k] - tm) * (t[k] - tm);
    }
    return sxy / sxx;
  };
  // 1. e-e equilibration (initial phase: the NRL rate assumes Maxwellian populations)
  if (run("ee")) {
    Config c = base(2e-5);
    c.species = {species("hot", -qe, me, 1e20, 2.0, 300), species("cold", -qe, me, 1e20, 1.0, 300)};
    c.coulomb.enabled = true;
    c.coulomb.log_lambda = lnL;
    c.dt = 5e-13;
    Simulation sim(c, MPI_COMM_WORLD);
    sim.initialise();
    const double nu = nrl(me, me, 1e20, temperature(sim, 0), temperature(sim, 1));
    double KE0 = sim.particles().kinetic_energy(0) + sim.particles().kinetic_energy(1);
    std::vector<double> t, f;
    const int N = 200;
    for (int s = 0; s <= N; ++s) {
      if (s % 20 == 0) {
        t.push_back(s * c.dt);
        f.push_back(temperature(sim, 0) - temperature(sim, 1));
      }
      if (s < N) sim.advance();
    }
    double rate = -log_slope(t, f);
    double KE = sim.particles().kinetic_energy(0) + sim.particles().kinetic_energy(1);
    log0("[coulomb  ] e-e equilibration: initial rate %.3e /s, Spitzer 2 nu = %.3e /s (ratio %.3f); "
         "kinetic energy drift %.1e\n", rate, 2 * nu, rate / (2 * nu), std::fabs(KE / KE0 - 1));
    CHECK(std::fabs(rate / (2 * nu) - 1) < 0.15, "e-e temperature equilibration rate");
    CHECK(std::fabs(KE / KE0 - 1) < 1e-9, "energy conservation (equal weights)");
  }
  // 2. e-i relaxation (all pairs: e-e and i-i collisions keep both species Maxwellian)
  if (run("ei")) {
    Config c = base(2e-5);
    double mi = 25 * me;
    c.species = {species("electron", -qe, me, 1e20, 2.0, 250), species("ion", qe, mi, 1e20, 0.5, 250)};
    c.coulomb.enabled = true;
    c.coulomb.log_lambda = lnL;
    c.dt = 2e-12;
    Simulation sim(c, MPI_COMM_WORLD);
    sim.initialise();
    double Te0 = temperature(sim, 0), Ti0 = temperature(sim, 1);
    std::vector<double> t, f;
    const int N = 300;
    for (int s = 0; s <= N; ++s) {
      if (s % 30 == 0) {
        t.push_back(s * c.dt);
        f.push_back(temperature(sim, 0) - temperature(sim, 1));
      }
      if (s < N) sim.advance();
    }
    double rate = -log_slope(t, f);
    double expected = nrl(me, mi, 1e20, Te0, Ti0) + nrl(mi, me, 1e20, Ti0, Te0);
    log0("[coulomb  ] e-i relaxation (m_i = 25 m_e): initial rate of Te - Ti %.3e /s, Spitzer %.3e /s (ratio %.3f)\n",
         rate, expected, rate / expected);
    CHECK(std::fabs(rate / expected - 1) < 0.15, "e-i relaxation rate");
  }
  // 3. ionisation rate
  if (run("ion")) {
    Config c = base(2e-3);
    c.species = {species("electron", -qe, me, 1e16, 10.0, 600), species("Ar+", qe, 39.95 * amu, 0, 0, 1)};
    c.species[1].load_initial = false;
    c.species[1].weight = 1;
    c.neutrals.enabled = true;
    c.neutrals.density = 1e20;
    c.neutrals.mass = 39.95 * amu;
    c.neutrals.temperature = 0.026;
    CollisionProcess ion;
    ion.species = "electron";
    ion.type = CollisionProcess::Ionization;
    ion.model = CollisionProcess::Step;
    ion.cross_section = 1e-20;
    ion.threshold = 15.0;
    ion.product = "Ar+";
    c.neutrals.processes = {ion};
    c.dt = 1e-10;
    Simulation sim(c, MPI_COMM_WORLD);
    sim.initialise();
    long long ne0 = sim.particles().count_global(0);
    const int N = 100;
    for (int s = 0; s < N; ++s) sim.advance();
    long long ni = sim.particles().count_global(1), ne = sim.particles().count_global(0);
    double T = 10.0, Ei = 15.0;
    double k = 1e-20 * std::sqrt(8 * T * qe / (pi * me)) * (1 + Ei / T) * std::exp(-Ei / T);
    double expected = ne0 * 1e20 * k * N * c.dt;
    log0("[ionize   ] ions created %lld, electrons %lld -> %lld; expected %.0f from k = %.3e m^3/s (ratio %.3f)\n", ni,
         ne0, ne, expected, k, ni / expected);
    CHECK(std::fabs(ni / expected - 1) < 0.1, "ionisation rate coefficient");
    CHECK(ne - ne0 == ni, "one secondary electron per ion");
  }
  // 4. cross-section models
  {
    CollisionProcess t;
    t.model = CollisionProcess::Table;
    t.table = {{0.0, 1e-20}, {10.0, 3e-20}, {20.0, 2e-20}};
    CollisionProcess l;
    l.model = CollisionProcess::Lotz;
    l.cross_section = 4.0e-18 * 1e-4 * 2.6;  // a = 4.5e-14 cm^2 eV^2 order of magnitude
    l.threshold = 15.76;
    bool ok = std::fabs(MonteCarloCollisions::sigma(t, 5.0) - 2e-20) < 1e-30 &&
              std::fabs(MonteCarloCollisions::sigma(t, 15.0) - 2.5e-20) < 1e-30 &&
              MonteCarloCollisions::sigma(t, 30.0) == 2e-20 && MonteCarloCollisions::sigma(l, 10.0) == 0.0 &&
              MonteCarloCollisions::sigma(l, 50.0) > 0;
    log0("[sigma    ] table interpolation and Lotz threshold behaviour: %s\n", ok ? "ok" : "wrong");
    CHECK(ok, "cross-section models");
  }
  log0(failures ? "test_collisions: %d FAILURE(S)\n" : "test_collisions: all checks passed\n", failures);
  MPI_Finalize();
  return failures ? 1 : 0;
}
