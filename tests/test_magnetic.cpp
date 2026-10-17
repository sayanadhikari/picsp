// Spatially varying magnetic fields:
//  1. current loop: exact (elliptic integral) field vs numerical Biot-Savart integration
//  2. dipole field on the equator and axis
//  3. magnetic mirror: adiabatic invariant mu conserved, trapping set by the loss cone
#include <cmath>
#include <cstdio>

#include "Simulation.hpp"
#include "fields/MagneticField.hpp"
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

static Vec3 biot_savart_loop(const MagnetSource& m, const Vec3& p) {
  const int N = 20000;
  const double mu0 = 1.25663706212e-6;
  Vec3 e1, e2, B{0, 0, 0};
  orthonormal_basis(m.axis, e1, e2);
  // orientation: e1 x e2 = axis
  for (int k = 0; k < N; ++k) {
    double t = twopi * (k + 0.5) / N;
    Vec3 x = m.center + (e1 * std::cos(t) + e2 * std::sin(t)) * m.radius;
    Vec3 dl = (e1 * -std::sin(t) + e2 * std::cos(t)) * (m.radius * twopi / N);
    Vec3 r = p - x;
    double rn = norm(r);
    B += cross(dl, r) * (mu0 * m.current / (4 * pi * rn * rn * rn));
  }
  return B;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  log0("test_magnetic\n");
  // 1. loop
  {
    MagnetSource m;
    m.type = MagnetSource::Loop;
    m.center = {0.1, -0.2, 0.05};
    m.axis = normalized(Vec3{0.3, 0.4, 1.0});
    m.radius = 0.5;
    m.current = 1000;
    double maxrel = 0;
    Vec3 pts[5] = {{0.1, -0.2, 0.05}, {0.4, 0.1, 0.2}, {1.0, 0.0, 0.0}, {-0.3, 0.3, 0.9}, {0.2, -0.6, -0.4}};
    for (auto& p : pts) {
      Vec3 a = magnet_field(m, p, 1e-9), b = biot_savart_loop(m, p);
      maxrel = std::max(maxrel, norm(a - b) / norm(b));
    }
    log0("[loop     ] max relative difference vs Biot-Savart: %.2e\n", maxrel);
    CHECK(maxrel < 1e-6, "loop field");
  }
  // 2. dipole
  {
    MagnetSource m;
    m.type = MagnetSource::Dipole;
    m.moment = {0, 0, 8e22};  // Earth-like
    double R = 6.371e6, mu0 = 1.25663706212e-6;
    Vec3 eq = magnet_field(m, {R, 0, 0}, 1), pole = magnet_field(m, {0, 0, R}, 1);
    double B0 = mu0 * 8e22 / (4 * pi * R * R * R);
    log0("[dipole   ] equator Bz = %.4e T (expect %.4e), pole Bz = %.4e T (expect %.4e)\n", eq.z, -B0, pole.z, 2 * B0);
    CHECK(std::fabs(eq.z / -B0 - 1) < 1e-12 && std::fabs(pole.z / (2 * B0) - 1) < 1e-12, "dipole field");
  }
  // 3. magnetic mirror with test particles
  {
    Config c;
    double L = 1.0;
    c.cells = {16, 16, 64};
    c.dx = {0.4 / 16, 0.4 / 16, L / 64};
    c.origin = {-0.2, -0.2, -0.5};
    for (auto& f : c.faces) {
      f.field = FieldBC::Dirichlet;
      f.particle = ParticleBC::Absorb;
    }
    c.electrostatic = false;
    MagnetSource m;
    m.type = MagnetSource::Mirror;
    m.axis = {0, 0, 1};
    m.radius = 0.25;
    m.separation = 0.8;
    m.current = 2e4;
    c.magnets = {m};
    SpeciesConfig e;
    e.name = "electron";
    e.charge = -qe;
    e.mass = me;
    e.weight = 1e-6;
    e.load_initial = false;
    e.inject = false;
    c.species = {e};
    c.output.dir = "test_output_mag";
    c.output.history_every = 0;
    c.output.status_every = 0;
    double Bmin = norm(total_magnetic_field(c, {0, 0, 0}, 1e-6));
    double Bmax = norm(total_magnetic_field(c, {0, 0, 0.4}, 1e-6));
    double wc = qe * Bmax / me;
    c.dt = 0.05 / wc;
    const double v = 1e7;
    double sin2_lc = Bmin / Bmax;
    double a_lc = std::asin(std::sqrt(sin2_lc));
    double angles[2] = {a_lc * 1.25, a_lc * 0.7};  // trapped, passing
    long long left[2];
    double muvar[2], zmx[2];
    for (int k = 0; k < 2; ++k) {
      Simulation sim(c, MPI_COMM_WORLD);
      sim.initialise();
      const auto& d = sim.decomposition();
      Vec3 p{0.01, 0.0, 0.0};
      Vec3 lo = d.local_lo(), hi = d.local_hi();
      if (p.x >= lo.x && p.x < hi.x && p.y >= lo.y && p.y < hi.y && p.z >= lo.z && p.z < hi.z)
        sim.particles().add_particle(0, p, {v * std::sin(angles[k]), 0, v * std::cos(angles[k])}, 1e-6);
      long long nsteps = (long long)(3 * 4 * 0.4 / v / c.dt);
      double mu_min = 1e300, mu_max = 0, zmax = 0;
      for (long long s = 0; s < nsteps; ++s) {
        sim.advance();
        auto& P = sim.particles().species[0].p;
        for (std::size_t i = 0; i < P.size(); ++i) {
          Vec3 x = P.pos(i), vv = P.vel(i);
          Vec3 B = total_magnetic_field(c, x, 1e-6);
          double b = norm(B);
          double vpar = dot(vv, B) / b;
          double mu = (dot(vv, vv) - vpar * vpar) / b;
          mu_min = std::min(mu_min, mu);
          mu_max = std::max(mu_max, mu);
          zmax = std::max(zmax, std::fabs(x.z));
        }
      }
      left[k] = sim.particles().count_global(0);
      double mmin = allreduce_min(mu_min, MPI_COMM_WORLD), mmax = allreduce_max(mu_max, MPI_COMM_WORLD);
      muvar[k] = (mmax - mmin) / mmax;
      zmx[k] = allreduce_max(zmax, MPI_COMM_WORLD);
    }
    double zm_theory = 0;  // mirror point: B(z_m) = B_min / sin^2(alpha)
    {
      double Bm = Bmin / std::pow(std::sin(angles[0]), 2);
      for (double z = 0; z < 0.4; z += 1e-4)
        if (norm(total_magnetic_field(c, {0.0, 0, z}, 1e-6)) >= Bm) {
          zm_theory = z;
          break;
        }
    }
    log0("[mirror   ] B_min = %.4f T, B_max = %.4f T, loss cone %.1f deg\n", Bmin, Bmax, a_lc * 180 / pi);
    log0("[mirror   ] trapped (%.1f deg): mu variation %.2f%%, mirror point %.3f m (theory %.3f m), retained %lld\n",
         angles[0] * 180 / pi, 100 * muvar[0], zmx[0], zm_theory, left[0]);
    log0("[mirror   ] passing (%.1f deg): retained %lld (escapes through the loss cone)\n", angles[1] * 180 / pi, left[1]);
    CHECK(left[0] == 1 && left[1] == 0, "trapping set by the loss cone");
    CHECK(muvar[0] < 0.05, "magnetic moment conserved");
    CHECK(std::fabs(zmx[0] - zm_theory) < 0.03, "mirror point");
  }
  log0(failures ? "test_magnetic: %d FAILURE(S)\n" : "test_magnetic: all checks passed\n", failures);
  MPI_Finalize();
  return failures ? 1 : 0;
}
