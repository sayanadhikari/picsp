// Poisson solver verification:
//  1. manufactured solution, Dirichlet walls (2nd-order convergence)
//  2. fully periodic (singular) problem
//  3. mixed Neumann/Dirichlet
//  4. concentric-sphere capacitor (embedded objects + CSG + cut cells vs staircase)
#include <cmath>
#include <cstdio>
#include <functional>

#include "fields/FieldSolver.hpp"
#include "parallel/Decomposition.hpp"
#include "util/Constants.hpp"
#include "util/Log.hpp"

using namespace picsp;
using constants::eps0;
using constants::pi;

static int failures = 0;
#define CHECK(cond, ...)                         \
  do {                                           \
    if (!(cond)) {                               \
      ++failures;                                \
      if (world_rank() == 0) {                   \
        std::printf("  FAIL: " __VA_ARGS__);     \
        std::printf("\n");                       \
      }                                          \
    }                                            \
  } while (0)

static Config base_config(int N, FieldBC fx, FieldBC fy, FieldBC fz) {
  Config c;
  c.cells = {N, N, N};
  c.dx = {1.0 / N, 1.0 / N, 1.0 / N};
  c.origin = {0, 0, 0};
  FieldBC b[3] = {fx, fy, fz};
  for (int a = 0; a < 3; ++a)
    for (int s = 0; s < 2; ++s) {
      c.faces[2 * a + s].field = b[a];
      c.faces[2 * a + s].particle = b[a] == FieldBC::Periodic ? ParticleBC::Periodic : ParticleBC::Absorb;
    }
  c.solver.tol = 1e-10;
  c.solver.max_iter = 200;
  return c;
}

// returns max error over owned nodes, number of iterations
static double run_manufactured(const Config& c, std::function<double(Vec3)> phi_ex, std::function<double(Vec3)> lap_ex,
                               int* iters) {
  Decomposition d = Decomposition::create(c, MPI_COMM_WORLD);
  GridGeom g = GridGeom::from(d);
  Scene scene;
  FieldSolver fs;
  fs.setup(c, g, scene);
  const int G = FieldSolver::kGhosts;
  const auto& n = fs.rho.size();
  for (int k = 0; k < n[2]; ++k)
    for (int j = 0; j < n[1]; ++j)
      for (int i = 0; i < n[0]; ++i) {
        Vec3 p{(g.c0[0] + i - G) * c.dx[0], (g.c0[1] + j - G) * c.dx[1], (g.c0[2] + k - G) * c.dx[2]};
        fs.rho(i, j, k) = -eps0 * lap_ex(p);
      }
  SolveStats st = fs.solve(0.0);
  *iters = st.iterations;
  double err = 0;
  for (int k = G; k < G + g.n1[2] - g.n0[2]; ++k)
    for (int j = G; j < G + g.n1[1] - g.n0[1]; ++j)
      for (int i = G; i < G + g.n1[0] - g.n0[0]; ++i) {
        Vec3 p{(g.c0[0] + i - G) * c.dx[0], (g.c0[1] + j - G) * c.dx[1], (g.c0[2] + k - G) * c.dx[2]};
        err = std::max(err, std::fabs(fs.phi(i, j, k) - phi_ex(p)));
      }
  MPI_Comm_free(&d.comm);
  return allreduce_max(err, MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  log0("test_poisson on %d rank(s)\n", size);

  // 1. Dirichlet manufactured solution
  {
    auto ex = [](Vec3 p) { return std::sin(pi * p.x) * std::sin(pi * p.y) * std::sin(pi * p.z); };
    auto lap = [&](Vec3 p) { return -3 * pi * pi * ex(p); };
    double e[2];
    int it[2];
    for (int r = 0; r < 2; ++r) {
      Config c = base_config(r == 0 ? 32 : 64, FieldBC::Dirichlet, FieldBC::Dirichlet, FieldBC::Dirichlet);
      e[r] = run_manufactured(c, ex, lap, &it[r]);
    }
    log0("[dirichlet] err(32)=%.3e err(64)=%.3e order=%.2f iters=%d,%d\n", e[0], e[1], std::log2(e[0] / e[1]), it[0], it[1]);
    CHECK(std::log2(e[0] / e[1]) > 1.8, "dirichlet convergence order");
    CHECK(it[1] < 30, "dirichlet MG-PCG iterations");
  }
  // 2. periodic
  {
    auto ex = [](Vec3 p) { return std::cos(2 * pi * p.x) * std::cos(2 * pi * p.y) + 0.5 * std::sin(2 * pi * p.z); };
    auto lap = [](Vec3 p) { return -8 * pi * pi * std::cos(2 * pi * p.x) * std::cos(2 * pi * p.y) - 2 * pi * pi * std::sin(2 * pi * p.z); };
    double e[2];
    int it[2];
    for (int r = 0; r < 2; ++r) {
      Config c = base_config(r == 0 ? 32 : 64, FieldBC::Periodic, FieldBC::Periodic, FieldBC::Periodic);
      e[r] = run_manufactured(c, ex, lap, &it[r]);
    }
    log0("[periodic ] err(32)=%.3e err(64)=%.3e order=%.2f iters=%d,%d\n", e[0], e[1], std::log2(e[0] / e[1]), it[0], it[1]);
    CHECK(std::log2(e[0] / e[1]) > 1.8, "periodic convergence order");
  }
  // 3. Neumann in x, Dirichlet y,z, periodic... (mixed)
  {
    auto ex = [](Vec3 p) { return std::cos(pi * p.x) * std::sin(pi * p.y) * std::sin(pi * p.z); };
    auto lap = [&](Vec3 p) { return -3 * pi * pi * ex(p); };
    double e[2];
    int it[2];
    for (int r = 0; r < 2; ++r) {
      Config c = base_config(r == 0 ? 32 : 64, FieldBC::Neumann, FieldBC::Dirichlet, FieldBC::Dirichlet);
      e[r] = run_manufactured(c, ex, lap, &it[r]);
    }
    log0("[neumann  ] err(32)=%.3e err(64)=%.3e order=%.2f iters=%d,%d\n", e[0], e[1], std::log2(e[0] / e[1]), it[0], it[1]);
    CHECK(std::log2(e[0] / e[1]) > 1.7, "neumann convergence order");
  }
  // 4. concentric spheres: inner sphere R1 at 1 V, grounded spherical cavity R2
  for (int cut = 0; cut < 2; ++cut) {
    double errs[2];
    for (int r = 0; r < 2; ++r) {
      int N = r == 0 ? 32 : 64;
      Config c = base_config(N, FieldBC::Dirichlet, FieldBC::Dirichlet, FieldBC::Dirichlet);
      c.origin = {-0.5, -0.5, -0.5};
      c.solver.cut_cell = cut;
      double R1 = 0.13, R2 = 0.43;
      ObjectConfig inner;
      inner.name = "inner";
      inner.bias = BiasMode::Fixed;
      inner.potential = 1.0;
      PartConfig s;
      s.shape = ShapeKind::Sphere;
      s.radius = R1;
      inner.parts.push_back(s);
      ObjectConfig outer;
      outer.name = "shell";
      outer.bias = BiasMode::Fixed;
      outer.potential = 0.0;
      PartConfig b;
      b.shape = ShapeKind::Box;
      b.size = {2, 2, 2};
      PartConfig hole;
      hole.shape = ShapeKind::Sphere;
      hole.radius = R2;
      hole.op = CsgOp::Subtract;
      outer.parts = {b, hole};
      c.objects = {inner, outer};
      Decomposition d = Decomposition::create(c, MPI_COMM_WORLD);
      GridGeom g = GridGeom::from(d);
      Scene scene;
      scene.build(c.objects);
      FieldSolver fs;
      fs.setup(c, g, scene);
      fs.rho.fill(0.0);
      fs.solve(0.0);
      double Cex = 4 * pi * eps0 / (1 / R1 - 1 / R2);
      double C = fs.capacitance()[0];
      // potential error between the spheres
      const int G = FieldSolver::kGhosts;
      double err = 0;
      for (int k = G; k < G + g.n1[2] - g.n0[2]; ++k)
        for (int j = G; j < G + g.n1[1] - g.n0[1]; ++j)
          for (int i = G; i < G + g.n1[0] - g.n0[0]; ++i) {
            Vec3 p{c.origin.x + (g.c0[0] + i - G) * c.dx[0], c.origin.y + (g.c0[1] + j - G) * c.dx[1],
                   c.origin.z + (g.c0[2] + k - G) * c.dx[2]};
            double rr = norm(p);
            if (rr < R1 + 0.02 || rr > R2 - 0.02) continue;
            double ex = (1 / rr - 1 / R2) / (1 / R1 - 1 / R2);
            err = std::max(err, std::fabs(fs.phi(i, j, k) - ex));
          }
      err = allreduce_max(err, MPI_COMM_WORLD);
      errs[r] = err;
      log0("[spheres %s N=%d] C=%.4f pF (exact %.4f, %+.2f%%), max|phi err|=%.3e\n", cut ? "cut-cell " : "staircase", N,
           C * 1e12, Cex * 1e12, 100 * (C / Cex - 1), err);
      if (cut && N == 64) {
        CHECK(std::fabs(C / Cex - 1) < 0.03, "cut-cell capacitance within 3%%");
        CHECK(err < 0.01, "cut-cell potential error < 1%% of bias");
      }
      MPI_Comm_free(&d.comm);
    }
    if (cut) {
      log0("[spheres cut-cell] potential error order = %.2f\n", std::log2(errs[0] / errs[1]));
      CHECK(std::log2(errs[0] / errs[1]) > 1.3, "cut-cell convergence order");
    }
  }
  log0(failures ? "test_poisson: %d FAILURE(S)\n" : "test_poisson: all checks passed\n", failures);
  MPI_Finalize();
  return failures ? 1 : 0;
}
