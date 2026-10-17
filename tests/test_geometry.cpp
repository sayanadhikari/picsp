// Geometry verification: primitives, CSG, STL import (winding numbers, overlapping
// shells), thin-feature segment detection and surface triangulation area.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <vector>

#include "Simulation.hpp"
#include "geometry/Geometry.hpp"
#include "util/Constants.hpp"
#include "util/Log.hpp"
#include "util/Random.hpp"

using namespace picsp;
using constants::pi;

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

static void write_stl(const std::string& path, const std::vector<Triangle>& tris) {
  std::ofstream f(path, std::ios::binary);
  char header[80] = "picsp3d test";
  f.write(header, 80);
  std::uint32_t n = tris.size();
  f.write(reinterpret_cast<char*>(&n), 4);
  for (auto& t : tris) {
    float buf[12] = {0, 0, 0, float(t.a.x), float(t.a.y), float(t.a.z), float(t.b.x), float(t.b.y), float(t.b.z),
                     float(t.c.x), float(t.c.y), float(t.c.z)};
    f.write(reinterpret_cast<char*>(buf), 48);
    std::uint16_t attr = 0;
    f.write(reinterpret_cast<char*>(&attr), 2);
  }
}

// icosphere with outward orientation
static std::vector<Triangle> icosphere(double R, int sub) {
  const double t = (1 + std::sqrt(5.0)) / 2;
  std::vector<Vec3> v = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t},
                         {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
  int f[20][3] = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4},
                  {11, 10, 2}, {10, 7, 6}, {7, 1, 8}, {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8},
                  {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
  std::vector<Triangle> tris;
  for (auto& q : f) tris.push_back({normalized(v[q[0]]), normalized(v[q[1]]), normalized(v[q[2]])});
  for (int s = 0; s < sub; ++s) {
    std::vector<Triangle> nt;
    for (auto& T : tris) {
      Vec3 ab = normalized((T.a + T.b) * 0.5), bc = normalized((T.b + T.c) * 0.5), ca = normalized((T.c + T.a) * 0.5);
      nt.push_back({T.a, ab, ca});
      nt.push_back({T.b, bc, ab});
      nt.push_back({T.c, ca, bc});
      nt.push_back({ab, bc, ca});
    }
    tris = nt;
  }
  for (auto& T : tris) {
    T.a = T.a * R;
    T.b = T.b * R;
    T.c = T.c * R;
  }
  return tris;
}

static std::vector<Triangle> cube(Vec3 lo, Vec3 hi) {
  Vec3 p[8];
  for (int c = 0; c < 8; ++c) p[c] = {(c & 1) ? hi.x : lo.x, (c & 2) ? hi.y : lo.y, (c & 4) ? hi.z : lo.z};
  int q[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
  std::vector<Triangle> t;
  for (auto& f : q) {
    t.push_back({p[f[0]], p[f[1]], p[f[2]]});
    t.push_back({p[f[0]], p[f[2]], p[f[3]]});
  }
  return t;
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  log0("test_geometry\n");
  Rng rng(7, 0);

  // 1. primitives: sign and distance at known points
  {
    PartConfig b;
    b.shape = ShapeKind::Box;
    b.size = {2, 4, 6};
    Part box(b);
    CHECK(std::fabs(box.sdf({2, 0, 0}, 1e300) - 1.0) < 1e-12, "box face distance");
    CHECK(std::fabs(box.sdf({0, 0, 0}, 1e300) + 1.0) < 1e-12, "box interior distance");
    PartConfig c;
    c.shape = ShapeKind::Cylinder;
    c.radius = 1;
    c.length = 4;
    c.axis = {0, 0, 1};
    Part cyl(c);
    CHECK(std::fabs(cyl.sdf({3, 0, 0}, 1e300) - 2.0) < 1e-12, "cylinder radial distance");
    CHECK(std::fabs(cyl.sdf({0, 0, 5}, 1e300) - 3.0) < 1e-12, "cylinder cap distance");
    PartConfig t;
    t.shape = ShapeKind::Torus;
    t.radius = 2;
    t.radius2 = 0.5;
    t.axis = {0, 0, 1};
    Part tor(t);
    CHECK(std::fabs(tor.sdf({2, 0, 0}, 1e300) + 0.5) < 1e-12, "torus tube centre");
    PartConfig k;
    k.shape = ShapeKind::Cone;
    k.radius = 1;
    k.radius2 = 0;
    k.length = 2;
    k.axis = {0, 0, 1};
    Part cone(k);
    CHECK(cone.sdf({0, 0, 0}, 1e300) < 0 && cone.sdf({0, 0, 1.5}, 1e300) > 0, "cone sign");
    log0("[primitives] box/cylinder/torus/cone distances ok\n");
  }
  // 2. STL icosphere vs analytic sphere
  {
    const double R = 0.3;
    write_stl("test_icosphere.stl", icosphere(R, 4));
    StlMesh m;
    CHECK(m.load("test_icosphere.stl", 1.0), "load binary STL");
    CHECK(m.is_closed(), "icosphere is watertight");
    int wrong = 0;
    double maxerr = 0;
    for (int n = 0; n < 20000; ++n) {
      Vec3 p{rng.uniform() - 0.5, rng.uniform() - 0.5, rng.uniform() - 0.5};
      double exact = norm(p) - R;
      double s = m.sdf(p);
      if (std::fabs(exact) > 0.01 && ((s < 0) != (exact < 0))) ++wrong;
      maxerr = std::max(maxerr, std::fabs(s - exact));
    }
    log0("[STL      ] icosphere (%zu tris): %d sign errors / 20000, max |sdf err| = %.2e (facet sag %.2e)\n", m.size(),
         wrong, maxerr, R * (1 - std::cos(0.5 * 0.0629 * 2)));
    CHECK(wrong == 0, "STL inside/outside classification");
    CHECK(maxerr < 2e-3, "STL distance accuracy");
  }
  // 3. overlapping closed shells (ray parity would fail in the overlap)
  {
    auto t1 = cube({-1, -1, -1}, {0.5, 1, 1});
    auto t2 = cube({-0.5, -0.5, -0.5}, {1.0, 0.5, 0.5});
    t1.insert(t1.end(), t2.begin(), t2.end());
    write_stl("test_overlap.stl", t1);
    StlMesh m;
    m.load("test_overlap.stl", 1.0);
    bool in_overlap = m.inside({0, 0, 0}), in_a = m.inside({-0.8, 0.8, 0}), in_b = m.inside({0.8, 0, 0});
    bool out = !m.inside({0.8, 0.8, 0});
    log0("[STL      ] overlapping shells: winding(overlap) = %.3f, inside tests %d%d%d%d, parity in overlap = %d\n",
         m.winding_number({0, 0, 0}), in_overlap, in_a, in_b, out, int(m.inside_parity({0, 0, 0})));
    CHECK(in_overlap && in_a && in_b && out, "winding-number inside test with overlapping shells");
  }
  // 4. CSG + thin plate segment detection
  {
    ObjectConfig o;
    o.name = "csg";
    PartConfig b;
    b.shape = ShapeKind::Box;
    b.size = {1, 1, 1};
    PartConfig hole;
    hole.shape = ShapeKind::Cylinder;
    hole.radius = 0.2;
    hole.length = 2;
    hole.axis = {0, 0, 1};
    hole.op = CsgOp::Subtract;
    PartConfig plate;
    plate.shape = ShapeKind::Plate;
    plate.size = {2, 2, 0.002};
    plate.center = {0, 0, 2};
    o.parts = {b, hole, plate};
    Scene sc;
    sc.build({o});
    CHECK(sc.sdf({0, 0, 0}) > 0, "CSG subtract hole is empty");
    CHECK(sc.sdf({0.4, 0, 0}) < 0, "CSG body solid");
    double t;
    int w;
    bool hit = sc.segment_hit({0.3, 0.3, 1.9}, {0.3, 0.3, 2.1}, 1e-6, &t, &w);
    log0("[CSG/thin ] hole empty, 2 mm plate hit by a 20 cm step: %d (t = %.4f, expected 0.4950)\n", int(hit), t);
    CHECK(hit && std::fabs(t - 0.495) < 1e-3, "thin plate segment hit");
  }
  // 5. surface triangulation of a sphere object: total area and outward normals
  {
    Config c;
    c.cells = {24, 24, 24};
    c.dx = {1.0 / 24, 1.0 / 24, 1.0 / 24};
    c.origin = {-0.5, -0.5, -0.5};
    for (auto& f : c.faces) {
      f.field = FieldBC::Dirichlet;
      f.particle = ParticleBC::Absorb;
    }
    SpeciesConfig e;
    e.name = "electron";
    e.charge = -constants::qe;
    e.mass = constants::me;
    e.weight = 1;
    e.load_initial = false;
    c.species = {e};
    c.dt = 1e-9;
    c.output.dir = "test_output";
    c.output.history_every = 0;
    c.output.status_every = 0;
    ObjectConfig o;
    o.name = "ball";
    o.bias = BiasMode::Fixed;
    PartConfig s;
    s.shape = ShapeKind::Sphere;
    s.radius = 0.2;
    o.parts = {s};
    c.objects = {o};
    Simulation sim(c, MPI_COMM_WORLD);
    sim.initialise();
    double A = 0, nout = 0;
    for (auto& p : sim.objects().panels()) {
      A += p.area;
      nout += dot(p.normal, normalized(p.centroid)) * p.area;
    }
    A = allreduce_sum(A, MPI_COMM_WORLD);
    nout = allreduce_sum(nout, MPI_COMM_WORLD);
    double Aex = 4 * pi * 0.04;
    log0("[surface  ] panel area %.5f vs 4 pi R^2 = %.5f (%+.2f%%), mean n.r_hat = %.4f\n", A, Aex, 100 * (A / Aex - 1), nout / A);
    CHECK(std::fabs(A / Aex - 1) < 0.03, "surface area");
    CHECK(nout / A > 0.99, "outward normals");
  }
  log0(failures ? "test_geometry: %d FAILURE(S)\n" : "test_geometry: all checks passed\n", failures);
  MPI_Finalize();
  return failures ? 1 : 0;
}
