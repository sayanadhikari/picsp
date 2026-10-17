#include "fields/EmbeddedBoundary.hpp"

#include <algorithm>
#include <cmath>

#include "util/Log.hpp"

namespace picsp {

namespace {
constexpr double kThetaMin = 0.1;  // below this a free node is snapped onto the surface

inline int wrap(int i, int N) { return ((i % N) + N) % N; }
}  // namespace

void EmbeddedBoundary::build(const Config& cfg, const GridGeom& g, int G, const Scene& scene) {
  has_objects_ = !scene.empty();
  sdf_.init(g, G);
  for (int a = 0; a < 3; ++a) {
    ncx_ = g.c1[0] - g.c0[0];
    ncy_ = g.c1[1] - g.c0[1];
    ncz_ = g.c1[2] - g.c0[2];
  }
  const auto& n = sdf_.size();
  node_obj_.assign(sdf_.total(), -1);
  cut_.clear();
  near_.clear();
  n_obj_nodes_ = 0;
  if (!has_objects_) {
    sdf_.fill(1e300);
    return;
  }

  // warn about objects crossing periodic boundaries and thin features
  const double hmin = std::min(g.dx[0], std::min(g.dx[1], g.dx[2]));
  const double hmax = std::max(g.dx[0], std::max(g.dx[1], g.dx[2]));
  for (std::size_t o = 0; o < scene.size(); ++o) {
    const auto& ob = scene.object(int(o));
    for (int a = 0; a < 3; ++a) {
      double lo = g.origin[a], hi = g.origin[a] + g.N[a] * g.dx[a];
      if (g.periodic[a] && (ob.bounds().lo[a] < lo || ob.bounds().hi[a] > hi))
        warn0("object '%s' crosses a periodic boundary; its periodic image is not represented\n", ob.name().c_str());
    }
    if (ob.min_feature() < 2 * hmin)
      log0(" note: object '%s' has features thinner than 2 cells (%.3g m); thin parts are\n"
           "       resolved as one-node-thick conductor sheets in the field solve\n",
           ob.name().c_str(), ob.min_feature());
  }

  // ---- padded working box (global node indices) ----
  const int M = 2;
  int lo[3], sz[3];
  for (int a = 0; a < 3; ++a) {
    lo[a] = g.c0[a] - G - M;
    sz[a] = n[a] + 2 * M;
  }
  const std::size_t total = std::size_t(sz[0]) * sz[1] * sz[2];
  auto pidx = [&](int i, int j, int k) { return std::size_t(i) + std::size_t(sz[0]) * (std::size_t(j) + std::size_t(sz[1]) * k); };
  auto coord = [&](int i, int j, int k) {
    int gi[3] = {lo[0] + i, lo[1] + j, lo[2] + k};
    Vec3 p;
    for (int a = 0; a < 3; ++a) {
      int q = g.periodic[a] ? wrap(gi[a], g.N[a]) : gi[a];
      p[a] = g.origin[a] + q * g.dx[a];
    }
    return p;
  };
  std::vector<double> ps(total);
  std::vector<std::int16_t> po(total, -1);
  const double cap = 3.0 * hmax;
#pragma omp parallel for collapse(2) schedule(dynamic, 4)
  for (int k = 0; k < sz[2]; ++k)
    for (int j = 0; j < sz[1]; ++j)
      for (int i = 0; i < sz[0]; ++i) {
        int w;
        double s = scene.sdf(coord(i, j, k), &w, cap);
        ps[pidx(i, j, k)] = s;
        if (s < 0) po[pidx(i, j, k)] = std::int16_t(w);
      }

  // ---- thin-feature pass: an edge between two outside nodes that crosses a surface
  //      marks the node nearer to the surface as an object node ----
  const double tol = 1e-3 * hmin;
  {
    std::vector<std::int16_t> add(total, -1);
#pragma omp parallel for collapse(2) schedule(dynamic, 4)
    for (int k = 0; k < sz[2]; ++k)
      for (int j = 0; j < sz[1]; ++j)
        for (int i = 0; i < sz[0]; ++i) {
          std::size_t p = pidx(i, j, k);
          if (po[p] >= 0 || ps[p] > hmax) continue;
          int ii[3] = {i, j, k};
          for (int a = 0; a < 3; ++a) {
            if (ii[a] + 1 >= sz[a]) continue;
            int jj[3] = {i, j, k};
            jj[a] += 1;
            std::size_t q = pidx(jj[0], jj[1], jj[2]);
            if (po[q] >= 0 || ps[q] > g.dx[a]) continue;
            double t;
            int w;
            if (scene.segment_hit(coord(i, j, k), coord(jj[0], jj[1], jj[2]), tol, &t, &w)) {
              std::size_t mark = ps[p] <= ps[q] ? p : q;
#pragma omp critical(eb_thin)
              add[mark] = std::int16_t(w);
            }
          }
        }
    for (std::size_t p = 0; p < total; ++p)
      if (add[p] >= 0 && po[p] < 0) po[p] = add[p];
  }

  // fraction along edge i->j (unit length) where object `o` surface is crossed
  auto edge_theta = [&](const Vec3& pi, const Vec3& pj, int o) -> double {
    const SolidObject& ob = scene.object(o);
    double si = ob.sdf(pi), sj = ob.sdf(pj);
    if (si <= 0) return 0.0;
    if (sj < 0) {
      double a = 0, b = 1;
      for (int it = 0; it < 40; ++it) {
        double m = 0.5 * (a + b);
        if (ob.sdf(pi + (pj - pi) * m) > 0) a = m;
        else b = m;
      }
      return 0.5 * (a + b);
    }
    double t;
    int w;
    if (scene.segment_hit(pi, pj, tol, &t, &w)) return t;
    return 1.0;
  };

  static const int dirs[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
  const bool cut_cell = cfg.solver.cut_cell;

  // ---- snapping pass: free nodes extremely close to a surface become object nodes ----
  if (cut_cell) {
    for (int pass = 0; pass < 2; ++pass) {
      std::vector<std::int16_t> snap(total, -1);
#pragma omp parallel for collapse(2) schedule(dynamic, 4)
      for (int k = 1; k < sz[2] - 1; ++k)
        for (int j = 1; j < sz[1] - 1; ++j)
          for (int i = 1; i < sz[0] - 1; ++i) {
            std::size_t p = pidx(i, j, k);
            if (po[p] >= 0 || ps[p] > hmax) continue;
            for (auto& d : dirs) {
              std::size_t q = pidx(i + d[0], j + d[1], k + d[2]);
              if (po[q] < 0) continue;
              double th = edge_theta(coord(i, j, k), coord(i + d[0], j + d[1], k + d[2]), po[q]);
              if (th < kThetaMin) {
                snap[p] = po[q];
                break;
              }
            }
          }
      for (std::size_t p = 0; p < total; ++p)
        if (snap[p] >= 0) po[p] = snap[p];
    }
  }

  // ---- copy into local field layout ----
  for (int k = 0; k < n[2]; ++k)
    for (int j = 0; j < n[1]; ++j)
      for (int i = 0; i < n[0]; ++i) {
        std::size_t p = pidx(i + M, j + M, k + M);
        sdf_(i, j, k) = ps[p];
        node_obj_[sdf_.idx(i, j, k)] = po[p];
      }

  // ---- cut-cell coefficients on owned free nodes ----
  int on[3] = {g.n1[0] - g.n0[0], g.n1[1] - g.n0[1], g.n1[2] - g.n0[2]};
  for (int k = G; k < G + on[2]; ++k)
    for (int j = G; j < G + on[1]; ++j)
      for (int i = G; i < G + on[0]; ++i) {
        std::size_t p = sdf_.idx(i, j, k);
        if (node_obj_[p] >= 0) {
          ++n_obj_nodes_;
          continue;
        }
        bool adj = false;
        for (auto& d : dirs)
          if (node_obj_[sdf_.idx(i + d[0], j + d[1], k + d[2])] >= 0) adj = true;
        if (!adj) continue;
        CutNode cn;
        cn.idx = p;
        for (int dd = 0; dd < 6; ++dd) {
          auto& d = dirs[dd];
          std::size_t q = sdf_.idx(i + d[0], j + d[1], k + d[2]);
          double th = 1.0;
          if (cut_cell && node_obj_[q] >= 0) {
            th = edge_theta(coord(i + M, j + M, k + M), coord(i + d[0] + M, j + d[1] + M, k + d[2] + M), node_obj_[q]);
            th = std::clamp(th, kThetaMin, 1.0);
          }
          cn.theta[dd] = th;
        }
        for (int a = 0; a < 3; ++a) {
          double tm = cn.theta[2 * a], tp = cn.theta[2 * a + 1], h2 = g.dx[a] * g.dx[a];
          cn.c[2 * a] = 2.0 / (h2 * tm * (tm + tp));
          cn.c[2 * a + 1] = 2.0 / (h2 * tp * (tm + tp));
        }
        cut_.push_back(cn);
      }

  // ---- near-surface cell flags (local cells plus one ghost layer) ----
  near_.assign(std::size_t(ncx_ + 2) * (ncy_ + 2) * (ncz_ + 2), 0);
  const double half_diag = 0.5 * std::sqrt(g.dx[0] * g.dx[0] + g.dx[1] * g.dx[1] + g.dx[2] * g.dx[2]);
#pragma omp parallel for collapse(2) schedule(dynamic, 4)
  for (int ck = -1; ck <= ncz_; ++ck)
    for (int cj = -1; cj <= ncy_; ++cj)
      for (int ci = -1; ci <= ncx_; ++ci) {
        int gc[3] = {g.c0[0] + ci, g.c0[1] + cj, g.c0[2] + ck};
        Vec3 c;
        for (int a = 0; a < 3; ++a) {
          int q = g.periodic[a] ? wrap(gc[a], g.N[a]) : gc[a];
          c[a] = g.origin[a] + (q + 0.5) * g.dx[a];
        }
        double s = scene.sdf(c, nullptr, cap);
        near_[std::size_t(ci + 1) + std::size_t(ncx_ + 2) * (std::size_t(cj + 1) + std::size_t(ncy_ + 2) * (ck + 1))] =
            s < 1.05 * half_diag + tol;
      }
}

}  // namespace picsp
