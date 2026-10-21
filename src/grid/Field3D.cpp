#include "grid/Field3D.hpp"

#include <algorithm>
#include <cmath>

#include "util/Log.hpp"

namespace picsp {

GridGeom GridGeom::from(const Decomposition& d) {
  GridGeom g;
  g.comm = d.comm;
  g.rank = d.rank;
  g.nbr = d.nbr;
  g.periodic = d.periodic;
  for (int a = 0; a < 3; ++a) {
    g.wlo[a] = d.wall_lo(a);
    g.whi[a] = d.wall_hi(a);
  }
  g.N = d.N;
  g.c0 = d.c0;
  g.c1 = d.c1;
  g.n0 = d.n0;
  g.n1 = d.n1;
  g.dx = d.dx;
  g.origin = d.origin;
  return g;
}

bool GridGeom::can_coarsen() const {
  int ok = 1;
  for (int a = 0; a < 3; ++a) {
    if (c0[a] % 2 || c1[a] % 2 || N[a] % 2) ok = 0;
    if ((c1[a] - c0[a]) / 2 < 2) ok = 0;
    if (N[a] / 2 < 2) ok = 0;
  }
  int all = 0;
  MPI_Allreduce(&ok, &all, 1, MPI_INT, MPI_MIN, comm);
  return all == 1;
}

GridGeom GridGeom::coarsen() const {
  GridGeom g = *this;
  for (int a = 0; a < 3; ++a) {
    g.N[a] = N[a] / 2;
    g.c0[a] = c0[a] / 2;
    g.c1[a] = c1[a] / 2;
    g.n0[a] = g.c0[a];
    g.n1[a] = g.c1[a] + (whi[a] ? 1 : 0);
    g.dx[a] = dx[a] * 2;
  }
  return g;
}

long long GridGeom::global_nodes() const {
  long long n = 1;
  for (int a = 0; a < 3; ++a) n *= (N[a] + (periodic[a] ? 0 : 1));
  return n;
}

Field3D::Field3D(const GridGeom& g, int ghosts) { init(g, ghosts); }

void Field3D::init(const GridGeom& g, int ghosts) {
  g_ = &g;
  G_ = ghosts;
  for (int a = 0; a < 3; ++a) n_[a] = (g.c1[a] - g.c0[a]) + 1 + 2 * G_;
  v_.assign(std::size_t(n_[0]) * n_[1] * n_[2], 0.0);
}

// threaded: these run several times per Poisson iteration
void Field3D::fill(double val) {
  const std::size_t N = v_.size();
  double* v = v_.data();
#pragma omp parallel for simd schedule(static)
  for (std::size_t p = 0; p < N; ++p) v[p] = val;
}

void Field3D::copy_from(const Field3D& o) {
  if (v_.size() != o.v_.size()) {
    v_ = o.v_;
    return;
  }
  const std::size_t N = v_.size();
  double* v = v_.data();
  const double* s = o.v_.data();
#pragma omp parallel for simd schedule(static)
  for (std::size_t p = 0; p < N; ++p) v[p] = s[p];
}

// Visit a slab [p0, p0+np) along axis a (full extent on the other axes) as contiguous
// x-rows: f(local index of row start, row length, running buffer offset). Rows touch disjoint
// memory, so large slabs are visited in parallel (small coarse-grid slabs stay serial).
namespace {
template <typename F>
inline void for_slab_rows(const std::array<int, 3>& n, int a, int p0, int np, F&& f) {
  int lo[3] = {0, 0, 0}, hi[3] = {n[0], n[1], n[2]};
  lo[a] = p0;
  hi[a] = p0 + np;
  const std::size_t len = std::size_t(hi[0] - lo[0]);
  const int nj = hi[1] - lo[1];
  const std::size_t total = len * std::size_t(nj) * std::size_t(hi[2] - lo[2]);
#pragma omp parallel for collapse(2) schedule(static) if (total > 16384)
  for (int k = lo[2]; k < hi[2]; ++k)
    for (int j = lo[1]; j < hi[1]; ++j) {
      const std::size_t c = (std::size_t(k - lo[2]) * nj + std::size_t(j - lo[1])) * len;
      f(std::size_t(lo[0]) + std::size_t(n[0]) * (std::size_t(j) + std::size_t(n[1]) * k), len, c);
    }
}
}  // namespace

void Field3D::exchange_axis(int a, bool add) {
  const GridGeom& g = *g_;
  const int own = g.c1[a] - g.c0[a];  // owned cells along a
  const int lo_n = g.nbr[a][0], hi_n = g.nbr[a][1];
  std::size_t plane = v_.size() / n_[a];
  // Two transfers per axis:
  //  (A) towards lo neighbour, (B) towards hi neighbour.
  // copy: A sends my planes [G, 2G] (G+1 planes) -> lo nbr high ghosts [own+G, own+2G]
  //       B sends my planes [own, own+G-1] (G planes) -> hi nbr low ghosts [0, G-1]
  // add : A sends my low ghosts [0, G-1] (G planes) -> lo nbr adds to [own, own+G-1]
  //       B sends my high ghosts [own+G, own+2G] (G+1) -> hi nbr adds to [G, 2G]
  int sA0, nA, sB0, nB, rA0, rB0;  // rA0: where I receive data coming from hi nbr (its A)
  if (!add) {
    sA0 = G_; nA = G_ + 1;
    sB0 = own; nB = G_;
    rA0 = own + G_;  // data from hi neighbour's A goes into my high ghosts
    rB0 = 0;         // data from lo neighbour's B goes into my low ghosts
  } else {
    sA0 = 0; nA = G_;
    sB0 = own + G_; nB = G_ + 1;
    rA0 = own;  // hi neighbour's low ghosts add into my top owned planes
    rB0 = G_;   // lo neighbour's high ghosts add into my bottom owned planes
  }
  for (int t = 0; t < 2; ++t) {
    sbuf_[t].resize(plane * (G_ + 1));
    rbuf_[t].resize(plane * (G_ + 1));
  }
  if (lo_n == MPI_PROC_NULL && hi_n == MPI_PROC_NULL) return;  // single rank, walls on both sides
  if (lo_n == g.rank && hi_n == g.rank) {
    // periodic axis owned by a single rank: direct in-place copy, no MPI
    int src[2] = {sA0, sB0}, dst[2] = {rA0, rB0}, cnt[2] = {nA, nB};
    for (int t = 0; t < 2; ++t) {
      const std::ptrdiff_t off = std::ptrdiff_t(dst[t] - src[t]) * stride(a);
      double* v = v_.data();
      for_slab_rows(n_, a, src[t], cnt[t], [&](std::size_t p, std::size_t len, std::size_t) {
        if (add)
          for (std::size_t q = 0; q < len; ++q) v[p + q + off] += v[p + q];
        else
          for (std::size_t q = 0; q < len; ++q) v[p + q + off] = v[p + q];
      });
    }
    return;
  }
  // A: send to lo, receive from hi
  for_slab_rows(n_, a, sA0, nA, [&](std::size_t p, std::size_t len, std::size_t c) { std::copy_n(&v_[p], len, &sbuf_[0][c]); });
  MPI_Sendrecv(sbuf_[0].data(), int(plane * nA), MPI_DOUBLE, lo_n, 100 + a, rbuf_[0].data(), int(plane * nA),
               MPI_DOUBLE, hi_n, 100 + a, g.comm, MPI_STATUS_IGNORE);
  // B: send to hi, receive from lo
  for_slab_rows(n_, a, sB0, nB, [&](std::size_t p, std::size_t len, std::size_t c) { std::copy_n(&v_[p], len, &sbuf_[1][c]); });
  MPI_Sendrecv(sbuf_[1].data(), int(plane * nB), MPI_DOUBLE, hi_n, 200 + a, rbuf_[1].data(), int(plane * nB),
               MPI_DOUBLE, lo_n, 200 + a, g.comm, MPI_STATUS_IGNORE);
  if (hi_n != MPI_PROC_NULL) {
    for_slab_rows(n_, a, rA0, nA, [&](std::size_t p, std::size_t len, std::size_t c) {
      if (add) for (std::size_t q = 0; q < len; ++q) v_[p + q] += rbuf_[0][c + q];
      else std::copy_n(&rbuf_[0][c], len, &v_[p]);
    });
  }
  if (lo_n != MPI_PROC_NULL) {
    for_slab_rows(n_, a, rB0, nB, [&](std::size_t p, std::size_t len, std::size_t c) {
      if (add) for (std::size_t q = 0; q < len; ++q) v_[p + q] += rbuf_[1][c + q];
      else std::copy_n(&rbuf_[1][c], len, &v_[p]);
    });
  }
}

void Field3D::halo_copy() {
  for (int a = 0; a < 3; ++a) exchange_axis(a, false);
}

void Field3D::halo_add() {
  for (int a = 0; a < 3; ++a) exchange_axis(a, true);
  halo_copy();
}

template <typename F>
static inline void for_owned(const Field3D& f, F&& fn) {
  const auto& g = f.geom();
  int G = f.ghosts();
  int ni = g.n1[0] - g.n0[0], nj = g.n1[1] - g.n0[1], nk = g.n1[2] - g.n0[2];
  for (int k = G; k < G + nk; ++k)
    for (int j = G; j < G + nj; ++j) {
      std::size_t base = f.idx(G, j, k);
      for (int i = 0; i < ni; ++i) fn(base + i);
    }
}

double dot_owned(const Field3D& a, const Field3D& b) {
  double s = 0;
  const auto& g = a.geom();
  int G = a.ghosts();
  int ni = g.n1[0] - g.n0[0], nj = g.n1[1] - g.n0[1], nk = g.n1[2] - g.n0[2];
#pragma omp parallel for collapse(2) reduction(+ : s) schedule(static)
  for (int k = G; k < G + nk; ++k)
    for (int j = G; j < G + nj; ++j) {
      std::size_t base = a.idx(G, j, k);
      const double* pa = a.data() + base;
      const double* pb = b.data() + base;
      double t = 0;
#pragma omp simd reduction(+ : t)
      for (int i = 0; i < ni; ++i) t += pa[i] * pb[i];
      s += t;
    }
  return allreduce_sum(s, g.comm);
}

double sum_owned(const Field3D& a) {
  double s = 0;
  for_owned(a, [&](std::size_t p) { s += a[p]; });
  return allreduce_sum(s, a.geom().comm);
}

double max_abs_owned(const Field3D& a) {
  double s = 0;
  for_owned(a, [&](std::size_t p) { s = std::max(s, std::fabs(a[p])); });
  return allreduce_max(s, a.geom().comm);
}

}  // namespace picsp
