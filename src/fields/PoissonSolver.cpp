#include "fields/PoissonSolver.hpp"

#include <algorithm>
#include <cmath>

#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

namespace {
inline bool neumann_lo(const std::array<FaceConfig, 6>& f, int a) { return f[2 * a].field == FieldBC::Neumann; }
inline bool neumann_hi(const std::array<FaceConfig, 6>& f, int a) { return f[2 * a + 1].field == FieldBC::Neumann; }
inline bool dirichlet_lo(const std::array<FaceConfig, 6>& f, int a) { return f[2 * a].field == FieldBC::Dirichlet; }
inline bool dirichlet_hi(const std::array<FaceConfig, 6>& f, int a) { return f[2 * a + 1].field == FieldBC::Dirichlet; }
}  // namespace

void PoissonSolver::setup(const Config& cfg, const GridGeom& g0, const EmbeddedBoundary& eb, int phi_ghosts) {
  sc_ = cfg.solver;
  faces_ = cfg.faces;
  phiG_ = phi_ghosts;
  L_.clear();
  GridGeom g = g0;
  int maxlev = std::max(1, sc_.mg_max_levels);
  std::string t = sc_.type;
  if (t == "jacobi-pcg" || t == "cg") precond_ = 0;
  for (int l = 0; l < maxlev; ++l) {
    L_.emplace_back();
    if (!(precond_ && g.can_coarsen()) || l == maxlev - 1) break;
    g = g.coarsen();
  }
  // geometry for each level (stored by value; Field3D keeps a pointer to it, so the
  // vector must not reallocate from here on)
  g = g0;
  for (std::size_t l = 0; l < L_.size(); ++l) {
    if (l > 0) g = g.coarsen();
    L_[l].geom = g;
  }
  for (std::size_t l = 0; l < L_.size(); ++l) {
    Level& L = L_[l];
    const GridGeom& gg = L.geom;
    L.u.init(gg, 1);
    L.f.init(gg, 1);
    L.r.init(gg, 1);
    const auto& n = L.u.size();
    for (int a = 0; a < 3; ++a) {
      L.c[a] = 1.0 / (gg.dx[a] * gg.dx[a]);
      L.own[a] = gg.n1[a] - gg.n0[a];
      L.wfac[a].assign(n[a], 1.0);
      for (int i = 0; i < n[a]; ++i) {
        int gi = gg.c0[a] + i - 1;
        if ((gg.wlo[a] && gi == 0 && neumann_lo(faces_, a)) || (gg.whi[a] && gi == gg.N[a] && neumann_hi(faces_, a)))
          L.wfac[a][i] = 0.5;
      }
    }
    L.iwfac0.resize(n[0]);
    for (int i = 0; i < n[0]; ++i) L.iwfac0[i] = 1.0 / L.wfac[0][i];
    L.fixed.assign(L.u.total(), 0);
    for (int k = 0; k < n[2]; ++k)
      for (int j = 0; j < n[1]; ++j)
        for (int i = 0; i < n[0]; ++i) {
          int gi[3] = {gg.c0[0] + i - 1, gg.c0[1] + j - 1, gg.c0[2] + k - 1};
          bool fx = false;
          for (int a = 0; a < 3; ++a) {
            if (gg.wlo[a] && gi[a] <= 0 && dirichlet_lo(faces_, a)) fx = true;
            if (gg.whi[a] && gi[a] >= gg.N[a] && dirichlet_hi(faces_, a)) fx = true;
          }
          if (!fx && eb.has_objects()) {
            // object node? (fine node at the same position)
            const Field3D& s = eb.sdf();
            int f = 1 << l;
            // global fine index of this node
            long long fgi[3] = {(long long)gi[0] * f, (long long)gi[1] * f, (long long)gi[2] * f};
            bool inside_fine_array = true;
            int fl[3];
            for (int a = 0; a < 3; ++a) {
              fl[a] = int(fgi[a] - g0.c0[a] + s.ghosts());
              if (fl[a] < 0 || fl[a] >= s.size()[a]) inside_fine_array = false;
            }
            if (inside_fine_array && eb.obj(s.idx(fl[0], fl[1], fl[2])) >= 0) fx = true;
          }
          L.fixed[L.u.idx(i, j, k)] = fx;
        }
  }
  // cut-cell coefficients on level 0 (remap indices from EB layout to G=1 layout)
  Level& L0 = L_[0];
  L0.cutid.assign(L0.u.total(), -1);
  L0.cut.clear();
  has_cut_ = false;
  if (eb.has_objects()) {
    const Field3D& s = eb.sdf();
    int Ge = s.ghosts();
    const auto& ne = s.size();
    for (const CutNode& cn : eb.cut_nodes()) {
      std::size_t p = cn.idx;
      int i = int(p % ne[0]), j = int((p / ne[0]) % ne[1]), k = int(p / (std::size_t(ne[0]) * ne[1]));
      CutNode c2 = cn;
      c2.idx = L0.u.idx(i - Ge + 1, j - Ge + 1, k - Ge + 1);
      if (L0.fixed[c2.idx]) continue;  // e.g. on a Dirichlet wall
      bool irregular = false;
      for (int d = 0; d < 6; ++d)
        if (cn.theta[d] < 1.0) irregular = true;
      if (!irregular) continue;
      L0.cutid[c2.idx] = int(L0.cut.size());
      L0.cut.push_back(c2);
      has_cut_ = true;
    }
  }
  int any_cut = has_cut_ ? 1 : 0, all_cut = 0;
  MPI_Allreduce(&any_cut, &all_cut, 1, MPI_INT, MPI_MAX, g0.comm);
  has_cut_ = all_cut;
  use_bicgstab_ = (t == "mg-bicgstab" || t == "bicgstab") || (has_cut_ && t == "auto");
  if (t == "mg" || t == "multigrid") use_bicgstab_ = false;

  // singular if no fixed node anywhere
  long long nfix = 0;
  for (auto v : L0.fixed) nfix += v;
  nfix = allreduce_sum(nfix, g0.comm);
  singular_ = (nfix == 0);

  for (Field3D* f : {&b_, &x_, &p_, &q_, &z_, &s_, &t_, &rhat_, &v_, &w_}) f->init(L0.geom, 1);
}

bool PoissonSolver::is_fixed_phi_layout(int i, int j, int k) const {
  const Level& L = L_[0];
  return L.fixed[L.u.idx(i - phiG_ + 1, j - phiG_ + 1, k - phiG_ + 1)];
}

double PoissonSolver::volume_weight_phi_layout(int i, int j, int k) const {
  const Level& L = L_[0];
  return L.wfac[0][i - phiG_ + 1] * L.wfac[1][j - phiG_ + 1] * L.wfac[2][k - phiG_ + 1];
}

double PoissonSolver::link_coef_phi_layout(int i, int j, int k, int dir) const {
  const Level& L = L_[0];
  std::size_t p = L.u.idx(i - phiG_ + 1, j - phiG_ + 1, k - phiG_ + 1);
  int cid = L.cutid[p];
  if (cid >= 0) return L.cut[cid].c[dir];
  return L.c[dir / 2];
}

void PoissonSolver::fill_ghosts(Level& L, Field3D& x) {
  x.halo_copy();
  const GridGeom& g = L.geom;
  const auto& n = x.size();
  for (int a = 0; a < 3; ++a) {
    int planes[2][2] = {{-1, -1}, {-1, -1}};  // {dst, src}
    if (g.wlo[a] && neumann_lo(faces_, a)) { planes[0][0] = 0; planes[0][1] = 2; }
    if (g.whi[a] && neumann_hi(faces_, a)) { planes[1][0] = L.own[a] + 1; planes[1][1] = L.own[a] - 1; }
    for (auto& pl : planes) {
      if (pl[0] < 0) continue;
      int lo[3] = {0, 0, 0}, hi[3] = {n[0], n[1], n[2]};
      lo[a] = 0;
      hi[a] = 1;
      for (int k = lo[2]; k < hi[2]; ++k)
        for (int j = lo[1]; j < hi[1]; ++j)
          for (int i = lo[0]; i < hi[0]; ++i) {
            int d[3] = {i, j, k}, s[3] = {i, j, k};
            d[a] = pl[0];
            s[a] = pl[1];
            x(d[0], d[1], d[2]) = x(s[0], s[1], s[2]);
          }
    }
  }
}

void PoissonSolver::zero_fixed(Level& L, Field3D& x) {
  const std::size_t N = x.total();
  double* X = x.data();
  const std::uint8_t* F = L.fixed.data();
#pragma omp parallel for simd schedule(static)
  for (std::size_t p = 0; p < N; ++p)
    if (F[p]) X[p] = 0.0;
}

void PoissonSolver::apply(Level& L, Field3D& x, Field3D& y) {
  fill_ghosts(L, x);
  const auto& n = x.size();
  const std::ptrdiff_t sy = n[0], sz = std::ptrdiff_t(n[0]) * n[1];
  const double cx = L.c[0], cy = L.c[1], cz = L.c[2], diag = 2 * (cx + cy + cz);
  const double* X = x.data();
  double* Y = y.data();
  const std::uint8_t* F = L.fixed.data();
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = 1; k <= L.own[2]; ++k)
    for (int j = 1; j <= L.own[1]; ++j) {
      const double wjk = L.wfac[1][j] * L.wfac[2][k];
      std::size_t base = x.idx(0, j, k);
#pragma omp simd
      for (int i = 1; i <= L.own[0]; ++i) {
        std::size_t p = base + i;
        double v = diag * X[p] - cx * (X[p - 1] + X[p + 1]) - cy * (X[p - sy] + X[p + sy]) - cz * (X[p - sz] + X[p + sz]);
        Y[p] = F[p] ? 0.0 : wjk * L.wfac[0][i] * v;
      }
    }
  if (!L.cut.empty()) {
    const auto& nn = n;
    for (const CutNode& cn : L.cut) {
      std::size_t p = cn.idx;
      int i = int(p % nn[0]), j = int((p / nn[0]) % nn[1]), k = int(p / (std::size_t(nn[0]) * nn[1]));
      double xp = X[p];
      double v = cn.c[0] * (xp - X[p - 1]) + cn.c[1] * (xp - X[p + 1]) + cn.c[2] * (xp - X[p - sy]) +
                 cn.c[3] * (xp - X[p + sy]) + cn.c[4] * (xp - X[p - sz]) + cn.c[5] * (xp - X[p + sz]);
      Y[p] = L.wfac[0][i] * L.wfac[1][j] * L.wfac[2][k] * v;
    }
  }
}

void PoissonSolver::residual(Level& L, Field3D& u, const Field3D& f, Field3D& r) {
  apply(L, u, r);
  const std::size_t N = r.total();
  double* R = r.data();
  const double* Fv = f.data();
#pragma omp parallel for simd schedule(static)
  for (std::size_t p = 0; p < N; ++p) R[p] = L.fixed[p] ? 0.0 : Fv[p] - R[p];
}

void PoissonSolver::smooth(Level& L, Field3D& u, const Field3D& f, int sweeps, bool reverse) {
  const auto& n = u.size();
  const std::ptrdiff_t sy = n[0], sz = std::ptrdiff_t(n[0]) * n[1];
  const double cx = L.c[0], cy = L.c[1], cz = L.c[2], inv_diag = 1.0 / (2 * (cx + cy + cz));
  const GridGeom& g = L.geom;
  double* U = u.data();
  const double* Fv = f.data();
  const std::uint8_t* Fx = L.fixed.data();
  const bool has_cut = !L.cut.empty();
  for (int s = 0; s < sweeps; ++s)
    for (int cc = 0; cc < 2; ++cc) {
      int color = reverse ? 1 - cc : cc;
      fill_ghosts(L, u);
#pragma omp parallel for collapse(2) schedule(static)
      for (int k = 1; k <= L.own[2]; ++k)
        for (int j = 1; j <= L.own[1]; ++j) {
          // global parity of node i is (c0x - 1 + i) + gy + gz; pick i0 with parity == color
          const int i0 = 1 + (((color - g.c0[0] - (g.c0[1] + j - 1) - (g.c0[2] + k - 1)) % 2 + 2) % 2);
          const double iwjk = 1.0 / (L.wfac[1][j] * L.wfac[2][k]);
          const double* iwx = L.iwfac0.data();
          std::size_t base = u.idx(0, j, k);
          for (int i = i0; i <= L.own[0]; i += 2) {
            std::size_t p = base + i;
            if (Fx[p]) continue;
            double rhs = Fv[p] * iwjk * iwx[i];
            if (has_cut && L.cutid[p] >= 0) {
              const CutNode& cn = L.cut[L.cutid[p]];
              double num = rhs + cn.c[0] * U[p - 1] + cn.c[1] * U[p + 1] + cn.c[2] * U[p - sy] + cn.c[3] * U[p + sy] +
                           cn.c[4] * U[p - sz] + cn.c[5] * U[p + sz];
              U[p] = num / (cn.c[0] + cn.c[1] + cn.c[2] + cn.c[3] + cn.c[4] + cn.c[5]);
            } else {
              U[p] = (rhs + cx * (U[p - 1] + U[p + 1]) + cy * (U[p - sy] + U[p + sy]) + cz * (U[p - sz] + U[p + sz])) *
                     inv_diag;
            }
          }
        }
    }
}

void PoissonSolver::restrict_to(Level& F, const Field3D& rf_in, Level& C, Field3D& fc) {
  Field3D& rf = const_cast<Field3D&>(rf_in);
  rf.halo_copy();
  const GridGeom& gf = F.geom;
  const GridGeom& gc = C.geom;
  // 1D tables: coarse local index -> three fine local indices with full-weighting weights
  // (weight 0 for fine nodes outside a non-periodic domain)
  std::vector<int> fi[3];
  std::vector<double> fw[3];
  for (int a = 0; a < 3; ++a) {
    fi[a].assign(3 * (C.own[a] + 2), 1);
    fw[a].assign(3 * (C.own[a] + 2), 0.0);
    for (int I = 1; I <= C.own[a]; ++I) {
      int G = gc.c0[a] + I - 1;
      for (int o = -1; o <= 1; ++o) {
        int f = 2 * G + o;
        bool inside = gf.periodic[a] || (f >= 0 && f <= gf.N[a]);
        fi[a][3 * I + o + 1] = inside ? f - gf.c0[a] + 1 : 1;
        fw[a][3 * I + o + 1] = inside ? (o == 0 ? 1.0 : 0.5) : 0.0;
      }
    }
  }
  const auto& nf = rf.size();
  const std::size_t sy = nf[0], sz = std::size_t(nf[0]) * nf[1];
  const double* R = rf.data();
#pragma omp parallel for collapse(2) schedule(static)
  for (int K = 1; K <= C.own[2]; ++K)
    for (int J = 1; J <= C.own[1]; ++J)
      for (int I = 1; I <= C.own[0]; ++I) {
        std::size_t pc = fc.idx(I, J, K);
        if (C.fixed[pc]) {
          fc[pc] = 0.0;
          continue;
        }
        double s = 0;
        for (int ok = 0; ok < 3; ++ok) {
          const double wz = fw[2][3 * K + ok];
          const std::size_t zo = sz * fi[2][3 * K + ok];
          for (int oj = 0; oj < 3; ++oj) {
            const double wzy = wz * fw[1][3 * J + oj];
            const double* row = R + zo + sy * fi[1][3 * J + oj];
            s += wzy * (fw[0][3 * I] * row[fi[0][3 * I]] + fw[0][3 * I + 1] * row[fi[0][3 * I + 1]] +
                        fw[0][3 * I + 2] * row[fi[0][3 * I + 2]]);
          }
        }
        fc[pc] = 0.125 * s;  // R = P^T / 2^d keeps the rediscretised coarse operator consistent
      }
}

void PoissonSolver::prolong_add(Level& C, Field3D& ec, Level& F, Field3D& uf) {
  ec.halo_copy();
  const GridGeom& gf = F.geom;
  const GridGeom& gc = C.geom;
  // 1D interpolation tables: fine local index -> two coarse local indices and weights
  int lo[3][2] = {}, nf[3];
  std::vector<int> ci[3][2];
  std::vector<double> cw[3][2];
  for (int a = 0; a < 3; ++a) {
    nf[a] = F.own[a];
    for (int t = 0; t < 2; ++t) {
      ci[a][t].resize(nf[a] + 2);
      cw[a][t].resize(nf[a] + 2);
    }
    for (int i = 1; i <= nf[a]; ++i) {
      int g = gf.c0[a] + i - 1;
      int c0 = (g - (g & 1)) / 2;
      ci[a][0][i] = c0 - gc.c0[a] + 1;
      ci[a][1][i] = ci[a][0][i] + ((g & 1) ? 1 : 0);
      cw[a][0][i] = (g & 1) ? 0.5 : 1.0;
      cw[a][1][i] = (g & 1) ? 0.5 : 0.0;
    }
    (void)lo;
  }
  const auto& nc = ec.size();
  const std::size_t sy = nc[0], sz = std::size_t(nc[0]) * nc[1];
  const double* E = ec.data();
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = 1; k <= F.own[2]; ++k)
    for (int j = 1; j <= F.own[1]; ++j) {
      const std::size_t z0 = sz * ci[2][0][k], z1 = sz * ci[2][1][k];
      const std::size_t y0 = sy * ci[1][0][j], y1 = sy * ci[1][1][j];
      const double wz0 = cw[2][0][k], wz1 = cw[2][1][k], wy0 = cw[1][0][j], wy1 = cw[1][1][j];
      const double w00 = wz0 * wy0, w01 = wz0 * wy1, w10 = wz1 * wy0, w11 = wz1 * wy1;
      const double* r00 = E + z0 + y0;
      const double* r01 = E + z0 + y1;
      const double* r10 = E + z1 + y0;
      const double* r11 = E + z1 + y1;
      std::size_t base = uf.idx(0, j, k);
      for (int i = 1; i <= F.own[0]; ++i) {
        std::size_t p = base + i;
        if (F.fixed[p]) continue;
        const int x0 = ci[0][0][i], x1 = ci[0][1][i];
        const double wx0 = cw[0][0][i], wx1 = cw[0][1][i];
        uf[p] += wx0 * (w00 * r00[x0] + w01 * r01[x0] + w10 * r10[x0] + w11 * r11[x0]) +
                 wx1 * (w00 * r00[x1] + w01 * r01[x1] + w10 * r10[x1] + w11 * r11[x1]);
      }
    }
}

void PoissonSolver::project_mean(Level& L, Field3D& b) {
  // remove the component along the null space (constant) for singular problems
  if (L.wsum == 0) {  // the weights are fixed: sum them once
    double sw = 0;
    for (int k = 1; k <= L.own[2]; ++k)
      for (int j = 1; j <= L.own[1]; ++j)
        for (int i = 1; i <= L.own[0]; ++i) sw += L.wfac[0][i] * L.wfac[1][j] * L.wfac[2][k];
    L.wsum = allreduce_sum(sw, L.geom.comm);
  }
  double sb = 0;
#pragma omp parallel for collapse(2) reduction(+ : sb) schedule(static)
  for (int k = 1; k <= L.own[2]; ++k)
    for (int j = 1; j <= L.own[1]; ++j)
      for (int i = 1; i <= L.own[0]; ++i) sb += b(i, j, k);
  const double m = allreduce_sum(sb, L.geom.comm) / L.wsum;
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = 1; k <= L.own[2]; ++k)
    for (int j = 1; j <= L.own[1]; ++j) {
      const double mw = m * L.wfac[1][j] * L.wfac[2][k];
      for (int i = 1; i <= L.own[0]; ++i) b(i, j, k) -= mw * L.wfac[0][i];
    }
}

void PoissonSolver::vcycle(int l) {
  Level& L = L_[l];
  L.u.fill(0.0);
  if (l == int(L_.size()) - 1) {
    if (singular_) project_mean(L, L.f);
    // sweeps sufficient for the coarsest grid (error decays ~ (1 - c/N^2) per sweep)
    int nmax = std::max(L.geom.N[0], std::max(L.geom.N[1], L.geom.N[2]));
    int sweeps = std::min(sc_.mg_coarse_sweeps, std::max(4, nmax * nmax / 4));
    smooth(L, L.u, L.f, sweeps, false);
    smooth(L, L.u, L.f, sweeps, true);
    return;
  }
  smooth(L, L.u, L.f, sc_.mg_pre, false);
  residual(L, L.u, L.f, L.r);
  restrict_to(L, L.r, L_[l + 1], L_[l + 1].f);
  vcycle(l + 1);
  prolong_add(L_[l + 1], L_[l + 1].u, L, L.u);
  smooth(L, L.u, L.f, sc_.mg_post, true);
}

void PoissonSolver::precondition(const Field3D& r, Field3D& z) {
  Level& L0 = L_[0];
  if (precond_ == 0) {  // Jacobi
    const double diag = 2 * (L0.c[0] + L0.c[1] + L0.c[2]);
    const std::size_t N = z.total();
#pragma omp parallel for simd schedule(static)
    for (std::size_t p = 0; p < N; ++p) z[p] = L0.fixed[p] ? 0.0 : r[p] / diag;
    for (const CutNode& cn : L0.cut) {
      double d = cn.c[0] + cn.c[1] + cn.c[2] + cn.c[3] + cn.c[4] + cn.c[5];
      z[cn.idx] = r[cn.idx] / d;
    }
    return;
  }
  L0.f.copy_from(r);
  vcycle(0);
  z.copy_from(L0.u);
  zero_fixed(L0, z);
}

SolveStats PoissonSolver::solve(Field3D& phi, const Field3D& rho, const Field3D& lift) {
  using constants::eps0;
  Level& L0 = L_[0];
  const auto& n = L0.u.size();
  const int d = phiG_ - 1;
  // x <- lift at fixed nodes, 0 elsewhere ; b <- W rho/eps0 - A(lift)
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = 0; k < n[2]; ++k)
    for (int j = 0; j < n[1]; ++j)
      for (int i = 0; i < n[0]; ++i) {
        std::size_t p = L0.u.idx(i, j, k);
        x_[p] = L0.fixed[p] ? lift(i + d, j + d, k + d) : 0.0;
      }
  apply(L0, x_, b_);
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = 1; k <= L0.own[2]; ++k)
    for (int j = 1; j <= L0.own[1]; ++j)
      for (int i = 1; i <= L0.own[0]; ++i) {
        std::size_t p = b_.idx(i, j, k);
        b_[p] = L0.fixed[p] ? 0.0
                            : L0.wfac[0][i] * L0.wfac[1][j] * L0.wfac[2][k] * rho(i + d, j + d, k + d) / eps0 - b_[p];
      }
  if (singular_) project_mean(L0, b_);
  // initial guess u = phi - lift (free nodes)
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = 0; k < n[2]; ++k)
    for (int j = 0; j < n[1]; ++j)
      for (int i = 0; i < n[0]; ++i) {
        std::size_t p = L0.u.idx(i, j, k);
        v_[p] = L0.fixed[p] ? 0.0 : phi(i + d, j + d, k + d);
      }
  SolveStats st;
  if (sc_.type == "mg" || sc_.type == "multigrid") st = solve_mg(v_, b_);
  else if (use_bicgstab_) st = solve_bicgstab(v_, b_);
  else {
    st = solve_pcg(v_, b_);
    if (!st.converged && has_cut_) {
      // CG on the (mildly non-symmetric) cut-cell operator stalled: continue with BiCGStab
      int it0 = st.iterations;
      st = solve_bicgstab(v_, b_);
      st.iterations += it0;
    }
  }
  if (singular_) {
    // fix the additive constant: zero mean potential
    double s = 0, c = 0;
#pragma omp parallel for collapse(2) reduction(+ : s, c) schedule(static)
    for (int k = 1; k <= L0.own[2]; ++k)
      for (int j = 1; j <= L0.own[1]; ++j)
        for (int i = 1; i <= L0.own[0]; ++i) {
          double w = L0.wfac[0][i] * L0.wfac[1][j] * L0.wfac[2][k];
          s += w * v_(i, j, k);
          c += w;
        }
    double m = allreduce_sum(s, L0.geom.comm) / allreduce_sum(c, L0.geom.comm);
    const std::size_t N = v_.total();
#pragma omp parallel for simd schedule(static)
    for (std::size_t p = 0; p < N; ++p) v_[p] -= m;
  }
  // phi = u + lift
  fill_ghosts(L0, v_);
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = 0; k < n[2]; ++k)
    for (int j = 0; j < n[1]; ++j)
      for (int i = 0; i < n[0]; ++i) {
        std::size_t p = L0.u.idx(i, j, k);
        phi(i + d, j + d, k + d) = L0.fixed[p] ? x_[p] : v_[p];
      }
  // ghost layers beyond the solver's single ghost layer
  phi.halo_copy();
  return st;
}

SolveStats PoissonSolver::solve_pcg(Field3D& u, const Field3D& b) {
  Level& L0 = L_[0];
  SolveStats st;
  double bnorm = std::sqrt(dot_owned(b, b));
  if (bnorm == 0) {
    u.fill(0.0);
    st.converged = true;
    return st;
  }
  Field3D& r = s_;
  residual(L0, u, b, r);
  if (singular_) project_mean(L0, r);
  double rn = std::sqrt(dot_owned(r, r));
  if (rn / bnorm < sc_.tol) {
    st.converged = true;
    st.residual = rn / bnorm;
    return st;
  }
  precondition(r, z_);
  p_.copy_from(z_);
  double rz = dot_owned(r, z_);
  for (int it = 1; it <= sc_.max_iter; ++it) {
    apply(L0, p_, q_);
    double pq = dot_owned(p_, q_);
    if (pq == 0) break;
    double alpha = rz / pq;
    const std::size_t N = u.total();
#pragma omp parallel for simd schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
      u[i] += alpha * p_[i];
      r[i] -= alpha * q_[i];
    }
    rn = std::sqrt(dot_owned(r, r));
    st.iterations = it;
    st.residual = rn / bnorm;
    if (st.residual < sc_.tol) {
      st.converged = true;
      break;
    }
    if (singular_) project_mean(L0, r);
    precondition(r, z_);
    double rz_new = dot_owned(r, z_);
    double beta = rz_new / rz;
    rz = rz_new;
#pragma omp parallel for simd schedule(static)
    for (std::size_t i = 0; i < N; ++i) p_[i] = z_[i] + beta * p_[i];
  }
  return st;
}

SolveStats PoissonSolver::solve_bicgstab(Field3D& u, const Field3D& b) {
  Level& L0 = L_[0];
  SolveStats st;
  double bnorm = std::sqrt(dot_owned(b, b));
  if (bnorm == 0) {
    u.fill(0.0);
    st.converged = true;
    return st;
  }
  const std::size_t N = u.total();
  Field3D& r = s_;
  residual(L0, u, b, r);
  if (singular_) project_mean(L0, r);
  rhat_.copy_from(r);
  double rho = 1, alpha = 1, omega = 1;
  Field3D& v = q_;
  Field3D& pp = p_;
  v.fill(0.0);
  pp.fill(0.0);
  Field3D& y = z_;
  Field3D& zz = x_;  // x_ holds lift; save and restore
  std::vector<double> xsave(x_.data(), x_.data() + N);
  Field3D& sv = t_;
  Field3D& tv = w_;
  double rn = std::sqrt(dot_owned(r, r));
  if (rn / bnorm < sc_.tol) {
    st.converged = true;
    st.residual = rn / bnorm;
    return st;
  }
  for (int it = 1; it <= sc_.max_iter; ++it) {
    double rho_new = dot_owned(rhat_, r);
    if (rho_new == 0) break;
    double beta = (rho_new / rho) * (alpha / omega);
    rho = rho_new;
#pragma omp parallel for simd schedule(static)
    for (std::size_t i = 0; i < N; ++i) pp[i] = r[i] + beta * (pp[i] - omega * v[i]);
    precondition(pp, y);
    apply(L0, y, v);
    double rv = dot_owned(rhat_, v);
    if (rv == 0) break;
    alpha = rho / rv;
#pragma omp parallel for simd schedule(static)
    for (std::size_t i = 0; i < N; ++i) sv[i] = r[i] - alpha * v[i];
    double sn = std::sqrt(dot_owned(sv, sv));
    if (sn / bnorm < sc_.tol) {
#pragma omp parallel for simd schedule(static)
      for (std::size_t i = 0; i < N; ++i) u[i] += alpha * y[i];
      st.iterations = it;
      st.residual = sn / bnorm;
      st.converged = true;
      break;
    }
    precondition(sv, zz);
    apply(L0, zz, tv);
    double tt = dot_owned(tv, tv);
    omega = tt > 0 ? dot_owned(tv, sv) / tt : 0.0;
#pragma omp parallel for simd schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
      u[i] += alpha * y[i] + omega * zz[i];
      r[i] = sv[i] - omega * tv[i];
    }
    if (singular_) project_mean(L0, r);
    rn = std::sqrt(dot_owned(r, r));
    st.iterations = it;
    st.residual = rn / bnorm;
    if (st.residual < sc_.tol) {
      st.converged = true;
      break;
    }
    if (omega == 0) break;
  }
  std::copy(xsave.begin(), xsave.end(), x_.data());
  return st;
}

SolveStats PoissonSolver::solve_mg(Field3D& u, const Field3D& b) {
  Level& L0 = L_[0];
  SolveStats st;
  double bnorm = std::sqrt(dot_owned(b, b));
  if (bnorm == 0) {
    u.fill(0.0);
    st.converged = true;
    return st;
  }
  Field3D& r = s_;
  for (int it = 1; it <= sc_.max_iter; ++it) {
    residual(L0, u, b, r);
    if (singular_) project_mean(L0, r);
    double rn = std::sqrt(dot_owned(r, r));
    st.residual = rn / bnorm;
    st.iterations = it - 1;
    if (st.residual < sc_.tol) {
      st.converged = true;
      return st;
    }
    precondition(r, z_);
    for (std::size_t i = 0; i < u.total(); ++i) u[i] += z_[i];
  }
  return st;
}

}  // namespace picsp
