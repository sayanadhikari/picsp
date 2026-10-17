#include "fields/FieldSolver.hpp"

#include <algorithm>
#include <cmath>
#include <map>

#include "fields/MagneticField.hpp"
#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

bool dense_solve(std::vector<double> A, std::vector<double> b, std::vector<double>& x, int n) {
  for (int c = 0; c < n; ++c) {
    int piv = c;
    for (int r = c + 1; r < n; ++r)
      if (std::fabs(A[r * n + c]) > std::fabs(A[piv * n + c])) piv = r;
    if (std::fabs(A[piv * n + c]) < 1e-300) return false;
    if (piv != c) {
      for (int k = 0; k < n; ++k) std::swap(A[c * n + k], A[piv * n + k]);
      std::swap(b[c], b[piv]);
    }
    for (int r = c + 1; r < n; ++r) {
      double f = A[r * n + c] / A[c * n + c];
      for (int k = c; k < n; ++k) A[r * n + k] -= f * A[c * n + k];
      b[r] -= f * b[c];
    }
  }
  x.assign(n, 0.0);
  for (int r = n - 1; r >= 0; --r) {
    double s = b[r];
    for (int k = r + 1; k < n; ++k) s -= A[r * n + k] * x[k];
    x[r] = s / A[r * n + r];
  }
  return true;
}

void FieldSolver::setup(const Config& cfg, const GridGeom& g, const Scene& scene) {
  cfg_ = cfg;
  g_ = &g;
  for (Field3D* f : {&phi, &rho, &Ex, &Ey, &Ez, &phi0_, &lift_}) f->init(g, kGhosts);

  // static magnetic field on the nodes (only when spatially varying)
  B_varies_ = !cfg.magnets.empty();
  if (B_varies_) {
    for (Field3D* f : {&Bx, &By, &Bz}) f->init(g, kGhosts);
    const auto& n = Bx.size();
    const double rmin = 0.5 * std::min(g.dx[0], std::min(g.dx[1], g.dx[2]));
    double bmin = 1e300, bmax = 0;
#pragma omp parallel for collapse(2) reduction(min : bmin) reduction(max : bmax) schedule(static)
    for (int k = 0; k < n[2]; ++k)
      for (int j = 0; j < n[1]; ++j)
        for (int i = 0; i < n[0]; ++i) {
          Vec3 p{g.origin.x + (g.c0[0] + i - kGhosts) * g.dx[0], g.origin.y + (g.c0[1] + j - kGhosts) * g.dx[1],
                 g.origin.z + (g.c0[2] + k - kGhosts) * g.dx[2]};
          Vec3 B = total_magnetic_field(cfg, p, rmin);
          std::size_t q = Bx.idx(i, j, k);
          Bx[q] = B.x;
          By[q] = B.y;
          Bz[q] = B.z;
          double b = norm(B);
          bmin = std::min(bmin, b);
          bmax = std::max(bmax, b);
        }
    B3.resize(3 * Bx.total());
    for (std::size_t q = 0; q < Bx.total(); ++q) {
      B3[3 * q] = Bx[q];
      B3[3 * q + 1] = By[q];
      B3[3 * q + 2] = Bz[q];
    }
    log0(" magnetic field: |B| in [%.4g, %.4g] T (spatially varying)\n", allreduce_min(bmin, g.comm),
         allreduce_max(bmax, g.comm));
  }

  eb_.build(cfg, g, kGhosts, scene);
  long long nobj = allreduce_sum(eb_.count_object_nodes(), g.comm);
  long long ncut = allreduce_sum(eb_.count_cut_nodes(), g.comm);
  if (!scene.empty()) log0(" embedded boundary: %lld object nodes, %lld cut-cell nodes\n", nobj, ncut);
  if (!scene.empty() && nobj == 0) warn0("no grid node lies inside any object (objects smaller than a cell?)\n");

  ps_.setup(cfg, g, eb_, kGhosts);
  log0(" Poisson solver: %s, %d multigrid level(s)%s%s\n",
       (cfg.solver.type == "mg" ? "multigrid" : (ps_.symmetric() ? "MG-PCG" : "MG-BiCGStab")), ps_.levels(),
       ps_.singular() ? ", singular (zero-mean)" : "", ncut > 0 && cfg.solver.cut_cell ? ", cut-cell boundaries" : "");

  // wall volume factors (half control volume on non-periodic boundary nodes)
  const auto& n = phi.size();
  for (int a = 0; a < 3; ++a) {
    volf_[a].assign(n[a], 1.0);
    for (int i = 0; i < n[a]; ++i) {
      int gi = g.c0[a] + i - kGhosts;
      if ((g.wlo[a] && gi == 0) || (g.whi[a] && gi == g.N[a])) volf_[a][i] = 0.5;
    }
  }

  // ---- conductors ----
  cond_.clear();
  obj_to_cond_.assign(cfg.objects.size(), -1);
  std::map<int, int> group_to_cond;
  for (std::size_t o = 0; o < cfg.objects.size(); ++o) {
    const ObjectConfig& oc = cfg.objects[o];
    int c;
    if (oc.group >= 0 && group_to_cond.count(oc.group)) {
      c = group_to_cond[oc.group];
      Conductor& C = cond_[c];
      C.name += "+" + oc.name;
      if (oc.bias == BiasMode::Fixed && !C.fixed) {
        C.fixed = true;
        C.V0 = oc.potential;
        C.ramp = oc.ramp;
      }
      C.Cext += oc.ext_capacitance;
      C.Iext += oc.ext_current;
    } else {
      c = int(cond_.size());
      Conductor C;
      C.name = oc.name;
      C.fixed = oc.bias == BiasMode::Fixed;
      C.V0 = oc.potential;
      C.ramp = oc.ramp;
      C.V = oc.potential;
      C.Cext = oc.ext_capacitance;
      C.Iext = oc.ext_current;
      cond_.push_back(C);
      if (oc.group >= 0) group_to_cond[oc.group] = c;
    }
    cond_[c].objects.push_back(int(o));
    obj_to_cond_[o] = c;
  }
  const int nc = int(cond_.size());

  // floating conductors need a potential reference
  bool any_dirichlet = false, any_floating = false, any_fixed = false;
  for (auto& f : cfg.faces) any_dirichlet |= f.field == FieldBC::Dirichlet;
  for (auto& c : cond_) (c.fixed ? any_fixed : any_floating) = true;
  if (any_floating && !any_dirichlet && !any_fixed)
    fatal("floating objects need a potential reference: use a Dirichlet outer boundary or a fixed-bias object");

  build_links();

  // ---- vacuum solutions psi_k and capacitance matrix ----
  psi_.clear();
  C_.assign(std::size_t(nc) * nc, 0.0);
  Field3D zero(g, kGhosts);
  for (int k = 0; k < nc; ++k) {
    psi_.emplace_back(g, kGhosts);
    lift_.fill(0.0);
    for (std::size_t p = 0; p < lift_.total(); ++p) {
      int o = eb_.obj(p);
      if (o >= 0 && obj_to_cond_[o] == k) lift_[p] = 1.0;
    }
    SolveStats st = ps_.solve(psi_[k], zero, lift_);
    if (!st.converged) warn0("vacuum solve for conductor %s: residual %.2e after %d it\n", cond_[k].name.c_str(), st.residual, st.iterations);
    for (int i = 0; i < nc; ++i) C_[i * nc + k] = allreduce_sum(conductor_charge(psi_[k], nullptr, i), g.comm);
  }
  if (nc > 0) {
    log0(" capacitance matrix [pF]:\n");
    for (int i = 0; i < nc; ++i) {
      log0("   %-12s", cond_[i].name.c_str());
      for (int k = 0; k < nc; ++k) log0(" %10.4f", C_[i * nc + k] * 1e12);
      log0("\n");
    }
  }
  charges_initialised_ = false;
  phi0_.fill(0.0);
}

void FieldSolver::build_links() {
  using constants::eps0;
  links_.clear();
  if (!eb_.has_objects()) return;
  const GridGeom& g = *g_;
  const int G = kGhosts;
  const double dV = g.dx[0] * g.dx[1] * g.dx[2];
  static const int dirs[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
  for (int k = G; k < G + g.n1[2] - g.n0[2]; ++k)
    for (int j = G; j < G + g.n1[1] - g.n0[1]; ++j)
      for (int i = G; i < G + g.n1[0] - g.n0[0]; ++i) {
        std::size_t p = phi.idx(i, j, k);
        if (ps_.is_fixed_phi_layout(i, j, k)) continue;
        double W = ps_.volume_weight_phi_layout(i, j, k);
        for (int d = 0; d < 6; ++d) {
          std::size_t q = phi.idx(i + dirs[d][0], j + dirs[d][1], k + dirs[d][2]);
          int o = eb_.obj(q);
          if (o < 0) continue;
          links_.push_back({p, q, eps0 * W * dV / (g.dx[d / 2] * g.dx[d / 2]), obj_to_cond_[o]});
        }
      }
  // cut-cell excess shares (only where the Shortley-Weller stencil is irregular)
  cut_shares_.clear();
  for (const CutNode& cn : eb_.cut_nodes()) {
    const auto& n = phi.size();
    int i = int(cn.idx % n[0]), j = int((cn.idx / n[0]) % n[1]), k = int(cn.idx / (std::size_t(n[0]) * n[1]));
    if (ps_.is_fixed_phi_layout(i, j, k)) continue;
    double W = ps_.volume_weight_phi_layout(i, j, k);
    double tot = 0;
    std::map<int, double> wk;
    for (int d = 0; d < 6; ++d) {
      std::size_t q = phi.idx(i + dirs[d][0], j + dirs[d][1], k + dirs[d][2]);
      int o = eb_.obj(q);
      if (o < 0) continue;
      double w = 1.0 / cn.theta[d];
      wk[obj_to_cond_[o]] += w;
      tot += w;
    }
    for (auto& kv : wk) cut_shares_.push_back({cn.idx, eps0 * W * dV, kv.first, kv.second / tot});
  }
}

double FieldSolver::conductor_charge(const Field3D& f, const Field3D* rho, int k) const {
  using constants::eps0;
  double q = 0;
  for (const Link& l : links_)
    if (l.cond == k) q += l.coef * (f[l.nb] - f[l.idx]);
  if (!cut_shares_.empty()) {
    const GridGeom& g = *g_;
    const double c[3] = {1 / (g.dx[0] * g.dx[0]), 1 / (g.dx[1] * g.dx[1]), 1 / (g.dx[2] * g.dx[2])};
    const std::ptrdiff_t st[3] = {1, f.stride(1), f.stride(2)};
    for (const CutShare& cs : cut_shares_) {
      if (cs.cond != k) continue;
      std::size_t p = cs.idx;
      double lap = 0;  // regular -lap(f)
      for (int a = 0; a < 3; ++a) lap += c[a] * (2 * f[p] - f[p - st[a]] - f[p + st[a]]);
      double excess = lap - (rho ? (*rho)[p] / eps0 : 0.0);
      q += cs.share * cs.wdv * excess;
    }
  }
  return q;
}

SolveStats FieldSolver::solve(double time) {
  const GridGeom& g = *g_;
  if (!cfg_.electrostatic) {  // test-particle mode: only external fields
    phi.fill(0.0);
    compute_E();
    SolveStats st;
    st.converged = true;
    return st;
  }
  const int nc = int(cond_.size());
  // lift for phi_0: Dirichlet walls at their values, conductors grounded
  lift_.fill(0.0);
  const auto& n = lift_.size();
  for (int a = 0; a < 3; ++a) {
    for (int side = 0; side < 2; ++side) {
      const FaceConfig& fc = cfg_.faces[2 * a + side];
      if (fc.field != FieldBC::Dirichlet) continue;
      if (side == 0 && !g.wlo[a]) continue;
      if (side == 1 && !g.whi[a]) continue;
      for (int k = 0; k < n[2]; ++k)
        for (int j = 0; j < n[1]; ++j)
          for (int i = 0; i < n[0]; ++i) {
            int ii[3] = {i, j, k};
            int gi = g.c0[a] + ii[a] - kGhosts;
            if ((side == 0 && gi <= 0) || (side == 1 && gi >= g.N[a])) lift_(i, j, k) = fc.value;
          }
    }
  }
  SolveStats st = ps_.solve(phi0_, rho, lift_);

  if (nc > 0) {
    std::vector<double> qind(nc, 0.0), qloc(nc, 0.0);
    for (int k = 0; k < nc; ++k) qloc[k] = conductor_charge(phi0_, &rho, k);
    MPI_Allreduce(qloc.data(), qind.data(), nc, MPI_DOUBLE, MPI_SUM, g.comm);
    for (auto& c : cond_)
      if (c.fixed) c.V = c.V0 + c.ramp * time;
    std::vector<int> fl;
    for (int k = 0; k < nc; ++k)
      if (!cond_[k].fixed) fl.push_back(k);
    const int nf = int(fl.size());
    if (!charges_initialised_) {
      // choose initial charges so that floating conductors start at their configured potential
      for (int a : fl) {
        double q = qind[a] + cond_[a].Cext * cond_[a].V0;
        for (int k = 0; k < nc; ++k) q += C_[a * nc + k] * (cond_[k].fixed ? cond_[k].V : cond_[k].V0);
        cond_[a].Q = q;
      }
      charges_initialised_ = true;
    }
    if (nf > 0) {
      std::vector<double> A(std::size_t(nf) * nf), b(nf), x;
      for (int r = 0; r < nf; ++r) {
        int a = fl[r];
        b[r] = cond_[a].Q - qind[a];
        for (int k = 0; k < nc; ++k)
          if (cond_[k].fixed) b[r] -= C_[a * nc + k] * cond_[k].V;
        for (int c2 = 0; c2 < nf; ++c2) A[r * nf + c2] = C_[a * nc + fl[c2]];
        A[r * nf + r] += cond_[a].Cext;
      }
      if (!dense_solve(A, b, x, nf)) fatal("singular capacitance matrix");
      for (int r = 0; r < nf; ++r) cond_[fl[r]].V = x[r];
    }
    // superpose
    const std::size_t N = phi.total();
    std::copy(phi0_.data(), phi0_.data() + N, phi.data());
    for (int k = 0; k < nc; ++k) {
      const double V = cond_[k].V;
      if (V == 0) continue;
      const double* ps = psi_[k].data();
      double* ph = phi.data();
#pragma omp parallel for simd schedule(static)
      for (std::size_t p = 0; p < N; ++p) ph[p] += V * ps[p];
    }
  } else {
    phi.copy_from(phi0_);
  }
  compute_E();
  return st;
}

void FieldSolver::compute_E() {
  const GridGeom& g = *g_;
  const int G = kGhosts;
  const auto& n = phi.size();
  const double ix2 = 0.5 / g.dx[0], iy2 = 0.5 / g.dx[1], iz2 = 0.5 / g.dx[2];
  const std::ptrdiff_t sy = n[0], sz = std::ptrdiff_t(n[0]) * n[1];
  const int on[3] = {g.n1[0] - g.n0[0], g.n1[1] - g.n0[1], g.n1[2] - g.n0[2]};
  const double* P = phi.data();
  bool dlo[3], dhi[3];
  for (int a = 0; a < 3; ++a) {
    dlo[a] = g.wlo[a] && cfg_.faces[2 * a].field == FieldBC::Dirichlet;
    dhi[a] = g.whi[a] && cfg_.faces[2 * a + 1].field == FieldBC::Dirichlet;
  }
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = G; k < G + on[2]; ++k)
    for (int j = G; j < G + on[1]; ++j)
      for (int i = G; i < G + on[0]; ++i) {
        std::size_t p = phi.idx(i, j, k);
        int ii[3] = {i, j, k};
        double e[3];
        const std::ptrdiff_t st[3] = {1, sy, sz};
        const double inv2[3] = {ix2, iy2, iz2};
        for (int a = 0; a < 3; ++a) {
          int gi = g.c0[a] + ii[a] - G;
          if (dlo[a] && gi == 0) e[a] = -(P[p + st[a]] - P[p]) * 2 * inv2[a];
          else if (dhi[a] && gi == g.N[a]) e[a] = -(P[p] - P[p - st[a]]) * 2 * inv2[a];
          else e[a] = -(P[p + st[a]] - P[p - st[a]]) * inv2[a];
        }
        Ex[p] = e[0] + cfg_.E0.x;
        Ey[p] = e[1] + cfg_.E0.y;
        Ez[p] = e[2] + cfg_.E0.z;
      }
  Ex.halo_copy();
  Ey.halo_copy();
  Ez.halo_copy();
  if (eb_.has_objects()) {
    Field3D* E[3] = {&Ex, &Ey, &Ez};
    const double E0[3] = {cfg_.E0.x, cfg_.E0.y, cfg_.E0.z};
    // (1) free nodes next to a surface: gradient of the quadratic through the boundary
    //     point(s) at fractional distance theta (consistent with the Shortley-Weller stencil)
    for (const CutNode& cn : eb_.cut_nodes()) {
      std::size_t p = cn.idx;
      for (int a = 0; a < 3; ++a) {
        double tm = cn.theta[2 * a], tp = cn.theta[2 * a + 1];
        if (tm >= 1.0 && tp >= 1.0) continue;
        const std::ptrdiff_t st = phi.stride(a);
        double pm = P[p - st], pp = P[p + st], pc = P[p];
        double d = (tm * tm * (pp - pc) + tp * tp * (pc - pm)) / (g.dx[a] * tm * tp * (tm + tp));
        (*E[a])[p] = -d + E0[a];
      }
    }
    for (auto* f : E) f->halo_copy();
    // (2) object nodes bordering free space: extrapolate E from the free neighbours so that
    //     particles in cut cells interpolate the correct near-surface field
    static const int dirs[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
    std::vector<std::pair<std::size_t, Vec3>> upd;
    for (int k = G; k < G + on[2]; ++k)
      for (int j = G; j < G + on[1]; ++j)
        for (int i = G; i < G + on[0]; ++i) {
          std::size_t p = phi.idx(i, j, k);
          if (eb_.obj(p) < 0) continue;
          Vec3 acc{0, 0, 0};
          int n = 0;
          for (auto& dd : dirs) {
            std::size_t q = phi.idx(i + dd[0], j + dd[1], k + dd[2]);
            if (eb_.obj(q) >= 0) continue;
            acc += Vec3{Ex[q], Ey[q], Ez[q]};
            ++n;
          }
          upd.emplace_back(p, n ? acc / double(n) : Vec3{E0[0], E0[1], E0[2]});
        }
    for (auto& u : upd) {
      Ex[u.first] = u.second.x;
      Ey[u.first] = u.second.y;
      Ez[u.first] = u.second.z;
    }
    for (auto* f : E) f->halo_copy();
  }
  // interleaved copy for cache-friendly particle gathers
  const std::size_t N = Ex.total();
  E3.resize(3 * N);
#pragma omp parallel for simd schedule(static)
  for (std::size_t q = 0; q < N; ++q) {
    E3[3 * q] = Ex[q];
    E3[3 * q + 1] = Ey[q];
    E3[3 * q + 2] = Ez[q];
  }
}

double FieldSolver::node_volume_factor(int i, int j, int k) const { return volf_[0][i] * volf_[1][j] * volf_[2][k]; }

double FieldSolver::field_energy() const {
  using constants::eps0;
  const GridGeom& g = *g_;
  const int G = kGhosts;
  const double dV = g.dx[0] * g.dx[1] * g.dx[2];
  double s = 0;
  for (int k = G; k < G + g.n1[2] - g.n0[2]; ++k)
    for (int j = G; j < G + g.n1[1] - g.n0[1]; ++j)
      for (int i = G; i < G + g.n1[0] - g.n0[0]; ++i) {
        std::size_t p = phi.idx(i, j, k);
        if (eb_.obj(p) >= 0) continue;
        double ex = Ex[p] - cfg_.E0.x, ey = Ey[p] - cfg_.E0.y, ez = Ez[p] - cfg_.E0.z;
        s += node_volume_factor(i, j, k) * (ex * ex + ey * ey + ez * ez);
      }
  return 0.5 * eps0 * dV * allreduce_sum(s, g.comm);
}

}  // namespace picsp
