#include "collisions/Coulomb.hpp"

#include <algorithm>
#include <cmath>
#include <set>

#include <omp.h>

#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

void CoulombCollisions::setup(const Config& cfg, ParticleSystem& ps, const Decomposition& d, std::vector<Rng>* rngs) {
  using namespace constants;
  enabled_ = cfg.coulomb.enabled;
  ps_ = &ps;
  d_ = &d;
  rng_ = rngs;
  if (!enabled_) return;
  every_ = cfg.coulomb.every;
  dt_ = cfg.dt * every_;
  lnL_ = cfg.coulomb.log_lambda;
  if (cfg.coulomb.auto_log) {
    if (cfg.ref_debye <= 0) fatal("coulomb_log = \"auto\" needs a reference species with density and temperature");
    // ln(Lambda) = ln(12 pi n lambda_D^3) of the reference plasma
    double n = 0;
    for (auto& s : cfg.species)
      if (s.charge < 0 && s.density > 0) {
        n = s.density;
        break;
      }
    lnL_ = std::log(12 * pi * n * std::pow(cfg.ref_debye, 3));
  }
  pairs_.clear();
  const int ns = int(ps.species.size());
  if (cfg.coulomb.pairs.empty()) {
    for (int a = 0; a < ns; ++a)
      for (int b = a; b < ns; ++b)
        if (ps.species[a].q != 0 && ps.species[b].q != 0) pairs_.emplace_back(a, b);
  } else {
    for (auto& pr : cfg.coulomb.pairs) {
      int a = ps.find_species(pr.first), b = ps.find_species(pr.second);
      if (a < 0 || b < 0) fatal("coulomb.pairs: unknown species " + (a < 0 ? pr.first : pr.second));
      pairs_.emplace_back(std::min(a, b), std::max(a, b));
    }
  }
  nx_ = d.c1[0] - d.c0[0];
  ny_ = d.c1[1] - d.c0[1];
  nz_ = d.c1[2] - d.c0[2];
  cellvol_ = d.dx[0] * d.dx[1] * d.dx[2];
  start_.assign(ns, {});
  idx_.assign(ns, {});
  log0(" Coulomb collisions: %zu species pair(s), ln(Lambda) = %.2f, every %d step(s)\n", pairs_.size(), lnL_, every_);
}

void CoulombCollisions::index_cells(int s) {
  const ParticleArray& p = ps_->species[s].p;
  const std::size_t np = p.size(), nc = std::size_t(nx_) * ny_ * nz_;
  const Vec3 lo = d_->local_lo();
  const double ix = 1 / d_->dx[0], iy = 1 / d_->dx[1], iz = 1 / d_->dx[2];
  std::vector<std::uint32_t> cell(np);
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < np; ++i) {
    int a = std::clamp(int((p.x[i] - lo.x) * ix), 0, nx_ - 1);
    int b = std::clamp(int((p.y[i] - lo.y) * iy), 0, ny_ - 1);
    int c = std::clamp(int((p.z[i] - lo.z) * iz), 0, nz_ - 1);
    cell[i] = std::uint32_t(a + std::size_t(nx_) * (b + std::size_t(ny_) * c));
  }
  auto& st = start_[s];
  auto& ix_ = idx_[s];
  st.assign(nc + 1, 0);
  for (std::size_t i = 0; i < np; ++i) st[cell[i] + 1]++;
  for (std::size_t c = 0; c < nc; ++c) st[c + 1] += st[c];
  ix_.resize(np);
  std::vector<std::size_t> fill(st.begin(), st.end() - 1);
  for (std::size_t i = 0; i < np; ++i) ix_[fill[cell[i]]++] = i;
}

namespace {
// Takizuka-Abe scattering of the relative velocity; returns delta u
inline Vec3 ta_delta(const Vec3& u, double var, Rng& rng) {
  double U = norm(u);
  if (U <= 0) return {0, 0, 0};
  double sint, omc;  // sin(theta), 1 - cos(theta)
  if (var > 1.0) {   // strong scattering: isotropic
    double ct = 2 * rng.uniform() - 1;
    sint = std::sqrt(std::max(0.0, 1 - ct * ct));
    omc = 1 - ct;
  } else {
    double del = std::sqrt(var) * rng.normal();
    double d2 = del * del;
    sint = 2 * del / (1 + d2);
    omc = 2 * d2 / (1 + d2);
  }
  double phi = constants::twopi * rng.uniform();
  double cp = std::cos(phi), sp = std::sin(phi);
  double up = std::sqrt(u.x * u.x + u.y * u.y);
  if (up < 1e-12 * U) return {U * sint * cp, U * sint * sp, -U * omc};
  return {u.x / up * u.z * sint * cp - u.y / up * U * sint * sp - u.x * omc,
          u.y / up * u.z * sint * cp + u.x / up * U * sint * sp - u.y * omc, -up * sint * cp - u.z * omc};
}
}  // namespace

void CoulombCollisions::collide(long long step) {
  using namespace constants;
  if (!enabled_ || step % every_ != 0) return;
  std::set<int> used;
  for (auto& pr : pairs_) {
    used.insert(pr.first);
    used.insert(pr.second);
  }
  for (int s : used) index_cells(s);
  const std::size_t nc = std::size_t(nx_) * ny_ * nz_;
  for (auto& pr : pairs_) {
    Species& A = ps_->species[pr.first];
    Species& B = ps_->species[pr.second];
    const double qq2 = A.q * A.q * B.q * B.q;
    const double mu = A.m * B.m / (A.m + B.m);
    const double fa = B.m / (A.m + B.m), fb = A.m / (A.m + B.m);
    const double pref = qq2 * lnL_ / (8 * pi * eps0 * eps0 * mu * mu);
    ParticleArray& pa = A.p;
    ParticleArray& pb = B.p;
    const auto& sa = start_[pr.first];
    const auto& sb = start_[pr.second];
    const auto& ia = idx_[pr.first];
    const auto& ib = idx_[pr.second];
    const bool intra = pr.first == pr.second;
#pragma omp parallel
    {
      Rng& rng = (*rng_)[omp_get_thread_num()];
      std::vector<std::size_t> la, lb;
      auto scatter = [&](std::size_t i, std::size_t j, double nfac, double dt) {
        Vec3 u{pa.vx[i] - pb.vx[j], pa.vy[i] - pb.vy[j], pa.vz[i] - pb.vz[j]};
        double U = norm(u);
        if (U <= 0) return;
        double var = pref * nfac * dt / (U * U * U);
        Vec3 du = ta_delta(u, var, rng);
        double wi = pa.w[i], wj = pb.w[j], wmax = std::max(wi, wj);
        if (wj >= wmax || rng.uniform() * wmax < wj) {
          pa.vx[i] += fa * du.x;
          pa.vy[i] += fa * du.y;
          pa.vz[i] += fa * du.z;
        }
        if (wi >= wmax || rng.uniform() * wmax < wi) {
          pb.vx[j] -= fb * du.x;
          pb.vy[j] -= fb * du.y;
          pb.vz[j] -= fb * du.z;
        }
      };
#pragma omp for schedule(dynamic, 64)
      for (std::size_t c = 0; c < nc; ++c) {
        if (intra) {
          std::size_t n = sa[c + 1] - sa[c];
          if (n < 2) continue;
          la.assign(ia.begin() + sa[c], ia.begin() + sa[c + 1]);
          for (std::size_t k = n - 1; k > 0; --k) std::swap(la[k], la[std::size_t(rng.uniform() * (k + 1))]);
          double wsum = 0, wmax = 0;
          for (auto i : la) {
            wsum += pa.w[i];
            wmax = std::max(wmax, pa.w[i]);
          }
          const double nfac = n * wmax / cellvol_;  // density with the Nanbu-Yonemura weight factor
          std::size_t k0 = 0;
          if (n % 2) {  // odd: first three particles collide pairwise with dt/2 (Takizuka-Abe)
            scatter(la[0], la[1], nfac, 0.5 * dt_);
            scatter(la[1], la[2], nfac, 0.5 * dt_);
            scatter(la[2], la[0], nfac, 0.5 * dt_);
            k0 = 3;
          }
          for (std::size_t k = k0; k + 1 < n; k += 2) scatter(la[k], la[k + 1], nfac, dt_);
          (void)wsum;
        } else {
          std::size_t na = sa[c + 1] - sa[c], nb = sb[c + 1] - sb[c];
          if (na == 0 || nb == 0) continue;
          la.assign(ia.begin() + sa[c], ia.begin() + sa[c + 1]);
          lb.assign(ib.begin() + sb[c], ib.begin() + sb[c + 1]);
          for (std::size_t k = na - 1; k > 0; --k) std::swap(la[k], la[std::size_t(rng.uniform() * (k + 1))]);
          for (std::size_t k = nb - 1; k > 0; --k) std::swap(lb[k], lb[std::size_t(rng.uniform() * (k + 1))]);
          double wmax = 0;
          for (auto i : la) wmax = std::max(wmax, pa.w[i]);
          for (auto j : lb) wmax = std::max(wmax, pb.w[j]);
          // every particle of the more numerous species collides once; partners of the
          // other species are reused cyclically, so the density factor is that species'
          const double nfac = std::min(na, nb) * wmax / cellvol_;
          if (na >= nb)
            for (std::size_t k = 0; k < na; ++k) scatter(la[k], lb[k % nb], nfac, dt_);
          else
            for (std::size_t k = 0; k < nb; ++k) scatter(la[k % na], lb[k], nfac, dt_);
        }
      }
    }
  }
}

}  // namespace picsp
