#include "collisions/MCC.hpp"

#include <algorithm>
#include <cmath>

#include <omp.h>

#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

double MonteCarloCollisions::sigma(const CollisionProcess& p, double E) {
  switch (p.model) {
    case CollisionProcess::Constant: return p.cross_section;
    case CollisionProcess::Step: return E >= p.threshold ? p.cross_section : 0.0;
    case CollisionProcess::Lotz:
      return E > p.threshold ? p.cross_section * std::log(E / p.threshold) / (E * p.threshold) : 0.0;
    case CollisionProcess::Table: {
      const auto& t = p.table;
      if (E <= t.front().first) return E < t.front().first ? 0.0 : t.front().second;
      if (E >= t.back().first) return t.back().second;
      auto it = std::upper_bound(t.begin(), t.end(), std::make_pair(E, -1e300));
      const auto& b = *it;
      const auto& a = *(it - 1);
      double f = (E - a.first) / (b.first - a.first);
      return a.second + f * (b.second - a.second);
    }
  }
  return 0.0;
}

void MonteCarloCollisions::setup(const Config& cfg, ParticleSystem& ps, const Decomposition& d, std::vector<Rng>* rngs) {
  nc_ = cfg.neutrals;
  ps_ = &ps;
  d_ = &d;
  rng_ = rngs;
  procs_.assign(ps.species.size(), {});
  enabled_ = false;
  if (!nc_.enabled || nc_.density <= 0) return;
  if (nc_.mass <= 0) fatal("collisions.neutral_mass_amu must be > 0");
  for (auto& p : nc_.processes) {
    int s = ps.find_species(p.species);
    if (s < 0) fatal("collision process refers to unknown species '" + p.species + "'");
    Proc pr;
    pr.cfg = p;
    if (p.type == CollisionProcess::Ionization) {
      pr.product = ps.find_species(p.product);
      if (pr.product < 0) fatal("ionization product species '" + p.product + "' not defined");
      pr.secondary = p.secondary.empty() ? s : ps.find_species(p.secondary);
      if (pr.secondary < 0) fatal("ionization secondary species '" + p.secondary + "' not defined");
    }
    procs_[s].push_back(pr);
    enabled_ = true;
  }
  vth_n_ = std::sqrt(nc_.temperature * constants::qe / nc_.mass);
}

void MonteCarloCollisions::collide(long long step) {
  using namespace constants;
  last_ = last_ion_ = 0;
  if (!enabled_) return;
  long long total = 0, total_ion = 0;
  const int nt = omp_get_max_threads();
  for (std::size_t si = 0; si < ps_->species.size(); ++si) {
    Species& sp = ps_->species[si];
    auto& pr = procs_[sp.id];
    if (pr.empty() || !sp.active(step)) continue;
    const double dt = ps_->dt_species(sp.id);
    ParticleArray& p = sp.p;
    const std::size_t np = p.size();
    const double m = sp.m, M = nc_.mass, mu = m * M / (m + M), mu_m = M / (m + M);
    // maximum relative speed (global) and the null-collision frequency max_g n sum(sigma) g
    double vmax = 0;
#pragma omp parallel for reduction(max : vmax) schedule(static)
    for (std::size_t i = 0; i < np; ++i)
      vmax = std::max(vmax, p.vx[i] * p.vx[i] + p.vy[i] * p.vy[i] + p.vz[i] * p.vz[i]);
    vmax = std::sqrt(allreduce_max(vmax, d_->comm)) + norm(nc_.drift) + 6 * vth_n_;
    double nu_max = 0;
    for (int k = 1; k <= 4000; ++k) {
      double g = vmax * k / 4000.0;
      double E = 0.5 * mu * g * g / qe;
      double s = 0;
      for (auto& q : pr) s += sigma(q.cfg, E);
      nu_max = std::max(nu_max, nc_.density * s * g);
    }
    nu_max *= 1.05;
    if (nu_max <= 0) continue;
    const double Pnull = 1.0 - std::exp(-nu_max * dt);
    std::vector<std::vector<double>> born_e(nt), born_i(nt);  // packed new particles per thread
    long long count = 0, nion = 0;
#pragma omp parallel reduction(+ : count, nion)
    {
      const int tid = omp_get_thread_num();
      Rng& rng = (*rng_)[tid];
      auto iso = [&]() {
        double ct = 2 * rng.uniform() - 1, st = std::sqrt(std::max(0.0, 1 - ct * ct));
        double ph = twopi * rng.uniform();
        return Vec3{st * std::cos(ph), st * std::sin(ph), ct};
      };
#pragma omp for schedule(static)
      for (std::size_t i = 0; i < np; ++i) {
        if (rng.uniform() >= Pnull) continue;
        Vec3 v{p.vx[i], p.vy[i], p.vz[i]};
        Vec3 vn = nc_.drift + Vec3{rng.normal(), rng.normal(), rng.normal()} * vth_n_;
        Vec3 g = v - vn;
        double gm = norm(g);
        double E = 0.5 * mu * gm * gm / qe;
        double R = rng.uniform() * nu_max;
        double acc = 0;
        for (auto& q : pr) {
          acc += nc_.density * sigma(q.cfg, E) * gm;
          if (R >= acc) continue;
          Vec3 vcm = (v * m + vn * M) / (m + M);
          switch (q.cfg.type) {
            case CollisionProcess::ChargeExchange: v = vn; break;
            case CollisionProcess::Elastic: v = vcm + iso() * (gm * mu_m); break;
            case CollisionProcess::Excitation: {
              double E2 = std::max(0.0, E - q.cfg.threshold);
              double g2 = std::sqrt(2 * E2 * qe / mu);
              v = vcm + iso() * (g2 * mu_m);
              break;
            }
            case CollisionProcess::Ionization: {
              // available energy shared randomly between primary and secondary electron;
              // the ion is born with the neutral's velocity
              double Ea = std::max(0.0, E - q.cfg.threshold);
              double f = rng.uniform();
              double g1 = std::sqrt(2 * (1 - f) * Ea * qe / mu), g2 = std::sqrt(2 * f * Ea * qe / mu);
              v = vn + iso() * g1;
              Vec3 v2 = vn + iso() * g2;
              double w = p.w[i];
              auto& be = born_e[tid];
              be.insert(be.end(), {p.x[i], p.y[i], p.z[i], v2.x, v2.y, v2.z, w, double(q.secondary)});
              auto& bi = born_i[tid];
              bi.insert(bi.end(), {p.x[i], p.y[i], p.z[i], vn.x, vn.y, vn.z, w, double(q.product)});
              ++nion;
              break;
            }
          }
          p.vx[i] = v.x;
          p.vy[i] = v.y;
          p.vz[i] = v.z;
          ++count;
          break;
        }
      }
    }
    for (int t = 0; t < nt; ++t) {
      for (auto* buf : {&born_e[t], &born_i[t]})
        for (std::size_t k = 0; k < buf->size(); k += 8) {
          const double* b = &(*buf)[k];
          ps_->add_particle(int(b[7]), {b[0], b[1], b[2]}, {b[3], b[4], b[5]}, b[6]);
        }
    }
    total += count;
    total_ion += nion;
  }
  last_ = allreduce_sum(total, d_->comm);
  last_ion_ = allreduce_sum(total_ion, d_->comm);
}

}  // namespace picsp
