#include "collisions/MCC.hpp"

#include <cmath>

#include <omp.h>

#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

void MonteCarloCollisions::setup(const Config& cfg, ParticleSystem& ps, const Decomposition& d, std::vector<Rng>* rngs) {
  nc_ = cfg.neutrals;
  ps_ = &ps;
  d_ = &d;
  rng_ = rngs;
  procs_.assign(ps.species.size(), {});
  enabled_ = false;
  if (!nc_.enabled || nc_.density <= 0) return;
  for (auto& p : nc_.processes) {
    int s = ps.find_species(p.species);
    if (s < 0) fatal("collision process refers to unknown species '" + p.species + "'");
    procs_[s].push_back({p.type, p.cross_section});
    enabled_ = true;
  }
  if (nc_.mass <= 0) fatal("collisions.neutral_mass_amu must be > 0");
  vth_n_ = std::sqrt(nc_.temperature * constants::qe / nc_.mass);
}

void MonteCarloCollisions::collide(long long step) {
  last_ = 0;
  if (!enabled_) return;
  long long total = 0;
  for (auto& sp : ps_->species) {
    auto& pr = procs_[sp.id];
    if (pr.empty() || !sp.active(step)) continue;
    const double dt = ps_->dt_species(sp.id);
    ParticleArray& p = sp.p;
    const std::size_t np = p.size();
    double sig_tot = 0;
    for (auto& q : pr) sig_tot += q.sigma;
    // maximum relative speed (global) for the null-collision frequency
    double vmax = 0;
#pragma omp parallel for reduction(max : vmax) schedule(static)
    for (std::size_t i = 0; i < np; ++i)
      vmax = std::max(vmax, p.vx[i] * p.vx[i] + p.vy[i] * p.vy[i] + p.vz[i] * p.vz[i]);
    vmax = std::sqrt(allreduce_max(vmax, d_->comm)) + norm(nc_.drift) + 6 * vth_n_;
    const double nu_max = nc_.density * sig_tot * vmax;
    const double Pnull = 1.0 - std::exp(-nu_max * dt);
    const double m = sp.m, M = nc_.mass, mu_m = M / (m + M);
    long long count = 0;
#pragma omp parallel reduction(+ : count)
    {
      Rng& rng = (*rng_)[omp_get_thread_num()];
#pragma omp for schedule(static)
      for (std::size_t i = 0; i < np; ++i) {
        if (rng.uniform() >= Pnull) continue;
        Vec3 v{p.vx[i], p.vy[i], p.vz[i]};
        Vec3 vn = nc_.drift + Vec3{rng.normal(), rng.normal(), rng.normal()} * vth_n_;
        Vec3 g = v - vn;
        double gm = norm(g);
        double R = rng.uniform() * nu_max;
        double acc = 0;
        for (auto& q : pr) {
          acc += nc_.density * q.sigma * gm;
          if (R < acc) {
            if (q.type == CollisionProcess::ChargeExchange) {
              v = vn;
            } else {
              // isotropic scattering in the centre-of-mass frame
              Vec3 vcm = (v * m + vn * M) / (m + M);
              double ct = 2 * rng.uniform() - 1, st = std::sqrt(std::max(0.0, 1 - ct * ct));
              double ph = constants::twopi * rng.uniform();
              Vec3 gp{gm * st * std::cos(ph), gm * st * std::sin(ph), gm * ct};
              v = vcm + gp * mu_m;
            }
            p.vx[i] = v.x;
            p.vy[i] = v.y;
            p.vz[i] = v.z;
            ++count;
            break;
          }
        }
      }
    }
    total += count;
  }
  last_ = allreduce_sum(total, d_->comm);
}

}  // namespace picsp
