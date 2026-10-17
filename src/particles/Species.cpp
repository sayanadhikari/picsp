#include "particles/Species.hpp"

#include <cmath>
#include <random>

#include "util/Constants.hpp"

namespace picsp {

void VelocitySampler::init(const SpeciesConfig& s, const Vec3& B0) {
  using constants::qe;
  dist_ = s.dist;
  drift_ = s.drift;
  vth_ = s.temperature > 0 ? std::sqrt(s.temperature * qe / s.mass) : 0.0;
  double tpar = s.t_par >= 0 ? s.t_par : s.temperature;
  double tperp = s.t_perp >= 0 ? s.t_perp : s.temperature;
  vpar_ = std::sqrt(std::max(tpar, 0.0) * qe / s.mass);
  vperp_ = std::sqrt(std::max(tperp, 0.0) * qe / s.mass);
  if (dist_ == VelocityDist::BiMaxwellian || dist_ == VelocityDist::Ring)
    vth_ = std::sqrt((vpar_ * vpar_ + 2 * vperp_ * vperp_) / 3.0);
  kappa_ = s.kappa;
  ring_ = s.ring_speed;
  if (dist_ == VelocityDist::Kappa) ksigma_ = vth_ * std::sqrt((2 * kappa_ - 3) / (2 * kappa_ - 1));
  if (s.aniso_axis) b_ = *s.aniso_axis;
  else if (norm(B0) > 0) b_ = normalized(B0);
  else b_ = {0, 0, 1};
  orthonormal_basis(b_, e1_, e2_);
}

Vec3 VelocitySampler::sample_thermal(Rng& rng) const {
  switch (dist_) {
    case VelocityDist::Cold: return {0, 0, 0};
    case VelocityDist::Maxwellian: return Vec3{rng.normal(), rng.normal(), rng.normal()} * vth_;
    case VelocityDist::BiMaxwellian:
      return b_ * (vpar_ * rng.normal()) + e1_ * (vperp_ * rng.normal()) + e2_ * (vperp_ * rng.normal());
    case VelocityDist::Kappa: {
      // multivariate Student-t with nu = 2 kappa - 1 degrees of freedom
      double nu = 2 * kappa_ - 1;
      std::gamma_distribution<double> gam(0.5 * nu, 2.0);
      double G = gam(rng);
      double s = ksigma_ * std::sqrt(nu / G);
      return Vec3{rng.normal(), rng.normal(), rng.normal()} * s;
    }
    case VelocityDist::Ring: {
      double vp = std::fabs(ring_ + vperp_ * rng.normal());
      double ph = constants::twopi * rng.uniform();
      return b_ * (vpar_ * rng.normal()) + e1_ * (vp * std::cos(ph)) + e2_ * (vp * std::sin(ph));
    }
  }
  return {0, 0, 0};
}

Vec3 VelocitySampler::sample(Rng& rng) const { return drift_ + sample_thermal(rng); }

Vec3 sample_emission_velocity(Rng& rng, const Vec3& n, double vth) {
  Vec3 e1, e2;
  orthonormal_basis(n, e1, e2);
  double vn = vth * std::sqrt(-2.0 * std::log(rng.uniform_open0()));
  return n * vn + e1 * (vth * rng.normal()) + e2 * (vth * rng.normal());
}

Species::Species(const SpeciesConfig& c, int i, const Vec3& B0) : cfg(c), id(i) {
  q = c.charge;
  m = c.mass;
  qm = q / m;
  sampler.init(c, B0);
}

}  // namespace picsp
