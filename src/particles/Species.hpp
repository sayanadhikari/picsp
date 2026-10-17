#pragma once
#include <string>
#include <vector>

#include "config/Config.hpp"
#include "particles/ParticleArray.hpp"
#include "util/Random.hpp"

namespace picsp {

// Velocity distribution sampler for a species (thread-safe given a per-thread Rng)
class VelocitySampler {
 public:
  void init(const SpeciesConfig& s, const Vec3& B0);
  Vec3 sample(Rng& rng) const;            // full velocity including drift
  Vec3 sample_thermal(Rng& rng) const;    // without drift
  double vth() const { return vth_; }     // sqrt(kT/m) (isotropic equivalent)

 private:
  VelocityDist dist_ = VelocityDist::Maxwellian;
  Vec3 drift_;
  double vth_ = 0, vpar_ = 0, vperp_ = 0, kappa_ = 4, ring_ = 0, ksigma_ = 0;
  Vec3 b_{0, 0, 1}, e1_{1, 0, 0}, e2_{0, 1, 0};
};

// Thermal emission from a surface: flux-weighted (half-)Maxwellian along normal n
Vec3 sample_emission_velocity(Rng& rng, const Vec3& n, double vth);

struct SpeciesStats {
  double absorbed_wall_charge = 0;  // charge lost through outer walls (this step)
  double injected = 0;              // macro particles injected (this step)
  long long migrated = 0;
};

class Species {
 public:
  Species(const SpeciesConfig& cfg, int id, const Vec3& B0);
  const std::string& name() const { return cfg.name; }

  SpeciesConfig cfg;
  int id;
  double q, m, qm;
  double w0 = 0;  // nominal macro weight
  VelocitySampler sampler;
  ParticleArray p;
  SpeciesStats stats;

  bool active(long long step) const { return step % cfg.subcycle == 0; }
};

}  // namespace picsp
