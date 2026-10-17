#pragma once
// Monte Carlo collisions with a uniform neutral background (null-collision method).
// Processes: elastic (isotropic in the centre-of-mass frame) and resonant charge
// exchange, with constant cross sections.
#include <vector>

#include "config/Config.hpp"
#include "parallel/Decomposition.hpp"
#include "particles/ParticleSystem.hpp"
#include "util/Random.hpp"

namespace picsp {

class MonteCarloCollisions {
 public:
  void setup(const Config& cfg, ParticleSystem& ps, const Decomposition& d, std::vector<Rng>* rngs);
  void collide(long long step);
  bool enabled() const { return enabled_; }
  long long collisions_last_step() const { return last_; }

 private:
  struct Proc {
    CollisionProcess::Type type;
    double sigma;
  };
  bool enabled_ = false;
  NeutralConfig nc_;
  ParticleSystem* ps_ = nullptr;
  const Decomposition* d_ = nullptr;
  std::vector<Rng>* rng_ = nullptr;
  std::vector<std::vector<Proc>> procs_;  // per species
  double vth_n_ = 0;
  long long last_ = 0;
};

}  // namespace picsp
