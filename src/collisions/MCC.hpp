#pragma once
// Monte Carlo collisions with a uniform neutral background (null-collision method).
// Processes: elastic (isotropic in the centre-of-mass frame), resonant charge exchange,
// excitation (inelastic energy loss) and electron-impact ionisation, which creates a
// secondary electron and an ion. Cross sections: constant, step above a threshold,
// Lotz-type ionisation formula, or tabulated sigma(E) (inline or from a two-column file).
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
  long long ionizations_last_step() const { return last_ion_; }

  // cross section [m^2] of a process at centre-of-mass energy E [eV]
  static double sigma(const CollisionProcess& p, double E);

 private:
  struct Proc {
    CollisionProcess cfg;
    int product = -1, secondary = -1;
  };
  bool enabled_ = false;
  NeutralConfig nc_;
  ParticleSystem* ps_ = nullptr;
  const Decomposition* d_ = nullptr;
  std::vector<Rng>* rng_ = nullptr;
  std::vector<std::vector<Proc>> procs_;  // per species
  double vth_n_ = 0;
  long long last_ = 0, last_ion_ = 0;
};

}  // namespace picsp
