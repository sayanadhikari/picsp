#pragma once
// Coulomb (small-angle, cumulative) collisions with the binary Monte Carlo method of
// Takizuka & Abe (J. Comput. Phys. 25, 205, 1977). Particles are paired randomly within
// each grid cell (intra- and inter-species); each pair scatters by an angle drawn from
//   <tan^2(theta/2)> = q_a^2 q_b^2 n lnL dt / (8 pi eps0^2 mu^2 u^3).
// Unequal macro-particle weights are handled with the Nanbu-Yonemura acceptance rule
// (each partner is updated with probability w_other / w_max), which keeps the correct
// collision rates for both species; momentum and energy are conserved on average
// (exactly for equal weights).
#include <string>
#include <vector>

#include "config/Config.hpp"
#include "parallel/Decomposition.hpp"
#include "particles/ParticleSystem.hpp"
#include "util/Random.hpp"

namespace picsp {

class CoulombCollisions {
 public:
  void setup(const Config& cfg, ParticleSystem& ps, const Decomposition& d, std::vector<Rng>* rngs);
  bool enabled() const { return enabled_; }
  void collide(long long step);
  double coulomb_log() const { return lnL_; }

 private:
  void index_cells(int s);
  bool enabled_ = false;
  int every_ = 1;
  double lnL_ = 10;
  double dt_ = 0;
  ParticleSystem* ps_ = nullptr;
  const Decomposition* d_ = nullptr;
  std::vector<Rng>* rng_ = nullptr;
  std::vector<std::pair<int, int>> pairs_;
  int nx_ = 0, ny_ = 0, nz_ = 0;
  double cellvol_ = 0;
  std::vector<std::vector<std::size_t>> start_, idx_;  // per species: CSR cell -> particle indices
};

}  // namespace picsp
