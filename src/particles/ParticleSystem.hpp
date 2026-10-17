#pragma once
// All particle operations on the local sub-domain: loading, Boris push with
// particle-object collision detection, outer boundaries, MPI migration, injection
// through open boundaries, charge/moment deposition and cell sorting.
#include <array>
#include <vector>

#include "config/Config.hpp"
#include "fields/FieldSolver.hpp"
#include "geometry/Geometry.hpp"
#include "grid/Field3D.hpp"
#include "parallel/Decomposition.hpp"
#include "particles/Species.hpp"
#include "util/Random.hpp"

namespace picsp {

struct Hit {
  std::size_t idx;  // particle index (before removal)
  int obj;          // object id
  Vec3 pos;         // impact point
  Vec3 vel;         // impact velocity
  double w;         // macro weight
};

struct FaceInjector {
  int axis, side;   // face
  int species;
  double flux;      // particles / (m^2 s)
  double vcap;      // rejection bound for the normal velocity
  double area;      // local face area
  // fast exact path for drifting Maxwellians: inverse CDF of the flux-weighted normal speed
  bool maxwellian = false;
  double vth = 0;
  Vec3 drift;
  std::vector<double> icdf;
};

class ParticleSystem {
 public:
  void setup(const Config& cfg, const Decomposition& d, const GridGeom& g, const Scene& scene, const FieldSolver& fs,
             std::vector<Rng>* rngs);

  std::vector<Species> species;

  void load_initial();
  // Advance velocities by -dt/2 at t = 0 (leapfrog start)
  void rewind_half();
  // Boris push of active species; collisions with objects are returned per species
  void push(long long step, std::vector<std::vector<Hit>>& hits);
  void remove_hits(int s, std::vector<Hit>& hits);
  void inject(long long step);
  void boundaries_and_migrate(long long step);
  void sort();
  void deposit_charge(Field3D& rho);
  // moments for species s: density [m^-3], mean velocity [m/s], temperature [eV]
  void deposit_moments(int s, Field3D& n, Field3D& ux, Field3D& uy, Field3D& uz, Field3D& T);

  // add particle (emission); position must be inside the global domain
  void add_particle(int s, const Vec3& pos, const Vec3& vel, double w) { species[s].p.push(pos, vel, w); }

  long long count_global(int s) const;
  double kinetic_energy(int s) const;  // global [J]
  double dt_species(int s) const { return dt_ * species[s].cfg.subcycle; }
  int find_species(const std::string& name) const;
  // true if p lies inside an object (cheap near-surface cell test first)
  bool inside_object(const Vec3& p) const;
  int default_electron_species() const;
  const std::vector<FaceInjector>& injectors() const { return inj_; }

 private:
  template <int O>
  void push_impl(Species& sp, double dt, std::vector<Hit>& hits);
  template <int O>
  void deposit_impl(const Species& sp, std::vector<double>& buf, int mode, double* const* out);
  void reduce_thread_buffers(Field3D& out);

  Config cfg_;
  const Decomposition* d_ = nullptr;
  const GridGeom* g_ = nullptr;
  const Scene* scene_ = nullptr;
  const FieldSolver* fs_ = nullptr;
  std::vector<Rng>* rng_ = nullptr;
  double dt_ = 0;
  int G_ = 2;
  std::array<double, 3> arr_lo_{}, inv_dx_{};  // field-array origin and 1/dx
  std::vector<FaceInjector> inj_;
  std::vector<std::vector<double>> tbuf_;     // per-thread deposition buffers
  std::vector<double> sendbuf_[2], recvbuf_[2];
};

}  // namespace picsp
