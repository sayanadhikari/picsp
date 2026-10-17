#pragma once
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "collisions/Coulomb.hpp"
#include "collisions/MCC.hpp"
#include "config/Config.hpp"
#include "fields/FieldSolver.hpp"
#include "geometry/Geometry.hpp"
#include "io/Vtk.hpp"
#include "objects/ObjectManager.hpp"
#include "parallel/Decomposition.hpp"
#include "particles/ParticleSystem.hpp"
#include "util/Random.hpp"
#include "util/Timer.hpp"

namespace picsp {

class Simulation {
 public:
  Simulation(const Config& cfg, MPI_Comm world);
  ~Simulation();
  void run();
  // single steps (for tests / embedding)
  void initialise();
  void advance();
  long long step() const { return step_; }
  double time() const { return time_; }

  FieldSolver& fields() { return fs_; }
  ParticleSystem& particles() { return ps_; }
  ObjectManager& objects() { return om_; }
  const Config& config() const { return cfg_; }
  const Decomposition& decomposition() const { return d_; }

 private:
  void write_fields();
  void write_particles();
  void write_surface();
  void history_line();
  void status_line();
  void checkpoint();
  void write_tracks();
  std::ofstream tracks_;
  bool restart(const std::string& dir);

  Config cfg_;
  Decomposition d_;
  GridGeom g_;
  Scene scene_;
  FieldSolver fs_;
  ParticleSystem ps_;
  ObjectManager om_;
  MonteCarloCollisions mcc_;
  CoulombCollisions coulomb_;
  VtkWriter vtk_;
  std::vector<Rng> rng_;
  Timers timers_;
  long long step_ = 0;
  double time_ = 0;
  bool initialised_ = false;
  SolveStats last_solve_;
  int solve_iters_acc_ = 0, solve_count_ = 0;
  // history accumulators (global)
  std::vector<double> hist_coll_, hist_emit_;
  double hist_time_ = 0;
  std::vector<double> hist_wall_, hist_inj_;
  std::ofstream hist_;
  double wall_start_ = 0, wall_last_status_ = 0;
  long long pushes_ = 0, pushes_last_status_ = 0;
  long long step_last_status_ = 0;
};

}  // namespace picsp
