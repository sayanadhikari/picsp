#pragma once
// Particle-surface interaction for embedded objects:
//   * absorption of particles -> conductor charging, per-panel current densities
//   * secondary electron emission (constant yield or Sternglass law, ion-induced)
//   * photoemission from sunlit, unshadowed surface panels
// The object surface is triangulated with marching tetrahedra on a refined sampling of
// the signed distance function ("panels"); panels are owned by the rank owning the cell.
#include <string>
#include <vector>

#include "config/Config.hpp"
#include "fields/FieldSolver.hpp"
#include "geometry/Geometry.hpp"
#include "particles/ParticleSystem.hpp"
#include "util/Random.hpp"

namespace picsp {

struct Panel {
  Vec3 a, b, c;
  Vec3 centroid, normal;
  double area = 0;
  int obj = -1;
  double sun_cos = 0;  // cos(incidence) if sunlit, 0 if dark or shadowed
};

class ObjectManager {
 public:
  void setup(const Config& cfg, const Decomposition& d, const GridGeom& g, const Scene& scene, FieldSolver& fs,
             ParticleSystem& ps, std::vector<Rng>* rngs);
  // absorb hits (removing the particles), emit secondaries
  void process_hits(long long step, std::vector<std::vector<Hit>>& hits);
  void photoemission(long long step);
  // add the step's external currents and finalise the per-step bookkeeping
  void end_step();
  // reset accumulated panel currents (after output)
  void reset_panel_accumulators();

  const std::vector<Panel>& panels() const { return panels_; }
  // accumulated collected charge per panel and species since last reset [C]
  const std::vector<double>& panel_charge() const { return panel_q_; }
  const std::vector<double>& panel_emitted() const { return panel_emit_; }
  double accum_time() const { return accum_time_; }
  int nspecies() const { return ns_; }

  // per-object (global, this step) collected charge by species and emitted charge
  const std::vector<double>& step_collected() const { return step_coll_; }  // [obj * ns + s]
  const std::vector<double>& step_emitted() const { return step_emit_; }    // [obj]
  long long global_panel_count() const;
  // Per-panel current density [A/m^2] of an accumulated per-panel charge array (stride
  // `ncomp`, component `comp`), averaged over all panels in the same grid cell so that
  // sliver triangles of the triangulation do not produce spurious peaks.
  std::vector<double> cell_averaged_density(const std::vector<double>& q, int ncomp, int comp) const;
  bool empty() const { return objs_ == 0; }

 private:
  void build_panels(const Scene& scene);
  int find_panel(const Vec3& p) const;
  double see_yield(const ObjectConfig& oc, int s, double energy_eV) const;
  void emit(int obj, int sp, const Vec3& pos, const Vec3& n, double vth, double w, int panel);

  Config cfg_;
  const Decomposition* d_ = nullptr;
  const GridGeom* g_ = nullptr;
  const Scene* scene_ = nullptr;
  FieldSolver* fs_ = nullptr;
  ParticleSystem* ps_ = nullptr;
  std::vector<Rng>* rng_ = nullptr;
  int objs_ = 0, ns_ = 0;
  std::vector<int> emit_species_;  // per object
  std::vector<Panel> panels_;
  std::vector<int> cell_start_, cell_panels_;  // CSR: owned cell -> panels
  std::vector<int> panel_cell_;               // owned cell of each panel
  std::vector<double> panel_q_, panel_emit_;
  double accum_time_ = 0;
  std::vector<double> step_coll_, step_emit_;
  std::vector<double> loc_coll_, loc_emit_;
};

}  // namespace picsp
