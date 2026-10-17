#pragma once
// Electrostatic field solver with embedded conductors.
//
// Conductors (one or more objects electrically connected via `group`) are either held
// at a prescribed (possibly time-ramped) bias, or floating. Floating potentials are
// obtained with the capacitance-matrix method:
//   phi = phi_0 + sum_k V_k psi_k
// where phi_0 solves Poisson with all conductors grounded and psi_k are the vacuum
// solutions with conductor k at 1 V. The conductor potentials V follow from charge
// conservation  Q_i = Q_ind,i[phi_0] + sum_k C_ik V_k  with C_ik = Q_i[psi_k].
// Only one Poisson solve per time step is required.
#include <string>
#include <vector>

#include "config/Config.hpp"
#include "fields/EmbeddedBoundary.hpp"
#include "fields/PoissonSolver.hpp"
#include "geometry/Geometry.hpp"
#include "grid/Field3D.hpp"

namespace picsp {

struct Conductor {
  std::string name;
  std::vector<int> objects;
  bool fixed = false;
  double V0 = 0, ramp = 0;  // fixed bias V(t) = V0 + ramp t
  double V = 0;             // current potential
  double Q = 0;             // net charge collected (floating) [C]
  double Cext = 0, Iext = 0;
  double Q_plasma_step = 0; // charge collected during the current step (diagnostic)
};

class FieldSolver {
 public:
  static constexpr int kGhosts = 2;

  void setup(const Config& cfg, const GridGeom& g, const Scene& scene);
  // rho must already contain the full (halo-summed) charge density
  SolveStats solve(double time);

  const EmbeddedBoundary& eb() const { return eb_; }
  const PoissonSolver& poisson() const { return ps_; }
  int object_conductor(int obj) const { return obj < 0 ? -1 : obj_to_cond_[obj]; }
  std::vector<Conductor>& conductors() { return cond_; }
  const std::vector<Conductor>& conductors() const { return cond_; }
  const std::vector<double>& capacitance() const { return C_; }

  double field_energy() const;  // 0.5 eps0 int E^2 dV (free nodes)
  void mark_charges_initialised() { charges_initialised_ = true; }
  Field3D& solver_state() { return phi0_; }  // warm-start potential (checkpointed)
  double node_volume_factor(int i, int j, int k) const;  // 1/2 per non-periodic wall face

  Field3D phi, rho, Ex, Ey, Ez;
  std::vector<double> E3;  // interleaved (Ex,Ey,Ez) per node for particle gathers
  // static external magnetic field on the nodes (only when spatially varying)
  bool B_varies() const { return B_varies_; }
  Field3D Bx, By, Bz;
  std::vector<double> B3;  // interleaved (Bx,By,Bz)

 private:
  void compute_E();
  // Local contribution to the charge on conductor k for potential f and space charge
  // rho (nullptr = vacuum). Conservative discrete Gauss law: regular-stencil flux out of
  // the conductor's object nodes plus the cut-cell excess of adjacent free nodes.
  double conductor_charge(const Field3D& f, const Field3D* rho, int k) const;
  void build_links();

  struct Link {
    std::size_t idx, nb;  // free node, object neighbour (phi layout)
    double coef;          // eps0 * W * dV / dx_a^2
    int cond;
  };
  struct CutShare {
    std::size_t idx;      // free cut node (phi layout)
    double wdv;           // eps0 * W * dV
    int cond;
    double share;         // fraction attributed to this conductor
  };
  std::vector<CutShare> cut_shares_;

  const GridGeom* g_ = nullptr;
  Config cfg_;
  EmbeddedBoundary eb_;
  PoissonSolver ps_;
  Field3D phi0_, lift_;
  std::vector<Field3D> psi_;
  std::vector<Conductor> cond_;
  std::vector<int> obj_to_cond_;
  std::vector<Link> links_;
  std::vector<double> C_;  // nc x nc capacitance matrix (row-major)
  bool charges_initialised_ = false;
  bool B_varies_ = false;
  std::array<std::vector<double>, 3> volf_;  // wall volume factors per axis
};

// Solve the dense system A x = b (n x n, row-major) with partial pivoting.
bool dense_solve(std::vector<double> A, std::vector<double> b, std::vector<double>& x, int n);

}  // namespace picsp
