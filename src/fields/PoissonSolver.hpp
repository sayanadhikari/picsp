#pragma once
// Distributed-memory solver for  -lap(phi) = rho/eps0  on the node grid with
//   * Dirichlet / Neumann / periodic outer boundaries,
//   * embedded conductors (Dirichlet object nodes),
//   * optional Shortley-Weller cut-cell stencils next to object surfaces.
// Method: Krylov (PCG for the symmetric staircase operator, BiCGStab when cut cells
// make it non-symmetric) preconditioned by a geometric multigrid V-cycle with
// red-black Gauss-Seidel smoothing, full-weighting restriction (= P^T/8) and
// trilinear prolongation. Neumann rows are scaled by the control-volume fraction so
// that the staircase operator stays symmetric.
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "config/Config.hpp"
#include "fields/EmbeddedBoundary.hpp"
#include "grid/Field3D.hpp"

namespace picsp {

struct SolveStats {
  int iterations = 0;
  double residual = 0;  // relative
  bool converged = false;
};

class PoissonSolver {
 public:
  void setup(const Config& cfg, const GridGeom& g0, const EmbeddedBoundary& eb, int phi_ghosts);

  // phi (G = phi_ghosts layout): on entry initial guess, on exit solution.
  // lift: values at fixed nodes (object potentials, Dirichlet walls), any value elsewhere.
  // rho: charge density [C/m^3].
  SolveStats solve(Field3D& phi, const Field3D& rho, const Field3D& lift);

  // Node classification in the phi (G = phi_ghosts) layout
  bool is_fixed_phi_layout(int i, int j, int k) const;
  // Control-volume weight (Neumann faces: 1/2 per face) for owned node, phi layout indices
  double volume_weight_phi_layout(int i, int j, int k) const;
  bool singular() const { return singular_; }
  int levels() const { return int(L_.size()); }
  bool symmetric() const { return !use_bicgstab_; }
  // coefficient of the link from owned free node (phi layout) towards neighbour dir (0..5)
  double link_coef_phi_layout(int i, int j, int k, int dir) const;

 private:
  struct Level {
    GridGeom geom;
    std::vector<std::uint8_t> fixed;   // per local node (G = 1)
    std::vector<int> cutid;            // index into cut list, -1 if regular (level 0)
    std::vector<CutNode> cut;          // with idx in this level's layout
    std::array<std::vector<double>, 3> wfac;  // Neumann weight factors per axis (local index)
    std::vector<double> iwfac0;                // 1 / wfac[0]
    std::array<double, 3> c;           // 1/dx^2
    Field3D u, f, r;                   // work arrays
    int own[3];
  };

  void fill_ghosts(Level& L, Field3D& x);
  void apply(Level& L, Field3D& x, Field3D& y);
  void residual(Level& L, Field3D& u, const Field3D& f, Field3D& r);
  void smooth(Level& L, Field3D& u, const Field3D& f, int sweeps, bool reverse);
  void restrict_to(Level& F, const Field3D& rf, Level& C, Field3D& fc);
  void prolong_add(Level& C, Field3D& ec, Level& F, Field3D& uf);
  void vcycle(int l);
  void precondition(const Field3D& r, Field3D& z);
  void project_mean(Level& L, Field3D& b);
  void zero_fixed(Level& L, Field3D& x);

  SolveStats solve_pcg(Field3D& u, const Field3D& b);
  SolveStats solve_bicgstab(Field3D& u, const Field3D& b);
  SolveStats solve_mg(Field3D& u, const Field3D& b);

  SolverConfig sc_;
  std::array<FaceConfig, 6> faces_;
  std::vector<Level> L_;
  int phiG_ = 2;
  bool singular_ = false;
  bool use_bicgstab_ = false;
  bool has_cut_ = false;
  int precond_ = 1;  // 0 none (Jacobi), 1 multigrid
  Field3D b_, x_, p_, q_, z_, s_, t_, rhat_, v_, w_;
};

}  // namespace picsp
