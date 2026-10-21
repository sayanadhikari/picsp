#pragma once
// Node-centred scalar field on the local sub-domain with G ghost layers.
// Storage is x-fastest (i + nx*(j + ny*k)), matching VTK ImageData ordering.
#include <array>
#include <cstddef>
#include <vector>

#include <mpi.h>

#include "parallel/Decomposition.hpp"

namespace picsp {

// Geometry of one (possibly coarsened) distributed grid level.
struct GridGeom {
  MPI_Comm comm = MPI_COMM_NULL;
  int rank = 0;
  std::array<std::array<int, 2>, 3> nbr{};
  std::array<bool, 3> periodic{}, wlo{}, whi{};
  std::array<int, 3> N{}, c0{}, c1{}, n0{}, n1{};
  std::array<double, 3> dx{};
  Vec3 origin;

  static GridGeom from(const Decomposition& d);
  bool can_coarsen() const;
  GridGeom coarsen() const;
  long long global_nodes() const;
};

class Field3D {
 public:
  Field3D() = default;
  Field3D(const GridGeom& g, int ghosts);
  void init(const GridGeom& g, int ghosts);

  const GridGeom& geom() const { return *g_; }
  int ghosts() const { return G_; }
  const std::array<int, 3>& size() const { return n_; }
  std::size_t total() const { return v_.size(); }

  // local index (0..n-1, ghosts included)
  inline std::size_t idx(int i, int j, int k) const { return std::size_t(i) + std::size_t(n_[0]) * (std::size_t(j) + std::size_t(n_[1]) * k); }
  inline double& operator()(int i, int j, int k) { return v_[idx(i, j, k)]; }
  inline double operator()(int i, int j, int k) const { return v_[idx(i, j, k)]; }
  inline double& operator[](std::size_t p) { return v_[p]; }
  inline double operator[](std::size_t p) const { return v_[p]; }
  double* data() { return v_.data(); }
  const double* data() const { return v_.data(); }

  // stride along axis
  inline std::ptrdiff_t stride(int a) const { return a == 0 ? 1 : (a == 1 ? n_[0] : std::ptrdiff_t(n_[0]) * n_[1]); }

  // local index of a global node index along axis a
  inline int loc(int a, int gidx) const { return gidx - g_->c0[a] + G_; }
  // owned local range [own_lo, own_hi)
  inline int own_lo(int /*a*/) const { return G_; }
  inline int own_hi(int a) const { return G_ + (g_->n1[a] - g_->n0[a]); }

  void fill(double val);
  void copy_from(const Field3D& o);

  // Ghost exchange. halo_copy: owners -> ghosts. halo_add: ghost contributions summed
  // into owners (deposition), followed by an implicit halo_copy.
  void halo_copy();
  void halo_add();

 private:
  void exchange_axis(int a, bool add);
  const GridGeom* g_ = nullptr;
  int G_ = 0;
  std::array<int, 3> n_{};
  std::vector<double> v_;
  std::vector<double> sbuf_[2], rbuf_[2];
};

// Owned-node global reductions
double dot_owned(const Field3D& a, const Field3D& b);
double sum_owned(const Field3D& a);
double max_abs_owned(const Field3D& a);

}  // namespace picsp
