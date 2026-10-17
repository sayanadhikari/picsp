#pragma once
// Cartesian MPI domain decomposition of a node-centred uniform grid.
//
// Global cells 0..N-1 along each axis, nodes 0..N. For periodic axes node N is the
// image of node 0. Each rank owns cells [c0, c1) and nodes [n0, n1) where n1 = c1,
// except the last rank of a non-periodic axis which also owns node N (n1 = c1 + 1).
#include <array>

#include <mpi.h>

#include "config/Config.hpp"
#include "util/Vec3.hpp"

namespace picsp {

struct Decomposition {
  MPI_Comm comm = MPI_COMM_NULL;  // Cartesian communicator
  int rank = 0, size = 1;
  std::array<int, 3> dims{1, 1, 1}, coords{0, 0, 0};
  std::array<bool, 3> periodic{false, false, false};
  std::array<std::array<int, 2>, 3> nbr{};  // neighbour rank [axis][lo/hi], MPI_PROC_NULL at walls
  std::array<int, 3> N{};                   // global cells
  std::array<int, 3> c0{}, c1{};            // owned cells
  std::array<int, 3> n0{}, n1{};            // owned nodes
  std::array<double, 3> dx{};
  Vec3 origin;
  // physical (non-periodic) boundary on this rank's lo/hi face
  bool wall_lo(int a) const { return !periodic[a] && coords[a] == 0; }
  bool wall_hi(int a) const { return !periodic[a] && coords[a] == dims[a] - 1; }

  Vec3 domain_lo() const { return origin; }
  Vec3 domain_hi() const { return {origin.x + N[0] * dx[0], origin.y + N[1] * dx[1], origin.z + N[2] * dx[2]}; }
  Vec3 local_lo() const { return {origin.x + c0[0] * dx[0], origin.y + c0[1] * dx[1], origin.z + c0[2] * dx[2]}; }
  Vec3 local_hi() const { return {origin.x + c1[0] * dx[0], origin.y + c1[1] * dx[1], origin.z + c1[2] * dx[2]}; }
  long long local_cells() const { return 1LL * (c1[0] - c0[0]) * (c1[1] - c0[1]) * (c1[2] - c0[2]); }

  static Decomposition create(const Config& cfg, MPI_Comm world);
};

// Global reductions
double allreduce_sum(double v, MPI_Comm comm);
double allreduce_max(double v, MPI_Comm comm);
double allreduce_min(double v, MPI_Comm comm);
long long allreduce_sum(long long v, MPI_Comm comm);

}  // namespace picsp
