#include "parallel/Decomposition.hpp"

#include <string>

#include "util/Log.hpp"

namespace picsp {

Decomposition Decomposition::create(const Config& cfg, MPI_Comm world) {
  Decomposition d;
  int wsize = 1;
  MPI_Comm_size(world, &wsize);
  int dims[3] = {cfg.decomposition[0], cfg.decomposition[1], cfg.decomposition[2]};
  if (dims[0] * dims[1] * dims[2] == 0) {
    // Automatic: let MPI factorise, but bias towards splitting the longest axes by
    // giving each axis a weight equal to its cell count.
    int best[3] = {0, 0, 0};
    double best_cost = 1e300;
    for (int px = 1; px <= wsize; ++px) {
      if (wsize % px) continue;
      for (int py = 1; py <= wsize / px; ++py) {
        if ((wsize / px) % py) continue;
        int pz = wsize / px / py;
        if ((dims[0] && dims[0] != px) || (dims[1] && dims[1] != py) || (dims[2] && dims[2] != pz)) continue;
        double lx = double(cfg.cells[0]) / px, ly = double(cfg.cells[1]) / py, lz = double(cfg.cells[2]) / pz;
        if (lx < 3 || ly < 3 || lz < 3) continue;
        double cost = lx * ly + ly * lz + lx * lz;  // halo surface
        if (cost < best_cost) {
          best_cost = cost;
          best[0] = px; best[1] = py; best[2] = pz;
        }
      }
    }
    if (best[0] == 0) fatal("cannot decompose grid over " + std::to_string(wsize) + " ranks (need >= 3 cells per rank)");
    for (int i = 0; i < 3; ++i) dims[i] = best[i];
  }
  if (dims[0] * dims[1] * dims[2] != wsize)
    fatal("parallel.decomposition does not match the number of MPI ranks");
  int periods[3];
  for (int a = 0; a < 3; ++a) {
    d.periodic[a] = cfg.faces[2 * a].field == FieldBC::Periodic;
    periods[a] = d.periodic[a] ? 1 : 0;
    d.dims[a] = dims[a];
  }
  MPI_Cart_create(world, 3, dims, periods, 0, &d.comm);
  MPI_Comm_rank(d.comm, &d.rank);
  MPI_Comm_size(d.comm, &d.size);
  int coords[3];
  MPI_Cart_coords(d.comm, d.rank, 3, coords);
  for (int a = 0; a < 3; ++a) {
    d.coords[a] = coords[a];
    MPI_Cart_shift(d.comm, a, 1, &d.nbr[a][0], &d.nbr[a][1]);
    d.N[a] = cfg.cells[a];
    d.dx[a] = cfg.dx[a];
    d.c0[a] = int((long long)coords[a] * d.N[a] / dims[a]);
    d.c1[a] = int((long long)(coords[a] + 1) * d.N[a] / dims[a]);
    if (d.c1[a] - d.c0[a] < 3) fatal("fewer than 3 cells per rank along an axis; use fewer ranks");
    d.n0[a] = d.c0[a];
    d.n1[a] = d.c1[a] + (d.wall_hi(a) ? 1 : 0);
  }
  d.origin = cfg.origin;
  return d;
}

double allreduce_sum(double v, MPI_Comm comm) {
  double r;
  MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_SUM, comm);
  return r;
}
double allreduce_max(double v, MPI_Comm comm) {
  double r;
  MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MAX, comm);
  return r;
}
double allreduce_min(double v, MPI_Comm comm) {
  double r;
  MPI_Allreduce(&v, &r, 1, MPI_DOUBLE, MPI_MIN, comm);
  return r;
}
long long allreduce_sum(long long v, MPI_Comm comm) {
  long long r;
  MPI_Allreduce(&v, &r, 1, MPI_LONG_LONG, MPI_SUM, comm);
  return r;
}

}  // namespace picsp
