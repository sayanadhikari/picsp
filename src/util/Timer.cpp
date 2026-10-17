#include "util/Timer.hpp"

#include <cstdio>

namespace picsp {

void Timers::report(MPI_Comm comm, double total) const {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  int size = 1;
  MPI_Comm_size(comm, &size);
  // All ranks have the same timer names (same code path), but be defensive:
  // use rank 0's list and broadcast names.
  std::vector<std::string> names;
  for (auto& kv : acc_) names.push_back(kv.first);
  int n = int(names.size());
  MPI_Bcast(&n, 1, MPI_INT, 0, comm);
  if (rank == 0) {
    std::printf("\n%-24s %12s %12s %12s %8s\n", "timer", "min [s]", "avg [s]", "max [s]", "% total");
    std::printf("-----------------------------------------------------------------------\n");
  }
  for (int i = 0; i < n; ++i) {
    char buf[64] = {0};
    if (rank == 0) std::snprintf(buf, sizeof(buf), "%s", names[i].c_str());
    MPI_Bcast(buf, 64, MPI_CHAR, 0, comm);
    double v = get(buf), vmin, vmax, vsum;
    MPI_Allreduce(&v, &vmin, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&v, &vmax, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&v, &vsum, 1, MPI_DOUBLE, MPI_SUM, comm);
    if (rank == 0)
      std::printf("%-24s %12.4f %12.4f %12.4f %7.1f%%\n", buf, vmin, vsum / size, vmax,
                  total > 0 ? 100.0 * vmax / total : 0.0);
  }
  if (rank == 0) std::printf("-----------------------------------------------------------------------\n");
}

}  // namespace picsp
