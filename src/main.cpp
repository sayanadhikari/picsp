// PICSP3D: parallel 3D3V electrostatic particle-in-cell code
// with embedded objects, background magnetic field and surface charging.
#include <cstdio>
#include <cstring>
#include <string>

#include <mpi.h>
#include <omp.h>

#include "Simulation.hpp"
#include "config/Config.hpp"
#include "util/Log.hpp"

using namespace picsp;

static void usage() {
  log0("usage: picsp3d <input.toml> [--restart <checkpoint_dir>] [--steps N] [--output DIR] [--check]\n"
       "  --restart DIR   continue from a checkpoint directory (same number of MPI ranks)\n"
       "  --steps N       override time.steps\n"
       "  --output DIR    override output.dir\n"
       "  --check         parse input, set up geometry and solver, then exit\n");
}

int main(int argc, char** argv) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  int rc = 0;
  {
    if (argc < 2 || !std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) {
      usage();
      MPI_Finalize();
      return argc < 2 ? 1 : 0;
    }
    std::string input = argv[1], restart, outdir;
    long long steps = -1;
    bool check = false;
    for (int i = 2; i < argc; ++i) {
      std::string a = argv[i];
      if (a == "--restart" && i + 1 < argc) restart = argv[++i];
      else if (a == "--steps" && i + 1 < argc) steps = std::stoll(argv[++i]);
      else if (a == "--output" && i + 1 < argc) outdir = argv[++i];
      else if (a == "--check") check = true;
      else {
        usage();
        fatal("unknown argument " + a);
      }
    }
    log0("PICSP3D %s  |  3D3V electrostatic PIC  |  MPI + OpenMP\n", PICSP_VERSION);
    Config cfg = Config::load(input);
    if (steps >= 0) cfg.steps = steps;
    if (!outdir.empty()) cfg.output.dir = outdir;
    cfg.restart_from = restart;
    cfg.print_summary();
    Simulation sim(cfg, MPI_COMM_WORLD);
    sim.initialise();
    if (!check) sim.run();
  }
  MPI_Finalize();
  return rc;
}
