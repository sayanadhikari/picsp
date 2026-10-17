#pragma once
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>

#include <mpi.h>

namespace picsp {

inline int world_rank() {
  int r = 0, init = 0;
  MPI_Initialized(&init);
  if (init) MPI_Comm_rank(MPI_COMM_WORLD, &r);
  return r;
}

// printf-style logging on rank 0 only
template <typename... Args>
inline void log0(const char* fmt, Args&&... args) {
  if (world_rank() == 0) {
    if constexpr (sizeof...(Args) == 0)
      std::fputs(fmt, stdout);
    else
      std::printf(fmt, std::forward<Args>(args)...);
    std::fflush(stdout);
  }
}

template <typename... Args>
inline void warn0(const char* fmt, Args&&... args) {
  if (world_rank() == 0) {
    std::fputs("\033[33mWARNING\033[0m ", stdout);
    if constexpr (sizeof...(Args) == 0)
      std::fputs(fmt, stdout);
    else
      std::printf(fmt, std::forward<Args>(args)...);
    std::fflush(stdout);
  }
}

[[noreturn]] inline void fatal(const std::string& msg) {
  std::fprintf(stderr, "\033[31mERROR\033[0m [rank %d] %s\n", world_rank(), msg.c_str());
  std::fflush(stderr);
  int init = 0;
  MPI_Initialized(&init);
  if (init) MPI_Abort(MPI_COMM_WORLD, 1);
  std::exit(1);
}

#define PICSP_REQUIRE(cond, msg)            \
  do {                                      \
    if (!(cond)) ::picsp::fatal(msg);       \
  } while (0)

}  // namespace picsp
