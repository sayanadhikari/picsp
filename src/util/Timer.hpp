#pragma once
#include <chrono>
#include <map>
#include <string>
#include <vector>

#include <mpi.h>

namespace picsp {

// Accumulating named wall-clock timers; report() prints min/avg/max over ranks.
class Timers {
 public:
  void start(const std::string& name) { t0_[name] = clock::now(); }
  void stop(const std::string& name) {
    auto it = t0_.find(name);
    if (it == t0_.end()) return;
    acc_[name] += std::chrono::duration<double>(clock::now() - it->second).count();
  }
  double get(const std::string& name) const {
    auto it = acc_.find(name);
    return it == acc_.end() ? 0.0 : it->second;
  }
  const std::map<std::string, double>& all() const { return acc_; }

  // Collective. Prints a table on rank 0.
  void report(MPI_Comm comm, double total) const;

 private:
  using clock = std::chrono::steady_clock;
  std::map<std::string, clock::time_point> t0_;
  std::map<std::string, double> acc_;
};

class ScopedTimer {
 public:
  ScopedTimer(Timers& t, const char* n) : t_(t), n_(n) { t_.start(n_); }
  ~ScopedTimer() { t_.stop(n_); }

 private:
  Timers& t_;
  std::string n_;
};

}  // namespace picsp
