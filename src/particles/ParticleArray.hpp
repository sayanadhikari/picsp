#pragma once
// Structure-of-arrays particle storage (cache- and SIMD-friendly).
#include <cstddef>
#include <vector>

#include "util/Vec3.hpp"

namespace picsp {

struct ParticleArray {
  std::vector<double> x, y, z, vx, vy, vz, w;

  std::size_t size() const { return x.size(); }
  void reserve(std::size_t n) {
    for (auto* a : arrays()) a->reserve(n);
  }
  void resize(std::size_t n) {
    for (auto* a : arrays()) a->resize(n);
  }
  void clear() {
    for (auto* a : arrays()) a->clear();
  }
  void push(const Vec3& p, const Vec3& v, double wt) {
    x.push_back(p.x); y.push_back(p.y); z.push_back(p.z);
    vx.push_back(v.x); vy.push_back(v.y); vz.push_back(v.z);
    w.push_back(wt);
  }
  // O(1) removal (order not preserved)
  void remove(std::size_t i) {
    std::size_t last = size() - 1;
    if (i != last) {
      x[i] = x[last]; y[i] = y[last]; z[i] = z[last];
      vx[i] = vx[last]; vy[i] = vy[last]; vz[i] = vz[last];
      w[i] = w[last];
    }
    x.pop_back(); y.pop_back(); z.pop_back();
    vx.pop_back(); vy.pop_back(); vz.pop_back();
    w.pop_back();
  }
  Vec3 pos(std::size_t i) const { return {x[i], y[i], z[i]}; }
  Vec3 vel(std::size_t i) const { return {vx[i], vy[i], vz[i]}; }

  // pack particle i into 7 doubles
  void pack(std::size_t i, double* buf) const {
    buf[0] = x[i]; buf[1] = y[i]; buf[2] = z[i];
    buf[3] = vx[i]; buf[4] = vy[i]; buf[5] = vz[i]; buf[6] = w[i];
  }
  void push_packed(const double* b) {
    x.push_back(b[0]); y.push_back(b[1]); z.push_back(b[2]);
    vx.push_back(b[3]); vy.push_back(b[4]); vz.push_back(b[5]); w.push_back(b[6]);
  }
  static constexpr int kPacked = 7;

  std::vector<std::vector<double>*> arrays() { return {&x, &y, &z, &vx, &vy, &vz, &w}; }
};

}  // namespace picsp
