#pragma once
// Fast, high-quality, splittable PRNG (xoshiro256**), one independent stream per
// (MPI rank, OpenMP thread). Satisfies UniformRandomBitGenerator so it can also be
// used with <random> distributions where convenience matters more than speed.
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

namespace picsp {

class Rng {
 public:
  using result_type = std::uint64_t;
  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

  explicit Rng(std::uint64_t seed = 0x9E3779B97F4A7C15ULL, std::uint64_t stream = 0) { reseed(seed, stream); }

  void reseed(std::uint64_t seed, std::uint64_t stream) {
    std::uint64_t z = seed ^ (0xD1B54A32D192ED03ULL * (stream + 1));
    for (auto& si : s_) si = splitmix64(z);
    has_spare_ = false;
  }

  inline result_type operator()() {
    const std::uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

  // Uniform in [0,1)
  inline double uniform() { return double((*this)() >> 11) * 0x1.0p-53; }
  // Uniform in (0,1]
  inline double uniform_open0() { return 1.0 - uniform(); }

  // Standard normal (polar Box-Muller with cached spare)
  inline double normal() {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    double u, v, s;
    do {
      u = 2.0 * uniform() - 1.0;
      v = 2.0 * uniform() - 1.0;
      s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    double f = std::sqrt(-2.0 * std::log(s) / s);
    spare_ = v * f;
    has_spare_ = true;
    return u * f;
  }

  // Full state for checkpointing: 4 x uint64 generator state + cached normal deviate
  static constexpr int kStateBytes = 48;
  void save(unsigned char* buf) const {
    std::memcpy(buf, s_, 32);
    std::memcpy(buf + 32, &spare_, 8);
    std::uint64_t f = has_spare_ ? 1 : 0;
    std::memcpy(buf + 40, &f, 8);
  }
  void load(const unsigned char* buf) {
    std::memcpy(s_, buf, 32);
    std::memcpy(&spare_, buf + 32, 8);
    std::uint64_t f;
    std::memcpy(&f, buf + 40, 8);
    has_spare_ = f != 0;
  }

 private:
  static inline std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  static std::uint64_t splitmix64(std::uint64_t& x) {
    std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }
  std::uint64_t s_[4];
  double spare_ = 0;
  bool has_spare_ = false;
};

// Radical inverse (van der Corput) for quasi-random "quiet" loading.
inline double radical_inverse(std::uint64_t i, unsigned base) {
  double inv = 1.0 / base, f = inv, r = 0.0;
  while (i > 0) {
    r += f * double(i % base);
    i /= base;
    f *= inv;
  }
  return r;
}

}  // namespace picsp
