#pragma once
// Particle shape functions (B-splines) on the node grid.
//   Order 1: cloud-in-cell (linear), 2 nodes per axis
//   Order 2: triangular-shaped cloud (quadratic), 3 nodes per axis
namespace picsp {

template <int Order>
struct Shape;

template <>
struct Shape<1> {
  static constexpr int n = 2;
  static inline void weights(double l, int& i0, double* w) {
    int i = int(l);
    double f = l - i;
    i0 = i;
    w[0] = 1.0 - f;
    w[1] = f;
  }
};

template <>
struct Shape<2> {
  static constexpr int n = 3;
  static inline void weights(double l, int& i0, double* w) {
    int i = int(l + 0.5);
    double d = l - i;
    i0 = i - 1;
    w[0] = 0.5 * (0.5 - d) * (0.5 - d);
    w[1] = 0.75 - d * d;
    w[2] = 0.5 * (0.5 + d) * (0.5 + d);
  }
};

}  // namespace picsp
