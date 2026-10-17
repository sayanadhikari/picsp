#pragma once
#include <array>
#include <cmath>

namespace picsp {

struct Vec3 {
  double x = 0, y = 0, z = 0;
  constexpr Vec3() = default;
  constexpr Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
  constexpr double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
  constexpr double& operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
  constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3 operator-() const { return {-x, -y, -z}; }
  constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  constexpr Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
  constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
  constexpr Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  constexpr Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};
constexpr Vec3 operator*(double s, const Vec3& v) { return v * s; }
constexpr double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalized(const Vec3& a) {
  double n = norm(a);
  return n > 0 ? a / n : Vec3{0, 0, 0};
}
inline Vec3 vmin(const Vec3& a, const Vec3& b) { return {std::fmin(a.x, b.x), std::fmin(a.y, b.y), std::fmin(a.z, b.z)}; }
inline Vec3 vmax(const Vec3& a, const Vec3& b) { return {std::fmax(a.x, b.x), std::fmax(a.y, b.y), std::fmax(a.z, b.z)}; }

// Build an orthonormal basis (e1, e2) perpendicular to unit vector n.
inline void orthonormal_basis(const Vec3& n, Vec3& e1, Vec3& e2) {
  Vec3 a = std::fabs(n.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
  e1 = normalized(cross(n, a));
  e2 = cross(n, e1);
}

// 3x3 rotation matrix (row-major)
struct Mat3 {
  std::array<double, 9> m{1, 0, 0, 0, 1, 0, 0, 0, 1};
  Vec3 operator*(const Vec3& v) const {
    return {m[0] * v.x + m[1] * v.y + m[2] * v.z, m[3] * v.x + m[4] * v.y + m[5] * v.z,
            m[6] * v.x + m[7] * v.y + m[8] * v.z};
  }
  Mat3 operator*(const Mat3& o) const {
    Mat3 r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double s = 0;
        for (int k = 0; k < 3; ++k) s += m[3 * i + k] * o.m[3 * k + j];
        r.m[3 * i + j] = s;
      }
    return r;
  }
  Mat3 transposed() const {
    Mat3 r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) r.m[3 * i + j] = m[3 * j + i];
    return r;
  }
  // Rotation from Euler angles in degrees, applied as Rz * Ry * Rx (x first).
  static Mat3 from_euler_deg(const Vec3& deg) {
    const double d2r = 3.14159265358979323846 / 180.0;
    double a = deg.x * d2r, b = deg.y * d2r, c = deg.z * d2r;
    Mat3 rx, ry, rz;
    rx.m = {1, 0, 0, 0, std::cos(a), -std::sin(a), 0, std::sin(a), std::cos(a)};
    ry.m = {std::cos(b), 0, std::sin(b), 0, 1, 0, -std::sin(b), 0, std::cos(b)};
    rz.m = {std::cos(c), -std::sin(c), 0, std::sin(c), std::cos(c), 0, 0, 0, 1};
    return rz * ry * rx;
  }
};

}  // namespace picsp
