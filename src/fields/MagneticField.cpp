#include "fields/MagneticField.hpp"

#include <cmath>

#include "util/Constants.hpp"

namespace picsp {

namespace {
constexpr double mu0 = 1.25663706212e-6;  // vacuum permeability [H/m]

// circular loop of radius a, current I, centred at c with unit normal n (right-hand rule)
Vec3 loop_field(const Vec3& c, const Vec3& n, double a, double I, const Vec3& p, double rmin) {
  Vec3 r = p - c;
  double z = dot(r, n);
  Vec3 rho_v = r - n * z;
  double rho = norm(rho_v);
  double ap2 = (a + rho) * (a + rho) + z * z;
  double am2 = std::max((a - rho) * (a - rho) + z * z, rmin * rmin);
  double k2 = std::min(4 * a * rho / ap2, 1.0 - 1e-12);
  double k = std::sqrt(k2);
  double K = std::comp_ellint_1(k), E = std::comp_ellint_2(k);
  double pref = mu0 * I / (2 * constants::pi * std::sqrt(ap2));
  double Bz = pref * (K + (a * a - rho * rho - z * z) / am2 * E);
  Vec3 B = n * Bz;
  if (rho > 1e-12 * a) {
    double Brho = pref * z / rho * (-K + (a * a + rho * rho + z * z) / am2 * E);
    B += rho_v * (Brho / rho);
  }
  return B;
}
}  // namespace

Vec3 magnet_field(const MagnetSource& m, const Vec3& p, double rmin) {
  switch (m.type) {
    case MagnetSource::Dipole: {
      Vec3 r = p - m.center;
      double rr = std::max(norm(r), rmin);
      Vec3 u = r / std::max(norm(r), 1e-300);
      if (norm(r) == 0) u = {0, 0, 1};
      return (u * (3 * dot(m.moment, u)) - m.moment) * (mu0 / (4 * constants::pi * rr * rr * rr));
    }
    case MagnetSource::Loop: return loop_field(m.center, m.axis, m.radius, m.current, p, rmin);
    case MagnetSource::Mirror:
      return loop_field(m.center - m.axis * (0.5 * m.separation), m.axis, m.radius, m.current, p, rmin) +
             loop_field(m.center + m.axis * (0.5 * m.separation), m.axis, m.radius, m.current, p, rmin);
    case MagnetSource::Solenoid: {
      int N = m.loops > 0 ? m.loops : 32;
      Vec3 B{0, 0, 0};
      for (int i = 0; i < N; ++i) {
        double s = (N == 1 ? 0.0 : (double(i) / (N - 1) - 0.5)) * m.length;
        B += loop_field(m.center + m.axis * s, m.axis, m.radius, m.current, p, rmin);
      }
      return B;
    }
    case MagnetSource::Wire: {
      Vec3 r = p - m.center;
      Vec3 perp = r - m.axis * dot(r, m.axis);
      double s = std::max(norm(perp), rmin);
      Vec3 phi = cross(m.axis, perp) / std::max(norm(perp), 1e-300);
      if (norm(perp) == 0) return {0, 0, 0};
      return phi * (mu0 * m.current / (2 * constants::pi * s));
    }
  }
  return {0, 0, 0};
}

Vec3 total_magnetic_field(const Config& cfg, const Vec3& p, double rmin) {
  Vec3 B = cfg.B0;
  for (auto& m : cfg.magnets) B += magnet_field(m, p, rmin);
  return B;
}

}  // namespace picsp
