#pragma once
// Static external magnetic field: uniform B0 plus analytic sources
//   dipole   B = mu0/(4 pi) [3 r^(m.r^) - m] / r^3
//   loop     circular current loop (exact, complete elliptic integrals)
//   mirror   two coaxial loops (magnetic bottle)
//   solenoid N coaxial loops spread over a length
//   wire     infinite straight line current
// Evaluated on the node grid once at start-up and interpolated to particles with the
// same shape functions as E.
#include <vector>

#include "config/Config.hpp"
#include "util/Vec3.hpp"

namespace picsp {

// field of one source at p; distances below rmin are clamped (removes singularities)
Vec3 magnet_field(const MagnetSource& m, const Vec3& p, double rmin);
// B0 + all sources
Vec3 total_magnetic_field(const Config& cfg, const Vec3& p, double rmin);

}  // namespace picsp
