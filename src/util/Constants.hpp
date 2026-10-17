#pragma once
// Physical constants (SI, CODATA 2018). All internal quantities in PICSP3D are SI.
namespace picsp::constants {
inline constexpr double pi = 3.14159265358979323846;
inline constexpr double twopi = 2.0 * pi;
inline constexpr double eps0 = 8.8541878128e-12;   // F/m
inline constexpr double qe = 1.602176634e-19;      // C (elementary charge)
inline constexpr double me = 9.1093837015e-31;     // kg
inline constexpr double amu = 1.66053906660e-27;   // kg
inline constexpr double kB = 1.380649e-23;         // J/K
inline constexpr double c = 299792458.0;           // m/s
inline constexpr double eV = qe;                   // J
inline constexpr double eV_to_K = qe / kB;         // 11604.518...
}  // namespace picsp::constants
