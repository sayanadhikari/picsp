#pragma once
// Input configuration for PICSP3D. Parsed from a TOML file (see docs/input_reference.md).
// Everything is converted to SI once here; the rest of the code is unit-free SI.
#include <array>
#include <optional>
#include <string>
#include <vector>

#include "util/Vec3.hpp"

namespace picsp {

enum class FieldBC { Periodic, Dirichlet, Neumann };
enum class ParticleBC { Periodic, Absorb, Reflect, Open };  // Open = absorb + inject ambient flux

struct FaceConfig {
  FieldBC field = FieldBC::Periodic;
  double value = 0.0;  // Dirichlet potential [V]
  ParticleBC particle = ParticleBC::Periodic;
};

enum class VelocityDist { Maxwellian, BiMaxwellian, Kappa, Ring, Cold };
enum class LoadKind { Random, Quiet };

struct Perturbation {
  int axis = 0;          // 0,1,2
  double amplitude = 0;  // relative density perturbation epsilon
  int mode = 1;          // wavenumber index over box length
  double velocity_amplitude = 0;  // optional velocity perturbation [m/s] along axis
};

struct SpeciesConfig {
  std::string name;
  double charge = 0;   // [C]
  double mass = 0;     // [kg]
  double density = 0;  // [m^-3] initial / ambient
  double temperature = 0;    // [eV] isotropic temperature
  double t_par = -1, t_perp = -1;  // [eV] bi-Maxwellian (relative to `aniso_axis`, default B direction)
  std::optional<Vec3> aniso_axis;
  Vec3 drift{0, 0, 0};       // [m/s]
  VelocityDist dist = VelocityDist::Maxwellian;
  double kappa = 4.0;
  double ring_speed = 0;     // [m/s] perpendicular ring speed
  double ppc = 8;            // macro particles per cell at density `density`
  double weight = 0;         // optional explicit macro weight (overrides ppc)
  LoadKind load = LoadKind::Random;
  bool load_initial = true;  // fill domain at t=0
  bool inject = true;        // inject at open faces
  int subcycle = 1;          // push every `subcycle` steps (heavy ions)
  std::optional<std::array<Vec3, 2>> region;  // optional loading box [lo, hi]
  Perturbation perturbation;
};

enum class ShapeKind { Sphere, Box, Cylinder, Cone, Torus, Ellipsoid, Plate, HalfSpace, STL };
enum class CsgOp { Union, Subtract, Intersect };

struct PartConfig {
  ShapeKind shape = ShapeKind::Sphere;
  CsgOp op = CsgOp::Union;
  Vec3 center{0, 0, 0};
  double radius = 0, radius2 = 0, length = 0;  // meaning depends on shape
  Vec3 size{0, 0, 0};                          // box/plate full extents, ellipsoid semi-axes
  Vec3 axis{0, 0, 1};                          // cylinder/cone/torus/half-space normal
  Vec3 rotate_deg{0, 0, 0};                    // extra rotation of this part
  double rounding = 0;                         // rounded edges for boxes
  std::string file;                            // STL path
  double scale = 1.0;                          // STL scale factor (e.g. 1e-3 for mm)
};

struct SeeConfig {
  enum Model { None, Constant, Sternglass } model = None;
  double yield = 0;        // constant yield
  double delta_max = 2.0;  // Sternglass peak yield
  double e_max = 400.0;    // [eV] energy of peak yield
  double temperature = 2.0;  // [eV] emitted secondary temperature
  double ion_yield = 0.0;    // ion-induced electron emission yield (constant)
};

struct PhotoConfig {
  bool enabled = false;
  double current_density = 0;  // [A/m^2] saturation photoemission current density at normal incidence
  double temperature = 1.5;    // [eV]
};

enum class BiasMode { Floating, Fixed };

struct ObjectConfig {
  std::string name;
  std::vector<PartConfig> parts;
  Vec3 translate{0, 0, 0};
  Vec3 rotate_deg{0, 0, 0};
  BiasMode bias = BiasMode::Floating;
  double potential = 0;      // [V] fixed bias, or initial potential for floating
  double ramp = 0;           // [V/s] linear bias ramp for fixed objects (probe sweeps)
  double ext_capacitance = 0;  // [F] extra capacitance to ground (floating)
  double ext_current = 0;      // [A] constant external current into object (floating)
  std::string emit_species;  // species name used for SEE/photo electrons
  SeeConfig see;
  PhotoConfig photo;
  int group = -1;            // objects with the same group >=0 are electrically connected
};

struct CollisionProcess {
  std::string species;
  enum Type { Elastic, ChargeExchange } type = Elastic;
  double cross_section = 0;  // [m^2], constant
};

struct NeutralConfig {
  bool enabled = false;
  double density = 0;       // [m^-3]
  double temperature = 0.025;  // [eV]
  double mass = 0;          // [kg]
  Vec3 drift{0, 0, 0};
  std::vector<CollisionProcess> processes;
};

struct OutputConfig {
  std::string dir = "output";
  int fields_every = 0;
  int particles_every = 0;
  int particle_stride = 100;
  int surface_every = 0;
  int history_every = 10;
  int checkpoint_every = 0;
  int status_every = 50;
  bool moments = true;     // per-species n, flux, T
  bool float32 = true;     // write field data as float32 to save space
};

struct SolverConfig {
  std::string type = "mg-pcg";  // mg-pcg | mg-bicgstab | sor | jacobi-pcg
  double tol = 1e-8;            // relative residual tolerance
  int max_iter = 500;
  bool cut_cell = true;         // Shortley-Weller boundary treatment near objects
  int mg_max_levels = 12;
  int mg_pre = 2, mg_post = 2;  // smoothing sweeps
  int mg_coarse_sweeps = 40;
  double sor_omega = 1.8;
};

struct Config {
  std::string input_file;
  std::string title = "picsp3d";
  unsigned long long seed = 12345;

  // grid
  std::array<int, 3> cells{32, 32, 32};
  std::array<double, 3> dx{0, 0, 0};
  Vec3 origin{0, 0, 0};
  std::array<int, 3> decomposition{0, 0, 0};  // 0 = automatic

  // time
  double dt = 0;
  long long steps = 0;

  // boundaries: index 2*axis + (0=lo, 1=hi)
  std::array<FaceConfig, 6> faces;

  // external fields
  Vec3 B0{0, 0, 0};
  Vec3 E0{0, 0, 0};
  Vec3 sun_dir{1, 0, 0};  // unit vector pointing toward the sun

  std::vector<SpeciesConfig> species;
  std::vector<ObjectConfig> objects;
  NeutralConfig neutrals;
  SolverConfig solver;
  OutputConfig output;

  int shape_order = 1;  // 1 = CIC (linear), 2 = TSC (quadratic)
  int sort_every = 20;
  int surface_refine = 2;  // sub-sampling of cells for the object surface triangulation
  std::string restart_from;  // checkpoint directory

  // derived reference quantities (from reference electron-like species)
  double ref_wpe = 0, ref_debye = 0;

  static Config load(const std::string& path);
  void print_summary() const;
};

const char* to_string(FieldBC);
const char* to_string(ParticleBC);

}  // namespace picsp
