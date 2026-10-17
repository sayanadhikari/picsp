#include "config/Config.hpp"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>

#define TOML_EXCEPTIONS 1
#include "toml.hpp"
#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

namespace {

using toml::node_view;
using toml::table;


template <typename NV>
double get_num(const NV& n, const std::string& where, double def) {
  if (!n) return def;
  if (auto v = n.template value<double>()) return *v;
  fatal("expected a number at '" + where + "'");
}
template <typename T>
std::optional<double> opt_num(const T& tbl, const char* key) {
  auto n = tbl[key];
  if (!n) return std::nullopt;
  if (auto v = n.template value<double>()) return *v;
  fatal(std::string("expected a number for key '") + key + "'");
}

template <typename NV>
Vec3 get_vec3(const NV& n, const std::string& where, Vec3 def) {
  if (!n) return def;
  auto arr = n.as_array();
  if (!arr || arr->size() != 3) fatal("expected an array of 3 numbers at '" + where + "'");
  Vec3 v;
  for (int i = 0; i < 3; ++i) {
    auto d = (*arr)[i].template value<double>();
    if (!d) fatal("expected numeric entries at '" + where + "'");
    v[i] = *d;
  }
  return v;
}

template <typename NV>
std::string get_str(const NV& n, const std::string& def) {
  if (!n) return def;
  if (auto v = n.template value<std::string>()) return *v;
  fatal("expected a string");
}

template <typename NV>
bool get_bool(const NV& n, bool def) {
  if (!n) return def;
  if (auto v = n.template value<bool>()) return *v;
  fatal("expected a boolean");
}

FieldBC parse_field_bc(const std::string& s) {
  if (s == "periodic") return FieldBC::Periodic;
  if (s == "dirichlet") return FieldBC::Dirichlet;
  if (s == "neumann") return FieldBC::Neumann;
  fatal("unknown field boundary '" + s + "' (periodic|dirichlet|neumann)");
}
ParticleBC parse_particle_bc(const std::string& s) {
  if (s == "periodic") return ParticleBC::Periodic;
  if (s == "absorb" || s == "absorbing" || s == "outflow") return ParticleBC::Absorb;
  if (s == "reflect" || s == "reflecting" || s == "specular") return ParticleBC::Reflect;
  if (s == "open" || s == "inject" || s == "inflow") return ParticleBC::Open;
  fatal("unknown particle boundary '" + s + "' (periodic|absorb|reflect|open)");
}

ShapeKind parse_shape(const std::string& s) {
  if (s == "sphere") return ShapeKind::Sphere;
  if (s == "box" || s == "cuboid") return ShapeKind::Box;
  if (s == "cylinder") return ShapeKind::Cylinder;
  if (s == "cone") return ShapeKind::Cone;
  if (s == "torus") return ShapeKind::Torus;
  if (s == "ellipsoid") return ShapeKind::Ellipsoid;
  if (s == "plate" || s == "panel") return ShapeKind::Plate;
  if (s == "halfspace" || s == "plane") return ShapeKind::HalfSpace;
  if (s == "stl" || s == "mesh") return ShapeKind::STL;
  fatal("unknown shape '" + s + "'");
}

void parse_face(const table& t, FaceConfig& f, const std::string& where) {
  if (auto n = t["field"]) f.field = parse_field_bc(get_str(n, ""));
  if (auto n = t["value"]) f.value = get_num(n, where + ".value", 0.0);
  if (auto n = t["particle"]) f.particle = parse_particle_bc(get_str(n, ""));
}

}  // namespace

const char* to_string(FieldBC b) {
  switch (b) {
    case FieldBC::Periodic: return "periodic";
    case FieldBC::Dirichlet: return "dirichlet";
    case FieldBC::Neumann: return "neumann";
  }
  return "?";
}
const char* to_string(ParticleBC b) {
  switch (b) {
    case ParticleBC::Periodic: return "periodic";
    case ParticleBC::Absorb: return "absorb";
    case ParticleBC::Reflect: return "reflect";
    case ParticleBC::Open: return "open";
  }
  return "?";
}

Config Config::load(const std::string& path) {
  using namespace constants;
  Config c;
  c.input_file = path;
  table root;
  try {
    root = toml::parse_file(path);
  } catch (const toml::parse_error& e) {
    std::ostringstream os;
    os << "failed to parse '" << path << "': " << e.description() << " at line " << e.source().begin.line;
    fatal(os.str());
  }
  const std::filesystem::path base = std::filesystem::path(path).parent_path();

  c.title = get_str(root["title"], c.title);
  c.seed = (unsigned long long)get_num(root["seed"], "seed", double(c.seed));

  // ---------------- species (parsed first: needed for normalised inputs) -------------
  if (auto arr = root["species"].as_array()) {
    for (auto& el : *arr) {
      auto* t = el.as_table();
      if (!t) fatal("[[species]] entries must be tables");
      SpeciesConfig s;
      const table& st = *t;
      s.name = get_str(st["name"], "species" + std::to_string(c.species.size()));
      std::string w = "species." + s.name;
      s.charge = get_num(st["charge"], w + ".charge", 0.0) * qe;  // in units of e
      if (auto v = opt_num(st, "mass_kg")) s.mass = *v;
      else if (auto v2 = opt_num(st, "mass_amu")) s.mass = *v2 * amu;
      else if (auto v3 = opt_num(st, "mass_me")) s.mass = *v3 * me;
      else if (auto v4 = opt_num(st, "mass")) s.mass = *v4 * amu;  // default unit: amu
      else fatal("species '" + s.name + "' needs mass_amu, mass_me or mass_kg");
      s.density = get_num(st["density"], w + ".density", 0.0);
      s.temperature = get_num(st["temperature"], w + ".temperature", 0.0);
      s.t_par = get_num(st["t_par"], w + ".t_par", -1.0);
      s.t_perp = get_num(st["t_perp"], w + ".t_perp", -1.0);
      if (st["aniso_axis"]) s.aniso_axis = normalized(get_vec3(st["aniso_axis"], w + ".aniso_axis", {0, 0, 1}));
      s.drift = get_vec3(st["drift"], w + ".drift", {0, 0, 0});
      std::string d = get_str(st["distribution"], "maxwellian");
      if (d == "maxwellian") s.dist = VelocityDist::Maxwellian;
      else if (d == "bimaxwellian" || d == "bi-maxwellian") s.dist = VelocityDist::BiMaxwellian;
      else if (d == "kappa") s.dist = VelocityDist::Kappa;
      else if (d == "ring") s.dist = VelocityDist::Ring;
      else if (d == "cold" || d == "beam") s.dist = VelocityDist::Cold;
      else fatal("unknown distribution '" + d + "' for species " + s.name);
      s.kappa = get_num(st["kappa"], w + ".kappa", 4.0);
      if (s.dist == VelocityDist::Kappa && s.kappa <= 1.5) fatal("kappa must be > 1.5");
      s.ring_speed = get_num(st["ring_speed"], w + ".ring_speed", 0.0);
      s.ppc = get_num(st["ppc"], w + ".ppc", 8.0);
      s.weight = get_num(st["weight"], w + ".weight", 0.0);
      std::string l = get_str(st["load"], "random");
      s.load = (l == "quiet" || l == "halton") ? LoadKind::Quiet : LoadKind::Random;
      s.load_initial = get_bool(st["load_initial"], s.density > 0);
      s.inject = get_bool(st["inject"], s.density > 0);
      s.subcycle = (int)get_num(st["subcycle"], w + ".subcycle", 1);
      if (s.subcycle < 1) s.subcycle = 1;
      if (auto r = st["region"].as_table()) {
        s.region = std::array<Vec3, 2>{get_vec3((*r)["min"], w + ".region.min", {}),
                                       get_vec3((*r)["max"], w + ".region.max", {})};
      }
      if (auto p = st["perturbation"].as_table()) {
        std::string ax = get_str((*p)["axis"], "x");
        s.perturbation.axis = ax == "x" ? 0 : (ax == "y" ? 1 : 2);
        s.perturbation.amplitude = get_num((*p)["amplitude"], w + ".perturbation.amplitude", 0.0);
        s.perturbation.mode = (int)get_num((*p)["mode"], w + ".perturbation.mode", 1);
        s.perturbation.velocity_amplitude = get_num((*p)["velocity_amplitude"], w, 0.0);
      }
      c.species.push_back(s);
    }
  }
  if (c.species.empty()) fatal("no [[species]] defined");

  // reference species: named, or first negatively charged species with n,T > 0
  {
    std::string refname = get_str(root["reference_species"], "");
    const SpeciesConfig* ref = nullptr;
    for (auto& s : c.species)
      if ((!refname.empty() && s.name == refname) ||
          (refname.empty() && !ref && s.charge < 0 && s.density > 0 && s.temperature > 0))
        ref = &s;
    if (!ref)
      for (auto& s : c.species)
        if (!ref && s.density > 0) ref = &s;
    if (ref && ref->density > 0) {
      c.ref_wpe = std::sqrt(ref->density * ref->charge * ref->charge / (eps0 * ref->mass));
      double T = ref->temperature > 0 ? ref->temperature : std::max(ref->t_par, ref->t_perp);
      if (T > 0) c.ref_debye = std::sqrt(eps0 * T * qe / (ref->density * ref->charge * ref->charge));
    }
  }

  // ---------------- grid ---------------------------------------------------------------
  if (auto g = root["grid"].as_table()) {
    const table& gt = *g;
    Vec3 cells = get_vec3(gt["cells"], "grid.cells", {32, 32, 32});
    for (int i = 0; i < 3; ++i) c.cells[i] = int(cells[i]);
    if (gt["dx"]) {
      if (gt["dx"].is_array()) {
        Vec3 d = get_vec3(gt["dx"], "grid.dx", {});
        c.dx = {d.x, d.y, d.z};
      } else {
        double d = get_num(gt["dx"], "grid.dx", 0);
        c.dx = {d, d, d};
      }
    } else if (gt["dx_debye"]) {
      if (c.ref_debye <= 0) fatal("grid.dx_debye requires a reference species with density and temperature");
      double d = get_num(gt["dx_debye"], "grid.dx_debye", 0) * c.ref_debye;
      c.dx = {d, d, d};
    } else if (gt["length"]) {
      Vec3 L = get_vec3(gt["length"], "grid.length", {});
      for (int i = 0; i < 3; ++i) c.dx[i] = L[i] / c.cells[i];
    } else {
      fatal("grid needs one of dx, dx_debye or length");
    }
    if (get_bool(gt["centered"], false))
      c.origin = Vec3{-0.5 * c.cells[0] * c.dx[0], -0.5 * c.cells[1] * c.dx[1], -0.5 * c.cells[2] * c.dx[2]};
    c.origin = get_vec3(gt["origin"], "grid.origin", c.origin);
  } else {
    fatal("missing [grid] section");
  }

  // ---------------- time ---------------------------------------------------------------
  if (auto t = root["time"].as_table()) {
    if ((*t)["dt"]) c.dt = get_num((*t)["dt"], "time.dt", 0);
    else if ((*t)["dt_wpe"]) {
      if (c.ref_wpe <= 0) fatal("time.dt_wpe requires a reference species with density");
      c.dt = get_num((*t)["dt_wpe"], "time.dt_wpe", 0) / c.ref_wpe;
    } else
      fatal("time needs dt or dt_wpe");
    if ((*t)["steps"]) c.steps = (long long)get_num((*t)["steps"], "time.steps", 0);
    else if ((*t)["duration"]) c.steps = (long long)std::ceil(get_num((*t)["duration"], "time.duration", 0) / c.dt);
    else if ((*t)["duration_wpe"]) c.steps = (long long)std::ceil(get_num((*t)["duration_wpe"], "", 0) / c.ref_wpe / c.dt);
  } else {
    fatal("missing [time] section");
  }

  // ---------------- boundaries ---------------------------------------------------------
  if (auto b = root["boundaries"].as_table()) {
    FaceConfig def;
    if (auto d = (*b)["default"].as_table()) parse_face(*d, def, "boundaries.default");
    for (auto& f : c.faces) f = def;
    const char* names[6] = {"x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi"};
    const char* axes[3] = {"x", "y", "z"};
    for (int a = 0; a < 3; ++a)
      if (auto ft = (*b)[axes[a]].as_table()) {  // both faces of an axis
        parse_face(*ft, c.faces[2 * a], std::string("boundaries.") + axes[a]);
        parse_face(*ft, c.faces[2 * a + 1], std::string("boundaries.") + axes[a]);
      }
    for (int i = 0; i < 6; ++i)
      if (auto ft = (*b)[names[i]].as_table()) parse_face(*ft, c.faces[i], std::string("boundaries.") + names[i]);
  }
  for (int a = 0; a < 3; ++a) {
    bool p0 = c.faces[2 * a].field == FieldBC::Periodic, p1 = c.faces[2 * a + 1].field == FieldBC::Periodic;
    bool q0 = c.faces[2 * a].particle == ParticleBC::Periodic, q1 = c.faces[2 * a + 1].particle == ParticleBC::Periodic;
    if (p0 != p1 || q0 != q1 || p0 != q0)
      fatal("periodic boundaries must be set on both faces of an axis, for both fields and particles");
  }

  // ---------------- fields -------------------------------------------------------------
  if (auto f = root["fields"].as_table()) {
    c.B0 = get_vec3((*f)["B0"], "fields.B0", {0, 0, 0});
    c.E0 = get_vec3((*f)["E0"], "fields.E0", {0, 0, 0});
    // motional electric field in the frame of an object moving with velocity v through B: E = -v x B
    if ((*f)["frame_velocity"]) {
      Vec3 v = get_vec3((*f)["frame_velocity"], "fields.frame_velocity", {});
      c.E0 += -1.0 * cross(v, c.B0);
    }
    c.electrostatic = get_bool((*f)["self_consistent"], true);
    if (auto arr = (*f)["magnet"].as_array()) {
      for (auto& el : *arr) {
        const table& mt = *el.as_table();
        MagnetSource m;
        std::string ty = get_str(mt["type"], "dipole");
        if (ty == "dipole") m.type = MagnetSource::Dipole;
        else if (ty == "loop" || ty == "coil") m.type = MagnetSource::Loop;
        else if (ty == "wire" || ty == "line") m.type = MagnetSource::Wire;
        else if (ty == "solenoid") m.type = MagnetSource::Solenoid;
        else if (ty == "mirror") m.type = MagnetSource::Mirror;
        else fatal("unknown magnet type '" + ty + "' (dipole|loop|wire|solenoid|mirror)");
        m.center = get_vec3(mt["center"], "fields.magnet.center", {0, 0, 0});
        m.axis = normalized(get_vec3(mt["axis"], "fields.magnet.axis", {0, 0, 1}));
        m.moment = get_vec3(mt["moment"], "fields.magnet.moment", {0, 0, 0});
        m.radius = get_num(mt["radius"], "fields.magnet.radius", 0.0);
        m.current = get_num(mt["current"], "fields.magnet.current", 0.0);
        m.length = get_num(mt["length"], "fields.magnet.length", 0.0);
        m.separation = get_num(mt["separation"], "fields.magnet.separation", 0.0);
        m.loops = (int)get_num(mt["loops"], "fields.magnet.loops", 0.0);
        if ((m.type == MagnetSource::Loop || m.type == MagnetSource::Mirror || m.type == MagnetSource::Solenoid) &&
            m.radius <= 0)
          fatal("loop / mirror / solenoid magnets need radius > 0");
        if (m.type == MagnetSource::Mirror && m.separation <= 0) fatal("mirror magnet needs separation > 0");
        if (m.type == MagnetSource::Solenoid && m.length <= 0) fatal("solenoid magnet needs length > 0");
        c.magnets.push_back(m);
      }
    }
  }
  if (auto s = root["sun"].as_table()) c.sun_dir = normalized(get_vec3((*s)["direction"], "sun.direction", {1, 0, 0}));

  // ---------------- objects ------------------------------------------------------------
  if (auto arr = root["objects"].as_array()) {
    for (auto& el : *arr) {
      const table& ot = *el.as_table();
      ObjectConfig o;
      o.name = get_str(ot["name"], "object" + std::to_string(c.objects.size()));
      std::string w = "objects." + o.name;
      std::string bias = get_str(ot["bias"], "floating");
      if (bias == "floating") o.bias = BiasMode::Floating;
      else if (bias == "fixed" || bias == "biased" || bias == "grounded") o.bias = BiasMode::Fixed;
      else fatal("objects.bias must be floating|fixed|grounded");
      o.potential = get_num(ot["potential"], w + ".potential", 0.0);
      o.ramp = get_num(ot["ramp"], w + ".ramp", 0.0);
      o.ext_capacitance = get_num(ot["ext_capacitance"], w + ".ext_capacitance", 0.0);
      o.ext_current = get_num(ot["ext_current"], w + ".ext_current", 0.0);
      o.group = (int)get_num(ot["group"], w + ".group", -1);
      o.translate = get_vec3(ot["translate"], w + ".translate", {0, 0, 0});
      o.rotate_deg = get_vec3(ot["rotate"], w + ".rotate", {0, 0, 0});
      o.emit_species = get_str(ot["emit_species"], "");
      if (auto see = ot["see"].as_table()) {
        std::string m = get_str((*see)["model"], "constant");
        o.see.model = m == "sternglass" ? SeeConfig::Sternglass : (m == "none" ? SeeConfig::None : SeeConfig::Constant);
        o.see.yield = get_num((*see)["yield"], w + ".see.yield", 0.0);
        o.see.delta_max = get_num((*see)["delta_max"], w + ".see.delta_max", 2.0);
        o.see.e_max = get_num((*see)["e_max"], w + ".see.e_max", 400.0);
        o.see.temperature = get_num((*see)["temperature"], w + ".see.temperature", 2.0);
        o.see.ion_yield = get_num((*see)["ion_yield"], w + ".see.ion_yield", 0.0);
      }
      if (auto ph = ot["photoemission"].as_table()) {
        o.photo.enabled = get_bool((*ph)["enabled"], true);
        o.photo.current_density = get_num((*ph)["current_density"], w + ".photoemission.current_density", 10e-6);
        o.photo.temperature = get_num((*ph)["temperature"], w + ".photoemission.temperature", 1.5);
      }
      auto parse_part = [&](const table& pt, PartConfig& p) {
        p.shape = parse_shape(get_str(pt["shape"], "sphere"));
        std::string op = get_str(pt["op"], "union");
        p.op = op == "subtract" ? CsgOp::Subtract : (op == "intersect" ? CsgOp::Intersect : CsgOp::Union);
        p.center = get_vec3(pt["center"], w + ".center", {0, 0, 0});
        p.radius = get_num(pt["radius"], w + ".radius", 0.0);
        p.radius2 = get_num(pt["radius2"], w + ".radius2", 0.0);
        p.length = get_num(pt["length"], w + ".length", 0.0);
        p.size = get_vec3(pt["size"], w + ".size", {0, 0, 0});
        p.axis = normalized(get_vec3(pt["axis"], w + ".axis", {0, 0, 1}));
        p.rotate_deg = get_vec3(pt["rotate"], w + ".rotate", {0, 0, 0});
        p.rounding = get_num(pt["rounding"], w + ".rounding", 0.0);
        p.scale = get_num(pt["scale"], w + ".scale", 1.0);
        p.file = get_str(pt["file"], "");
        if (!p.file.empty() && std::filesystem::path(p.file).is_relative()) p.file = (base / p.file).string();
        if (p.shape == ShapeKind::STL && p.file.empty()) fatal(w + ": STL part needs 'file'");
      };
      if (auto parts = ot["parts"].as_array()) {
        for (auto& pe : *parts) {
          PartConfig p;
          parse_part(*pe.as_table(), p);
          o.parts.push_back(p);
        }
      } else {
        PartConfig p;  // single-part shorthand: shape keys directly on the object
        parse_part(ot, p);
        o.parts.push_back(p);
      }
      c.objects.push_back(o);
    }
  }

  // ---------------- collisions ---------------------------------------------------------
  if (auto col = root["collisions"].as_table()) {
    auto& n = c.neutrals;
    n.enabled = get_bool((*col)["enabled"], true);
    n.density = get_num((*col)["neutral_density"], "collisions.neutral_density", 0.0);
    n.temperature = get_num((*col)["neutral_temperature"], "collisions.neutral_temperature", 0.025);
    n.mass = get_num((*col)["neutral_mass_amu"], "collisions.neutral_mass_amu", 1.0) * amu;
    n.drift = get_vec3((*col)["neutral_drift"], "collisions.neutral_drift", {0, 0, 0});
    if (auto pr = (*col)["process"].as_array())
      for (auto& pe : *pr) {
        const table& pt = *pe.as_table();
        CollisionProcess p;
        p.species = get_str(pt["species"], "");
        std::string ty = get_str(pt["type"], "elastic");
        if (ty == "elastic") p.type = CollisionProcess::Elastic;
        else if (ty == "cx" || ty == "charge_exchange") p.type = CollisionProcess::ChargeExchange;
        else if (ty == "excitation") p.type = CollisionProcess::Excitation;
        else if (ty == "ionization" || ty == "ionisation") p.type = CollisionProcess::Ionization;
        else fatal("unknown collision type '" + ty + "' (elastic|charge_exchange|excitation|ionization)");
        p.cross_section = get_num(pt["cross_section"], "collisions.process.cross_section", 0.0);
        p.threshold = get_num(pt["threshold"], "collisions.process.threshold", 0.0);
        p.product = get_str(pt["product"], "");
        p.secondary = get_str(pt["secondary"], "");
        std::string model = get_str(pt["model"], p.threshold > 0 ? "step" : "constant");
        if (model == "constant") p.model = CollisionProcess::Constant;
        else if (model == "step") p.model = CollisionProcess::Step;
        else if (model == "lotz") p.model = CollisionProcess::Lotz;
        else if (model == "table") p.model = CollisionProcess::Table;
        else fatal("unknown cross-section model '" + model + "' (constant|step|lotz|table)");
        if (p.model == CollisionProcess::Table) {
          if (auto tb = pt["table"].as_array()) {
            for (auto& row : *tb) {
              auto* r = row.as_array();
              if (!r || r->size() != 2) fatal("collision table rows must be [energy_eV, sigma_m2]");
              p.table.emplace_back(*(*r)[0].value<double>(), *(*r)[1].value<double>());
            }
          } else if (pt["file"]) {
            std::string fn = get_str(pt["file"], "");
            if (std::filesystem::path(fn).is_relative()) fn = (base / fn).string();
            std::ifstream tf(fn);
            if (!tf) fatal("cannot read cross-section file " + fn);
            std::string line;
            while (std::getline(tf, line)) {
              auto h = line.find_first_of("#%");
              if (h != std::string::npos) line = line.substr(0, h);
              std::istringstream ls(line);
              double e, sg;
              if (ls >> e >> sg) p.table.emplace_back(e, sg);
            }
          }
          if (p.table.size() < 2) fatal("cross-section table needs at least two rows (table = [[E, sigma], ...] or file)");
          std::sort(p.table.begin(), p.table.end());
        }
        if ((p.type == CollisionProcess::Excitation || p.type == CollisionProcess::Ionization) && p.threshold <= 0)
          fatal("excitation / ionization processes need threshold > 0 (eV)");
        if (p.type == CollisionProcess::Ionization && p.product.empty())
          fatal("ionization process needs 'product' (the ion species created)");
        n.processes.push_back(p);
      }
  }

  // ---------------- Coulomb collisions ---------------------------------------------------
  if (auto cc = root["coulomb"].as_table()) {
    auto& C = c.coulomb;
    C.enabled = get_bool((*cc)["enabled"], true);
    if ((*cc)["coulomb_log"]) {
      if ((*cc)["coulomb_log"].is_string()) C.auto_log = get_str((*cc)["coulomb_log"], "") == "auto";
      else C.log_lambda = get_num((*cc)["coulomb_log"], "coulomb.coulomb_log", 10.0);
    }
    C.every = std::max(1, (int)get_num((*cc)["every"], "coulomb.every", 1));
    if (auto pr = (*cc)["pairs"].as_array())
      for (auto& el : *pr) {
        auto* r = el.as_array();
        if (!r || r->size() != 2) fatal("coulomb.pairs entries must be [species_a, species_b]");
        C.pairs.emplace_back(*(*r)[0].value<std::string>(), *(*r)[1].value<std::string>());
      }
  }

  // ---------------- solver -------------------------------------------------------------
  if (auto s = root["solver"].as_table()) {
    auto& sc = c.solver;
    sc.type = get_str((*s)["type"], sc.type);
    sc.tol = get_num((*s)["tolerance"], "solver.tolerance", sc.tol);
    sc.max_iter = (int)get_num((*s)["max_iter"], "solver.max_iter", sc.max_iter);
    sc.cut_cell = get_bool((*s)["cut_cell"], sc.cut_cell);
    sc.mg_max_levels = (int)get_num((*s)["mg_levels"], "solver.mg_levels", sc.mg_max_levels);
    sc.mg_pre = (int)get_num((*s)["mg_pre"], "solver.mg_pre", sc.mg_pre);
    sc.mg_post = (int)get_num((*s)["mg_post"], "solver.mg_post", sc.mg_post);
    sc.mg_coarse_sweeps = (int)get_num((*s)["mg_coarse_sweeps"], "solver.mg_coarse_sweeps", sc.mg_coarse_sweeps);
    sc.sor_omega = get_num((*s)["sor_omega"], "solver.sor_omega", sc.sor_omega);
  }

  // ---------------- numerics / parallel -------------------------------------------------
  if (auto p = root["numerics"].as_table()) {
    c.shape_order = (int)get_num((*p)["shape_order"], "numerics.shape_order", 1);
    c.sort_every = (int)get_num((*p)["sort_every"], "numerics.sort_every", 20);
    c.surface_refine = (int)get_num((*p)["surface_refine"], "numerics.surface_refine", 2);
  }
  if (c.shape_order != 1 && c.shape_order != 2) fatal("numerics.shape_order must be 1 (CIC) or 2 (TSC)");
  if (auto p = root["parallel"].as_table()) {
    Vec3 d = get_vec3((*p)["decomposition"], "parallel.decomposition", {0, 0, 0});
    c.decomposition = {int(d.x), int(d.y), int(d.z)};
  }

  // ---------------- output -------------------------------------------------------------
  if (auto o = root["output"].as_table()) {
    auto& oc = c.output;
    oc.dir = get_str((*o)["dir"], oc.dir);
    if (std::filesystem::path(oc.dir).is_relative() && get_bool((*o)["relative_to_input"], false))
      oc.dir = (base / oc.dir).string();
    oc.fields_every = (int)get_num((*o)["fields_every"], "", oc.fields_every);
    oc.particles_every = (int)get_num((*o)["particles_every"], "", oc.particles_every);
    oc.particle_stride = (int)get_num((*o)["particle_stride"], "", oc.particle_stride);
    oc.surface_every = (int)get_num((*o)["surface_every"], "", oc.surface_every);
    oc.history_every = (int)get_num((*o)["history_every"], "", oc.history_every);
    oc.checkpoint_every = (int)get_num((*o)["checkpoint_every"], "", oc.checkpoint_every);
    oc.status_every = (int)get_num((*o)["status_every"], "", oc.status_every);
    oc.track_every = (int)get_num((*o)["track_every"], "", oc.track_every);
    oc.track_count = (int)get_num((*o)["track_count"], "", oc.track_count);
    oc.moments = get_bool((*o)["moments"], oc.moments);
    oc.float32 = get_bool((*o)["float32"], oc.float32);
    if (oc.particle_stride < 1) oc.particle_stride = 1;
  }
  return c;
}

void Config::print_summary() const {
  using namespace constants;
  log0("================================================================\n");
  log0(" PICSP3D  -  %s\n", title.c_str());
  log0("================================================================\n");
  log0(" grid      : %d x %d x %d cells, dx = (%.4g, %.4g, %.4g) m\n", cells[0], cells[1], cells[2], dx[0], dx[1], dx[2]);
  log0(" domain    : [%.4g, %.4g] x [%.4g, %.4g] x [%.4g, %.4g] m\n", origin.x, origin.x + cells[0] * dx[0], origin.y,
       origin.y + cells[1] * dx[1], origin.z, origin.z + cells[2] * dx[2]);
  log0(" time      : dt = %.4g s, %lld steps (T = %.4g s)\n", dt, steps, dt * steps);
  if (ref_wpe > 0) log0(" reference : w_p = %.4g rad/s (dt*w_p = %.3f)\n", ref_wpe, dt * ref_wpe);
  if (ref_debye > 0) log0("             lambda_D = %.4g m (dx/lambda_D = %.3f)\n", ref_debye, dx[0] / ref_debye);
  log0(" B0        : (%.3g, %.3g, %.3g) T,  E0 = (%.3g, %.3g, %.3g) V/m\n", B0.x, B0.y, B0.z, E0.x, E0.y, E0.z);
  if (!magnets.empty()) log0(" magnets   : %zu spatially varying field source(s)\n", magnets.size());
  if (!electrostatic) log0(" fields    : self-consistent field solve OFF (test-particle mode)\n");
  if (coulomb.enabled)
    log0(" coulomb   : binary collisions every %d step(s), ln(Lambda) = %s\n", coulomb.every,
         coulomb.auto_log ? "auto" : std::to_string(coulomb.log_lambda).c_str());
  const char* fn[6] = {"x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi"};
  for (int i = 0; i < 6; ++i)
    log0(" %-9s : field=%-9s (%.3g V)  particles=%s\n", fn[i], to_string(faces[i].field), faces[i].value,
         to_string(faces[i].particle));
  for (auto& s : species) {
    double B = norm(B0);
    double wc = std::fabs(s.charge) * B / s.mass;
    log0(" species %-10s q=%+.3g e  m=%.4g amu  n=%.3g m^-3  T=%.3g eV  ppc=%.3g%s\n", s.name.c_str(), s.charge / qe,
         s.mass / amu, s.density, s.temperature, s.ppc, s.subcycle > 1 ? "  (subcycled)" : "");
    if (wc > 0 && s.temperature > 0) {
      double vth = std::sqrt(std::max(s.temperature, 1e-30) * qe / s.mass);
      log0("            w_c = %.4g rad/s (dt*w_c = %.3f), r_L = %.4g m\n", wc, wc * dt * s.subcycle, vth / wc);
      if (wc * dt * s.subcycle > 0.5) warn0("species %s: dt*w_c > 0.5, gyration poorly resolved\n", s.name.c_str());
    }
  }
  for (auto& o : objects)
    log0(" object %-10s %zu part(s), %s%s\n", o.name.c_str(), o.parts.size(),
         o.bias == BiasMode::Floating ? "floating" : "fixed bias", o.photo.enabled ? ", photoemission" : "");
  if (neutrals.enabled && !neutrals.processes.empty())
    log0(" neutrals  : n = %.3g m^-3, T = %.3g eV, %zu process(es)\n", neutrals.density, neutrals.temperature,
         neutrals.processes.size());
  log0(" solver    : %s, tol = %.1e, cut-cell = %s, shape order = %d\n", solver.type.c_str(), solver.tol,
       solver.cut_cell ? "on" : "off", shape_order);
  if (ref_debye > 0 && dx[0] > 3.4 * ref_debye) warn0("dx > 3.4 lambda_D: expect finite-grid heating\n");
  if (ref_wpe > 0 && dt * ref_wpe > 0.25) warn0("dt*w_p > 0.25: leapfrog may be inaccurate (unstable above 2)\n");
  log0("================================================================\n");
}

}  // namespace picsp
