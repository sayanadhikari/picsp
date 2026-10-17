#include "geometry/Geometry.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>

#include "util/Log.hpp"

namespace picsp {

namespace {

inline double length2d(double a, double b) { return std::sqrt(a * a + b * b); }

double sd_box(const Vec3& p, const Vec3& half, double r) {
  Vec3 q{std::fabs(p.x) - half.x + r, std::fabs(p.y) - half.y + r, std::fabs(p.z) - half.z + r};
  Vec3 qm{std::max(q.x, 0.0), std::max(q.y, 0.0), std::max(q.z, 0.0)};
  return norm(qm) + std::min(std::max(q.x, std::max(q.y, q.z)), 0.0) - r;
}

// axis-aligned decomposition: h = coordinate along axis, rho = distance from axis
inline void axial(const Vec3& p, const Vec3& axis, double& h, double& rho) {
  h = dot(p, axis);
  Vec3 perp = p - axis * h;
  rho = norm(perp);
}

double sd_cylinder(const Vec3& p, const Vec3& axis, double r, double L) {
  double h, rho;
  axial(p, axis, h, rho);
  double dr = rho - r, dh = std::fabs(h) - 0.5 * L;
  return std::min(std::max(dr, dh), 0.0) + length2d(std::max(dr, 0.0), std::max(dh, 0.0));
}

// Capped cone (frustum): radius r1 at h = -L/2, r2 at h = +L/2 (I. Quilez)
double sd_cone(const Vec3& p, const Vec3& axis, double r1, double r2, double L) {
  double qy, qx;
  axial(p, axis, qy, qx);
  double hh = 0.5 * L;
  double k1x = r2, k1y = hh;
  double k2x = r2 - r1, k2y = 2.0 * hh;
  double cax = qx - std::min(qx, qy < 0.0 ? r1 : r2), cay = std::fabs(qy) - hh;
  double t = ((k1x - qx) * k2x + (k1y - qy) * k2y) / (k2x * k2x + k2y * k2y);
  t = std::clamp(t, 0.0, 1.0);
  double cbx = qx - k1x + k2x * t, cby = qy - k1y + k2y * t;
  double s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
  return s * std::sqrt(std::min(cax * cax + cay * cay, cbx * cbx + cby * cby));
}

double sd_torus(const Vec3& p, const Vec3& axis, double R, double r) {
  double h, rho;
  axial(p, axis, h, rho);
  return length2d(rho - R, h) - r;
}

double sd_ellipsoid(const Vec3& p, const Vec3& r) {
  Vec3 a{p.x / r.x, p.y / r.y, p.z / r.z};
  Vec3 b{p.x / (r.x * r.x), p.y / (r.y * r.y), p.z / (r.z * r.z)};
  double k0 = norm(a), k1 = norm(b);
  if (k1 == 0) return -std::min(r.x, std::min(r.y, r.z));
  return k0 * (k0 - 1.0) / k1;
}

// STL meshes are cached so that multiple objects/parts referencing the same file
// share one BVH.
std::shared_ptr<StlMesh> load_mesh_cached(const std::string& file, double scale) {
  static std::map<std::pair<std::string, double>, std::shared_ptr<StlMesh>> cache;
  static std::mutex mtx;
  std::lock_guard<std::mutex> lock(mtx);
  auto key = std::make_pair(file, scale);
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;
  auto m = std::make_shared<StlMesh>();
  if (!m->load(file, scale)) fatal("cannot read STL file '" + file + "'");
  if (!m->is_closed())
    log0(" note: STL '%s' has open edges or overlapping shells; using generalised winding numbers\n", file.c_str());
  log0(" loaded STL %s: %zu triangles, bbox [%.3g %.3g %.3g]-[%.3g %.3g %.3g] m\n", file.c_str(), m->size(),
       m->bounds().lo.x, m->bounds().lo.y, m->bounds().lo.z, m->bounds().hi.x, m->bounds().hi.y, m->bounds().hi.z);
  cache[key] = m;
  return m;
}

Aabb transform_box(const Aabb& b, const Mat3& R, const Vec3& t) {
  Aabb out;
  for (int c = 0; c < 8; ++c) {
    Vec3 p{(c & 1) ? b.hi.x : b.lo.x, (c & 2) ? b.hi.y : b.lo.y, (c & 4) ? b.hi.z : b.lo.z};
    out.expand(R * p + t);
  }
  return out;
}

}  // namespace

Part::Part(const PartConfig& cfg) : cfg_(cfg) {
  R_ = Mat3::from_euler_deg(cfg.rotate_deg);
  Rt_ = R_.transposed();
  if (cfg.shape == ShapeKind::STL) mesh_ = load_mesh_cached(cfg.file, cfg.scale);
  auto bad = [&](const char* what) { fatal(std::string("object part: ") + what); };
  switch (cfg.shape) {
    case ShapeKind::Sphere: if (cfg.radius <= 0) bad("sphere needs radius > 0"); break;
    case ShapeKind::Box:
    case ShapeKind::Plate:
    case ShapeKind::Ellipsoid:
      if (cfg.size.x <= 0 || cfg.size.y <= 0 || cfg.size.z <= 0) bad("box/plate/ellipsoid need size = [sx, sy, sz] > 0");
      break;
    case ShapeKind::Cylinder: if (cfg.radius <= 0 || cfg.length <= 0) bad("cylinder needs radius and length"); break;
    case ShapeKind::Cone: if (cfg.length <= 0) bad("cone needs length, radius (base) and radius2 (top)"); break;
    case ShapeKind::Torus: if (cfg.radius <= 0 || cfg.radius2 <= 0) bad("torus needs radius (major) and radius2 (minor)"); break;
    default: break;
  }
}

double Part::sdf(const Vec3& pw, double cap) const {
  Vec3 p = Rt_ * (pw - cfg_.center);
  switch (cfg_.shape) {
    case ShapeKind::Sphere: return norm(p) - cfg_.radius;
    case ShapeKind::Box:
    case ShapeKind::Plate: return sd_box(p, cfg_.size * 0.5, std::min(cfg_.rounding, 0.5 * std::min(cfg_.size.x, std::min(cfg_.size.y, cfg_.size.z))));
    case ShapeKind::Cylinder: return sd_cylinder(p, cfg_.axis, cfg_.radius, cfg_.length);
    case ShapeKind::Cone: return sd_cone(p, cfg_.axis, cfg_.radius, cfg_.radius2, cfg_.length);
    case ShapeKind::Torus: return sd_torus(p, cfg_.axis, cfg_.radius, cfg_.radius2);
    case ShapeKind::Ellipsoid: return sd_ellipsoid(p, cfg_.size * 0.5);
    case ShapeKind::HalfSpace: return dot(p, cfg_.axis);
    case ShapeKind::STL: return mesh_->sdf(p, cap);
  }
  return 1e300;
}

double Part::udf(const Vec3& pw, double cap) const {
  if (cfg_.shape == ShapeKind::STL) return mesh_->distance(Rt_ * (pw - cfg_.center), cap);
  return std::fabs(sdf(pw, cap));
}

Aabb Part::bounds() const {
  Aabb b;
  auto sym = [&](double r) {
    b.lo = {-r, -r, -r};
    b.hi = {r, r, r};
  };
  switch (cfg_.shape) {
    case ShapeKind::Sphere: sym(cfg_.radius); break;
    case ShapeKind::Box:
    case ShapeKind::Plate:
    case ShapeKind::Ellipsoid:
      b.lo = cfg_.size * -0.5;
      b.hi = cfg_.size * 0.5;
      break;
    case ShapeKind::Cylinder: sym(std::sqrt(cfg_.radius * cfg_.radius + 0.25 * cfg_.length * cfg_.length)); break;
    case ShapeKind::Cone: sym(std::sqrt(std::pow(std::max(cfg_.radius, cfg_.radius2), 2) + 0.25 * cfg_.length * cfg_.length)); break;
    case ShapeKind::Torus: sym(cfg_.radius + cfg_.radius2); break;
    case ShapeKind::HalfSpace: sym(1e30); break;
    case ShapeKind::STL: b = mesh_->bounds(); break;
  }
  return transform_box(b, R_, cfg_.center);
}

double Part::min_feature() const {
  switch (cfg_.shape) {
    case ShapeKind::Sphere: return 2 * cfg_.radius;
    case ShapeKind::Box:
    case ShapeKind::Plate:
    case ShapeKind::Ellipsoid: return std::min(cfg_.size.x, std::min(cfg_.size.y, cfg_.size.z));
    case ShapeKind::Cylinder: return std::min(2 * cfg_.radius, cfg_.length);
    case ShapeKind::Cone: return std::min(2 * std::max(cfg_.radius, cfg_.radius2), cfg_.length);
    case ShapeKind::Torus: return 2 * cfg_.radius2;
    default: return 1e300;
  }
}

SolidObject::SolidObject(const ObjectConfig& cfg, int id) : cfg_(cfg), id_(id) {
  R_ = Mat3::from_euler_deg(cfg.rotate_deg);
  Rt_ = R_.transposed();
  if (cfg.parts.empty()) fatal("object '" + cfg.name + "' has no parts");
  if (cfg.parts[0].op != CsgOp::Union) fatal("object '" + cfg.name + "': first part must be a union");
  for (auto& pc : cfg.parts) parts_.emplace_back(pc);
  Aabb b;
  for (auto& p : parts_)
    if (p.op() == CsgOp::Union) b.expand(p.bounds());
  bounds_ = transform_box(b, R_, cfg.translate);
}

double SolidObject::sdf(const Vec3& pw, double cap) const {
  Vec3 p = Rt_ * (pw - cfg_.translate);
  double d = parts_[0].sdf(p, cap);
  for (std::size_t i = 1; i < parts_.size(); ++i) {
    double di = parts_[i].sdf(p, cap);
    switch (parts_[i].op()) {
      case CsgOp::Union: d = std::min(d, di); break;
      case CsgOp::Subtract: d = std::max(d, -di); break;
      case CsgOp::Intersect: d = std::max(d, di); break;
    }
  }
  return d;
}

double SolidObject::distance_lb(const Vec3& pw, double cap) const {
  Vec3 p = Rt_ * (pw - cfg_.translate);
  double d = cap;
  for (auto& part : parts_) d = std::min(d, part.udf(p, d));
  return d;
}

Vec3 SolidObject::gradient(const Vec3& p, double h) const {
  Vec3 g;
  for (int a = 0; a < 3; ++a) {
    Vec3 e{0, 0, 0};
    e[a] = h;
    g[a] = sdf(p + e) - sdf(p - e);
  }
  return g;
}

double SolidObject::min_feature() const {
  double m = 1e300;
  for (auto& p : parts_)
    if (p.op() == CsgOp::Union) m = std::min(m, p.min_feature());
  return m;
}

void Scene::build(const std::vector<ObjectConfig>& objs) {
  objs_.clear();
  for (std::size_t i = 0; i < objs.size(); ++i) objs_.emplace_back(objs[i], int(i));
}

double Scene::sdf(const Vec3& p, int* which, double cap) const {
  double best = 1e300;
  int w = -1;
  for (auto& o : objs_) {
    double bd = o.bounds().dist2(p);
    if (bd > best * best && best > 0) continue;  // cannot beat current best
    double d = o.sdf(p, cap);
    if (d < best) {
      best = d;
      w = o.id();
    }
  }
  if (which) *which = w;
  return best;
}

double Scene::distance_lb(const Vec3& p, double cap, int* which) const {
  double best = cap;
  int w = -1;
  for (auto& o : objs_) {
    if (o.bounds().dist2(p) >= best * best) continue;
    double d = o.distance_lb(p, best);
    if (d < best) {
      best = d;
      w = o.id();
    }
  }
  if (which) *which = w;
  return best;
}

bool Scene::segment_hit(const Vec3& p0, const Vec3& p1, double tol, double* t_hit, int* which) const {
  Vec3 d = p1 - p0;
  const double L = norm(d);
  if (L <= 0) return false;
  const Vec3 u = d / L;
  double t = 0;
  int w;
  // Sphere tracing with a cheap distance lower bound (valid since p0 is outside). A
  // near-zero bound may belong to a part surface that is not on the CSG boundary (e.g. a
  // subtracted tool in free space): confirm with the exact signed distance, which is
  // itself a safe step size.
  for (int it = 0; it < 96; ++it) {
    Vec3 q = p0 + u * t;
    double s = distance_lb(q, L - t + 2 * tol, &w);
    if (s <= tol) {
      double st = sdf(q, &w);
      if (st <= tol) {
        *t_hit = std::min(t / L, 1.0);
        *which = w;
        return true;
      }
      s = st;
    }
    t += std::max(s, 0.25 * tol);
    if (t > L) return false;
  }
  // grazing trajectory that did not converge: decide by the end point
  if (sdf(p1, &w) <= 0) {
    *t_hit = 1.0;
    *which = w;
    return true;
  }
  return false;
}

bool Scene::ray_blocked(const Vec3& p, const Vec3& d, const Vec3& lo, const Vec3& hi, double tol) const {
  double t = 4 * tol;
  for (int it = 0; it < 4096; ++it) {
    Vec3 q = p + d * t;
    for (int a = 0; a < 3; ++a)
      if (q[a] < lo[a] - tol || q[a] > hi[a] + tol) return false;
    double s = sdf(q);
    if (s < 0.5 * tol) return true;
    t += std::max(s, 0.5 * tol);
  }
  return false;
}

Vec3 Scene::project_to_surface(const Vec3& p, int id, double h) const {
  const SolidObject& o = objs_[id];
  Vec3 q = p;
  for (int it = 0; it < 4; ++it) {
    double s = o.sdf(q);
    Vec3 g = normalized(o.gradient(q, 0.25 * h));
    if (norm(g) == 0) break;
    q = q - g * s;
    if (std::fabs(s) < 1e-3 * h) break;
  }
  return q;
}

}  // namespace picsp
