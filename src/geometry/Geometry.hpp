#pragma once
// Solid geometry described by signed distance functions (negative inside).
// Objects are CSG combinations of analytic primitives and/or STL meshes, each with
// their own rigid transform. Every MPI rank holds the full (global) geometry, so
// any point can be classified locally without communication.
#include <memory>
#include <string>
#include <vector>

#include "config/Config.hpp"
#include "geometry/StlMesh.hpp"
#include "util/Vec3.hpp"

namespace picsp {

class Part {
 public:
  explicit Part(const PartConfig& cfg);
  // signed distance in object coordinates
  double sdf(const Vec3& p, double cap) const;
  // unsigned distance to this part's surface (no inside test), capped
  double udf(const Vec3& p, double cap) const;
  Aabb bounds() const;  // in object coordinates (conservative)
  CsgOp op() const { return cfg_.op; }
  double min_feature() const;  // smallest characteristic thickness (for warnings)
  const PartConfig& config() const { return cfg_; }

 private:
  PartConfig cfg_;
  Mat3 R_, Rt_;  // part rotation and its transpose
  std::shared_ptr<StlMesh> mesh_;
};

class SolidObject {
 public:
  SolidObject(const ObjectConfig& cfg, int id);
  double sdf(const Vec3& p, double cap = 1e300) const;
  // lower bound of the distance to the object's surface (min over part surfaces), capped
  double distance_lb(const Vec3& p, double cap) const;
  Vec3 gradient(const Vec3& p, double h) const;  // outward normal direction (unnormalised)
  const Aabb& bounds() const { return bounds_; }  // world coordinates (conservative)
  const std::string& name() const { return cfg_.name; }
  int id() const { return id_; }
  const ObjectConfig& config() const { return cfg_; }
  double min_feature() const;

 private:
  ObjectConfig cfg_;
  int id_;
  Mat3 R_, Rt_;
  std::vector<Part> parts_;
  Aabb bounds_;
};

// Collection of all solid objects
class Scene {
 public:
  void build(const std::vector<ObjectConfig>& objs);
  bool empty() const { return objs_.empty(); }
  std::size_t size() const { return objs_.size(); }
  const SolidObject& object(int i) const { return objs_[i]; }

  // Minimum SDF over all objects; `which` receives the nearest object id (-1 if none)
  double sdf(const Vec3& p, int* which = nullptr, double cap = 1e300) const;

  // Lower bound of the distance to any object surface (cheap: no inside tests), capped
  double distance_lb(const Vec3& p, double cap, int* which = nullptr) const;

  // Sphere-trace the segment p0->p1 (p0 outside all objects) for the first surface crossing.
  // Returns true and fills hit fraction t in [0,1] and object id on a hit.
  bool segment_hit(const Vec3& p0, const Vec3& p1, double tol, double* t_hit, int* which) const;

  // Is the ray from p in direction d (unit) blocked by any object before leaving
  // the box [lo, hi]? (used for photoemission shadowing)
  bool ray_blocked(const Vec3& p, const Vec3& d, const Vec3& lo, const Vec3& hi, double tol) const;

  // Project a point onto the nearest surface of object `id`
  Vec3 project_to_surface(const Vec3& p, int id, double h) const;

 private:
  std::vector<SolidObject> objs_;
};

}  // namespace picsp
