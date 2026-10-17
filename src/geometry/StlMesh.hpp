#pragma once
// Triangle mesh (STL ASCII/binary) with an AABB bounding-volume hierarchy for fast
// closest-point distance and ray-parity inside/outside queries.
#include <array>
#include <string>
#include <vector>

#include "util/Vec3.hpp"

namespace picsp {

struct Aabb {
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  void expand(const Vec3& p) { lo = vmin(lo, p); hi = vmax(hi, p); }
  void expand(const Aabb& b) { lo = vmin(lo, b.lo); hi = vmax(hi, b.hi); }
  bool valid() const { return lo.x <= hi.x; }
  Vec3 center() const { return (lo + hi) * 0.5; }
  // squared distance from point to box (0 if inside)
  double dist2(const Vec3& p) const {
    double d = 0;
    for (int a = 0; a < 3; ++a) {
      double v = p[a] < lo[a] ? lo[a] - p[a] : (p[a] > hi[a] ? p[a] - hi[a] : 0.0);
      d += v * v;
    }
    return d;
  }
};

struct Triangle {
  Vec3 a, b, c;
};

class StlMesh {
 public:
  bool load(const std::string& path, double scale = 1.0);
  void build_bvh();
  std::size_t size() const { return tris_.size(); }
  const Aabb& bounds() const { return bounds_; }
  const std::vector<Triangle>& triangles() const { return tris_; }

  // unsigned distance to the surface, optionally capped (search pruned beyond cap)
  double distance(const Vec3& p, double cap = 1e300) const;
  // inside test: |generalised winding number| > 1/2 (robust to overlapping shells,
  // small gaps and inconsistent orientation); evaluated with a BVH dipole far-field
  bool inside(const Vec3& p) const;
  double winding_number(const Vec3& p) const;
  bool inside_parity(const Vec3& p) const;  // ray-parity fallback
  // signed distance (negative inside); |value| capped at `cap`
  double sdf(const Vec3& p, double cap = 1e300) const;
  bool is_closed() const { return closed_; }

 private:
  struct Node {
    Aabb box;
    int left = -1, right = -1;  // children; leaf if left < 0
    int first = 0, count = 0;   // triangle range for leaves
    Vec3 an{0, 0, 0};           // sum of area-weighted normals (dipole moment)
    Vec3 c{0, 0, 0};            // area centroid
    double r = 0;               // radius of the node's triangles about c
  };
  void build_moments(int node);
  int build(int first, int count, int depth);
  int ray_crossings(const Vec3& o, const Vec3& d) const;
  std::vector<Triangle> tris_;
  std::vector<Node> nodes_;
  Aabb bounds_;
  bool closed_ = true;
};

// Closest point on triangle (Ericson, Real-Time Collision Detection 5.1.5)
Vec3 closest_point_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c);

}  // namespace picsp
