#include "geometry/StlMesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <tuple>

#include "util/Log.hpp"

namespace picsp {

Vec3 closest_point_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
  Vec3 ab = b - a, ac = c - a, ap = p - a;
  double d1 = dot(ab, ap), d2 = dot(ac, ap);
  if (d1 <= 0 && d2 <= 0) return a;
  Vec3 bp = p - b;
  double d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) return b;
  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) return a + ab * (d1 / (d1 - d3));
  Vec3 cp = p - c;
  double d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) return c;
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) return a + ac * (d2 / (d2 - d6));
  double va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
  double denom = 1.0 / (va + vb + vc);
  double v = vb * denom, w = vc * denom;
  return a + ab * v + ac * w;
}

bool StlMesh::load(const std::string& path, double scale) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return false;
  std::string content((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  tris_.clear();
  bool ascii = content.size() >= 5 && content.compare(0, 5, "solid") == 0;
  if (ascii) {
    // Some binary files also start with "solid": verify there is a "facet" keyword.
    if (content.find("facet") == std::string::npos) ascii = false;
  }
  if (content.size() >= 84) {
    std::uint32_t n;
    std::memcpy(&n, content.data() + 80, 4);
    if (84 + std::size_t(n) * 50 == content.size()) ascii = false;
  }
  if (ascii) {
    std::istringstream is(content);
    std::string tok;
    std::vector<Vec3> verts;
    while (is >> tok) {
      if (tok == "vertex") {
        Vec3 v;
        is >> v.x >> v.y >> v.z;
        verts.push_back(v * scale);
        if (verts.size() == 3) {
          tris_.push_back({verts[0], verts[1], verts[2]});
          verts.clear();
        }
      } else if (tok == "endloop") {
        verts.clear();
      }
    }
  } else {
    if (content.size() < 84) return false;
    std::uint32_t n;
    std::memcpy(&n, content.data() + 80, 4);
    if (content.size() < 84 + std::size_t(n) * 50) return false;
    tris_.reserve(n);
    for (std::uint32_t t = 0; t < n; ++t) {
      const char* p = content.data() + 84 + std::size_t(t) * 50 + 12;  // skip normal
      float v[9];
      std::memcpy(v, p, sizeof(v));
      tris_.push_back({Vec3{v[0], v[1], v[2]} * scale, Vec3{v[3], v[4], v[5]} * scale, Vec3{v[6], v[7], v[8]} * scale});
    }
  }
  // remove degenerate triangles
  tris_.erase(std::remove_if(tris_.begin(), tris_.end(),
                             [](const Triangle& t) { return norm(cross(t.b - t.a, t.c - t.a)) < 1e-300; }),
              tris_.end());
  // watertightness check: every edge must be shared by an even number of triangles
  {
    std::map<std::tuple<long long, long long, long long, long long, long long, long long>, int> edges;
    double q = 1e-9 * std::max(1.0, scale);
    auto key = [&](const Vec3& a, const Vec3& b) {
      long long ax = llround(a.x / q), ay = llround(a.y / q), az = llround(a.z / q);
      long long bx = llround(b.x / q), by = llround(b.y / q), bz = llround(b.z / q);
      if (std::tie(ax, ay, az) > std::tie(bx, by, bz)) {
        std::swap(ax, bx); std::swap(ay, by); std::swap(az, bz);
      }
      return std::make_tuple(ax, ay, az, bx, by, bz);
    };
    if (tris_.size() < 200000) {
      for (auto& t : tris_) {
        edges[key(t.a, t.b)]++;
        edges[key(t.b, t.c)]++;
        edges[key(t.c, t.a)]++;
      }
      closed_ = true;
      for (auto& kv : edges)
        if (kv.second % 2) { closed_ = false; break; }
    }
  }
  build_bvh();
  return !tris_.empty();
}

void StlMesh::build_bvh() {
  nodes_.clear();
  bounds_ = Aabb{};
  for (auto& t : tris_) {
    bounds_.expand(t.a);
    bounds_.expand(t.b);
    bounds_.expand(t.c);
  }
  if (!tris_.empty()) {
    nodes_.reserve(2 * tris_.size() / 4 + 1);
    build(0, int(tris_.size()), 0);
    build_moments(0);
  }
}

void StlMesh::build_moments(int id) {
  Node& n = nodes_[id];
  Vec3 an{0, 0, 0}, ac{0, 0, 0};
  double A = 0;
  for (int i = n.first; i < n.first + n.count; ++i) {
    const Triangle& t = tris_[i];
    Vec3 cr = cross(t.b - t.a, t.c - t.a);
    double a = 0.5 * norm(cr);
    an += cr * 0.5;
    ac += (t.a + t.b + t.c) * (a / 3.0);
    A += a;
  }
  n.an = an;
  n.c = A > 0 ? ac / A : n.box.center();
  double r = 0;
  for (int i = n.first; i < n.first + n.count; ++i)
    for (const Vec3* v : {&tris_[i].a, &tris_[i].b, &tris_[i].c}) r = std::max(r, norm(*v - n.c));
  n.r = r;
  if (n.left >= 0) {
    build_moments(n.left);
    build_moments(n.right);
  }
}

double StlMesh::winding_number(const Vec3& p) const {
  if (nodes_.empty()) return 0.0;
  constexpr double beta = 2.3;  // far-field acceptance (distance > beta * radius)
  double omega = 0;
  int stack[128];
  int sp = 0;
  stack[sp++] = 0;
  while (sp) {
    const Node& n = nodes_[stack[--sp]];
    Vec3 d = n.c - p;
    double dist = norm(d);
    if (dist > beta * n.r && n.count > 4) {
      omega += dot(n.an, d) / (dist * dist * dist);  // dipole approximation
      continue;
    }
    if (n.left < 0) {
      for (int i = n.first; i < n.first + n.count; ++i) {
        const Triangle& t = tris_[i];
        Vec3 a = t.a - p, b = t.b - p, c = t.c - p;
        double la = norm(a), lb = norm(b), lc = norm(c);
        double num = dot(a, cross(b, c));
        double den = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la;
        omega += 2.0 * std::atan2(num, den);  // Van Oosterom & Strackee
      }
    } else {
      stack[sp++] = n.left;
      stack[sp++] = n.right;
    }
  }
  return omega / (4.0 * 3.14159265358979323846);
}

int StlMesh::build(int first, int count, int depth) {
  Node node;
  node.first = first;
  node.count = count;
  Aabb cb;  // centroid bounds
  for (int i = first; i < first + count; ++i) {
    node.box.expand(tris_[i].a);
    node.box.expand(tris_[i].b);
    node.box.expand(tris_[i].c);
    cb.expand((tris_[i].a + tris_[i].b + tris_[i].c) / 3.0);
  }
  int id = int(nodes_.size());
  nodes_.push_back(node);
  if (count <= 4 || depth > 40) return id;
  Vec3 ext = cb.hi - cb.lo;
  int axis = ext.x > ext.y ? (ext.x > ext.z ? 0 : 2) : (ext.y > ext.z ? 1 : 2);
  if (ext[axis] <= 0) return id;
  int mid = first + count / 2;
  std::nth_element(tris_.begin() + first, tris_.begin() + mid, tris_.begin() + first + count,
                   [axis](const Triangle& u, const Triangle& v) {
                     return (u.a[axis] + u.b[axis] + u.c[axis]) < (v.a[axis] + v.b[axis] + v.c[axis]);
                   });
  int l = build(first, mid - first, depth + 1);
  int r = build(mid, first + count - mid, depth + 1);
  nodes_[id].left = l;
  nodes_[id].right = r;
  return id;
}

double StlMesh::distance(const Vec3& p, double cap) const {
  if (nodes_.empty()) return cap;
  double best2 = cap * cap;
  int stack[128];
  int sp = 0;
  stack[sp++] = 0;
  while (sp) {
    const Node& n = nodes_[stack[--sp]];
    if (n.box.dist2(p) >= best2) continue;
    if (n.left < 0) {
      for (int i = n.first; i < n.first + n.count; ++i) {
        const Triangle& t = tris_[i];
        Vec3 q = closest_point_triangle(p, t.a, t.b, t.c);
        Vec3 d = q - p;
        double d2 = dot(d, d);
        if (d2 < best2) best2 = d2;
      }
    } else {
      // visit nearer child last (popped first)
      double dl = nodes_[n.left].box.dist2(p), dr = nodes_[n.right].box.dist2(p);
      if (dl < dr) {
        stack[sp++] = n.right;
        stack[sp++] = n.left;
      } else {
        stack[sp++] = n.left;
        stack[sp++] = n.right;
      }
    }
  }
  return std::sqrt(best2);
}

namespace {
inline bool ray_box(const Vec3& o, const Vec3& inv, const Aabb& b) {
  double tmin = 0, tmax = 1e300;
  for (int a = 0; a < 3; ++a) {
    double t1 = (b.lo[a] - o[a]) * inv[a], t2 = (b.hi[a] - o[a]) * inv[a];
    if (t1 > t2) std::swap(t1, t2);
    tmin = std::max(tmin, t1);
    tmax = std::min(tmax, t2);
  }
  return tmin <= tmax;
}
// Moller-Trumbore, counts hits with t > 0
inline bool ray_tri(const Vec3& o, const Vec3& d, const Triangle& t) {
  Vec3 e1 = t.b - t.a, e2 = t.c - t.a;
  Vec3 pv = cross(d, e2);
  double det = dot(e1, pv);
  if (std::fabs(det) < 1e-300) return false;
  double inv = 1.0 / det;
  Vec3 tv = o - t.a;
  double u = dot(tv, pv) * inv;
  if (u < 0 || u > 1) return false;
  Vec3 qv = cross(tv, e1);
  double v = dot(d, qv) * inv;
  if (v < 0 || u + v > 1) return false;
  return dot(e2, qv) * inv > 0;
}
}  // namespace

int StlMesh::ray_crossings(const Vec3& o, const Vec3& d) const {
  Vec3 inv{1.0 / d.x, 1.0 / d.y, 1.0 / d.z};
  int hits = 0;
  int stack[128];
  int sp = 0;
  stack[sp++] = 0;
  while (sp) {
    const Node& n = nodes_[stack[--sp]];
    if (!ray_box(o, inv, n.box)) continue;
    if (n.left < 0) {
      for (int i = n.first; i < n.first + n.count; ++i) hits += ray_tri(o, d, tris_[i]);
    } else {
      stack[sp++] = n.left;
      stack[sp++] = n.right;
    }
  }
  return hits;
}

bool StlMesh::inside(const Vec3& p) const {
  if (nodes_.empty() || bounds_.dist2(p) > 0) return false;
  return std::fabs(winding_number(p)) > 0.5;
}

bool StlMesh::inside_parity(const Vec3& p) const {
  if (nodes_.empty() || bounds_.dist2(p) > 0) return false;
  // Three generic (non axis-aligned) directions; majority vote guards against rays
  // that graze edges/vertices of the tessellation.
  static const Vec3 dirs[3] = {normalized(Vec3{0.5773, 0.5779, 0.5767}), normalized(Vec3{-0.6012, 0.3071, 0.7377}),
                               normalized(Vec3{0.2113, -0.8931, 0.3971})};
  int votes = 0;
  for (auto& d : dirs) votes += ray_crossings(p, d) & 1;
  return votes >= 2;
}

double StlMesh::sdf(const Vec3& p, double cap) const {
  double d = bounds_.dist2(p);
  if (d > cap * cap) return std::sqrt(d);  // far outside the bounding box
  double dist = distance(p, cap);
  return inside(p) ? -dist : dist;
}

}  // namespace picsp
