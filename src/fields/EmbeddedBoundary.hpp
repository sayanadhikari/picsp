#pragma once
// Grid representation of the embedded objects on the local sub-domain:
//  * node classification (free / object node, owning object id)
//  * exact node signed distances
//  * Shortley-Weller cut-cell coefficients for free nodes next to a surface
//  * "near surface" cell flags used to restrict particle-object collision tests
#include <cstdint>
#include <vector>

#include "config/Config.hpp"
#include "geometry/Geometry.hpp"
#include "grid/Field3D.hpp"

namespace picsp {

struct CutNode {
  std::size_t idx;    // local field index
  double c[6];        // coefficients for neighbours -x,+x,-y,+y,-z,+z
  double theta[6];    // fractional distances (1 = regular neighbour)
};

class EmbeddedBoundary {
 public:
  void build(const Config& cfg, const GridGeom& g, int ghosts, const Scene& scene);

  bool has_objects() const { return has_objects_; }
  // node classification in local field layout (-1 = free)
  int obj(std::size_t idx) const { return node_obj_.empty() ? -1 : node_obj_[idx]; }
  const std::vector<std::int16_t>& node_obj() const { return node_obj_; }
  const Field3D& sdf() const { return sdf_; }
  const std::vector<CutNode>& cut_nodes() const { return cut_; }

  // cell flags: local cell (ci, cj, ck) relative to c0 may be in [-1, nc]
  inline bool cell_near(int ci, int cj, int ck) const {
    if (near_.empty()) return false;
    return near_[std::size_t(ci + 1) + std::size_t(ncx_ + 2) * (std::size_t(cj + 1) + std::size_t(ncy_ + 2) * (ck + 1))];
  }
  long long count_object_nodes() const { return n_obj_nodes_; }
  long long count_cut_nodes() const { return (long long)cut_.size(); }

 private:
  bool has_objects_ = false;
  std::vector<std::int16_t> node_obj_;
  Field3D sdf_;
  std::vector<CutNode> cut_;
  std::vector<std::uint8_t> near_;
  int ncx_ = 0, ncy_ = 0, ncz_ = 0;
  long long n_obj_nodes_ = 0;
};

}  // namespace picsp
