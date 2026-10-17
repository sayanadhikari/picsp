#pragma once
// Parallel VTK XML output (no external dependencies).
//   fields   : one .vti piece per rank + .pvti master  (ImageData, point data)
//   particles: one .vtp piece per rank + .pvtp master  (PolyData vertices)
//   surfaces : one .vtp piece per rank + .pvtp master  (PolyData triangles, cell data)
// plus .pvd collections (time series) that ParaView / pyvista open as animations.
// Binary data are written raw-appended (UInt64 headers) for speed and compactness.
#include <map>
#include <string>
#include <vector>

#include "grid/Field3D.hpp"
#include "parallel/Decomposition.hpp"
#include "util/Vec3.hpp"

namespace picsp {

struct NamedScalar {
  std::string name;
  const Field3D* f;
};
struct NamedVector {
  std::string name;
  const Field3D *x, *y, *z;
};

// generic per-point / per-cell array for poly data
struct PolyArray {
  std::string name;
  int ncomp = 1;
  std::vector<float> data;  // size n * ncomp
};

class VtkWriter {
 public:
  void init(const std::string& dir, const Decomposition& d, bool float32);

  void write_fields(long long step, double time, const std::vector<NamedScalar>& scalars,
                    const std::vector<NamedVector>& vectors);
  // points: xyz triplets
  void write_points(const std::string& series, long long step, double time, const std::vector<float>& xyz,
                    const std::vector<PolyArray>& point_data);
  // triangles: 9 floats per triangle; arrays are cell data (one value per triangle)
  void write_triangles(const std::string& series, long long step, double time, const std::vector<float>& tri,
                       const std::vector<PolyArray>& cell_data);

 private:
  void add_to_collection(const std::string& series, double time, const std::string& file);
  std::string dir_;
  const Decomposition* d_ = nullptr;
  bool f32_ = true;
  std::map<std::string, std::vector<std::pair<double, std::string>>> collections_;
};

}  // namespace picsp
