#include "io/Vtk.hpp"

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "util/Log.hpp"

namespace picsp {

namespace fs = std::filesystem;

namespace {

std::string stepname(const std::string& prefix, long long step) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%s_%08lld", prefix.c_str(), step);
  return buf;
}

// Collects appended binary blocks and the XML that references them
struct Appender {
  std::string blob;
  std::uint64_t offset = 0;
  template <typename T>
  std::uint64_t add(const T* data, std::size_t n) {
    std::uint64_t off = offset;
    std::uint64_t bytes = n * sizeof(T);
    blob.append(reinterpret_cast<const char*>(&bytes), sizeof(bytes));
    blob.append(reinterpret_cast<const char*>(data), bytes);
    offset += sizeof(bytes) + bytes;
    return off;
  }
};

void write_file(const std::string& path, const std::string& xml_head, const Appender& app, const std::string& xml_tail) {
  std::ofstream f(path, std::ios::binary);
  if (!f) fatal("cannot write " + path);
  f << xml_head;
  f << "  <AppendedData encoding=\"raw\">\n   _";
  f.write(app.blob.data(), std::streamsize(app.blob.size()));
  f << "\n  </AppendedData>\n" << xml_tail;
}

}  // namespace

void VtkWriter::init(const std::string& dir, const Decomposition& d, bool float32) {
  dir_ = dir;
  d_ = &d;
  f32_ = float32;
  if (d.rank == 0) {
    fs::create_directories(dir_);
    fs::create_directories(dir_ + "/pieces");
  }
  MPI_Barrier(d.comm);
}

void VtkWriter::add_to_collection(const std::string& series, double time, const std::string& file) {
  if (d_->rank != 0) return;
  auto& c = collections_[series];
  c.emplace_back(time, file);
  std::ofstream f(dir_ + "/" + series + ".pvd");
  f << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n  <Collection>\n";
  for (auto& e : c) f << "    <DataSet timestep=\"" << e.first << "\" part=\"0\" file=\"" << e.second << "\"/>\n";
  f << "  </Collection>\n</VTKFile>\n";
}

void VtkWriter::write_fields(long long step, double time, const std::vector<NamedScalar>& scalars,
                             const std::vector<NamedVector>& vectors) {
  const Decomposition& d = *d_;
  const std::string base = stepname("fields", step);
  const char* ftype = f32_ ? "Float32" : "Float64";
  // piece extent: owned nodes plus the shared node at c1 (always present in the array)
  int ext[6];
  for (int a = 0; a < 3; ++a) {
    ext[2 * a] = d.c0[a];
    ext[2 * a + 1] = d.c1[a];
  }
  auto extent_str = [](const int* e) {
    std::ostringstream os;
    os << e[0] << " " << e[1] << " " << e[2] << " " << e[3] << " " << e[4] << " " << e[5];
    return os.str();
  };
  const std::size_t np = std::size_t(ext[1] - ext[0] + 1) * (ext[3] - ext[2] + 1) * (ext[5] - ext[4] + 1);
  Appender app;
  std::ostringstream pd;
  auto extract = [&](const Field3D& f, std::vector<double>& out) {
    out.resize(np);
    const int G = f.ghosts();
    std::size_t c = 0;
    for (int k = ext[4]; k <= ext[5]; ++k)
      for (int j = ext[2]; j <= ext[3]; ++j)
        for (int i = ext[0]; i <= ext[1]; ++i) out[c++] = f(i - d.c0[0] + G, j - d.c0[1] + G, k - d.c0[2] + G);
  };
  std::vector<double> tmp;
  for (auto& s : scalars) {
    extract(*s.f, tmp);
    std::uint64_t off;
    if (f32_) {
      std::vector<float> v(tmp.begin(), tmp.end());
      off = app.add(v.data(), v.size());
    } else
      off = app.add(tmp.data(), tmp.size());
    pd << "      <DataArray type=\"" << ftype << "\" Name=\"" << s.name << "\" format=\"appended\" offset=\"" << off
       << "\"/>\n";
  }
  for (auto& v : vectors) {
    std::vector<double> X, Y, Z;
    extract(*v.x, X);
    extract(*v.y, Y);
    extract(*v.z, Z);
    std::uint64_t off;
    if (f32_) {
      std::vector<float> w(3 * np);
      for (std::size_t i = 0; i < np; ++i) {
        w[3 * i] = float(X[i]);
        w[3 * i + 1] = float(Y[i]);
        w[3 * i + 2] = float(Z[i]);
      }
      off = app.add(w.data(), w.size());
    } else {
      std::vector<double> w(3 * np);
      for (std::size_t i = 0; i < np; ++i) {
        w[3 * i] = X[i];
        w[3 * i + 1] = Y[i];
        w[3 * i + 2] = Z[i];
      }
      off = app.add(w.data(), w.size());
    }
    pd << "      <DataArray type=\"" << ftype << "\" Name=\"" << v.name
       << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << off << "\"/>\n";
  }
  int whole[6] = {0, d.N[0], 0, d.N[1], 0, d.N[2]};
  std::ostringstream head;
  head << "<?xml version=\"1.0\"?>\n<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" "
          "header_type=\"UInt64\">\n";
  head << "  <ImageData WholeExtent=\"" << extent_str(whole) << "\" Origin=\"" << d.origin.x << " " << d.origin.y << " "
       << d.origin.z << "\" Spacing=\"" << d.dx[0] << " " << d.dx[1] << " " << d.dx[2] << "\">\n";
  head << "    <FieldData>\n      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">"
       << time << "</DataArray>\n    </FieldData>\n";
  head << "    <Piece Extent=\"" << extent_str(ext) << "\">\n     <PointData>\n" << pd.str()
       << "     </PointData>\n     <CellData/>\n    </Piece>\n  </ImageData>\n";
  char piece[64];
  std::snprintf(piece, sizeof(piece), "_r%05d.vti", d.rank);
  std::string piece_rel = "pieces/" + base + piece;
  write_file(dir_ + "/" + piece_rel, head.str(), app, "</VTKFile>\n");

  // master file: gather extents
  std::vector<int> all(6 * d.size);
  MPI_Gather(ext, 6, MPI_INT, all.data(), 6, MPI_INT, 0, d.comm);
  if (d.rank == 0) {
    std::ofstream f(dir_ + "/" + base + ".pvti");
    f << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PImageData\" version=\"1.0\" byte_order=\"LittleEndian\" "
         "header_type=\"UInt64\">\n";
    f << "  <PImageData WholeExtent=\"" << extent_str(whole) << "\" GhostLevel=\"0\" Origin=\"" << d.origin.x << " "
      << d.origin.y << " " << d.origin.z << "\" Spacing=\"" << d.dx[0] << " " << d.dx[1] << " " << d.dx[2] << "\">\n";
    f << "    <PPointData>\n";
    for (auto& s : scalars) f << "      <PDataArray type=\"" << ftype << "\" Name=\"" << s.name << "\"/>\n";
    for (auto& v : vectors)
      f << "      <PDataArray type=\"" << ftype << "\" Name=\"" << v.name << "\" NumberOfComponents=\"3\"/>\n";
    f << "    </PPointData>\n";
    for (int r = 0; r < d.size; ++r) {
      char pc[64];
      std::snprintf(pc, sizeof(pc), "_r%05d.vti", r);
      f << "    <Piece Extent=\"" << extent_str(&all[6 * r]) << "\" Source=\"pieces/" << base << pc << "\"/>\n";
    }
    f << "  </PImageData>\n</VTKFile>\n";
  }
  add_to_collection("fields", time, base + ".pvti");
}

namespace {
std::string poly_master(const std::string& base, int nranks, const std::vector<PolyArray>& arrays, bool cell) {
  std::ostringstream f;
  f << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
       "header_type=\"UInt64\">\n  <PPolyData GhostLevel=\"0\">\n";
  f << (cell ? "    <PCellData>\n" : "    <PPointData>\n");
  for (auto& a : arrays)
    f << "      <PDataArray type=\"Float32\" Name=\"" << a.name << "\" NumberOfComponents=\"" << a.ncomp << "\"/>\n";
  f << (cell ? "    </PCellData>\n" : "    </PPointData>\n");
  f << "    <PPoints>\n      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n    </PPoints>\n";
  for (int r = 0; r < nranks; ++r) {
    char pc[64];
    std::snprintf(pc, sizeof(pc), "_r%05d.vtp", r);
    f << "    <Piece Source=\"pieces/" << base << pc << "\"/>\n";
  }
  f << "  </PPolyData>\n</VTKFile>\n";
  return f.str();
}
}  // namespace

void VtkWriter::write_points(const std::string& series, long long step, double time, const std::vector<float>& xyz,
                             const std::vector<PolyArray>& point_data) {
  const Decomposition& d = *d_;
  const std::string base = stepname(series, step);
  const std::size_t n = xyz.size() / 3;
  Appender app;
  std::ostringstream x;
  x << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
       "header_type=\"UInt64\">\n  <PolyData>\n";
  x << "    <Piece NumberOfPoints=\"" << n << "\" NumberOfVerts=\"" << n
    << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
  x << "     <PointData>\n";
  for (auto& a : point_data) {
    auto off = app.add(a.data.data(), a.data.size());
    x << "      <DataArray type=\"Float32\" Name=\"" << a.name << "\" NumberOfComponents=\"" << a.ncomp
      << "\" format=\"appended\" offset=\"" << off << "\"/>\n";
  }
  x << "     </PointData>\n";
  auto offp = app.add(xyz.data(), xyz.size());
  x << "     <Points>\n      <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offp
    << "\"/>\n     </Points>\n";
  std::vector<std::int32_t> conn(n), offs(n);
  for (std::size_t i = 0; i < n; ++i) {
    conn[i] = std::int32_t(i);
    offs[i] = std::int32_t(i + 1);
  }
  auto oc = app.add(conn.data(), conn.size());
  auto oo = app.add(offs.data(), offs.size());
  x << "     <Verts>\n      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << oc
    << "\"/>\n      <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << oo
    << "\"/>\n     </Verts>\n";
  x << "    </Piece>\n  </PolyData>\n";
  char pc[64];
  std::snprintf(pc, sizeof(pc), "_r%05d.vtp", d.rank);
  write_file(dir_ + "/pieces/" + base + pc, x.str(), app, "</VTKFile>\n");
  if (d.rank == 0) {
    std::ofstream f(dir_ + "/" + base + ".pvtp");
    f << poly_master(base, d.size, point_data, false);
  }
  add_to_collection(series, time, base + ".pvtp");
}

void VtkWriter::write_triangles(const std::string& series, long long step, double time, const std::vector<float>& tri,
                                const std::vector<PolyArray>& cell_data) {
  const Decomposition& d = *d_;
  const std::string base = stepname(series, step);
  const std::size_t nt = tri.size() / 9;
  Appender app;
  std::ostringstream x;
  x << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" "
       "header_type=\"UInt64\">\n  <PolyData>\n";
  x << "    <Piece NumberOfPoints=\"" << 3 * nt
    << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"" << nt << "\">\n";
  x << "     <CellData>\n";
  for (auto& a : cell_data) {
    auto off = app.add(a.data.data(), a.data.size());
    x << "      <DataArray type=\"Float32\" Name=\"" << a.name << "\" NumberOfComponents=\"" << a.ncomp
      << "\" format=\"appended\" offset=\"" << off << "\"/>\n";
  }
  x << "     </CellData>\n";
  auto offp = app.add(tri.data(), tri.size());
  x << "     <Points>\n      <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offp
    << "\"/>\n     </Points>\n";
  std::vector<std::int32_t> conn(3 * nt), offs(nt);
  for (std::size_t i = 0; i < 3 * nt; ++i) conn[i] = std::int32_t(i);
  for (std::size_t i = 0; i < nt; ++i) offs[i] = std::int32_t(3 * (i + 1));
  auto oc = app.add(conn.data(), conn.size());
  auto oo = app.add(offs.data(), offs.size());
  x << "     <Polys>\n      <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << oc
    << "\"/>\n      <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << oo
    << "\"/>\n     </Polys>\n";
  x << "    </Piece>\n  </PolyData>\n";
  char pc[64];
  std::snprintf(pc, sizeof(pc), "_r%05d.vtp", d.rank);
  write_file(dir_ + "/pieces/" + base + pc, x.str(), app, "</VTKFile>\n");
  if (d.rank == 0) {
    std::ofstream f(dir_ + "/" + base + ".pvtp");
    f << poly_master(base, d.size, cell_data, true);
  }
  add_to_collection(series, time, base + ".pvtp");
}

}  // namespace picsp
