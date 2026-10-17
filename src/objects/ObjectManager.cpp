#include "objects/ObjectManager.hpp"

#include <algorithm>
#include <cmath>

#include <omp.h>

#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

void ObjectManager::setup(const Config& cfg, const Decomposition& d, const GridGeom& g, const Scene& scene,
                          FieldSolver& fs, ParticleSystem& ps, std::vector<Rng>* rngs) {
  cfg_ = cfg;
  d_ = &d;
  g_ = &g;
  scene_ = &scene;
  fs_ = &fs;
  ps_ = &ps;
  rng_ = rngs;
  objs_ = int(cfg.objects.size());
  ns_ = int(ps.species.size());
  emit_species_.assign(objs_, -1);
  for (int o = 0; o < objs_; ++o) {
    const ObjectConfig& oc = cfg.objects[o];
    int s = oc.emit_species.empty() ? ps.default_electron_species() : ps.find_species(oc.emit_species);
    bool needs = oc.see.model != SeeConfig::None || oc.photo.enabled;
    if (needs && s < 0) fatal("object '" + oc.name + "' emits electrons but no electron species '" + oc.emit_species + "' exists");
    emit_species_[o] = s;
  }
  step_coll_.assign(std::size_t(objs_) * ns_, 0.0);
  loc_coll_ = step_coll_;
  step_emit_.assign(objs_, 0.0);
  loc_emit_ = step_emit_;
  if (objs_ == 0) return;
  build_panels(scene);
  long long np = global_panel_count();
  log0(" object surface: %lld panels\n", np);
  // sunlit panels for photoemission (shadowing by any object within the domain)
  bool any_photo = false;
  for (auto& oc : cfg.objects) any_photo |= oc.photo.enabled;
  if (any_photo) {
    Vec3 lo = d.domain_lo(), hi = d.domain_hi();
    double tol = 0.05 * std::min(g.dx[0], std::min(g.dx[1], g.dx[2]));
    double lit_area = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : lit_area)
    for (std::size_t i = 0; i < panels_.size(); ++i) {
      Panel& p = panels_[i];
      if (!cfg.objects[p.obj].photo.enabled) continue;
      double c = dot(p.normal, cfg.sun_dir);
      if (c <= 0) continue;
      if (scene.ray_blocked(p.centroid + p.normal * (2 * tol), cfg.sun_dir, lo, hi, tol)) continue;
      p.sun_cos = c;
      lit_area += c * p.area;
    }
    lit_area = allreduce_sum(lit_area, d.comm);
    log0(" photoemission: projected sunlit area %.4g m^2\n", lit_area);
  }
  reset_panel_accumulators();
}

long long ObjectManager::global_panel_count() const { return allreduce_sum((long long)panels_.size(), d_->comm); }

void ObjectManager::reset_panel_accumulators() {
  panel_q_.assign(panels_.size() * ns_, 0.0);
  panel_emit_.assign(panels_.size(), 0.0);
  accum_time_ = 0;
}

void ObjectManager::build_panels(const Scene& scene) {
  const Decomposition& d = *d_;
  const GridGeom& g = *g_;
  const EmbeddedBoundary& eb = fs_->eb();
  const int nx = d.c1[0] - d.c0[0], ny = d.c1[1] - d.c0[1], nz = d.c1[2] - d.c0[2];
  double hmax = std::max(g.dx[0], std::max(g.dx[1], g.dx[2]));
  const int r = std::clamp(cfg_.surface_refine, 1, 8);
  // Extract the iso-surface sdf = delta (slightly outside the body) so that features
  // thinner than the sampling spacing are still captured, then project the vertices
  // back onto the true surface (sdf = 0).
  const double delta = 0.55 * hmax / r;
  static const int tets[6][4] = {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
  const int nthreads = omp_get_max_threads();
  std::vector<std::vector<Panel>> tp(nthreads);
  std::vector<std::vector<int>> tcell(nthreads);
  const double cap = 3 * hmax;
#pragma omp parallel for collapse(2) schedule(dynamic, 2)
  for (int ck = 0; ck < nz; ++ck)
    for (int cj = 0; cj < ny; ++cj) {
      int tid = omp_get_thread_num();
      std::vector<double> sv((r + 1) * (r + 1) * (r + 1));
      std::vector<Vec3> pv(sv.size());
      for (int ci = 0; ci < nx; ++ci) {
        if (!eb.cell_near(ci, cj, ck)) continue;
        Vec3 c0{g.origin.x + (d.c0[0] + ci) * g.dx[0], g.origin.y + (d.c0[1] + cj) * g.dx[1],
                g.origin.z + (d.c0[2] + ck) * g.dx[2]};
        for (int k = 0; k <= r; ++k)
          for (int j = 0; j <= r; ++j)
            for (int i = 0; i <= r; ++i) {
              Vec3 p = c0 + Vec3{g.dx[0] * i / r, g.dx[1] * j / r, g.dx[2] * k / r};
              int q = i + (r + 1) * (j + (r + 1) * k);
              pv[q] = p;
              sv[q] = scene.sdf(p, nullptr, cap);
            }
        for (int k = 0; k < r; ++k)
          for (int j = 0; j < r; ++j)
            for (int i = 0; i < r; ++i) {
              int corner[8];
              for (int c = 0; c < 8; ++c)
                corner[c] = (i + (c & 1)) + (r + 1) * ((j + ((c >> 1) & 1)) + (r + 1) * (k + ((c >> 2) & 1)));
              bool anyin = false, anyout = false;
              for (int c = 0; c < 8; ++c) (sv[corner[c]] < delta ? anyin : anyout) = true;
              if (!(anyin && anyout)) continue;
              for (auto& t : tets) {
                int in[4], out[4], ni = 0, no = 0;
                for (int v = 0; v < 4; ++v) {
                  int q = corner[t[v]];
                  if (sv[q] < delta) in[ni++] = q;
                  else out[no++] = q;
                }
                if (ni == 0 || no == 0) continue;
                auto ip = [&](int a, int b) {
                  double sa = sv[a] - delta, sb = sv[b] - delta;
                  double f = sa / (sa - sb);
                  Vec3 q = pv[a] + (pv[b] - pv[a]) * f;
                  // project onto the zero level set (two Newton steps along the sdf gradient)
                  for (int it = 0; it < 2; ++it) {
                    const double e = 0.05 * hmax / r;
                    Vec3 gr{scene.sdf(q + Vec3{e, 0, 0}) - scene.sdf(q - Vec3{e, 0, 0}),
                            scene.sdf(q + Vec3{0, e, 0}) - scene.sdf(q - Vec3{0, e, 0}),
                            scene.sdf(q + Vec3{0, 0, e}) - scene.sdf(q - Vec3{0, 0, e})};
                    double gn = norm(gr);
                    if (gn <= 0) break;
                    q = q - gr * (scene.sdf(q) / gn);
                  }
                  return q;
                };
                Vec3 gin{0, 0, 0}, gout{0, 0, 0};
                for (int v = 0; v < ni; ++v) gin += pv[in[v]];
                for (int v = 0; v < no; ++v) gout += pv[out[v]];
                Vec3 outward = gout / no - gin / ni;
                auto add = [&](Vec3 A, Vec3 B, Vec3 C) {
                  Vec3 nrm = cross(B - A, C - A);
                  double ar = 0.5 * norm(nrm);
                  if (ar <= 0) return;
                  if (dot(nrm, outward) < 0) {
                    std::swap(B, C);
                    nrm = -nrm;
                  }
                  Panel P;
                  P.a = A;
                  P.b = B;
                  P.c = C;
                  P.centroid = (A + B + C) / 3.0;
                  P.normal = nrm / (2 * ar);
                  P.area = ar;
                  scene.sdf(P.centroid, &P.obj, cap);
                  tp[tid].push_back(P);
                  tcell[tid].push_back(ci + nx * (cj + ny * ck));
                };
                if (ni == 1) add(ip(in[0], out[0]), ip(in[0], out[1]), ip(in[0], out[2]));
                else if (ni == 3) add(ip(out[0], in[0]), ip(out[0], in[1]), ip(out[0], in[2]));
                else {
                  Vec3 p0 = ip(in[0], out[0]), p1 = ip(in[0], out[1]), p2 = ip(in[1], out[1]), p3 = ip(in[1], out[0]);
                  add(p0, p1, p2);
                  add(p0, p2, p3);
                }
              }
            }
      }
    }
  panels_.clear();
  std::vector<int> pcell;
  for (int t = 0; t < nthreads; ++t) {
    panels_.insert(panels_.end(), tp[t].begin(), tp[t].end());
    pcell.insert(pcell.end(), tcell[t].begin(), tcell[t].end());
  }
  // CSR cell -> panels
  const std::size_t ncell = std::size_t(nx) * ny * nz;
  cell_start_.assign(ncell + 1, 0);
  for (int c : pcell) cell_start_[c + 1]++;
  for (std::size_t c = 0; c < ncell; ++c) cell_start_[c + 1] += cell_start_[c];
  cell_panels_.assign(panels_.size(), 0);
  std::vector<int> fill(cell_start_.begin(), cell_start_.end() - 1);
  for (std::size_t p = 0; p < pcell.size(); ++p) cell_panels_[fill[pcell[p]]++] = int(p);
  panel_cell_ = pcell;
}

std::vector<double> ObjectManager::cell_averaged_density(const std::vector<double>& q, int ncomp, int comp) const {
  std::vector<double> out(panels_.size(), 0.0);
  const double T = accum_time_ > 0 ? accum_time_ : 1.0;
  const std::size_t ncell = cell_start_.empty() ? 0 : cell_start_.size() - 1;
  for (std::size_t c = 0; c < ncell; ++c) {
    // group panels by object within the cell
    for (int a = cell_start_[c]; a < cell_start_[c + 1]; ++a) {
      int pa = cell_panels_[a];
      double sq = 0, sa = 0;
      for (int b = cell_start_[c]; b < cell_start_[c + 1]; ++b) {
        int pb = cell_panels_[b];
        if (panels_[pb].obj != panels_[pa].obj) continue;
        sq += q[std::size_t(pb) * ncomp + comp];
        sa += panels_[pb].area;
      }
      out[pa] = sa > 0 ? sq / (sa * T) : 0.0;
    }
  }
  return out;
}

int ObjectManager::find_panel(const Vec3& p) const {
  if (panels_.empty()) return -1;
  const Decomposition& d = *d_;
  const GridGeom& g = *g_;
  const int nx = d.c1[0] - d.c0[0], ny = d.c1[1] - d.c0[1], nz = d.c1[2] - d.c0[2];
  int ci = int(std::floor((p.x - g.origin.x) / g.dx[0])) - d.c0[0];
  int cj = int(std::floor((p.y - g.origin.y) / g.dx[1])) - d.c0[1];
  int ck = int(std::floor((p.z - g.origin.z) / g.dx[2])) - d.c0[2];
  int best = -1;
  double bd = 1e300;
  for (int dk = -1; dk <= 1; ++dk)
    for (int dj = -1; dj <= 1; ++dj)
      for (int di = -1; di <= 1; ++di) {
        int a = ci + di, b = cj + dj, c = ck + dk;
        if (a < 0 || b < 0 || c < 0 || a >= nx || b >= ny || c >= nz) continue;
        std::size_t cell = a + std::size_t(nx) * (b + std::size_t(ny) * c);
        for (int q = cell_start_[cell]; q < cell_start_[cell + 1]; ++q) {
          int pi = cell_panels_[q];
          Vec3 dd = panels_[pi].centroid - p;
          double d2 = dot(dd, dd);
          if (d2 < bd) {
            bd = d2;
            best = pi;
          }
        }
      }
  return best;
}

double ObjectManager::see_yield(const ObjectConfig& oc, int s, double E) const {
  const Species& sp = ps_->species[s];
  if (sp.q > 0) return oc.see.ion_yield;
  if (sp.q < 0 && sp.m < 100 * constants::me) {
    if (oc.see.model == SeeConfig::Constant) return oc.see.yield;
    if (oc.see.model == SeeConfig::Sternglass) {
      double x = E / oc.see.e_max;
      return 7.4 * oc.see.delta_max * x * std::exp(-2.0 * std::sqrt(x));
    }
  }
  return 0.0;
}

void ObjectManager::emit(int obj, int s, const Vec3& pos, const Vec3& n, double vth, double w, int panel) {
  Rng& rng = (*rng_)[0];
  const Species& sp = ps_->species[s];
  Vec3 v = sample_emission_velocity(rng, n, vth);
  double h = std::min(g_->dx[0], std::min(g_->dx[1], g_->dx[2]));
  Vec3 p = pos + n * (0.02 * h) + v * (ps_->dt_species(s) * rng.uniform());
  // keep inside the global domain
  Vec3 lo = d_->domain_lo(), hi = d_->domain_hi();
  for (int a = 0; a < 3; ++a)
    if (p[a] < lo[a] || p[a] >= hi[a]) return;
  if (ps_->inside_object(p)) return;
  ps_->add_particle(s, p, v, w);
  double q = sp.q * w;  // charge leaving the object
  int cnd = fs_->object_conductor(obj);
  fs_->conductors()[cnd].Q -= q;
  loc_emit_[obj] -= q;
  if (panel >= 0) panel_emit_[panel] -= q;
}

void ObjectManager::process_hits(long long step, std::vector<std::vector<Hit>>& hits) {
  using constants::qe;
  auto& conds = fs_->conductors();
  const double hmin = std::min(g_->dx[0], std::min(g_->dx[1], g_->dx[2]));
  for (int s = 0; s < ns_; ++s) {
    if (hits[s].empty()) continue;
    const Species& sp = ps_->species[s];
    // copy hit data first: removal reorders the particle array
    std::vector<Hit> hs = hits[s];
    ps_->remove_hits(s, hits[s]);
    for (const Hit& h : hs) {
      double q = sp.q * h.w;
      int cnd = fs_->object_conductor(h.obj);
      conds[cnd].Q += q;
      loc_coll_[std::size_t(h.obj) * ns_ + s] += q;
      int panel = find_panel(h.pos);
      if (panel >= 0) panel_q_[std::size_t(panel) * ns_ + s] += q;
      const ObjectConfig& oc = cfg_.objects[h.obj];
      if (oc.see.model == SeeConfig::None && oc.see.ion_yield <= 0) continue;
      double E = 0.5 * sp.m * dot(h.vel, h.vel) / qe;
      double y = see_yield(oc, s, E);
      if (y <= 0) continue;
      int es = emit_species_[h.obj];
      double w_emit = ps_->species[es].w0;
      double expected = y * h.w / w_emit;
      int n = int(expected);
      if ((*rng_)[0].uniform() < expected - n) ++n;
      if (n == 0) continue;
      Vec3 nrm = normalized(scene_->object(h.obj).gradient(h.pos, 0.25 * hmin));
      double vth = std::sqrt(oc.see.temperature * qe / ps_->species[es].m);
      for (int k = 0; k < n; ++k) emit(h.obj, es, h.pos, nrm, vth, w_emit, panel);
    }
  }
}

void ObjectManager::photoemission(long long step) {
  using constants::qe;
  if (panels_.empty()) return;
  Rng& rng = (*rng_)[0];
  for (std::size_t i = 0; i < panels_.size(); ++i) {
    const Panel& P = panels_[i];
    if (P.sun_cos <= 0) continue;
    const ObjectConfig& oc = cfg_.objects[P.obj];
    int es = emit_species_[P.obj];
    if (!ps_->species[es].active(step)) continue;
    double w = ps_->species[es].w0;
    double expected = oc.photo.current_density * P.sun_cos * P.area * ps_->dt_species(es) / (qe * w);
    int n = int(expected);
    if (rng.uniform() < expected - n) ++n;
    double vth = std::sqrt(oc.photo.temperature * qe / ps_->species[es].m);
    for (int k = 0; k < n; ++k) {
      double u = rng.uniform(), v = rng.uniform();
      if (u + v > 1) {
        u = 1 - u;
        v = 1 - v;
      }
      Vec3 pos = P.a + (P.b - P.a) * u + (P.c - P.a) * v;
      emit(P.obj, es, pos, P.normal, vth, w, int(i));
    }
  }
}

void ObjectManager::end_step() {
  auto& conds = fs_->conductors();
  for (auto& c : conds) c.Q += c.Iext * cfg_.dt;
  // conductor charges must be identical on all ranks: reduce this step's increments
  if (objs_ > 0) {
    MPI_Allreduce(loc_coll_.data(), step_coll_.data(), int(loc_coll_.size()), MPI_DOUBLE, MPI_SUM, d_->comm);
    MPI_Allreduce(loc_emit_.data(), step_emit_.data(), int(loc_emit_.size()), MPI_DOUBLE, MPI_SUM, d_->comm);
    // Q was incremented locally; replace the local increments by the global sums
    for (int o = 0; o < objs_; ++o) {
      int c = fs_->object_conductor(o);
      double loc = loc_emit_[o], glob = step_emit_[o];
      for (int s = 0; s < ns_; ++s) {
        loc += loc_coll_[std::size_t(o) * ns_ + s];
        glob += step_coll_[std::size_t(o) * ns_ + s];
      }
      conds[c].Q += glob - loc;
    }
    std::fill(loc_coll_.begin(), loc_coll_.end(), 0.0);
    std::fill(loc_emit_.begin(), loc_emit_.end(), 0.0);
  }
  accum_time_ += cfg_.dt;
}

}  // namespace picsp
