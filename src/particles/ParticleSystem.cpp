#include "particles/ParticleSystem.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

#include <omp.h>

#include "particles/Shape.hpp"
#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

void ParticleSystem::setup(const Config& cfg, const Decomposition& d, const GridGeom& g, const Scene& scene,
                           const FieldSolver& fs, std::vector<Rng>* rngs) {
  cfg_ = cfg;
  d_ = &d;
  g_ = &g;
  scene_ = &scene;
  fs_ = &fs;
  rng_ = rngs;
  dt_ = cfg.dt;
  G_ = FieldSolver::kGhosts;
  for (int a = 0; a < 3; ++a) {
    arr_lo_[a] = g.origin[a] + (g.c0[a] - G_) * g.dx[a];
    inv_dx_[a] = 1.0 / g.dx[a];
  }
  species.clear();
  for (std::size_t s = 0; s < cfg.species.size(); ++s) species.emplace_back(cfg.species[s], int(s), cfg.B0);

  // nominal macro weights
  const double dV = g.dx[0] * g.dx[1] * g.dx[2];
  double ref_w = 0;
  for (auto& sp : species) {
    if (sp.cfg.weight > 0) sp.w0 = sp.cfg.weight;
    else if (sp.cfg.density > 0) sp.w0 = sp.cfg.density * dV / sp.cfg.ppc;
    if (ref_w == 0 && sp.w0 > 0 && sp.q < 0) ref_w = sp.w0;
  }
  if (ref_w == 0)
    for (auto& sp : species)
      if (ref_w == 0 && sp.w0 > 0) ref_w = sp.w0;
  for (auto& sp : species)
    if (sp.w0 <= 0) sp.w0 = ref_w > 0 ? ref_w : 1.0;

  // injectors on open faces touching this rank
  inj_.clear();
  for (int a = 0; a < 3; ++a)
    for (int side = 0; side < 2; ++side) {
      if (cfg.faces[2 * a + side].particle != ParticleBC::Open) continue;
      bool mine = side == 0 ? d.wall_lo(a) : d.wall_hi(a);
      for (auto& sp : species) {
        if (!sp.cfg.inject || sp.cfg.density <= 0) continue;
        // Monte-Carlo estimate of the one-sided flux (identical on all ranks)
        Rng r(cfg.seed + 7777, 1000 + 10 * sp.id + 2 * a + side);
        Vec3 nin{0, 0, 0};
        nin[a] = side == 0 ? 1.0 : -1.0;
        const int Ns = 200000;
        double sum = 0, vmax = 0;
        for (int k = 0; k < Ns; ++k) {
          double vn = dot(sp.sampler.sample(r), nin);
          if (vn > 0) {
            sum += vn;
            vmax = std::max(vmax, vn);
          }
        }
        FaceInjector fi;
        fi.axis = a;
        fi.side = side;
        fi.species = sp.id;
        fi.flux = sp.cfg.density * sum / Ns;
        fi.vcap = vmax;
        if (sp.cfg.dist == VelocityDist::Maxwellian && sp.sampler.vth() > 0) {
          // exact: flux-weighted drifting Maxwellian, f(v_n) ~ v_n exp(-(v_n - u)^2 / 2 vth^2), v_n > 0
          const double vt = sp.sampler.vth(), u = dot(sp.cfg.drift, nin);
          fi.maxwellian = true;
          fi.vth = vt;
          fi.drift = sp.cfg.drift;
          fi.flux = sp.cfg.density * (vt / std::sqrt(constants::twopi) * std::exp(-u * u / (2 * vt * vt)) +
                                      0.5 * u * (1 + std::erf(u / (std::sqrt(2.0) * vt))));
          const int M = 20000, Q = 4096;
          const double vmx = std::max(u, 0.0) + 10 * vt;
          std::vector<double> cdf(M + 1, 0.0);
          auto pdf = [&](double v) { return v * std::exp(-(v - u) * (v - u) / (2 * vt * vt)); };
          for (int k = 1; k <= M; ++k) {
            double v0 = vmx * (k - 1) / M, v1 = vmx * k / M;
            cdf[k] = cdf[k - 1] + 0.5 * (pdf(v0) + pdf(v1)) * (v1 - v0);
          }
          fi.icdf.resize(Q + 1);
          int k = 0;
          for (int q = 0; q <= Q; ++q) {
            double target = cdf[M] * q / Q;
            while (k < M && cdf[k + 1] < target) ++k;
            double c0 = cdf[k], c1 = cdf[std::min(k + 1, M)];
            double f = c1 > c0 ? (target - c0) / (c1 - c0) : 0.0;
            fi.icdf[q] = vmx * (k + std::clamp(f, 0.0, 1.0)) / M;
          }
        }
        int b = (a + 1) % 3, c = (a + 2) % 3;
        fi.area = (d.c1[b] - d.c0[b]) * g.dx[b] * (d.c1[c] - d.c0[c]) * g.dx[c];
        if (mine && fi.flux > 0) inj_.push_back(fi);
      }
    }

  int nt = omp_get_max_threads();
  tbuf_.assign(nt, std::vector<double>());
}

bool ParticleSystem::inside_object(const Vec3& p) const {
  if (scene_->empty()) return false;
  const Decomposition& d = *d_;
  int c[3];
  for (int a = 0; a < 3; ++a) c[a] = int(std::floor((p[a] - g_->origin[a]) * inv_dx_[a])) - d.c0[a];
  bool in_range = true;
  for (int a = 0; a < 3; ++a)
    if (c[a] < -1 || c[a] > d.c1[a] - d.c0[a]) in_range = false;
  if (in_range && !fs_->eb().cell_near(c[0], c[1], c[2])) return false;
  return scene_->sdf(p) <= 0;
}

int ParticleSystem::find_species(const std::string& name) const {
  for (auto& s : species)
    if (s.name() == name) return s.id;
  return -1;
}

int ParticleSystem::default_electron_species() const {
  for (auto& s : species)
    if (s.q < 0 && std::fabs(s.m - constants::me) < 0.01 * constants::me) return s.id;
  for (auto& s : species)
    if (s.q < 0) return s.id;
  return -1;
}

// ------------------------------------------------------------------------------------
// Loading
// ------------------------------------------------------------------------------------
void ParticleSystem::load_initial() {
  const Decomposition& d = *d_;
  Vec3 llo = d.local_lo(), lhi = d.local_hi();
  Vec3 dlo = d.domain_lo(), dhi = d.domain_hi();
  const EmbeddedBoundary& eb = fs_->eb();
  for (auto& sp : species) {
    if (!sp.cfg.load_initial || sp.cfg.density <= 0) continue;
    Vec3 lo = llo, hi = lhi;
    if (sp.cfg.region) {
      lo = vmax(lo, (*sp.cfg.region)[0]);
      hi = vmin(hi, (*sp.cfg.region)[1]);
    }
    double vol = 1;
    for (int a = 0; a < 3; ++a) vol *= std::max(0.0, hi[a] - lo[a]);
    if (vol <= 0) continue;
    const double dV = d.dx[0] * d.dx[1] * d.dx[2];
    double expected = sp.cfg.density * vol / sp.w0;
    long long npart = (long long)std::floor(expected);
    Rng& rng0 = (*rng_)[0];
    if (rng0.uniform() < expected - npart) ++npart;
    (void)dV;
    // perturbation wavenumber
    const auto& pert = sp.cfg.perturbation;
    int pa = pert.axis;
    double Lp = dhi[pa] - dlo[pa];
    double kp = constants::twopi * pert.mode / Lp;
    sp.p.reserve(std::size_t(npart * 1.2) + 16);
    // quasi-random sequence offset so ranks do not share points
    unsigned long long qoff = 1 + (unsigned long long)d.rank * (unsigned long long)(npart + 1000);
    std::vector<Vec3> pos(npart), vel(npart);
    std::vector<char> ok(npart, 1);
#pragma omp parallel
    {
      Rng& rng = (*rng_)[omp_get_thread_num()];
#pragma omp for schedule(static)
      for (long long n = 0; n < npart; ++n) {
        Vec3 p;
        if (sp.cfg.load == LoadKind::Quiet) {
          unsigned long long q = qoff + n;
          p = {lo.x + (hi.x - lo.x) * radical_inverse(q, 2), lo.y + (hi.y - lo.y) * radical_inverse(q, 3),
               lo.z + (hi.z - lo.z) * radical_inverse(q, 5)};
        } else {
          p = {lo.x + (hi.x - lo.x) * rng.uniform(), lo.y + (hi.y - lo.y) * rng.uniform(),
               lo.z + (hi.z - lo.z) * rng.uniform()};
        }
        Vec3 v = sp.sampler.sample(rng);
        if (pert.amplitude != 0) {
          double xa = p[pa] - dlo[pa];
          p[pa] -= pert.amplitude / kp * std::sin(kp * xa);
        }
        if (pert.velocity_amplitude != 0) v[pa] += pert.velocity_amplitude * std::sin(kp * (p[pa] - dlo[pa]));
        // keep inside the global domain (periodic wrap / clamp)
        for (int a = 0; a < 3; ++a) {
          double L = dhi[a] - dlo[a];
          if (d.periodic[a]) {
            while (p[a] < dlo[a]) p[a] += L;
            while (p[a] >= dhi[a]) p[a] -= L;
          } else {
            p[a] = std::clamp(p[a], dlo[a], std::nextafter(dhi[a], dlo[a]));
          }
        }
        if (!scene_->empty()) {
          int ci = int(std::floor((p.x - llo.x) * inv_dx_[0])), cj = int(std::floor((p.y - llo.y) * inv_dx_[1])),
              ck = int(std::floor((p.z - llo.z) * inv_dx_[2]));
          bool near = true;
          if (ci >= -1 && cj >= -1 && ck >= -1 && ci <= d.c1[0] - d.c0[0] && cj <= d.c1[1] - d.c0[1] && ck <= d.c1[2] - d.c0[2])
            near = eb.cell_near(ci, cj, ck);
          if (near && scene_->sdf(p) <= 0) ok[n] = 0;
        }
        pos[n] = p;
        vel[n] = v;
      }
    }
    for (long long n = 0; n < npart; ++n)
      if (ok[n]) sp.p.push(pos[n], vel[n], sp.w0);
  }
  // perturbations may move particles across rank boundaries
  boundaries_and_migrate(0);
  for (auto& sp : species) {
    long long n = count_global(sp.id);
    if (n > 0) log0(" loaded %-10s %12lld macro-particles (weight %.4g)\n", sp.name().c_str(), n, sp.w0);
  }
}

long long ParticleSystem::count_global(int s) const {
  return allreduce_sum((long long)species[s].p.size(), d_->comm);
}

double ParticleSystem::kinetic_energy(int s) const {
  const ParticleArray& p = species[s].p;
  double e = 0;
  const std::size_t n = p.size();
#pragma omp parallel for reduction(+ : e) schedule(static)
  for (std::size_t i = 0; i < n; ++i) e += p.w[i] * (p.vx[i] * p.vx[i] + p.vy[i] * p.vy[i] + p.vz[i] * p.vz[i]);
  return 0.5 * species[s].m * allreduce_sum(e, d_->comm);
}

// ------------------------------------------------------------------------------------
// Push (Boris) + object collision detection
// ------------------------------------------------------------------------------------
template <int O>
void ParticleSystem::push_impl(Species& sp, double dt, std::vector<Hit>& hits) {
  using S = Shape<O>;
  ParticleArray& p = sp.p;
  const std::size_t np = p.size();
  const auto& n = fs_->Ex.size();
  const std::size_t sy = n[0], sz = std::size_t(n[0]) * n[1];
  const double* E3 = fs_->E3.data();
  const double ox = arr_lo_[0], oy = arr_lo_[1], oz = arr_lo_[2];
  const double ix = inv_dx_[0], iy = inv_dx_[1], iz = inv_dx_[2];
  const double h = 0.5 * sp.qm * dt;
  const double tx = h * cfg_.B0.x, ty = h * cfg_.B0.y, tz = h * cfg_.B0.z;
  const double t2 = tx * tx + ty * ty + tz * tz;
  const double sf = 2.0 / (1.0 + t2);
  const double sx = sf * tx, syy = sf * ty, szz = sf * tz;
  double* X = p.x.data();
  double* Y = p.y.data();
  double* Z = p.z.data();
  double* VX = p.vx.data();
  double* VY = p.vy.data();
  double* VZ = p.vz.data();
  const bool objects = !scene_->empty();
  const EmbeddedBoundary& eb = fs_->eb();
  const Decomposition& d = *d_;
  const int ncx = d.c1[0] - d.c0[0], ncy = d.c1[1] - d.c0[1], ncz = d.c1[2] - d.c0[2];
  const int G = G_;
  // is the cell containing array-local coordinate (lx,ly,lz) near an object surface?
  auto near = [&](double lx, double ly, double lz) {
    int ci = int(lx) - G, cj = int(ly) - G, ck = int(lz) - G;  // lx > 0 in array coordinates
    if (ci < -1 || cj < -1 || ck < -1 || ci > ncx || cj > ncy || ck > ncz) return true;
    return eb.cell_near(ci, cj, ck);
  };
  const int nt = omp_get_max_threads();
  std::vector<std::vector<std::size_t>> cand(nt);

#pragma omp parallel
  {
    auto& mycand = cand[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (std::size_t i = 0; i < np; ++i) {
      const double lx = (X[i] - ox) * ix, ly = (Y[i] - oy) * iy, lz = (Z[i] - oz) * iz;
      int i0, j0, k0;
      double wx[S::n], wy[S::n], wz[S::n];
      S::weights(lx, i0, wx);
      S::weights(ly, j0, wy);
      S::weights(lz, k0, wz);
      double ex = 0, ey = 0, ez = 0;
#pragma GCC unroll 3
      for (int c = 0; c < S::n; ++c)
#pragma GCC unroll 3
        for (int b = 0; b < S::n; ++b) {
          const double* e = E3 + 3 * (std::size_t(i0) + sy * std::size_t(j0 + b) + sz * std::size_t(k0 + c));
          const double wbc = wy[b] * wz[c];
#pragma GCC unroll 3
          for (int a = 0; a < S::n; ++a) {
            const double ww = wx[a] * wbc;
            ex += ww * e[3 * a];
            ey += ww * e[3 * a + 1];
            ez += ww * e[3 * a + 2];
          }
        }
      double vx = VX[i] + h * ex, vy = VY[i] + h * ey, vz = VZ[i] + h * ez;
      // Boris rotation (t = 0 reduces to the identity)
      double px = vx + (vy * tz - vz * ty);
      double py = vy + (vz * tx - vx * tz);
      double pz = vz + (vx * ty - vy * tx);
      vx += py * szz - pz * syy;
      vy += pz * sx - px * szz;
      vz += px * syy - py * sx;
      vx += h * ex;
      vy += h * ey;
      vz += h * ez;
      VX[i] = vx;
      VY[i] = vy;
      VZ[i] = vz;
      X[i] += vx * dt;
      Y[i] += vy * dt;
      Z[i] += vz * dt;
      if (objects && (near(lx, ly, lz) || near(lx + vx * dt * ix, ly + vy * dt * iy, lz + vz * dt * iz)))
        mycand.push_back(i);
    }
  }
  if (!objects) return;
  // exact segment / surface intersection for particles near a surface
  const double tol = 1e-4 * std::min(g_->dx[0], std::min(g_->dx[1], g_->dx[2]));
  std::vector<std::size_t> all;
  for (auto& c : cand) all.insert(all.end(), c.begin(), c.end());
  std::vector<std::vector<Hit>> th(nt);
#pragma omp parallel
  {
    auto& my = th[omp_get_thread_num()];
#pragma omp for schedule(dynamic, 64)
    for (std::size_t k = 0; k < all.size(); ++k) {
      std::size_t i = all[k];
      Vec3 b{X[i], Y[i], Z[i]};
      Vec3 a{b.x - VX[i] * dt, b.y - VY[i] * dt, b.z - VZ[i] * dt};
      double t;
      int o;
      if (scene_->segment_hit(a, b, tol, &t, &o)) my.push_back({i, o, a + (b - a) * t, {VX[i], VY[i], VZ[i]}, p.w[i]});
    }
  }
  for (auto& v : th) hits.insert(hits.end(), v.begin(), v.end());
}

void ParticleSystem::push(long long step, std::vector<std::vector<Hit>>& hits) {
  hits.assign(species.size(), {});
  for (auto& sp : species) {
    if (!sp.active(step)) continue;
    double dt = dt_ * sp.cfg.subcycle;
    if (cfg_.shape_order == 2) push_impl<2>(sp, dt, hits[sp.id]);
    else push_impl<1>(sp, dt, hits[sp.id]);
  }
}

void ParticleSystem::remove_hits(int s, std::vector<Hit>& hits) {
  std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.idx > b.idx; });
  for (auto& h : hits) species[s].p.remove(h.idx);
}

void ParticleSystem::rewind_half() {
  // v(-dt/2) = v(0) - (q/m)(E + v x B) dt/2, using a Boris step with -dt/2
  for (auto& sp : species) {
    double dt = -0.5 * dt_ * sp.cfg.subcycle;
    std::vector<double> x = sp.p.x, y = sp.p.y, z = sp.p.z;
    std::vector<Hit> dummy;
    const Scene* sc = scene_;
    Scene empty;
    scene_ = &empty;  // no collision detection for the rewind
    if (cfg_.shape_order == 2) push_impl<2>(sp, dt, dummy);
    else push_impl<1>(sp, dt, dummy);
    scene_ = sc;
    sp.p.x = std::move(x);
    sp.p.y = std::move(y);
    sp.p.z = std::move(z);
  }
}

// ------------------------------------------------------------------------------------
// Injection through open boundaries
// ------------------------------------------------------------------------------------
void ParticleSystem::inject(long long step) {
  const Decomposition& d = *d_;
  Vec3 llo = d.local_lo(), lhi = d.local_hi();
  Rng& rng = (*rng_)[0];
  for (auto& fi : inj_) {
    Species& sp = species[fi.species];
    sp.stats.injected = 0;
    if (!sp.active(step)) continue;
    double dt = dt_ * sp.cfg.subcycle;
    double expected = fi.flux * fi.area * dt / sp.w0;
    long long n = (long long)std::floor(expected);
    if (rng.uniform() < expected - n) ++n;
    Vec3 nin{0, 0, 0};
    nin[fi.axis] = fi.side == 0 ? 1.0 : -1.0;
    double wall = fi.side == 0 ? llo[fi.axis] : lhi[fi.axis];
    for (long long k = 0; k < n; ++k) {
      Vec3 v;
      if (fi.maxwellian) {
        double q = rng.uniform() * (fi.icdf.size() - 1);
        int qi = int(q);
        double f = q - qi;
        double vn = fi.icdf[qi] * (1 - f) + fi.icdf[std::min<std::size_t>(qi + 1, fi.icdf.size() - 1)] * f;
        v = fi.drift + Vec3{rng.normal(), rng.normal(), rng.normal()} * fi.vth;
        v[fi.axis] = vn * nin[fi.axis];
      } else {
        for (int tries = 0; tries < 10000; ++tries) {
          v = sp.sampler.sample(rng);
          double vn = dot(v, nin);
          if (vn <= 0) continue;
          if (vn >= fi.vcap || rng.uniform() * fi.vcap < vn) break;
        }
      }
      Vec3 p;
      for (int a = 0; a < 3; ++a) p[a] = llo[a] + (lhi[a] - llo[a]) * rng.uniform();
      p[fi.axis] = wall;
      p += v * (dt * rng.uniform());
      // nudge off the wall
      if (fi.side == 0) p[fi.axis] = std::max(p[fi.axis], wall + 1e-9 * g_->dx[fi.axis]);
      else p[fi.axis] = std::min(p[fi.axis], wall - 1e-9 * g_->dx[fi.axis]);
      if (inside_object(p)) continue;
      sp.p.push(p, v, sp.w0);
      sp.stats.injected += 1;
    }
  }
}

// ------------------------------------------------------------------------------------
// Outer boundaries and MPI migration (dimension by dimension)
// ------------------------------------------------------------------------------------
void ParticleSystem::boundaries_and_migrate(long long step) {
  const Decomposition& d = *d_;
  const Vec3 llo = d.local_lo(), lhi = d.local_hi();
  const Vec3 dlo = d.domain_lo(), dhi = d.domain_hi();
  const int K = ParticleArray::kPacked;
  for (auto& sp : species) sp.stats.absorbed_wall_charge = 0;
  const int nt = omp_get_max_threads();
  std::vector<std::vector<std::size_t>> tl(nt);
  for (auto& sp : species) {
    if (!sp.active(step) && step > 0) continue;
    ParticleArray& p = sp.p;
    // (1) one parallel pass: indices of particles outside the local box (any axis)
    const std::size_t np = p.size();
    const double *X = p.x.data(), *Y = p.y.data(), *Z = p.z.data();
#pragma omp parallel
    {
      auto& my = tl[omp_get_thread_num()];
      my.clear();
#pragma omp for schedule(static) nowait
      for (std::size_t i = 0; i < np; ++i)
        if (X[i] < llo.x || X[i] >= lhi.x || Y[i] < llo.y || Y[i] >= lhi.y || Z[i] < llo.z || Z[i] >= lhi.z)
          my.push_back(i);
    }
    std::vector<std::size_t> out;
    for (auto& v : tl) out.insert(out.end(), v.begin(), v.end());
    // (2) move all leavers into a packed buffer (descending: swap-removal keeps indices valid)
    std::sort(out.begin(), out.end(), std::greater<std::size_t>());
    std::vector<double> buf(out.size() * K);
    for (std::size_t k = 0; k < out.size(); ++k) {
      p.pack(out[k], &buf[k * K]);
      p.remove(out[k]);
    }
    // (3) per axis: walls, periodic wrap and neighbour exchange on the buffer only
    double lost = 0;
    for (int a = 0; a < 3; ++a) {
      const double L = dhi[a] - dlo[a];
      const ParticleBC bc[2] = {cfg_.faces[2 * a].particle, cfg_.faces[2 * a + 1].particle};
      sendbuf_[0].clear();
      sendbuf_[1].clear();
      std::vector<double> keep;
      keep.reserve(buf.size());
      for (std::size_t k = 0; k < buf.size(); k += K) {
        double* q = &buf[k];
        double xa = q[a];
        int side = xa < llo[a] ? 0 : (xa >= lhi[a] ? 1 : -1);
        if (side < 0) {
          keep.insert(keep.end(), q, q + K);
          continue;
        }
        bool wall = side == 0 ? d.wall_lo(a) : d.wall_hi(a);
        if (wall) {
          if (bc[side] == ParticleBC::Reflect) {
            double r = side == 0 ? 2 * dlo[a] - xa : 2 * dhi[a] - xa;
            q[a] = std::clamp(r, dlo[a], std::nextafter(dhi[a], dlo[a]));
            q[3 + a] = -q[3 + a];
            keep.insert(keep.end(), q, q + K);
          } else {
            lost += q[6];
          }
          continue;
        }
        if (d.periodic[a]) {
          if (side == 0 && d.coords[a] == 0) q[a] += L;
          if (side == 1 && d.coords[a] == d.dims[a] - 1) q[a] -= L;
        }
        if (d.nbr[a][side] == d.rank) {
          keep.insert(keep.end(), q, q + K);
          continue;
        }
        sendbuf_[side].insert(sendbuf_[side].end(), q, q + K);
      }
      buf.swap(keep);
      if (d.dims[a] == 1) continue;  // no neighbours along this axis
      for (int dir = 0; dir < 2; ++dir) {
        int to = d.nbr[a][dir], from = d.nbr[a][1 - dir];
        long long ns = (long long)sendbuf_[dir].size(), nr = 0;
        MPI_Sendrecv(&ns, 1, MPI_LONG_LONG, to, 300 + a, &nr, 1, MPI_LONG_LONG, from, 300 + a, d.comm, MPI_STATUS_IGNORE);
        if (from == MPI_PROC_NULL) nr = 0;
        recvbuf_[dir].resize(nr);
        MPI_Sendrecv(sendbuf_[dir].data(), int(ns), MPI_DOUBLE, to, 400 + a, recvbuf_[dir].data(), int(nr), MPI_DOUBLE,
                     from, 400 + a, d.comm, MPI_STATUS_IGNORE);
        sp.stats.migrated += nr / K;
        buf.insert(buf.end(), recvbuf_[dir].begin(), recvbuf_[dir].end());
      }
    }
    // (4) whatever remains is inside the local box: append back
    for (std::size_t k = 0; k < buf.size(); k += K) p.push_packed(&buf[k]);
    sp.stats.absorbed_wall_charge += lost * sp.q;
  }
}

// ------------------------------------------------------------------------------------
// Counting sort by cell (improves cache locality of gather/deposit)
// ------------------------------------------------------------------------------------
void ParticleSystem::sort() {
  const Decomposition& d = *d_;
  const int nx = d.c1[0] - d.c0[0], ny = d.c1[1] - d.c0[1], nz = d.c1[2] - d.c0[2];
  const std::size_t ncell = std::size_t(nx) * ny * nz;
  Vec3 llo = d.local_lo();
  for (auto& sp : species) {
    ParticleArray& p = sp.p;
    const std::size_t np = p.size();
    if (np < 2) continue;
    std::vector<std::uint32_t> cell(np);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < np; ++i) {
      int ci = std::clamp(int((p.x[i] - llo.x) * inv_dx_[0]), 0, nx - 1);
      int cj = std::clamp(int((p.y[i] - llo.y) * inv_dx_[1]), 0, ny - 1);
      int ck = std::clamp(int((p.z[i] - llo.z) * inv_dx_[2]), 0, nz - 1);
      cell[i] = std::uint32_t(ci + std::size_t(nx) * (cj + std::size_t(ny) * ck));
    }
    std::vector<std::size_t> start(ncell + 1, 0);
    for (std::size_t i = 0; i < np; ++i) start[cell[i] + 1]++;
    std::partial_sum(start.begin(), start.end(), start.begin());
    std::vector<std::size_t> perm(np);
    for (std::size_t i = 0; i < np; ++i) perm[start[cell[i]]++] = i;
    std::vector<double> tmp(np);
    for (auto* arr : p.arrays()) {
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < np; ++i) tmp[i] = (*arr)[perm[i]];
      arr->swap(tmp);
    }
  }
}

// ------------------------------------------------------------------------------------
// Deposition
// ------------------------------------------------------------------------------------
// mode: 0 = q w, 1 = w, 2..4 = w v_{x,y,z}, 5 = w |v|^2
template <int O>
void ParticleSystem::deposit_impl(const Species& sp, std::vector<double>& buf, int mode, double* const* out) {
  using S = Shape<O>;
  const ParticleArray& p = sp.p;
  const std::size_t np = p.size();
  const auto& n = fs_->rho.size();
  const std::size_t sy = n[0], sz = std::size_t(n[0]) * n[1];
  const double ox = arr_lo_[0], oy = arr_lo_[1], oz = arr_lo_[2];
  const double ix = inv_dx_[0], iy = inv_dx_[1], iz = inv_dx_[2];
  (void)buf;
#pragma omp parallel
  {
    double* B = out[omp_get_thread_num()];
#pragma omp for schedule(static)
    for (std::size_t i = 0; i < np; ++i) {
      double val;
      switch (mode) {
        case 0: val = sp.q * p.w[i]; break;
        case 1: val = p.w[i]; break;
        case 2: val = p.w[i] * p.vx[i]; break;
        case 3: val = p.w[i] * p.vy[i]; break;
        case 4: val = p.w[i] * p.vz[i]; break;
        default: val = p.w[i] * (p.vx[i] * p.vx[i] + p.vy[i] * p.vy[i] + p.vz[i] * p.vz[i]); break;
      }
      int i0, j0, k0;
      double wx[S::n], wy[S::n], wz[S::n];
      S::weights((p.x[i] - ox) * ix, i0, wx);
      S::weights((p.y[i] - oy) * iy, j0, wy);
      S::weights((p.z[i] - oz) * iz, k0, wz);
#pragma GCC unroll 3
      for (int c = 0; c < S::n; ++c)
#pragma GCC unroll 3
        for (int b = 0; b < S::n; ++b) {
          const std::size_t base = std::size_t(i0) + sy * std::size_t(j0 + b) + sz * std::size_t(k0 + c);
          const double vbc = val * wy[b] * wz[c];
#pragma GCC unroll 3
          for (int a = 0; a < S::n; ++a) B[base + a] += vbc * wx[a];
        }
    }
  }
}

void ParticleSystem::reduce_thread_buffers(Field3D& out) {
  const std::size_t N = out.total();
  const int nt = int(tbuf_.size());
  double* o = out.data();
#pragma omp parallel for schedule(static)
  for (std::size_t p = 0; p < N; ++p) {
    double s = 0;
    for (int t = 0; t < nt; ++t) s += tbuf_[t][p];
    o[p] = s;
  }
}

namespace {
void scale_by_volume(Field3D& f, const FieldSolver& fs, const GridGeom& g) {
  const double inv = 1.0 / (g.dx[0] * g.dx[1] * g.dx[2]);
  const auto& n = f.size();
#pragma omp parallel for collapse(2) schedule(static)
  for (int k = 0; k < n[2]; ++k)
    for (int j = 0; j < n[1]; ++j)
      for (int i = 0; i < n[0]; ++i) f(i, j, k) *= inv / fs.node_volume_factor(i, j, k);
}
}  // namespace

void ParticleSystem::deposit_charge(Field3D& rho) {
  const std::size_t N = rho.total();
  std::vector<double*> ptrs(tbuf_.size());
#pragma omp parallel
  {
    int t = omp_get_thread_num();
    tbuf_[t].assign(N, 0.0);
    ptrs[t] = tbuf_[t].data();
  }
  for (auto& sp : species) {
    if (cfg_.shape_order == 2) deposit_impl<2>(sp, tbuf_[0], 0, ptrs.data());
    else deposit_impl<1>(sp, tbuf_[0], 0, ptrs.data());
  }
  reduce_thread_buffers(rho);
  rho.halo_add();
  scale_by_volume(rho, *fs_, *g_);
}

void ParticleSystem::deposit_moments(int s, Field3D& nf, Field3D& ux, Field3D& uy, Field3D& uz, Field3D& T) {
  const std::size_t N = nf.total();
  std::vector<double*> ptrs(tbuf_.size());
  Field3D* outs[6] = {nullptr, &nf, &ux, &uy, &uz, &T};
  for (int mode = 1; mode <= 5; ++mode) {
#pragma omp parallel
    {
      int t = omp_get_thread_num();
      tbuf_[t].assign(N, 0.0);
      ptrs[t] = tbuf_[t].data();
    }
    if (cfg_.shape_order == 2) deposit_impl<2>(species[s], tbuf_[0], mode, ptrs.data());
    else deposit_impl<1>(species[s], tbuf_[0], mode, ptrs.data());
    reduce_thread_buffers(*outs[mode]);
    outs[mode]->halo_add();
  }
  const double m = species[s].m;
  for (std::size_t p = 0; p < N; ++p) {
    double w = nf[p];
    if (w > 0) {
      double vx = ux[p] / w, vy = uy[p] / w, vz = uz[p] / w;
      double v2 = T[p] / w;
      ux[p] = vx;
      uy[p] = vy;
      uz[p] = vz;
      T[p] = std::max(0.0, m * (v2 - (vx * vx + vy * vy + vz * vz)) / (3.0 * constants::qe));
    } else {
      ux[p] = uy[p] = uz[p] = T[p] = 0.0;
    }
  }
  scale_by_volume(nf, *fs_, *g_);
}

}  // namespace picsp
