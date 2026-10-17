#include "Simulation.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <sstream>

#include <omp.h>
#include <sched.h>

#include "util/Constants.hpp"
#include "util/Log.hpp"

namespace picsp {

namespace fsys = std::filesystem;

Simulation::Simulation(const Config& cfg, MPI_Comm world) : cfg_(cfg) {
  d_ = Decomposition::create(cfg_, world);
  g_ = GridGeom::from(d_);
}

Simulation::~Simulation() {
  if (d_.comm != MPI_COMM_NULL) MPI_Comm_free(&d_.comm);
}

void Simulation::initialise() {
  wall_start_ = MPI_Wtime();
  int nt = omp_get_max_threads();
  log0(" MPI ranks: %d (decomposition %d x %d x %d), OpenMP threads/rank: %d\n", d_.size, d_.dims[0], d_.dims[1],
       d_.dims[2], nt);
  {
    // warn if the MPI launcher pinned this rank to fewer cores than OpenMP threads
    cpu_set_t set;
    CPU_ZERO(&set);
    int avail = sched_getaffinity(0, sizeof(set), &set) == 0 ? CPU_COUNT(&set) : nt;
    int bad = avail < nt ? 1 : 0, anybad = 0;
    MPI_Allreduce(&bad, &anybad, 1, MPI_INT, MPI_MAX, d_.comm);
    if (anybad)
      warn0("a rank is bound to %d core(s) but runs %d OpenMP threads; launch with e.g.\n"
            "        mpirun --map-by ppr:<ranks per node>:node:pe=%d  (or --bind-to none)\n",
            avail, nt, nt);
  }
  rng_.clear();
  for (int t = 0; t < nt; ++t) rng_.emplace_back(cfg_.seed, std::uint64_t(d_.rank) * 4096 + t);

  timers_.start("setup");
  scene_.build(cfg_.objects);
  fs_.setup(cfg_, g_, scene_);
  ps_.setup(cfg_, d_, g_, scene_, fs_, &rng_);
  om_.setup(cfg_, d_, g_, scene_, fs_, ps_, &rng_);
  mcc_.setup(cfg_, ps_, d_, &rng_);
  vtk_.init(cfg_.output.dir, d_, cfg_.output.float32);

  bool restarted = false;
  if (!cfg_.restart_from.empty()) restarted = restart(cfg_.restart_from);
  if (!restarted) {
    ps_.load_initial();
    ps_.deposit_charge(fs_.rho);
    last_solve_ = fs_.solve(0.0);
    ps_.rewind_half();
  } else {
    ps_.deposit_charge(fs_.rho);
    last_solve_ = fs_.solve(time_);
  }
  timers_.stop("setup");

  const int no = int(cfg_.objects.size()), ns = int(ps_.species.size());
  hist_coll_.assign(std::size_t(no) * ns, 0.0);
  hist_emit_.assign(no, 0.0);
  hist_inj_.assign(ns, 0.0);
  hist_wall_.assign(ns, 0.0);
  hist_time_ = 0;

  if (d_.rank == 0 && cfg_.output.history_every > 0) {
    std::string path = cfg_.output.dir + "/history.csv";
    bool append = restarted && fsys::exists(path);
    hist_.open(path, append ? std::ios::app : std::ios::trunc);
    if (!append) {
      hist_ << "step,time";
      for (auto& s : ps_.species) hist_ << ",N_" << s.name() << ",KE_" << s.name();
      hist_ << ",E_field,E_total";
      for (auto& c : fs_.conductors()) hist_ << ",V_" << c.name << ",Q_" << c.name;
      for (int o = 0; o < no; ++o) {
        for (auto& s : ps_.species) hist_ << ",I_" << cfg_.objects[o].name << "_" << s.name();
        hist_ << ",Iemit_" << cfg_.objects[o].name;
      }
      for (auto& s : ps_.species) hist_ << ",Iwall_" << s.name();
      hist_ << ",solver_iters,collisions,wall_s\n";
    }
  }
  // write the initial state
  if (!restarted) {
    if (cfg_.output.fields_every > 0) write_fields();
    if (cfg_.output.particles_every > 0) write_particles();
    if (cfg_.output.surface_every > 0 && !om_.empty()) write_surface();
    if (cfg_.output.history_every > 0) history_line();
  }
  wall_last_status_ = MPI_Wtime();
  step_last_status_ = step_;
  initialised_ = true;
}

void Simulation::advance() {
  std::vector<std::vector<Hit>> hits;
  {
    ScopedTimer t(timers_, "push");
    for (auto& sp : ps_.species)
      if (sp.active(step_)) pushes_ += (long long)sp.p.size();
    ps_.push(step_, hits);
  }
  {
    ScopedTimer t(timers_, "objects");
    om_.process_hits(step_, hits);
    om_.photoemission(step_);
  }
  {
    ScopedTimer t(timers_, "inject+migrate");
    ps_.inject(step_);
    ps_.boundaries_and_migrate(step_);
  }
  if (mcc_.enabled()) {
    ScopedTimer t(timers_, "collisions");
    mcc_.collide(step_);
  }
  if (cfg_.sort_every > 0 && step_ % cfg_.sort_every == 0) {
    ScopedTimer t(timers_, "sort");
    ps_.sort();
  }
  {
    ScopedTimer t(timers_, "deposit");
    ps_.deposit_charge(fs_.rho);
  }
  {
    ScopedTimer t(timers_, "objects");
    om_.end_step();
  }
  ++step_;
  time_ = step_ * cfg_.dt;
  {
    ScopedTimer t(timers_, "field solve");
    last_solve_ = fs_.solve(time_);
    solve_iters_acc_ += last_solve_.iterations;
    ++solve_count_;
  }
  // history accumulators
  const int no = int(cfg_.objects.size()), ns = int(ps_.species.size());
  for (std::size_t i = 0; i < hist_coll_.size(); ++i) hist_coll_[i] += om_.step_collected()[i];
  for (int o = 0; o < no; ++o) hist_emit_[o] += om_.step_emitted()[o];
  for (int s = 0; s < ns; ++s) hist_wall_[s] += ps_.species[s].stats.absorbed_wall_charge;
  hist_time_ += cfg_.dt;

  ScopedTimer t(timers_, "output");
  const auto& oc = cfg_.output;
  if (oc.history_every > 0 && step_ % oc.history_every == 0) history_line();
  if (oc.fields_every > 0 && step_ % oc.fields_every == 0) write_fields();
  if (oc.particles_every > 0 && step_ % oc.particles_every == 0) write_particles();
  if (oc.surface_every > 0 && step_ % oc.surface_every == 0 && !om_.empty()) write_surface();
  if (oc.checkpoint_every > 0 && step_ % oc.checkpoint_every == 0) checkpoint();
  if (oc.status_every > 0 && step_ % oc.status_every == 0) status_line();
}

void Simulation::run() {
  if (!initialised_) initialise();
  log0("---------------------------------------------------------------- running %lld steps\n", cfg_.steps - step_);
  double t0 = MPI_Wtime();
  long long p0 = pushes_;
  while (step_ < cfg_.steps) advance();
  if (cfg_.output.checkpoint_every > 0 && step_ % cfg_.output.checkpoint_every != 0) checkpoint();
  double el = MPI_Wtime() - t0;
  long long pushes = allreduce_sum(pushes_ - p0, d_.comm);
  log0("----------------------------------------------------------------\n");
  log0(" finished %lld steps in %.2f s  (%.3f ms/step, %.3g particle-steps/s, %.1f ns/particle-step/core)\n",
       cfg_.steps, el, cfg_.steps > 0 ? 1e3 * el / std::max<long long>(1, cfg_.steps) : 0.0, pushes / std::max(el, 1e-30),
       pushes > 0 ? 1e9 * el * d_.size * omp_get_max_threads() / pushes : 0.0);
  log0(" average Poisson iterations per step: %.1f\n", solve_count_ ? double(solve_iters_acc_) / solve_count_ : 0.0);
  timers_.report(d_.comm, MPI_Wtime() - wall_start_);
}

void Simulation::status_line() {
  double now = MPI_Wtime();
  long long n = 0;
  for (auto& sp : ps_.species) n += (long long)sp.p.size();
  n = allreduce_sum(n, d_.comm);
  double ms = 1e3 * (now - wall_last_status_) / std::max<long long>(1, step_ - step_last_status_);
  std::ostringstream os;
  os << std::setprecision(4);
  for (auto& c : fs_.conductors()) os << "  V(" << c.name << ")=" << c.V;
  log0(" step %8lld  t=%.4e s  N=%-11lld  solver %2d it (%.1e)  %7.2f ms/step%s\n", step_, time_, n,
       last_solve_.iterations, last_solve_.residual, ms, os.str().c_str());
  if (!last_solve_.converged) warn0("Poisson solver did not converge (residual %.2e)\n", last_solve_.residual);
  wall_last_status_ = now;
  step_last_status_ = step_;
}

void Simulation::history_line() {
  const int no = int(cfg_.objects.size()), ns = int(ps_.species.size());
  std::vector<long long> N(ns);
  std::vector<double> KE(ns);
  for (int s = 0; s < ns; ++s) {
    N[s] = ps_.count_global(s);
    KE[s] = ps_.kinetic_energy(s);
  }
  double EF = fs_.field_energy();
  if (d_.rank == 0 && hist_.is_open()) {
    double tot = EF;
    hist_ << step_ << "," << std::setprecision(10) << time_;
    hist_ << std::setprecision(8);
    for (int s = 0; s < ns; ++s) {
      hist_ << "," << N[s] << "," << KE[s];
      tot += KE[s];
    }
    hist_ << "," << EF << "," << tot;
    for (auto& c : fs_.conductors()) hist_ << "," << c.V << "," << c.Q;
    double T = hist_time_ > 0 ? hist_time_ : 1.0;
    for (int o = 0; o < no; ++o) {
      for (int s = 0; s < ns; ++s) hist_ << "," << hist_coll_[std::size_t(o) * ns + s] / T;
      hist_ << "," << hist_emit_[o] / T;
    }
    for (int s = 0; s < ns; ++s) hist_ << "," << hist_wall_[s] / T;
    hist_ << "," << last_solve_.iterations << "," << mcc_.collisions_last_step() << "," << (MPI_Wtime() - wall_start_)
          << "\n";
    hist_.flush();
  }
  std::fill(hist_coll_.begin(), hist_coll_.end(), 0.0);
  std::fill(hist_emit_.begin(), hist_emit_.end(), 0.0);
  std::fill(hist_wall_.begin(), hist_wall_.end(), 0.0);
  hist_time_ = 0;
}

void Simulation::write_fields() {
  std::vector<NamedScalar> sc = {{"phi", &fs_.phi}, {"rho", &fs_.rho}};
  std::vector<NamedVector> vec = {{"E", &fs_.Ex, &fs_.Ey, &fs_.Ez}};
  Field3D obj;
  if (!scene_.empty()) {
    obj.init(g_, FieldSolver::kGhosts);
    for (std::size_t p = 0; p < obj.total(); ++p) obj[p] = fs_.eb().obj(p);
    sc.push_back({"object_id", &obj});
    sc.push_back({"sdf", &fs_.eb().sdf()});  // signed distance to the nearest surface [m] (capped)
  }
  std::vector<Field3D> mom;
  if (cfg_.output.moments) {
    const int ns = int(ps_.species.size());
    mom.resize(5 * ns);
    for (auto& f : mom) f.init(g_, FieldSolver::kGhosts);
    for (int s = 0; s < ns; ++s) {
      Field3D* m = &mom[5 * s];
      ps_.deposit_moments(s, m[0], m[1], m[2], m[3], m[4]);
    }
    for (int s = 0; s < ns; ++s) {
      const std::string& n = ps_.species[s].name();
      sc.push_back({"n_" + n, &mom[5 * s]});
      sc.push_back({"T_" + n, &mom[5 * s + 4]});
      vec.push_back({"u_" + n, &mom[5 * s + 1], &mom[5 * s + 2], &mom[5 * s + 3]});
    }
  }
  vtk_.write_fields(step_, time_, sc, vec);
}

void Simulation::write_particles() {
  const int stride = cfg_.output.particle_stride;
  std::vector<float> xyz;
  PolyArray vel{"velocity", 3, {}}, sid{"species", 1, {}}, wt{"weight", 1, {}}, ke{"energy_eV", 1, {}};
  for (auto& sp : ps_.species) {
    const ParticleArray& p = sp.p;
    for (std::size_t i = 0; i < p.size(); i += stride) {
      xyz.insert(xyz.end(), {float(p.x[i]), float(p.y[i]), float(p.z[i])});
      vel.data.insert(vel.data.end(), {float(p.vx[i]), float(p.vy[i]), float(p.vz[i])});
      sid.data.push_back(float(sp.id));
      wt.data.push_back(float(p.w[i]));
      double v2 = p.vx[i] * p.vx[i] + p.vy[i] * p.vy[i] + p.vz[i] * p.vz[i];
      ke.data.push_back(float(0.5 * sp.m * v2 / constants::qe));
    }
  }
  vtk_.write_points("particles", step_, time_, xyz, {vel, sid, wt, ke});
}

void Simulation::write_surface() {
  const auto& P = om_.panels();
  const int ns = om_.nspecies();
  std::vector<float> tri;
  tri.reserve(P.size() * 9);
  PolyArray oid{"object_id", 1, {}}, pot{"potential", 1, {}}, nrm{"normal", 3, {}}, sun{"sun_cos", 1, {}},
      jnet{"J_net", 1, {}}, jem{"J_emit", 1, {}};
  std::vector<PolyArray> js(ns);
  for (int s = 0; s < ns; ++s) js[s] = {"J_" + ps_.species[s].name(), 1, {}};
  std::vector<std::vector<double>> jsp(ns);
  for (int s = 0; s < ns; ++s) jsp[s] = om_.cell_averaged_density(om_.panel_charge(), ns, s);
  std::vector<double> jemit = om_.cell_averaged_density(om_.panel_emitted(), 1, 0);
  for (std::size_t i = 0; i < P.size(); ++i) {
    const Panel& p = P[i];
    for (const Vec3* v : {&p.a, &p.b, &p.c}) tri.insert(tri.end(), {float(v->x), float(v->y), float(v->z)});
    oid.data.push_back(float(p.obj));
    pot.data.push_back(float(fs_.conductors()[fs_.object_conductor(p.obj)].V));
    nrm.data.insert(nrm.data.end(), {float(p.normal.x), float(p.normal.y), float(p.normal.z)});
    sun.data.push_back(float(p.sun_cos));
    double net = 0;
    for (int s = 0; s < ns; ++s) {
      js[s].data.push_back(float(jsp[s][i]));
      net += jsp[s][i];
    }
    jem.data.push_back(float(jemit[i]));
    jnet.data.push_back(float(net + jemit[i]));
  }
  std::vector<PolyArray> arrays = {oid, pot, nrm, sun, jnet, jem};
  for (auto& a : js) arrays.push_back(a);
  vtk_.write_triangles("surface", step_, time_, tri, arrays);
  om_.reset_panel_accumulators();
}

// ------------------------------------------------------------------------------------
// Checkpoint / restart: one binary file per rank (same rank count required)
// ------------------------------------------------------------------------------------
namespace {
constexpr std::uint64_t kMagic = 0x5049435350334421ULL;  // "PICSP3D!"
}

void Simulation::checkpoint() {
  ScopedTimer t(timers_, "checkpoint");
  std::string dir = cfg_.output.dir + "/checkpoint";
  if (d_.rank == 0) fsys::create_directories(dir);
  MPI_Barrier(d_.comm);
  char name[64];
  std::snprintf(name, sizeof(name), "/rank_%05d.bin", d_.rank);
  std::string tmp = dir + name + ".tmp";
  {
    std::ofstream f(tmp, std::ios::binary);
    auto put = [&](const void* p, std::size_t n) { f.write(reinterpret_cast<const char*>(p), std::streamsize(n)); };
    std::uint64_t magic = kMagic;
    put(&magic, 8);
    int nr = d_.size;
    put(&nr, 4);
    put(&step_, 8);
    put(&time_, 8);
    int nt = int(rng_.size());
    put(&nt, 4);
    for (auto& r : rng_) {
      unsigned char b[Rng::kStateBytes];
      r.save(b);
      put(b, Rng::kStateBytes);
    }
    int nc = int(fs_.conductors().size());
    put(&nc, 4);
    for (auto& c : fs_.conductors()) {
      put(&c.Q, 8);
      put(&c.V, 8);
    }
    std::uint64_t nf = fs_.solver_state().total();
    put(&nf, 8);
    put(fs_.solver_state().data(), nf * 8);
    int ns = int(ps_.species.size());
    put(&ns, 4);
    for (auto& sp : ps_.species) {
      std::uint64_t n = sp.p.size();
      put(&n, 8);
      for (auto* a : sp.p.arrays()) put(a->data(), n * 8);
    }
  }
  MPI_Barrier(d_.comm);
  fsys::rename(tmp, dir + name);
  if (d_.rank == 0) {
    std::ofstream m(dir + "/info.txt");
    m << "step " << step_ << "\ntime " << std::setprecision(17) << time_ << "\nranks " << d_.size << "\n";
  }
  log0(" checkpoint written at step %lld\n", step_);
}

bool Simulation::restart(const std::string& dir) {
  char name[64];
  std::snprintf(name, sizeof(name), "/rank_%05d.bin", d_.rank);
  std::ifstream f(dir + name, std::ios::binary);
  if (!f) fatal("cannot open checkpoint " + dir + name);
  auto get = [&](void* p, std::size_t n) { f.read(reinterpret_cast<char*>(p), std::streamsize(n)); };
  std::uint64_t magic;
  get(&magic, 8);
  if (magic != kMagic) fatal("not a PICSP3D checkpoint: " + dir);
  int nr;
  get(&nr, 4);
  if (nr != d_.size) fatal("checkpoint was written with a different number of MPI ranks");
  get(&step_, 8);
  get(&time_, 8);
  int nt;
  get(&nt, 4);
  for (int t = 0; t < nt; ++t) {
    unsigned char b[Rng::kStateBytes];
    get(b, Rng::kStateBytes);
    if (t < int(rng_.size())) rng_[t].load(b);
  }
  if (nt != int(rng_.size())) warn0("checkpoint used %d OpenMP threads, now %d: random streams differ\n", nt, int(rng_.size()));
  int nc;
  get(&nc, 4);
  if (nc != int(fs_.conductors().size())) fatal("checkpoint conductor count mismatch");
  for (auto& c : fs_.conductors()) {
    get(&c.Q, 8);
    get(&c.V, 8);
  }
  std::uint64_t nf;
  get(&nf, 8);
  if (nf != fs_.solver_state().total()) fatal("checkpoint grid mismatch");
  get(fs_.solver_state().data(), nf * 8);
  int ns;
  get(&ns, 4);
  if (ns != int(ps_.species.size())) fatal("checkpoint species count mismatch");
  for (auto& sp : ps_.species) {
    std::uint64_t n;
    get(&n, 8);
    sp.p.resize(n);
    for (auto* a : sp.p.arrays()) get(a->data(), n * 8);
  }
  if (!f) fatal("truncated checkpoint " + dir + name);
  // floating charges were restored: do not re-initialise them on the first solve
  fs_.mark_charges_initialised();
  log0(" restarted from %s at step %lld (t = %.4e s)\n", dir.c_str(), step_, time_);
  return true;
}

}  // namespace picsp
