#!/usr/bin/env bash
# Install the build dependencies for PICSP3D, then configure, build and (optionally) test it.
#
#   ./scripts/install.sh               # packages (uses sudo) + build in ./build
#   ./scripts/install.sh --no-deps     # build only (clusters, or packages already present)
#   ./scripts/install.sh --test        # also run the verification suite
#   ./scripts/install.sh --python      # also install the picsp_post package (pip --user)
#
# Supported: Ubuntu >= 22.04, Debian >= 12, RHEL / Rocky / AlmaLinux 8 and 9, Fedora,
# openSUSE Leap 15.6 / Tumbleweed, Arch Linux. Other systems: install a C++20 compiler
# (GCC >= 8 or Clang >= 11), CMake >= 3.18 and an MPI library, then use --no-deps.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$ROOT/build"
DEPS=1 TEST=0 PY=0
for a in "$@"; do
  case "$a" in
    --no-deps) DEPS=0 ;;
    --test) TEST=1 ;;
    --python) PY=1 ;;
    --build-dir=*) BUILD="${a#*=}" ;;
    -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
    *) echo "unknown option: $a" >&2; exit 2 ;;
  esac
done

SUDO=""
if [ "$(id -u)" -ne 0 ] && [ "$DEPS" = 1 ]; then SUDO="sudo"; fi

. /etc/os-release 2>/dev/null || { ID=unknown; VERSION_ID=0; }
family="$ID ${ID_LIKE:-}"
major="${VERSION_ID%%.*}"
echo ">> detected: ${PRETTY_NAME:-$ID}"

# Extra environment for the build, filled in per distribution:
#   MPI_BIN   directory with mpicxx / mpirun when it is not on PATH
#   CXX_SEL   compiler to use when the default one is too old for C++20
MPI_BIN="" CXX_SEL=""

case "$family" in
  *debian*|*ubuntu*)
    if [ "$DEPS" = 1 ]; then
      $SUDO apt-get update -q
      $SUDO env DEBIAN_FRONTEND=noninteractive apt-get install -y -q --no-install-recommends \
        build-essential cmake ninja-build libopenmpi-dev openmpi-bin python3 python3-pip ca-certificates
    fi
    ;;
  *rhel*|*centos*|*rocky*|*almalinux*|*fedora*)
    if [ "$DEPS" = 1 ]; then
      py="python3 python3-pip"
      # RHEL 8 ships Python 3.6; the post-processing package needs >= 3.9
      [ "$ID" != fedora ] && [ "$major" -lt 9 ] && py="python3.11 python3.11-pip"
      $SUDO dnf install -y -q gcc-c++ make cmake openmpi-devel $py
    fi
    MPI_BIN=/usr/lib64/openmpi/bin   # Open MPI is not on PATH without 'module load mpi/openmpi-x86_64'
    ;;
  *suse*)
    if [ "$DEPS" = 1 ]; then
      if [ "$ID" = opensuse-leap ]; then
        $SUDO zypper -n install gcc13-c++ cmake make openmpi4-devel python311 python311-pip
      else
        $SUDO zypper -n install gcc-c++ cmake make openmpi4-devel python3 python3-pip
      fi
    fi
    [ "$ID" = opensuse-leap ] && CXX_SEL=g++-13   # Leap's default GCC 7 lacks C++20
    MPI_BIN=/usr/lib64/mpi/gcc/openmpi4/bin
    ;;
  *arch*)
    [ "$DEPS" = 1 ] && $SUDO pacman -Sy --noconfirm --needed base-devel cmake ninja openmpi python python-pip
    ;;
  *)
    if [ "$DEPS" = 1 ]; then
      echo "!! unknown distribution '$ID': install a C++20 compiler, CMake >= 3.18 and MPI, then rerun with --no-deps" >&2
      exit 1
    fi
    ;;
esac

# Environment modules (clusters, RHEL) may already provide MPI; otherwise use the distribution path
if ! command -v mpicxx >/dev/null 2>&1 && [ -n "$MPI_BIN" ] && [ -x "$MPI_BIN/mpicxx" ]; then
  export PATH="$MPI_BIN:$PATH"
fi
command -v mpicxx >/dev/null 2>&1 || { echo "!! mpicxx not found; install an MPI library or load its module" >&2; exit 1; }

cmake_args=(-S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release
            -DMPI_CXX_COMPILER="$(command -v mpicxx)" -DMPIEXEC_EXECUTABLE="$(command -v mpirun || command -v mpiexec)")
[ -n "$CXX_SEL" ] && cmake_args+=(-DCMAKE_CXX_COMPILER="$CXX_SEL")
command -v ninja >/dev/null 2>&1 && cmake_args+=(-G Ninja)
[ "$TEST" = 0 ] && cmake_args+=(-DPICSP_TESTS=OFF)

echo ">> mpicxx: $(command -v mpicxx)"
cmake "${cmake_args[@]}"
cmake --build "$BUILD" -j "$(nproc 2>/dev/null || echo 4)"

if [ "$TEST" = 1 ]; then
  # containers often run as root, which Open MPI refuses without these
  [ "$(id -u)" -eq 0 ] && export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
  export OMPI_MCA_rmaps_base_oversubscribe=1 PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
  ctest --test-dir "$BUILD" --output-on-failure --timeout 1800
fi

if [ "$PY" = 1 ]; then
  py=python3
  command -v python3.11 >/dev/null 2>&1 && py=python3.11
  "$py" -m pip install --user -e "$ROOT/python"
fi

echo
echo ">> built $BUILD/picsp3d"
if [ -n "$MPI_BIN" ] && [ "${PATH#"$MPI_BIN"}" != "$PATH" ]; then
  echo ">> MPI is not on your PATH; before running:"
  echo "     export PATH=$MPI_BIN:\$PATH        # RHEL/Fedora alternative: module load mpi/openmpi-x86_64"
fi
echo ">> run:  OMP_NUM_THREADS=4 mpirun -np 2 $BUILD/picsp3d examples/langmuir.toml"
