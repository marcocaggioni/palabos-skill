#!/usr/bin/env bash
# Fetch the Palabos C++ library (there is no PyPI package for it) and
# prove the toolchain by building one of Palabos' own showcase examples.
# Usage:
#   bash scripts/install-palabos.sh              # clone to ~/palabos
#   PALABOS_ROOT=/somewhere bash scripts/install-palabos.sh
#   bash scripts/install-palabos.sh --build-example   # also build cavity2d
# After this: export PALABOS_ROOT=~/palabos
# Then build this repo's Couette case:
#   cd examples/couette && mkdir -p build && cd build
#   cmake .. -DPALABOS_ROOT="$PALABOS_ROOT" && make
set -e
ROOT="${PALABOS_ROOT:-$HOME/palabos}"

if [ -d "$ROOT/.git" ]; then
  echo "Palabos already cloned at $ROOT — pulling latest."
  git -C "$ROOT" pull --ff-only
else
  echo "Cloning Palabos into $ROOT ..."
  git clone --depth 1 https://gitlab.com/unigespc/palabos.git "$ROOT"
fi

if [ "$1" = "--build-example" ]; then
  echo "Building Palabos' own cavity2d example to prove the toolchain ..."
  cd "$ROOT/examples/showCases/cavity2d/build"
  cmake .. && make -j"$(nproc)"
  echo "cavity2d built OK."
fi

echo "Palabos is at $ROOT"
echo "Add to your shell: export PALABOS_ROOT=$ROOT"
