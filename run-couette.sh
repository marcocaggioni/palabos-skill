#!/usr/bin/env bash
# Couette example — dry-run by default (lattice check + build check only).
#   ./run-couette.sh            dry run
#   ./run-couette.sh --run      run the sweep (needs a built ./couette binary)
#   ./run-couette.sh --post     post-process runs/ into flow_curve.csv/png
set -e
cd "$(dirname "$0")/examples/couette"
if command -v uv >/dev/null 2>&1; then PY="uv run python"; else PY="python3"; fi
$PY lattice.py
if [ "$1" = "--post" ] 2>/dev/null; then $PY couette_post.py; exit 0; fi
if [ "$1" = "--run" ] 2>/dev/null; then $PY couette_sweep.py --run; $PY couette_post.py; exit 0; fi
$PY couette_sweep.py   # dry run: prints the plan, runs nothing
echo "Dry run only. Pass --run to execute the solver sweep."
