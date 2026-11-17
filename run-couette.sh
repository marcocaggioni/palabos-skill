#!/usr/bin/env bash
# Couette example — dry-run by default (lattice check + build check only).
#   ./run-couette.sh            dry run
#   ./run-couette.sh --run      run the sweep (needs a built ./couette binary)
#   ./run-couette.sh --post     post-process runs/ into flow_curve.csv/png
set -e
cd "$(dirname "$0")/examples/couette"
python3 lattice.py
if [ "$1" = "--post" ] 2>/dev/null; then python3 couette_post.py; exit 0; fi
if [ "$1" = "--run" ] 2>/dev/null; then python3 couette_sweep.py --run; python3 couette_post.py; exit 0; fi
python3 couette_sweep.py   # dry run: prints the plan, runs nothing
echo "Dry run only. Pass --run to execute the solver sweep."
