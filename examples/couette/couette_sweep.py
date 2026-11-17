#!/usr/bin/env python3
"""Sweep the Couette case over the shear rates in case.yaml.

For each gamma: build the lattice parameters, run the compiled Palabos
binary (./couette, built in examples/couette/build/), and collect
runs/gamma_<value>/torque.csv. Default is a dry run: print the plan
and the lattice parameters, run nothing.

Usage:
    python3 couette_sweep.py            # dry run
    python3 couette_sweep.py --run      # actually run the sweep
    python3 couette_sweep.py --run --mpi 4
"""
import argparse
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from lattice import lattice_params, load_case  # noqa: E402

HERE = Path(__file__).parent
BIN = HERE / "build" / "couette"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="actually run the solver")
    ap.add_argument("--mpi", type=int, default=1)
    args = ap.parse_args()

    case = load_case()
    lat = lattice_params(case)
    if not lat["mach_ok"]:
        print("FAIL: lattice check failed (u_lb too high) — fix case.yaml first")
        return 1
    print(f"lattice: dx={lat['dx_m']*1e3:.3f} mm dt={lat['dt_s']:.3e} s "
          f"steps={lat['n_steps']:,} u_lb,max={lat['u_lb_max']:.4f}")

    gammas = case["sweep"]["shear_rate_1_s"]
    if not args.run:
        print(f"DRY RUN: would run {len(gammas)} points: {gammas} 1/s")
        print("Pass --run to execute. Binary expected at", BIN)
        return 0
    if not BIN.exists():
        print(f"Binary not found: {BIN} — build it first (see README.md)")
        return 1

    for gamma in gammas:
        out = HERE / "runs" / f"gamma_{gamma}"
        out.mkdir(parents=True, exist_ok=True)
        cmd = (["mpirun", "-np", str(args.mpi)] if args.mpi > 1 else []) + [
            str(BIN), "--gamma", str(gamma), "--out", str(out)]
        print("running:", " ".join(cmd))
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print(f"run failed for gamma={gamma}")
            return r.returncode
    print("sweep done — post-process with: python3 couette_post.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
