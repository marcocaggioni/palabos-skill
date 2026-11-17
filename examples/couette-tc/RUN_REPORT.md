# Run report — Palabos TC Couette example (real run, 2026-11-17)

First end-to-end TC run in this repo. TC validation **PASSED**
(max analytic-torque error 3.80%, gate 10%). Nothing in this file
is synthetic: every number below comes from the Palabos runs in
`runs/`, via `couette_post.py`.

## Environment

Same as the Newtonian run (`examples/couette/RUN_REPORT.md`):
Palabos master at commit `4127697e90169bbef982295f1d1c933cf6e90caa`
(2026-04-03), checkout at `~/palabos`, `PALABOS_ROOT=~/palabos`.
Machine: this VM, 4 CPUs, ~7 GB RAM. Serial build, no MPI
(`ENABLE_MPI=OFF`). g++ 13.3, CMake 3.28, C++20. System python3
(numpy, matplotlib, pyyaml) for the scripts.

## Case (from case.yaml)

2D annulus cross-section, Ri = 0.020 m, Ro = 0.025 m, L = 0.060 m.
TC fluid: sigma_y = 70 Pa, gamma_dot_c = 100 1/s, eta_bg = 0.05 Pa s,
rho = 1000 kg/m3. Regularization: eta = min(eta(gamma_dot),
eta_max), eta_max = eta(1e-3 1/s) = 70221.4 Pa s. Sweep
gamma = [0.1, 0.3, 1, 3, 10, 30, 100] 1/s (Newtonian
inner-wall-rate label convention). 10 cells across the gap
(dx = 0.5 mm), omega_target = 1.0 at the label viscosity, 24,000
steps/point, quarter-sine ramp over the first 20% of the steps,
steady window last 20%, initial condition rest (u = 0).

### Per-point lattice (from lattice.py; the C++ case recomputes the same dt from gamma)

| gamma (1/s) | eta_label (Pa s) | dt (s) | run time (s) | omega_plug | u_lb |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 722.19 | 5.77e-8 | 1.38e-3 | 0.02036 | 0.00000 |
| 0.3 | 246.16 | 1.69e-7 | 4.06e-3 | 0.00699 | 0.00000 |
| 1 | 77.05 | 5.41e-7 | 1.30e-2 | 0.00219 | 0.00000 |
| 3 | 27.42 | 1.52e-6 | 3.65e-2 | 0.00078 | 0.00003 |
| 10 | 9.26 | 4.50e-6 | 1.08e-1 | 0.00026 | 0.00032 |
| 30 | 3.66 | 1.14e-5 | 2.73e-1 | 0.00010 | 0.00246 |
| 100 | 1.45 | 2.87e-5 | 6.90e-1 | 0.00004 | 0.02069 |

(Correction: an early brief wrote eta(100) = 70.75 Pa s, conflating
a stress term with a viscosity. The correct value is 1.45 Pa s —
eta = sigma_y/gamma + sigma_y/sqrt(gamma*gamma_dot_c) + eta_bg —
and it is what lattice.py, the C++ case, and every run below use.
A dt solved from 70.75 would have pushed omega toward ~2.0. The
eta_max cap, 70221.4 Pa s, was correct as stated.)

## How the TC rheology is implemented

`TCBGKdynamics` in `couette.cpp`, modelled line-for-line on
Palabos' own `BGKCarreauDynamics`
(`src/complexDynamics/carreauDynamics.hh`): in `collide()`,
PiNeq from the populations via
`momentTemplates::compute_rhoBar_j_PiNeq`, local shear rate from
PiNeq by Chapman-Enskog,
gamma_dot_lb = omega*sqrt(2*|PiNeq|^2)/(2*rho*cs^2) (the factor
matches Palabos' own Carreau alpha = lambda^2*gamma_dot^2),
local eta from the TC law with the cap, local
omega = 1/(3*nu_lb+0.5), iterated to a fixed point (10 iterations
max, 1e-10 relative tolerance — the map is a contraction, factor
~0.24 near the inner wall at gamma = 10), then a standard BGK
collision with that omega. No per-cell state is stored: omega is
recomputed from the cell's own populations at every collide.
Everything else (staircase annulus, `VelocityBounceBack` moving
inner wall with `setUwall()` ramp, outer `BounceBack`,
momentum-exchange torque, `rho*dx^4/dt^2` conversion) is identical
to the Newtonian case — see `examples/couette/couette.cpp` and its
`RUN_REPORT.md`.

One fix during this run: the first build warned that
`gammaDotReg` was unused (etaMax was hard-coded). Fixed by
computing etaMax from gammaDotReg at the top of `main()`, so the
code cannot drift from `case.yaml` / `lattice.py` (which checks
eta_max against eta(gamma_dot_reg) and refuses to run on mismatch).

## Build

`cmake .. -DPALABOS_ROOT=$HOME/palabos && make -j2` in
`examples/couette-tc/build/` — clean, one warning (fixed as
above), ~7 min, mostly rebuilding `libpalabos.a` from the Palabos
sources. Binary `build/couette`.

## Gamma = 10 check (first run, before the sweep)

- Simulated steady torque (momentum exchange): **0.0164296 N·m**
  (per-length 0.273827 N·m/m)
- Analytic TC torque: 0.0162797 N·m
- **Relative error: 0.92%** (gate: 10%)
- Torque from velocity profile: 0.0152631 N·m (6.2% low — see
  "Still unverified / known limitations")
- omega seen in fluid cells: min 0.000264, max 1.381 (target 1.0
  at the label viscosity; low end is the capped plug)
- Torque time series: settled by step ~4800 (mean 0.27384,
  sd 0.00005); dead flat after (sd 0.00000 in the steady window)
- Velocity profile vs analytic: +7% at the innermost bin rising
  to +30% mid-gap and more near the outer wall — see the profile
  discussion below; the torque (the validated quantity) is
  unaffected
- Run time: ~72 s for 24,000 steps (the PiNeq + fixed-point
  collide costs ~12x the Newtonian BGK collide)

## Full sweep (`couette_sweep.py --run`, then `couette_post.py`)

`couette_post.py` output — **validation PASSED** (max analytic
error 3.80%, tolerance 10%):

| gamma (1/s) | tau (Pa) | tau_an (Pa) | torque (N·m) | analytic err | R_y (mm) | plug? | omega range |
|---:|---:|---:|---:|---:|---:|---|---:|
| 0.1 | 80.33 | 77.39 | 1.2114e-02 | 3.80% | 21.03 | yes | 0.020–1.61 |
| 0.3 | 83.13 | 80.80 | 1.2536e-02 | 2.89% | 21.49 | yes | 0.007–1.60 |
| 1 | 88.27 | 86.46 | 1.3311e-02 | 2.10% | 22.23 | yes | 0.002–1.57 |
| 3 | 95.86 | 94.37 | 1.4455e-02 | 1.57% | 23.22 | yes | 0.0008–1.50 |
| 10 | 108.95 | 107.96 | 1.6430e-02 | 0.92% | 24.84 | yes | 0.0003–1.38 |
| 30 | 128.97 | 128.08 | 1.9448e-02 | 0.69% | 25.00 | no | 0.0001–1.27 |
| 100 | 168.44 | 167.34 | 2.5400e-02 | 0.66% | 25.00 | no | 0.00004–1.18 |

The error is positive at every point and falls monotonically with
gamma — systematic, not noise. It is largest where the yielded
zone is thinnest (R_y = 21.03 mm at gamma = 0.1 leaves a yielded
layer ~1 mm wide, ~2 cells) and smallest where the gap is fully
yielded. The Newtonian run showed the same sign and character
(+1.43% constant, staircase-wall bias: the staircase makes the
effective inner radius very slightly larger than nominal). The TC
error is the same bias amplified at low gamma by the thin yielded
layer. Total solver time for the sweep: ~8 min (54–86 s/point).
Torque steadiness: sd <= 0.007% in the steady window at every
point, and the mid-run (40–60%) mean equals the steady mean to 6
digits — every point was settled long before the averaging window.

Note on determinism: gamma = 0.1 was run twice (once in the
first sweep attempt, once alone to capture its diagnostics) and
gave the identical torque 0.0121142 N·m both times.

## Videos

Rendered from the field snapshots in `runs/` by `make_video.py`
(matplotlib + ffmpeg, H.264, 600x600, 10 fps):

- `couette-tc-gamma10.mp4` (2.5 s, 48 KB) — velocity field vs time
  at gamma = 10 1/s, 25 snapshots (every 1000 steps). Spin-up to
  steady rotation; the bob renders dark (solid-body nodes carry no
  fluid velocity), the sheared annulus hugs the inner wall and is
  bright, the plug near the outer wall is dark/slow.
- `couette-tc-sweep.mp4` (7.0 s, 73 KB) — steady-state field at
  each sweep shear rate, 1 s per point, gamma = 0.1 -> 100 1/s.
- Snapshots: `runs/gamma_*/field_*.csv`, 25 per run, 51x51 grid
  (every 2nd node), columns `iX,iY,ux_m_s,uy_m_s`. One frame of
  each video was visually checked (steady state, gamma = 10):
  clean rotating annulus, tangential arrows, speed falling from
  inner to outer wall, dark plug at the outer wall.

## Files left in the repo (not pushed)

- `couette.cpp`, `CMakeLists.txt` (working, as run here)
- `runs/gamma_*/torque.csv` (241 samples per run: time,
  torque per unit length), `profile.csv`, `params.json`
- `runs/gamma_*/field_*.csv` (25 snapshots per run)
- `flow_curve.csv`, `flow_curve.png`, `lattice.json`
- `make_video.py`, `couette-tc-gamma10.mp4`, `couette-tc-sweep.mp4`
- this file

Note: `.gitignore` excludes `runs/` and `build/` (and
`flow_curve.csv/png`, `lattice.json`) for `couette-tc/` as for
`couette/`; the push script uploads the deliverables via the
Contents API regardless.

## Still unverified / known limitations

- **Velocity profile vs analytic (gamma = 10).** The simulated
  azimuthal velocity is higher than the analytic TC profile:
  +7% at the innermost bin, +15–20% mid-gap, more near the outer
  wall. This was checked three ways and is not a measurement
  artifact: (1) bin-averaging the analytic profile over each bin
  changes it by <1%; (2) the raw field snapshots along a ray
  agree with the binned profile; (3) applying the same
  finite-difference shear-rate estimator to the analytic profile
  recovers the analytic shear rate to +1.5–4% (mid-gap), while
  the simulated shear rate exceeds it by +15% (22 mm), +34%
  (23 mm). At the same time the torque matches to 0.92% and the
  local (tau, gamma_dot) state sits on the TC curve to ~1% in
  stress. The resolution: the TC law tau(gamma_dot) is FLAT near
  the yield stress (tau ~ sigma_y + small corrections), so a
  ~1–2% stress error — the size of the staircase-wall bias —
  maps to a 10–30% shear-rate error where gamma_dot is small.
  The likely source is the Chapman-Enskog PiNeq -> gamma_dot
  relation at this resolution: gamma_dot varies by a factor ~3
  over 1 mm (~2 cells) in the outer gap, so the higher-order
  (Burnett-level) terms in PiNeq are not small there, and the
  near-wall cells are contaminated by bounce-back. The torque,
  a stress/integral quantity, is robust to this; the local
  shear-rate field is not. A 20-cells/gap run at gamma = 10
  would test this directly (not run yet — ~8x the cost).
  Practical consequence: use the torque (flow curve) from this
  case, not the velocity profile, and do not read a yield radius
  off the simulated profile.
- **Profile-derived torque ~5–6.5% low** at every point
  (0.01113 vs 0.01167 N·m at gamma = 0.1; 0.02358 vs 0.02523 at
  gamma = 100): the staircase wall displaces the effective wall
  and biases the near-wall velocity gradient, as in the
  Newtonian case (~9% low there). The momentum-exchange torque
  is the primary measurement; the profile is a cross-check only.
- **First torque sample** in each `torque.csv` (t = 0) is an
  unphysical transient (the equilibrium rest initial condition
  read before any step). `couette_post.py` only averages the
  last 20%, so it is harmless, but do not plot the raw first
  sample as a result.
- **Resolution**: 10 cells across the gap. The 3.80% max error
  should shrink with finer resolution; not tested (a 20-cell
  run is ~8x the cost in 2D).
- **3D**: 2D cross-section only; end effects are not modelled.
- **Serial only**: no MPI in this VM; the build has `ENABLE_MPI`
  off. Multi-node / GPU runs are untested.
- **Sweep interruption**: the first `couette_sweep.py --run`
  invocation was killed (SIGTERM) after the gamma = 0.1 point
  when the agent turn that launched it ended; the remaining
  points were then run one by one in the foreground
  (0.3, 1, 3, 30, 100) plus a repeat of 0.1 for its diagnostics.
  All seven `torque.csv` files are complete (241 samples each)
  and post-processing passed on the first attempt.
