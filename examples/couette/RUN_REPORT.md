# Run report — Palabos Couette example (real run, 2026-11-17)

First end-to-end run of the Couette example on a real Palabos install.
Newtonian validation **PASSED**. Nothing in this file is synthetic:
every number below comes from the Palabos runs in `runs/`.

## Environment

- Palabos: cloned 2026-11-17 from https://gitlab.com/unigespc/palabos.git,
  `--depth 1`, master at commit
  `4127697e90169bbef982295f1d1c933cf6e90caa` (2026-04-03, "Merge branch
  'weno3_new' into 'master'"). Checkout at `~/palabos` (80 MB),
  `PALABOS_ROOT=~/palabos`.
- First clone attempt failed with a TLS error (`curl 56 GnuTLS recv
  error (-110)`); a plain retry succeeded.
- Machine: this VM, 4 CPUs, ~7.7 GB RAM. Serial build, no MPI
  (`ENABLE_MPI=OFF` — no MPI installed here). g++ 13.3, CMake 3.28,
  C++20. `libtbb-dev` 2021.11.0 was installed via apt (Palabos'
  examples require TBB); nothing else was installed.
- Toolchain proved first with Palabos' own `cavity2d` example
  (built clean, binary `examples/showCases/cavity2d/cavity2d`).
- Python via the repo's uv venv (`uv run python ...` from the repo
  root); `uv sync` already set up in the repo.

## Case (unchanged from case.yaml)

2D annulus cross-section, Ri = 0.020 m, Ro = 0.025 m, L = 0.060 m,
Newtonian eta = 1.0 Pa·s, rho = 1000 kg/m³, sweep
γ̇ = [0.1, 0.3, 1, 3, 10, 30, 100] s⁻¹, run time 1.0 s/point,
quarter-sine ramp 0.2 s, steady window last 20%, 10 cells across the
gap, omega = 1.0. Lattice (from `lattice.py`, verified):
dx = 0.5 mm, dt = 4.167e-5 s, 24,000 steps/point, u_lb,max = 0.030.
Grid 102×102, 312 inner-wall nodes.

## How the moving wall and torque are implemented

All API names below were read from the cloned Palabos sources, not
invented. The C++ case is `couette.cpp` (D2Q9, BGK, staircase annulus:
fluid nodes Ri < r < Ro; outer solid = plain `BounceBack`).

### Moving inner wall — `VelocityBounceBack` + `setUwall()`

- `VelocityBounceBack` dynamics: `src/core/dynamics.h:765`
  (class), `setUwall()` at `src/core/dynamics.h:786`.
  Its `collide()` (`src/core/dynamics.hh:1165`) does the standard
  moving-wall bounce-back: swap the opposite populations and add the
  correction `2·rho·w_i·(c_i·u_w)/cs²`.
- Each inner-wall node gets its own `VelocityBounceBack(rhoWall=1,
  uWall)` with the local wall velocity `u_wall = Omega × r`
  (lattice units). The wall velocity is updated every step during the
  ramp via `setUwall()` on the dynamics object fetched from the
  lattice cell (`lattice.get(iX,iY).getDynamics()`), with
  `Omega(t) = Omega_final·sin(pi/2 · t/0.2 s)`, then held constant.
- The initial condition is the analytic Couette profile, so the run
  starts close to steady state.

### Torque — moving-wall momentum exchange

- Generalisation of Palabos' own `MomentumExchangeBounceBack`
  (`src/boundaryCondition/bounceBackModels.h:61`, collide at
  `src/boundaryCondition/bounceBackModels.hh:112`): per inner-wall
  link, `F_wall = -c_i·(2·f_in + C_i)`, with
  `C_i = 2·rho·w_i·(c_i·u_w)/cs²` and `f_in` the population that just
  streamed in from the fluid. For `u_w = 0` this reduces exactly to
  Palabos' formula (force = `-2·c_i·f_in`). The fluid-direction
  lists are built the same way Palabos builds them in
  `initializeMomentumExchange`
  (`src/boundaryCondition/bounceBackModels2D.hh`: a direction counts
  if the neighbour in that direction is a fluid node). Reference
  example: `examples/codesByTopic/bounceBack/computeDrag.cpp`.
- Torque = Σ over wall nodes of `r × F`. Sampling is done right
  after `collideAndStream()` — populations are then post-streaming,
  i.e. the state Palabos' own momentum-exchange collide reads
  (`BlockLattice2D::collideAndStream` collides first, then streams:
  `src/atomicBlock/blockLattice2D.hh:249`).
- Sign convention: the raw sum is the torque the fluid exerts on
  the bob (opposes rotation, negative for Omega > 0). The code
  returns the driving torque (negated) so it is positive and
  directly comparable with the analytic torque and `couette_post.py`.
- **Unit conversion (the subtle one):** this is a 2D run, so the
  torque is per unit length and has units of force. The lattice
  torque is multiplied by `rho·dx⁴/dt²` (the force conversion in
  `lattice.py`), NOT `rho·dx⁵/dt²` (that would be a 3D torque).
  Also: Palabos stores populations as deviations from the rest
  weights (`f_i − w_i`); the missing `w_i` part is a pure pressure
  force (radial, no torque), so the deviation form is correct here.
- Independent cross-check: torque from the simulated azimuthal
  velocity profile — shear rate at the inner wall from the two
  innermost fluid bins, `tau = −eta·r·d(omega)/dr`, times
  `2·pi·Ri²·L`. Written to `profile.csv` (radius, simulated and
  analytic `v_theta`).

### Bugs found and fixed during this run

1. **Dangling dynamics pointer (segfault after step 0).**
   `defineDynamics()` clones the dynamics into the cell and deletes
   the object passed in
   (`InstantiateDotDynamicsFunctional2D::process`,
   `src/dataProcessors/dataInitializerFunctional2D.hh:491`:
   `lattice.attributeDynamics(dot.x, dot.y, dynamics->clone())`).
   The first version stored the passed `VelocityBounceBack*` and
   called `setUwall()` on it — use-after-free. Fixed by fetching
   the live dynamics from the lattice cell each step.
   (Per-cell `defineDynamics` overload:
   `src/dataProcessors/dataInitializerWrapper2D.h:219`.)
2. **Wrong torque conversion (100× error, then sign).** First
   version used `rho·dx⁵/dt²` and got a torque ~2000× too small;
   dimensional analysis against the stress (`tau_lb·2·pi·Ri_lb²`)
   showed the 2D per-length torque needs `rho·dx⁴/dt²`. The raw
   sign was also negative (fluid drag opposes rotation); negated
   to the driving-torque convention. After both fixes: 1.43% error.
3. **Binary location.** The first `CMakeLists.txt` put the binary
   in the source dir, but `couette_sweep.py` expects
   `build/couette`. Fixed in `CMakeLists.txt` (no
   `CMAKE_RUNTIME_OUTPUT_DIRECTORY` override). The old
   `CMakeLists.txt` also referenced a non-existent
   `${PALABOS_ROOT}/cmake/Palabos.cmake`; the current one follows
   the cavity2d pattern (glob Palabos sources into a static lib).

## Gamma = 10 check (first run, before the sweep)

- Simulated steady torque (momentum exchange): **0.00152945 N·m**
  (per-length 0.0254909 N·m/m)
- Analytic torque: 0.00150796 N·m
- **Relative error: 1.42%** (gate: 10%)
- Torque from velocity profile: 0.00137 N·m (9.1% low — see
  "Still unverified / known limitations")
- Velocity profile vs analytic: 0.1–0.5% error at every radius
  except the outermost bin (−5.6%, staircase outer wall)
- Run time: ~6 s for 24,000 steps

## Full sweep (2026-11-17, `couette_sweep.py --run`, 43 s total)

`couette_post.py` output — **validation PASSED** (max analytic error
1.43%, tolerance 10%):

| gamma (1/s) | tau (Pa) | eta_app (Pa·s) | torque (N·m) | analytic err |
|---:|---:|---:|---:|---:|
| 0.1 | 0.1014 | 1.0143 | 1.529e-05 | 1.43% |
| 0.3 | 0.3043 | 1.0142 | 4.588e-05 | 1.42% |
| 1 | 1.0143 | 1.0143 | 1.529e-04 | 1.43% |
| 3 | 3.0427 | 1.0142 | 4.588e-04 | 1.42% |
| 10 | 10.1425 | 1.0143 | 1.529e-03 | 1.43% |
| 30 | 30.4275 | 1.0142 | 4.588e-03 | 1.42% |
| 100 | 101.4251 | 1.0143 | 1.529e-02 | 1.43% |

eta_app = 1.0143 Pa·s at every point (true value 1.0). The error is
essentially constant across the sweep, i.e. systematic, not noise:
consistent with the staircase wall making the effective inner radius
very slightly larger than the nominal Ri = 0.020 m. The flow curve
is a straight line of slope 1 on a log-log plot (Newtonian, as it
should be). Files: `flow_curve.csv`, `flow_curve.png`.

## Videos

Rendered from the field snapshots in `runs/` by `make_video.py`
(matplotlib + ffmpeg, H.264, 600×600, 10 fps):

- `couette-gamma10.mp4` (2.5 s, 59 KB) — velocity field vs time at
  γ̇ = 10 s⁻¹, 25 snapshots (every 1000 steps). Spin-up to steady
  rotation; the bob renders dark (solid-body nodes carry no fluid
  velocity), the sheared annulus is bright, speed highest at the
  inner wall.
- `couette-sweep.mp4` (7.0 s, 81 KB) — steady-state field at each
  sweep shear rate, 1 s per point, γ̇ = 0.1 → 100 s⁻¹.
- Snapshots: `runs/gamma_*/field_*.csv`, 25 per run, 51×51 grid
  (every 2nd node), columns `iX,iY,ux_m_s,uy_m_s`. One frame was
  visually checked (steady state, γ̇ = 10): clean rotating annulus,
  tangential arrows, speed falling from inner to outer wall.

## Files left in the repo (not pushed)

- `couette.cpp`, `CMakeLists.txt` (working, as run here)
- `runs/gamma_*/torque.csv` (241 samples per run: time,
  torque per unit length) and `profile.csv`
- `runs/gamma_*/field_*.csv` (25 snapshots per run)
- `flow_curve.csv`, `flow_curve.png`, `lattice.json`
- `make_video.py`, `couette-gamma10.mp4`, `couette-sweep.mp4`
- this file

Note: `.gitignore` currently excludes `runs/` and
`flow_curve.csv/png` (and `lattice.json`), so a plain `git add -A`
would NOT include the run outputs. The videos are not ignored.
Decide before pushing whether the run outputs should be committed
or published another way (e.g. a GitHub release asset).

## Still unverified / known limitations

- **Non-Newtonian fluids**: no model implemented yet — Newtonian
  only. The documented next step is a local relaxation frequency
  from the shear rate (Carreau-Yasuda pattern; see
  `docs/palabos-guide.md`).
- **3D**: 2D cross-section only; end effects are not modelled.
- **Profile-derived torque ~9% low** at γ̇ = 10 (0.00137 vs
  0.00151 N·m): the staircase wall displaces the effective wall
  and biases the near-wall velocity gradient. The momentum-exchange
  torque (1.4%) is the primary measurement; the profile is a
  cross-check only.
- **First torque sample** in each `torque.csv` (t = 0) is an
  unphysical transient (up to −0.49 N·m/m) — it is the equilibrium
  initial condition read before any step. `couette_post.py` only
  averages the last 20%, so it is harmless, but do not plot the
  raw first sample as a result.
- **Resolution**: 10 cells across the gap. The 1.43% systematic
  error should shrink with finer resolution; not tested (a 20-cell
  run is ~8× the cost in 2D).
- **Serial only**: no MPI in this VM; the build has `ENABLE_MPI`
  off. Multi-node / GPU runs are untested.
