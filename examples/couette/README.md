# Couette example — Palabos

Fluid in a Couette rheometer: inner cylinder rotates, outer cylinder
fixed, torque on the bob gives the flow curve. Same geometry and
fluid as the M-Star example (`marcocaggioni/mstar-skill`), so the two
solvers can be compared point for point.

## What is here

| File | What it does |
|---|---|
| `case.yaml` | The contract: geometry, fluid, sweep, lattice. Confirm before running. |
| `lattice.py` | Physical → lattice parameters (dx, dt, omega, Mach check). **Tested here.** |
| `couette.cpp` + `CMakeLists.txt` | Palabos case (2D annulus). **Template — not yet compiled.** See below. |
| `couette_sweep.py` | Runs the sweep, one run per shear rate. Dry-run by default. |
| `couette_post.py` | Torque → flow curve + Newtonian validation gate. **Tested on synthetic torque.** |

## Defaults

R_i = 0.020 m, R_o = 0.025 m, L = 0.060 m, Newtonian η = 1.0 Pa·s,
ρ = 1000 kg/m³, γ̇ = 0.1 → 100 s⁻¹ log-spaced, 10 cells across the
5 mm gap, ω = 1.0, 1 s per point, steady window = last 20%.
2D cross-section first (torque per unit length); 3D is a follow-up.

## Run it

```bash
bash ../../scripts/check-setup.sh
python3 lattice.py            # dx, dt, omega, Mach check — must PASS
mkdir build && cd build && cmake .. -DPALABOS_ROOT=<path> && make
cd .. && python3 couette_sweep.py --run
python3 couette_post.py       # flow_curve.csv + flow_curve.png, validation gate
```

## Honest status (v0.1 — run for real, 2026-11-17)

- `lattice.py`: verified — dx 0.5 mm, dt 4.17e-5 s, omega = 1.0,
  u_lb,max = 0.030 at shear rate 100 1/s, 24,000 steps per point.
- `couette.cpp`: **compiled and run** against Palabos master
  (commit 4127697, 2026-04-03). Moving inner wall =
  `VelocityBounceBack` per node with a `setUwall()` ramp; torque =
  moving-wall momentum exchange at inner-wall links (see
  `RUN_REPORT.md` for the API references and the two bugs found
  on the way: a dangling dynamics pointer, and a 2D torque
  conversion that must use rho*dx^4/dt^2, not dx^5).
- **Newtonian validation PASSED on a real run**: full sweep
  shear rate 0.1 -> 100 1/s, eta_app = 1.0143 Pa.s at every point
  (true 1.0), analytic torque error 1.43% at every point — see
  `flow_curve.csv`, `flow_curve.png`, `RUN_REPORT.md`.
- Videos: `couette-gamma10.mp4` (spin-up at shear rate 10 1/s),
  `couette-sweep.mp4` (steady state across the sweep) — rendered
  from the field snapshots in `runs/` by `make_video.py`.
- Still unverified: non-Newtonian fluids (no model implemented
  yet — Newtonian only), 3D (2D cross-section here), and the
  profile-derived torque, which reads ~9% low (staircase-wall
  effect on the near-wall velocity gradient; the momentum-exchange
  torque is the primary measurement).

## Analytic check

γ̇_i = 2·Ω·R_o²/(R_o²−R_i²), τ_i = T/(2π·R_i²·L),
T_analytic = 4π·L·η·Ω·R_i²·R_o²/(R_o²−R_i²).
At γ̇ = 10 s⁻¹: Ω = 1.8 rad/s (17.2 rpm), T = 1.508e-3 N·m, τ = 10 Pa.
