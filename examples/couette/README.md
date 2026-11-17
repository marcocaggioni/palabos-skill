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

## Honest status (v0.1)

- `lattice.py`: verified in this VM — dx 0.5 mm, dt 4.17e-5 s,
  ω = 1.0, u_lb,max = 0.030 at γ̇ = 100 s⁻¹, 24,000 steps per point.
- `couette_post.py`: verified on synthetic torque data — recovers
  η_app ≈ 1.0 Pa·s, analytic error ~0%, validation PASSED.
- `couette.cpp`: **not yet compiled.** It follows the Palabos
  showcase pattern, but the moving-wall boundary and the torque
  extraction (`computeInnerTorque`) must be pinned against the
  installed Palabos before any run is trusted. The `TODO (verify)`
  markers say exactly where. Until the Newtonian validation passes
  on a real run, do not use this for a non-Newtonian flow curve.

## Analytic check

γ̇_i = 2·Ω·R_o²/(R_o²−R_i²), τ_i = T/(2π·R_i²·L),
T_analytic = 4π·L·η·Ω·R_i²·R_o²/(R_o²−R_i²).
At γ̇ = 10 s⁻¹: Ω = 1.8 rad/s (17.2 rpm), T = 1.508e-3 N·m, τ = 10 Pa.
