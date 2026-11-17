# Couette example, three-component (TC) fluid — Palabos

Fluid in a Couette rheometer: inner cylinder rotates, outer cylinder
fixed, torque on the bob gives the flow curve. Same geometry and
sweep as the Newtonian example (`examples/couette/`), so the two can
be compared point for point — but the fluid is the TC law instead of
Newtonian:

```
tau = sigma_y + sigma_y*sqrt(gamma_dot/gamma_dot_c) + eta_bg*gamma_dot
sigma_y = 70 Pa, gamma_dot_c = 100 1/s, eta_bg = 0.05 Pa s (50 cP)
```

## What is here

| File | What it does |
|---|---|
| `case.yaml` | The contract: geometry, TC fluid + regularization, sweep, lattice. Confirm before running. |
| `lattice.py` | Physical → lattice parameters, PER POINT (dt is solved from the TC viscosity at the label rate). **Tested here.** |
| `couette.cpp` + `CMakeLists.txt` | Palabos case (2D annulus) with `TCBGKdynamics`: BGK with a local relaxation frequency from the TC law, modelled on Palabos' own `BGKCarreauDynamics`. **Compiled and run.** |
| `couette_sweep.py` | Runs the sweep, one run per shear rate. Dry-run by default. |
| `couette_post.py` | Torque → flow curve + TC validation gate (torque vs analytic TC torque). **Tested here.** |
| `make_video.py` | Renders `couette-tc-gamma10.mp4` and `couette-tc-sweep.mp4` from the field snapshots in `runs/`. |

## Defaults

R_i = 0.020 m, R_o = 0.025 m, L = 0.060 m, ρ = 1000 kg/m³,
γ̇ = 0.1 → 100 s⁻¹ log-spaced (labels use the Newtonian
inner-wall-rate convention, Ω = γ̇(R_o²−R_i²)/(2R_o²)), 10 cells
across the 5 mm gap, ω_target = 1.0 at the label viscosity,
24,000 steps per point, quarter-sine ramp over the first 20% of the
steps, steady window = last 20%. 2D cross-section (torque per unit
length); 3D is a follow-up. Initial condition: rest (u = 0) — the
Newtonian analytic profile is not the TC profile, so the run spins
up from rest under the ramp.

## Regularization (required)

The TC apparent viscosity
η(γ̇) = σ_y/γ̇ + σ_y/√(γ̇·γ̇_c) + η_bg diverges as γ̇ → 0, so the
unyielded plug would need infinite viscosity. The run caps it:
η = min(η(γ̇), η_max), η_max = η(10⁻³ s⁻¹) = 7.0221e4 Pa·s. The plug
then creeps instead of being exactly rigid; the torque is set by the
yielded region and is insensitive to the exact cap (the cap only has
to be ≫ the yielded viscosity, ~1.5–722 Pa·s).

## Lattice — dt is per point

η at the label rate spans 722.19 Pa·s (γ̇ = 0.1) down to 1.45 Pa·s
(γ̇ = 100), so one dt cannot serve the whole sweep: dt is solved
per point from the label viscosity at ω_target = 1.0
(dt = ν_lb·dx²/(η_label/ρ), ν_lb = 1/6). The step count (24,000),
not the physical run time, is held fixed — 24,000 steps is ~40
viscous diffusion times at the label viscosity at every point.

| γ̇ (1/s) | η_label (Pa·s) | dt (s) | run time (s) | ω in plug | u_lb |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 722.19 | 5.77e-8 | 1.38e-3 | 0.0204 | 0.00000 |
| 0.3 | 246.16 | 1.69e-7 | 4.06e-3 | 0.0070 | 0.00000 |
| 1 | 77.05 | 5.41e-7 | 1.30e-2 | 0.0022 | 0.00000 |
| 3 | 27.42 | 1.52e-6 | 3.65e-2 | 0.0008 | 0.00003 |
| 10 | 9.26 | 4.50e-6 | 1.08e-1 | 0.0003 | 0.00032 |
| 30 | 3.66 | 1.14e-5 | 2.73e-1 | 0.0001 | 0.00246 |
| 100 | 1.45 | 2.87e-5 | 6.90e-1 | 0.00004 | 0.02069 |

Local ω in the run spans ~1.6 in the yielded zone near the bob
down to ω(η_max) in the capped plug (see the table). The low end is
the plug behaving as a near-solid, which is the intent, not an
instability. u_lb,max = 0.021 (γ̇ = 100), well under 0.1.
(Correction, 2026-11-17: an early brief wrote η(100) = 70.75 Pa·s,
conflating a stress term with a viscosity; the correct value is
1.45 Pa·s and all runs use it. A dt solved from 70.75 would have
pushed ω toward ~2.0.)

## Run it

```bash
bash ../../scripts/check-setup.sh
python3 lattice.py            # per-point dx, dt, omega, Mach check — must PASS
mkdir build && cd build && cmake .. -DPALABOS_ROOT=<path> && make
cd .. && python3 couette_sweep.py --run
python3 couette_post.py       # flow_curve.csv + flow_curve.png, validation gate
python3 make_video.py         # couette-tc-gamma10.mp4, couette-tc-sweep.mp4
```

## Honest status (v0.1 — run for real, 2026-11-17)

- `lattice.py`: verified — per-point table above; `lattice.json`
  written next to `case.yaml`.
- `couette.cpp`: **compiled and run** against Palabos master
  (commit 4127697, 2026-04-03), the same checkout as the Newtonian
  example. TC rheology = `TCBGKdynamics`, a local relaxation
  frequency from the local shear rate via Π^neq (Chapman-Enskog),
  fixed-point iterated exactly like Palabos' own
  `BGKCarreauDynamics` (see the header comment in `couette.cpp`
  for the derivation and the references).
- **TC validation PASSED on a real run**: full sweep γ̇ = 0.1 →
  100 s⁻¹, max analytic-torque error **3.80%** (gate 10%), at
  γ̇ = 0.1; errors fall monotonically with γ̇ (2.89, 2.10, 1.57,
  0.92, 0.69, 0.66%) — a systematic positive bias, largest where
  the yielded zone is thinnest (~2 cells at γ̇ = 0.1), consistent
  with the staircase-wall bias the Newtonian run showed (+1.43%).
  See `flow_curve.csv`, `flow_curve.png`, `RUN_REPORT.md`.
- The plug is real: at γ̇ ≤ 10 part of the gap is unyielded and
  attached to the outer wall (R_y = 21.03 mm at γ̇ = 0.1 up to
  24.84 mm at γ̇ = 10); γ̇ ≥ 30 is fully yielded. Validation is
  torque vs analytic TC torque at the same Ω — not η_app vs the
  TC law at the label rate (the inner-wall rate exceeds the label
  rate).
- Videos: `couette-tc-gamma10.mp4` (spin-up at γ̇ = 10 s⁻¹),
  `couette-tc-sweep.mp4` (steady state across the sweep) —
  rendered from the field snapshots in `runs/` by `make_video.py`.
  The sheared annulus hugs the bob; the plug near the outer wall
  is dark/slow.
- Still unverified: 3D (2D cross-section here), finer resolution
  (the profile check below suggests 10 cells/gap is coarse for
  this fluid — see `RUN_REPORT.md`), and the profile-derived
  torque, which reads ~5–6.5% low (staircase-wall effect on the
  near-wall velocity gradient, as in the Newtonian case where it
  read ~9% low; the momentum-exchange torque is the primary
  measurement).

## Analytic check

τ(r) = T/(2π·L·r²); invert the TC law for γ̇(τ) (quadratic in
√γ̇; γ̇ = 0 where τ < σ_y); Ω = ∫_{R_i}^{R_y} γ̇ dr/r with
R_y = min(R_o, √(T/(2π·L·σ_y))); solve for T given Ω by
bisection. Reference torques (N·m): γ̇ = 0.1 → 1.167063e-02;
0.3 → 1.218414e-02; 1 → 1.303766e-02; 3 → 1.423130e-02;
10 → 1.627965e-02; 30 → 1.931375e-02; 100 → 2.523402e-02.
At γ̇ = 10 s⁻¹: Ω = 1.8 rad/s (17.2 rpm), T = 1.628e-2 N·m,
τ_i = 108.0 Pa, R_y = 24.84 mm.
