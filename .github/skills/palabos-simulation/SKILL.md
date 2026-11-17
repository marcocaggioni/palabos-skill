---
name: palabos-simulation
description: "Run an open-source Palabos lattice Boltzmann simulation from an English discussion: elicit the geometry, fluid and operating conditions in plain language, turn them into a checked case file, build and run the Palabos case from the command line, and turn the wall torque back into a physical result (e.g. a Couette flow curve). Use when the user wants a free, no-license LBM simulation — a Couette rheometer, a channel, a stirred tank — or asks for the open-source counterpart to an M-Star run."
---

# Palabos Simulation

Palabos (University of Geneva) is an open-source lattice Boltzmann
framework for general-purpose CFD, AGPLv3, native interface in C++,
with Python and Java interfaces also available. It has "practically
no external dependencies (only Posix and MPI)" and is therefore easy
to deploy. A GPU-capable version is in pre-release. This skill is
the open-source counterpart of the `mstar-simulation` skill in the
private `mstar-skill` repo: same discussion → case → run → report
workflow, same Couette example, different solver — no license, no
GUI, no `.msb` file.

## The workflow, in order

```
1. Discuss (English)     -> what fluid, what geometry, what do you want to measure?
2. Write case.yaml       -> every number in one file, with units. Show it. Get a yes.
3. Lattice parameters    -> lattice.py: dx, dt, omega, lattice velocity. Check Mach < 0.1.
4. Build the case        -> C++ case in examples/<name>/, cmake + make. No GUI.
5. Sweep / run           -> one run per shear rate, torque CSV per run.
6. Post-process          -> torque -> stress, omega -> shear rate, flow curve + validation.
7. Report                -> plot + CSV + what was assumed, checked, and unverified.
```

**Never skip step 2.** Same rule as in the M-Star skill: a simulation
with an unconfirmed shear-rate range, fluid model, or geometry is a
wasted run. If the user has not named a number, ask — do not invent
one and bury it in code. Defaults in the examples are labelled
`default:` and must be confirmed.

**Running is free here** (no license seat), but it is still not free:
a 3D run at fine resolution costs real CPU time. Default to a cheap
first run (coarse lattice, short time, one shear rate) and only then
sweep.

## 1. Discuss — what to elicit

Same checklist as the M-Star skill, plus one Palabos-specific item:

### Geometry
- Inner (rotating) radius `R_i`, outer (fixed) radius `R_o`, immersed
  length `L`. Metres.
- **First version is a 2D cross-section** (annulus in the plane,
  torque per unit length) unless the user asks for 3D. The flow in an
  ideal Couette cell has no axial variation except end effects, so
  2D validates the physics at a fraction of the cost; 3D is the
  follow-up, with end effects included. Say which you are running.

### Fluid — pick one rheology, with numbers
Palabos has **no rheology dropdown**. Newtonian is built in (constant
relaxation frequency). Non-Newtonian means making the relaxation
frequency a function of the local shear rate — Carreau-Yasuda is the
documented pattern in the Palabos literature (apparent viscosity
from the local strain rate, computed from the distribution
functions). Herschel–Bulkley needs the same treatment plus a yield
stress, with a viscosity cap at low shear rate (the yield stress
makes the apparent viscosity diverge as γ̇ → 0 — cap it, and say
what cap you used). If the user names a real material, do not invent
its parameters: ask for a measured flow curve, a rheofit fit, or a
literature value, or start Newtonian as a validation case.

### Operating conditions — the sweep
- Shear-rate range and spacing, e.g. 0.1 → 100 s⁻¹, log-spaced.
- Lattice resolution: cells across the gap (default 10, i.e.
  dx = gap / 10). More cells = more accurate torque, much more cost.
- Run time per point and the steady window for averaging torque.
  Viscous diffusion time is gap²/ν — for the default case that is
  ~0.025 s, so 1 s per point is generous.

### What is being measured
- Couette: torque on the rotating bob → shear stress, rotation rate
  → shear rate. The flow curve. Name the file it comes from
  (`torque.csv`, written by the C++ case) before running.

## 2. Write case.yaml

One file per simulation, in `examples/<name>/case.yaml`. It is the
contract: the scripts read it, they do not re-derive numbers from
chat. Template: `examples/couette/case.yaml`. Show it to the user
and get a yes before building.

## 3. Lattice parameters — the step M-Star hides

Palabos works in **lattice units**. Converting physical → lattice
is the agent's job, and getting it wrong is the most common failure.
`examples/couette/lattice.py` does it; the rules are:

```
dx = gap / N_gap                     lattice spacing (m)
nu_lb = nu * dt / dx^2               lattice kinematic viscosity
omega = 1 / (3*nu_lb + 0.5)         relaxation frequency (D2Q9/D3Q19, cs^2 = 1/3)
u_lb  = u_phys * dt / dx             lattice velocity — keep max u_lb < 0.1 (low Mach)
rho_lb = 1.0                         by convention
```

Choose `dt` by fixing `omega` (default 1.0, so nu_lb = 1/6) and
solving `dt = nu_lb * dx^2 / nu`. Then check `u_lb` at the highest
shear rate. If `u_lb >= 0.1`, reduce `dt` (more steps) or coarsen
`dx` — and say which you did. Never run with `omega` outside
roughly [0.8, 1.9] without saying so: near 2 the scheme goes
unstable, near 0.5 it is over-damped.

Force/torque conversion: a lattice force converts with
`rho * dx^4 / dt^2` (per the Palabos unit-conversion pattern).
The C++ case writes torque already converted to N·m per unit
length; `couette_post.py` checks the units by comparing with the
analytic torque (step 6).

## 4. Build the case

Palabos is a C++ **library**, not an application: you write a small
program against it and compile it. The repo layout follows the
Palabos showcase pattern:

```bash
git clone https://gitlab.com/unigespc/palabos.git
cd palabos/examples/showCases/cavity2d/build
cmake ..; make        # the documented first build: compile an example
```

For this repo's case, from the repo root:

```bash
bash scripts/check-setup.sh     # cmake, g++, MPI, Palabos checkout
./run-couette.sh                # dry-run: lattice params + build only
```

The C++ case (`examples/couette/couette.cpp`) follows the showcase
pattern: `MultiBlockLattice`, BGK dynamics, `BounceBack` on solid
nodes, velocity boundary conditions on the walls
(`setVelocityConditionOnBlockBoundaries`), `collideAndStream` loop,
VTK output. The moving inner wall uses a velocity bounce-back with
local wall velocity `u_wall = Omega x r`. **Status: the C++ case
has been compiled and run against Palabos master (commit 4127697,
2026-04-03) and the Newtonian validation in §6 passes (1.43%
analytic torque error across the sweep — see
`examples/couette/RUN_REPORT.md`).** The moving wall is
`VelocityBounceBack` per node with `setUwall()` updates during the
ramp; torque is moving-wall momentum exchange at the inner-wall
links. The API names were checked against the installed Palabos'
own examples (`cavity2d`, Poiseuille, `computeDrag`) and sources.

Non-Newtonian: implement it as a local relaxation frequency,
`omega(x, t)` from the local shear rate via the chosen model
(Carreau-Yasuda first; Herschel–Bulkley second, with a viscosity
cap). Keep the Newtonian path as the default and the validation
gate — never debug a non-Newtonian model on an unvalidated geometry.

## 5. Run

```bash
./run-couette.sh --run           # one run per shear rate, in examples/couette/runs/
# or, per point:
./couette --gamma 10 --out runs/gamma_10
mpirun -np 4 ./couette --gamma 10 --out runs/gamma_10   # MPI, if built with it
```

Each run writes `torque.csv` (time, torque per unit length) and
`profile.csv` (radius, time-averaged azimuthal velocity). If a run
is slow, reduce the resolution or the run time first — do not
silently change the geometry to make it fast.

## 6. Post-process — torque to flow curve (Couette)

Same conversion as the M-Star skill, per unit length in 2D
(multiply by `L` for the 3D torque):

```
shear rate at the inner wall:  gamma_i = 2*Omega*Ro^2 / (Ro^2 - Ri^2)
shear stress at the inner wall: tau_i  = T / (2*pi*Ri^2*L)   (T = torque * L in 2D)
apparent viscosity:             eta_app = tau_i / gamma_i
```

Average torque over the steady window named in `case.yaml`, report
its standard deviation as the steadiness check.

**Validation gate (mandatory for the first run of a new geometry).**
Run a Newtonian fluid of known eta first and compare with the
analytic Couette torque:

```
T_analytic = 4*pi*L*eta*Omega*Ri^2*Ro^2 / (Ro^2 - Ri^2)
```

If the simulated torque differs by more than ~10%, do **not**
proceed to a non-Newtonian flow curve. Stop and diagnose: lattice
resolution across the gap, `omega` / `dt`, run time, moving-wall
boundary, torque extraction. Record the validation number in the
report. `couette_post.py` enforces this gate and exits non-zero on
failure.

Then each sweep point gives one (gamma_i, tau_i) point. Plot tau vs
gamma (log-log) and eta_app vs gamma, write `flow_curve.csv`, and —
if the user wants parameters — fit the curve with rheofit
(Herschel–Bulkley, Carreau, …) in a separate step.

## 7. Report

Deliver, per run: `case.yaml` (as run), `lattice.json` (dx, dt,
omega, max u_lb), `flow_curve.csv`, `flow_curve.png`, and a short
note stating the geometry (2D or 3D), fluid, sweep, validation
result, steady-state criterion, and anything unverified. If a run
failed, report the compiler / solver log tail — do not silently drop
the point from the curve.

## Starting a session with GitHub Copilot CLI

Same setup as `mstar-skill`: Copilot CLI discovers project skills in
`.github/skills/`, and this skill is mirrored there. From the repo
root:

```bash
git clone https://github.com/marcocaggioni/palabos-skill.git
cd palabos-skill
bash scripts/check-setup.sh
copilot
```

Then `/skills list` (you should see `palabos-simulation`), and:

```
/palabos-simulation run the Couette example in examples/couette/case.yaml
```

This skill deliberately has **no `allowed-tools` frontmatter**:
Copilot must ask before running shell commands, and in particular
before a long solver run. If the skill was added during a session,
run `/skills reload` first.

## What this skill does NOT do

- It does not install Palabos for you beyond `git clone` + `cmake`:
  you need a C++ compiler, CMake, and MPI on the machine.
- It does not claim the C++ case is verified until it compiles
  against an installed Palabos and the Newtonian validation passes.
- It does not replace M-Star for tanks, impellers, multiphase, or
  particle-laden flows — those are out of scope for v0.1. This repo
  is the Couette flow curve, done openly, as the first proof that
  the open route works.

## Files in this repo

- `skills/palabos-simulation/SKILL.md` — this file. Mirrored
  byte-for-byte at `.github/skills/palabos-simulation/SKILL.md`
  (the copy GitHub Copilot CLI discovers). **Keep the two copies
  in sync** — edit both in the same commit.
- `docs/palabos-guide.md` — install, build, lattice units,
  non-Newtonian notes, limits.
- `examples/couette/` — the worked Couette example: `case.yaml`,
  `lattice.py`, `couette.cpp` + `CMakeLists.txt`,
  `couette_sweep.py`, `couette_post.py`, `README.md`.
- Palabos itself: https://palabos.unige.ch / https://gitlab.com/unigespc/palabos (AGPLv3)
