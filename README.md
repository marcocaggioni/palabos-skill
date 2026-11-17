# palabos-skill

Skill teaching an agent to run open-source Palabos lattice Boltzmann
simulations from an English discussion — starting with a Couette
rheometer flow curve. **Public repo. No license needed.**

Palabos (University of Geneva) is a free LBM framework for
general-purpose CFD, AGPLv3, C++ with MPI, "practically no external
dependencies (only Posix and MPI)". This is the open-source
counterpart of the private `mstar-skill` repo: same workflow, same
Couette example, different solver — no license seat, no GUI, no
`.msb` file. Instead you write a small C++ case and compile it.

## Quickstart

```bash
git clone https://github.com/marcocaggioni/palabos-skill.git
cd palabos-skill
uv sync              # Python side only (numpy, matplotlib, pyyaml)
bash scripts/install-palabos.sh   # Palabos C++ — no PyPI package exists for it
bash scripts/check-setup.sh
copilot
```

First prompt:

```
/palabos-simulation run the Couette example in examples/couette/case.yaml
```

Or one go: `copilot -i "Use the /palabos-simulation skill to run the Couette example in examples/couette/case.yaml"`. Copilot must be started from the repo root; if the skill was added mid-session, run `/skills reload`. See `QUICKSTART.md`.

## The workflow

1. **Discuss** in English — geometry, fluid + rheology, sweep, what is measured.
2. **Write `case.yaml`**, show it, get a yes. Never skip this, never invent a number.
3. **Lattice parameters** — `lattice.py` converts physical → lattice units and checks the low-Mach condition. This is the step M-Star hides.
4. **Build** the C++ case (`cmake` + `make`). Dry-run by default.
5. **Sweep / run** — one run per shear rate, torque CSV per run.
6. **Post-process** torque → flow curve, with a mandatory Newtonian validation gate (analytic Couette torque, ~10%).
7. **Report** — plot + CSV + assumptions + validation + unverified items.

## The Couette example

Fluid in a Couette rheometer, inner cylinder rotating, outer fixed,
torque on the bob → flow curve. Defaults (labelled, to confirm):
R_i = 0.020 m, R_o = 0.025 m, L = 0.060 m, Newtonian η = 1.0 Pa·s,
ρ = 1000 kg/m³, γ̇ = 0.1 → 100 s⁻¹ log-spaced, 10 cells across the
gap. 2D cross-section first (torque per unit length); 3D is a
follow-up. Same numbers as the M-Star example, so the solvers can
be compared point for point.

## Results — it ran

<video src="https://raw.githubusercontent.com/marcocaggioni/palabos-skill/main/examples/couette/couette-sweep.mp4" controls></video>

Steady state at each shear rate, γ̇ = 0.1 → 100 s⁻¹. The spin-up at γ̇ = 10 s⁻¹ is in [couette-gamma10.mp4](examples/couette/couette-gamma10.mp4).

![Flow curve](examples/couette/flow_curve.png)

Newtonian validation **PASSED**: max analytic-torque error **1.43%**, η_app = 1.0143 Pa·s (true 1.0), constant across the sweep — a systematic staircase-wall bias. Full numbers in [examples/couette/RUN_REPORT.md](examples/couette/RUN_REPORT.md).

## How this compares

Same Couette case, same validation, five solvers. Full write-up: [docs/comparison.md](docs/comparison.md).

| Solver | Method | Repo | Validation (max error) | η_app | Sweep video |
|---|---|---|---|---|---|
| M-Star | LBM, GPU, licensed | [mstar-skill](https://github.com/marcocaggioni/mstar-skill) (private) | not yet run | — | — |
| Palabos | LBM, C++/MPI, AGPLv3 | [palabos-skill](https://github.com/marcocaggioni/palabos-skill) | 1.43% | 1.0143 | [video](https://github.com/marcocaggioni/palabos-skill/blob/main/examples/couette/couette-sweep.mp4) |
| XLB | LBM, Python/JAX | [xlb-skill](https://github.com/marcocaggioni/xlb-skill) | 1.16% | 0.988–1.000 | [video](https://github.com/marcocaggioni/xlb-skill/blob/main/examples/couette/couette-sweep.mp4) |
| OpenFOAM | Finite-volume, GPL | [openfoam-skill](https://github.com/marcocaggioni/openfoam-skill) | 0.03% | 1.0003 | [video](https://github.com/marcocaggioni/openfoam-skill/blob/main/examples/couette/couette-sweep.mp4) |
| OpenLB | LBM, C++/CUDA, GPLv2 | [openlb-skill](https://github.com/marcocaggioni/openlb-skill) | running | — | — |

## Status — v0.1

- `lattice.py` and `couette_post.py`: **tested in this VM** (lattice
  check passes; post-processing recovers η ≈ 1.0 Pa·s on synthetic
  torque, validation PASSED).
- `couette.cpp`: **compiled and run** against Palabos master
  (commit 4127697, 2026-04-03) on 2026-11-17. Moving wall =
  `VelocityBounceBack` + `setUwall()` ramp; torque = moving-wall
  momentum exchange. Full sweep validated: eta_app = 1.0143 Pa.s
  (true 1.0), analytic torque error 1.43% at every point — see
  `examples/couette/RUN_REPORT.md`, `flow_curve.csv/png`, and the
  videos `couette-gamma10.mp4` / `couette-sweep.mp4`.
- Still open: non-Newtonian fluids and 3D (2D cross-section here).

## Files

- `skills/palabos-simulation/SKILL.md` — the skill, mirrored at `.github/skills/palabos-simulation/SKILL.md`
- `docs/palabos-guide.md` — install, build, lattice units, non-Newtonian notes, limits
- `examples/couette/` — the worked example
- `scripts/check-setup.sh` — readiness check (copilot, cmake, g++, MPI, Palabos)

Palabos: https://palabos.unige.ch — https://gitlab.com/unigespc/palabos (AGPLv3).
This repo's own code: MIT.
