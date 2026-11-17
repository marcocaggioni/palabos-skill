# Palabos — short guide for this repo

Sources: https://palabos.unige.ch, https://gitlab.com/unigespc/palabos,
Palabos user guide (ships with the repo). Anything not in those
sources is marked as unverified in this repo — do not treat it as fact.

## What Palabos is

Open-source lattice Boltzmann framework from the University of
Geneva, AGPLv3, native interface in C++, with Python and Java
interfaces. "Practically no external dependencies (only Posix and
MPI)". A GPU-capable version is in pre-release. It is a **library**,
not an application: you write a C++ program against it and compile
it. There is no GUI, no case file format, no built-in post-processor
— output is VTK (ParaView) plus whatever your program writes.

## Install / first build

```bash
git clone https://gitlab.com/unigespc/palabos.git
cd palabos/examples/showCases/cavity2d/build
cmake ..; make
```

Requires a modern C++ compiler (gcc/clang), CMake, make, MPI for
parallel runs. CMake options include `ENABLE_MPI`, `BUILD_HDF5`.

## Lattice units — the one thing to get right

Palabos works in lattice units internally. Convert physical → lattice
before setting anything, and convert forces/torques back afterwards:

```
dx = characteristic length / resolution
dt from omega:  omega = 1 / (3*nu_lb + 0.5),  nu_lb = nu * dt / dx^2
u_lb = u_phys * dt / dx          keep < 0.1
rho_lb = 1.0 (by convention)
force conversion = rho * dx^4 / dt^2
```

`examples/couette/lattice.py` does this for the Couette case.
Choosing `omega` near 1 is a good default; near 2 the scheme goes
unstable, near 0.5 it is over-damped.

## Boundaries

- No-slip walls: `BounceBack` dynamics on solid nodes (the pattern
  used in the showcase Poiseuille example).
- Velocity boundaries: `setVelocityConditionOnBlockBoundaries`
  (also in the Poiseuille example).
- Curved / moving walls: Palabos ships off-lattice (curved)
  boundary conditions; a staircase bounce-back annulus is the v0.1
  approximation in this repo, and its torque error is measured by
  the Newtonian validation gate, not assumed away.

## Non-Newtonian fluids

There is no rheology menu. The documented pattern (in the
Palabos literature, e.g. Carreau-Yasuda foaming / cement-paste
studies) is to make the relaxation frequency a local function of
the shear rate, with the shear-rate tensor computed locally from
the distribution functions. Herschel–Bulkley needs the same plus
a yield stress and a viscosity cap at low shear rate. Implement
Newtonian first, validate it, then add the model — never the
other way round.

## Limits (v0.1)

- No shipped Taylor-Couette example in Palabos — the Couette case
  here is written from the showcase pattern, not copied.
- Post-processing is your own code (this repo's `couette_post.py`)
  or ParaView on the VTK output.
- GPU support is pre-release; the v0.1 path is CPU + MPI.
- Palabos is AGPLv3: fine for a public repo, check before using it
  inside anything you would not publish.
