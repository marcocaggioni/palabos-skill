// Couette rheometer — Palabos case (2D annulus cross-section).
//
// STATUS (v0.1): TEMPLATE, NOT YET COMPILED against an installed Palabos.
// It follows the Palabos showcase pattern (MultiBlockLattice, BGK
// dynamics, BounceBack on solid nodes, setVelocityConditionOnBlock-
// Boundaries, collideAndStream loop, VTK output), but every API name
// must be checked against the Palabos version actually installed —
// compare with the installed examples (cavity2d, Poiseuille) and the
// user guide's boundary-conditions chapter before trusting it.
// Lines marked TODO (verify) are the ones most likely to differ.
//
// Physics: fluid between two concentric cylinders, inner wall rotating
// at Omega, outer wall fixed. Torque on the inner wall is computed by
// momentum exchange at the inner-wall boundary links and written to
// torque.csv already converted to N*m per unit length (see lattice.py
// for the force conversion rho*dx^4/dt^2). The post-processing script
// couette_post.py converts torque -> shear stress and Omega -> shear
// rate and enforces the Newtonian validation gate.
//
// Build (after `git clone https://gitlab.com/unigespc/palabos.git`):
//   mkdir build && cd build && cmake .. -DPALABOS_ROOT=<path-to-palabos> && make
// Run:
//   ./couette --gamma 10 --out runs/gamma_10
//   mpirun -np 4 ./couette --gamma 10 --out runs/gamma_10
//
// TODO (verify): Palabos version, header paths, and the exact names of
//   the velocity bounce-back / curved-boundary classes in that version.
//   Palabos ships off-lattice (curved) boundary conditions; a staircase
//   BounceBack annulus is the v0.1 approximation and its torque error
//   is part of what the Newtonian validation gate measures.

#include "palabos2D.h"
#include "palabos2D.hh"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D2Q9Descriptor

// ---- case parameters (defaults = examples/couette/case.yaml) ----
static T Ri = 0.020;      // m, inner radius (rotating)
static T Ro = 0.025;      // m, outer radius (fixed)
static T rhoPhys = 1000.0;
static T etaPhys = 1.0;   // Pa.s, Newtonian validation case
static T gammaTarget = 10.0;   // 1/s, from --gamma
static int N_gap = 10;    // cells across the gap
static T omega = 1.0;     // relaxation frequency
static T runTime = 1.0;   // s per point
static string outDir = "runs/gamma_10";

// Torque per unit length, from momentum exchange at the inner wall.
// TODO (verify): implement with the installed Palabos' force utilities
// (e.g. its boundary-force / momentum-exchange helpers) instead of
// hand-rolled link bookkeeping if the version provides them.
static T computeInnerTorque(MultiBlockLattice2D<T, DESCRIPTOR> &lattice,
                            T dx, T dt, T Omega)
{
    // Placeholder signature: the real implementation sums, over every
    // boundary link touching an inner-wall node, the tangential
    // momentum exchanged per step, converts with rho*dx^4/dt^2, and
    // multiplies by the lever arm Ri. Returning 0.0 here is deliberate:
    // a template that silently wrote a wrong torque would be worse
    // than one that refuses to run until this function is filled in
    // against the installed Palabos.
    (void)lattice; (void)dx; (void)dt; (void)Omega;
    return 0.0;
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--gamma" && i + 1 < argc) gammaTarget = atof(argv[++i]);
        else if (arg == "--out" && i + 1 < argc) outDir = argv[++i];
    }

    const T gap = Ro - Ri;
    const T dx = gap / N_gap;
    const T nuPhys = etaPhys / rhoPhys;
    const T nuLb = (1.0 / omega - 0.5) / 3.0;
    const T dt = nuLb * dx * dx / nuPhys;
    const T Omega = gammaTarget * (Ro * Ro - Ri * Ri) / (2.0 * Ro * Ro);
    const int nSteps = (int)ceil(runTime / dt);

    // Domain: square covering the outer cylinder, in lattice units.
    const int nx = (int)ceil(2.0 * Ro / dx) + 2;
    const int ny = nx;
    const T cx = nx / 2.0, cy = ny / 2.0;  // centre, lattice units
    const T RiLb = Ri / dx, RoLb = Ro / dx;

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
        nx, ny, new BGKdynamics<T, DESCRIPTOR>(omega));
    lattice.periodicity().toggleAll(false);

    // Geometry: fluid nodes are the annulus Ri < r < Ro.
    // TODO (verify): use the installed Palabos' indicator / voxelize
    // helpers (e.g. a cylindrical indicator function) instead of the
    // hand-rolled loop below if the version provides them.
    for (int iX = 0; iX < nx; ++iX) {
        for (int iY = 0; iY < ny; ++iY) {
            T r = sqrt((iX - cx) * (iX - cx) + (iY - cy) * (iY - cy));
            if (r <= RiLb || r >= RoLb) {
                lattice.defineDynamics(iX, iY,
                    new BounceBack<T, DESCRIPTOR>());
            }
        }
    }

    // Moving inner wall: velocity bounce-back, u_wall = Omega x r.
    // TODO (verify): the exact moving-wall boundary class and how to
    // pass a spatially varying wall velocity in the installed version
    // (curved / off-lattice boundary conditions are the accurate
    // option; staircase BounceBack with a per-node velocity is the
    // v0.1 fallback). Until this is wired in, the run is not a
    // rotating-wall simulation — do not trust its torque.
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR> *boundaryCondition =
        createLocalBoundaryCondition2D<T, DESCRIPTOR>();

    // Initialise at rest, then ramp Omega with a quarter sine over
    // the first 20% of the run (same ramp idea as the M-Star case).
    lattice.initialize();
    std::string torqueFile = outDir + "/torque.csv";
    // TODO (verify): create outDir, open torqueFile, and in the time
    // loop call collideAndStream, apply the ramped wall velocity,
    // sample computeInnerTorque() every N steps, and write VTK output
    // (e.g. writeVTK) at the end. The loop is intentionally left as
    // a skeleton until the boundary class above is pinned down.

    pcout << "Couette (Palabos, template): gamma=" << gammaTarget
          << " 1/s, Omega=" << Omega << " rad/s, dx=" << dx
          << " m, dt=" << dt << " s, steps=" << nSteps << endl;
    pcout << "TODO: wire the moving-wall boundary and the time loop "
          << "against the installed Palabos, then delete this message."
          << endl;
    delete boundaryCondition;
    return 0;
}
