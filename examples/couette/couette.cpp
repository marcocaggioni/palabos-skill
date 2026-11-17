// Couette rheometer — Palabos case (2D annulus cross-section).
//
// WORKING case, built and validated against the analytic Couette torque
// (see RUN_REPORT.md). Palabos version/commit used: see RUN_REPORT.md.
//
// Physics: fluid between two concentric cylinders, inner wall rotating
// at Omega, outer wall fixed. Torque on the inner wall is computed by
// momentum exchange at the inner-wall links, converted to N*m per unit
// length. The post-processing script couette_post.py converts torque ->
// shear stress and Omega -> shear rate and enforces the Newtonian
// validation gate.
//
// Method (all API names verified against the cloned Palabos sources):
//  - D2Q9, BGK dynamics (palabos2D.h, as in examples/showCases/poiseuille).
//  - Staircase annulus: fluid nodes Ri < r < Ro; solid nodes elsewhere.
//  - Outer wall: plain BounceBack dynamics (src/core/dynamics.h).
//  - Inner wall: VelocityBounceBack(rhoWall=1, uWall) per node
//    (src/core/dynamics.h, class VelocityBounceBack). Its collide() does
//    the standard moving-wall bounce-back: swap + correction
//    2*rho*w_i*(c_i.u_w)/cs^2. Wall velocity is updated every step via
//    setUwall() during the quarter-sine ramp, then held constant.
//  - Torque: moving-wall generalisation of Palabos' own
//    MomentumExchangeBounceBack (src/boundaryCondition/bounceBackModels.hh):
//    per inner-wall link, F_wall = -c_i * (2*f_in + C_i), with
//    C_i = 2*rho*w_i*(c_i.u_w)/cs^2 and f_in the population that just
//    streamed in from the fluid. For u_w = 0 this reduces exactly to
//    Palabos' formula (force = -2*c_i*f_in). Torque = sum over wall
//    nodes of r x F. NOTE on units: Palabos stores populations as
//    deviations from the rest weights (f_i - w_i); the missing w_i
//    part is a pure pressure force (radial, no torque), so the
//    deviation form is correct here. Conversion: this is a 2D run, so
//    the torque is PER UNIT LENGTH and has units of force: multiply
//    the lattice torque by rho*dx^4/dt^2 (the force conversion from
//    lattice.py), NOT rho*dx^5/dt^2 (that would be a 3D torque).
//    Verified at gamma=10: lattice torque ~0.70 -> 0.0251 N*m/m,
//    analytic 0.02513 N*m/m (see RUN_REPORT.md).
//    A second, independent torque is computed from the simulated
//    azimuthal velocity profile (shear rate at the inner wall) and both
//    are printed for cross-checking.
//  - Initial condition: analytic Couette profile, so the run starts
//    close to steady state; the ramp then only removes the start-up shock.
//
// Build (serial, no MPI):
//   mkdir -p build && cd build && cmake .. -DPALABOS_ROOT=$HOME/palabos && make -j4
// Run:
//   ./couette --gamma 10 --out runs/gamma_10
//
// Defaults below are examples/couette/case.yaml (do not change the case
// by editing numbers here — change case.yaml and keep this file in sync).

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "latticeBoltzmann/indexTemplates.h"
#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
using namespace plb::descriptors;
using namespace indexTemplates;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

// ---- case parameters (case.yaml) ----
static T Ri = 0.020;    // m, inner radius (rotating)
static T Ro = 0.025;    // m, outer radius (fixed)
static T Lphys = 0.060; // m, wetted length (only used by couette_post.py)
static T rhoPhys = 1000.0;
static T etaPhys = 1.0; // Pa.s, Newtonian validation case
static T gammaTarget = 10.0; // 1/s, from --gamma
static int N_gap = 10;  // cells across the gap
static T omega = 1.0;   // relaxation frequency
static T runTime = 1.0; // s per point
static T rampTime = 0.2; // s, quarter-sine ramp of Omega
static string outDir = "runs/gamma_10";

// ---- derived lattice parameters (same formulas as lattice.py) ----
static T dx, dt, nuLb, Omega, OmegaLb, torqueConv;
static plint nx, ny, nSteps, rampSteps;
static T cx, cy; // centre, lattice units
static T RiLb, RoLb;

struct WallNode {
    plint iX, iY;
    T xRel, yRel; // position relative to centre, lattice units
    vector<plint> fluidDirs; // iPop with c_iPop pointing to a fluid neighbour
    // NOTE: do NOT store a VelocityBounceBack* here. defineDynamics()
    // clones the dynamics into the cell and deletes the object passed
    // in, so a stored pointer would dangle. Fetch the live dynamics
    // from the lattice cell (lattice.get(iX,iY).getDynamics()) instead.
};

static vector<WallNode> innerWall;

// Analytic Couette azimuthal velocity (m/s) at radius r (m).
static T couetteVtheta(T r)
{
    T A = Omega * Ri * Ri / (Ro * Ro - Ri * Ri);
    return A * (Ro * Ro / r - r);
}

// Initial-condition functional: analytic profile in fluid, rest elsewhere.
struct AnalyticCouette {
    void operator()(plint iX, plint iY, T &rho, Array<T, 2> &u) const
    {
        rho = (T)1;
        T x = (T)iX - cx, y = (T)iY - cy;
        T rLb = std::sqrt(x * x + y * y);
        T r = rLb * dx;
        if (rLb > RiLb && rLb < RoLb && r > (T)0) {
            T vTheta = couetteVtheta(r); // m/s
            T uLb = vTheta * dt / dx;    // lattice units
            u[0] = -uLb * y / rLb;
            u[1] = uLb * x / rLb;
        } else {
            u[0] = (T)0;
            u[1] = (T)0;
        }
    }
};

static bool isFluid(plint iX, plint iY)
{
    T x = (T)iX - cx, y = (T)iY - cy;
    T r = std::sqrt(x * x + y * y);
    return r > RiLb && r < RoLb;
}

// Wall velocity (lattice units) at a wall node, for a given Omega_lb.
static Array<T, 2> wallVelocity(T xRel, T yRel, T OmegaLb_)
{
    Array<T, 2> u;
    u[0] = -OmegaLb_ * yRel;
    u[1] = OmegaLb_ * xRel;
    return u;
}

static T rampFactorCurrent = (T)1; // 0..1, set each step
static const T PI = std::acos((T)-1);

// Torque on the inner wall (N*m per unit length) from momentum exchange,
// sampled right after collideAndStream (populations just streamed in).
static T computeTorque(MultiBlockLattice2D<T, DESCRIPTOR> &lattice)
{
    T torqueLb = (T)0;
    for (auto &wn : innerWall) {
        Cell<T, DESCRIPTOR> &cell = lattice.get(wn.iX, wn.iY);
        Array<T, 2> uW = wallVelocity(wn.xRel, wn.yRel, OmegaLb * rampFactorCurrent);
        Array<T, 2> F;
        F.resetToZero();
        for (plint iPop : wn.fluidDirs) {
            plint iOpp = opposite<DESCRIPTOR<T> >(iPop);
            T fIn = cell[iOpp];
            T dot = (T)0;
            for (int iD = 0; iD < 2; ++iD) dot += DESCRIPTOR<T>::c[iPop][iD] * uW[iD];
            T C = (T)2 * (T)1.0 * DESCRIPTOR<T>::t[iPop] * DESCRIPTOR<T>::invCs2 * dot;
            for (int iD = 0; iD < 2; ++iD) {
                F[iD] -= DESCRIPTOR<T>::c[iPop][iD] * ((T)2 * fIn + C);
            }
        }
        torqueLb += wn.xRel * F[1] - wn.yRel * F[0];
    }
    // Sign convention: the raw sum is the torque the fluid exerts ON the
    // bob, which opposes the rotation (negative for Omega > 0). Return
    // the driving torque instead — the torque that must be applied to
    // the bob to sustain the rotation — so it is positive and directly
    // comparable with the analytic Couette torque and with couette_post.py.
    return -torqueLb * torqueConv;
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--gamma" && i + 1 < argc) gammaTarget = atof(argv[++i]);
        else if (arg == "--out" && i + 1 < argc) outDir = argv[++i];
    }
    (void)Lphys;

    // Lattice parameters (identical to lattice.py).
    const T gap = Ro - Ri;
    dx = gap / N_gap;
    const T nuPhys = etaPhys / rhoPhys;
    nuLb = (1.0 / omega - 0.5) / 3.0;
    dt = nuLb * dx * dx / nuPhys;
    Omega = gammaTarget * (Ro * Ro - Ri * Ri) / (2.0 * Ro * Ro);
    OmegaLb = Omega * dt;
    nSteps = (plint)std::ceil(runTime / dt);
    rampSteps = (plint)std::ceil(rampTime / dt);
    torqueConv = rhoPhys * std::pow(dx, 4) / (dt * dt); // lattice torque -> N*m/m (2D: force-like)
    nx = (plint)std::ceil(2.0 * Ro / dx) + 2;
    ny = nx;
    cx = nx / 2.0;
    cy = ny / 2.0;
    RiLb = Ri / dx;
    RoLb = Ro / dx;

    pcout << "Couette (Palabos): gamma=" << gammaTarget << " 1/s, Omega=" << Omega
          << " rad/s, dx=" << dx << " m, dt=" << dt << " s, omega=" << omega
          << ", steps=" << nSteps << ", grid=" << nx << "x" << ny << endl;

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
        nx, ny, new BGKdynamics<T, DESCRIPTOR>(omega));
    lattice.periodicity().toggleAll(false);

    // Geometry: inner solid gets VelocityBounceBack, outer solid BounceBack.
    for (plint iX = 0; iX < nx; ++iX) {
        for (plint iY = 0; iY < ny; ++iY) {
            T x = (T)iX - cx, y = (T)iY - cy;
            T r = std::sqrt(x * x + y * y);
            if (r <= RiLb) {
                Array<T, 2> u0 = wallVelocity(x, y, (T)0);
                auto *dyn = new VelocityBounceBack<T, DESCRIPTOR>((T)1.0, u0);
                defineDynamics(lattice, iX, iY, dyn);
                // Record wall nodes (inner solid with a fluid neighbour).
                vector<plint> fluidDirs;
                for (plint iPop = 1; iPop < DESCRIPTOR<T>::q; ++iPop) {
                    plint nX = iX + DESCRIPTOR<T>::c[iPop][0];
                    plint nY = iY + DESCRIPTOR<T>::c[iPop][1];
                    if (nX >= 0 && nX < nx && nY >= 0 && nY < ny && isFluid(nX, nY)) {
                        fluidDirs.push_back(iPop);
                    }
                }
                if (!fluidDirs.empty()) {
                    innerWall.push_back({iX, iY, x, y, fluidDirs});
                }
            } else if (r >= RoLb) {
                defineDynamics(lattice, iX, iY, new BounceBack<T, DESCRIPTOR>());
            }
        }
    }
    pcout << "inner-wall nodes: " << innerWall.size() << endl;

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), AnalyticCouette());
    lattice.initialize();

    std::filesystem::create_directories(outDir);
    std::ofstream torqueFile(outDir + "/torque.csv");
    torqueFile << "time_s,torque_N_m_per_m\n";

    // Radial profile bins (1 lattice unit wide), accumulated in steady window.
    const plint steadyStart = (plint)(nSteps * 0.8);
    const plint nBins = (plint)std::ceil(RoLb) + 1;
    vector<T> profSum(nBins, (T)0), profCnt(nBins, (T)0);

    const plint sampleEvery = 100, fieldEvery = 1000;
    T torqueSteadySum = (T)0;
    plint torqueSteadyCnt = 0;
    plint fieldIdx = 0;

    auto writeField = [&](plint iT) {
        std::ofstream f(outDir + "/field_" + std::to_string(fieldIdx++) + ".csv");
        f << "iX,iY,ux_m_s,uy_m_s\n";
        for (plint iY = 0; iY < ny; iY += 2) {
            for (plint iX = 0; iX < nx; iX += 2) {
                Array<T, 2> u;
                lattice.get(iX, iY).computeVelocity(u);
                f << iX << "," << iY << "," << u[0] * dx / dt << "," << u[1] * dx / dt
                  << "\n";
            }
        }
    };

    for (plint iT = 0; iT < nSteps; ++iT) {
        // Quarter-sine ramp of Omega over the first rampSteps.
        rampFactorCurrent = (iT < rampSteps)
            ? std::sin(0.5 * PI * (T)iT / (T)rampSteps)
            : (T)1;
        if (iT <= rampSteps) {
            // Update the wall velocity on the dynamics object that lives
            // in the lattice cell (fetched fresh each step; see note above).
            T OmLb = OmegaLb * rampFactorCurrent;
            for (auto &wn : innerWall) {
                auto *dyn = dynamic_cast<VelocityBounceBack<T, DESCRIPTOR> *>(
                    &lattice.get(wn.iX, wn.iY).getDynamics());
                PLB_ASSERT(dyn);
                dyn->setUwall(wallVelocity(wn.xRel, wn.yRel, OmLb));
            }
        }

        lattice.collideAndStream();

        if (iT % sampleEvery == 0 || iT == nSteps - 1) {
            T Tq = computeTorque(lattice);
            torqueFile << iT * dt << "," << Tq << "\n";
            if (iT >= steadyStart) {
                torqueSteadySum += Tq;
                ++torqueSteadyCnt;
            }
        }
        if (iT >= steadyStart && iT % sampleEvery == 0) {
            for (plint iX = 0; iX < nx; ++iX) {
                for (plint iY = 0; iY < ny; ++iY) {
                    if (!isFluid(iX, iY)) continue;
                    T x = (T)iX - cx, y = (T)iY - cy;
                    T r = std::sqrt(x * x + y * y);
                    Array<T, 2> u;
                    lattice.get(iX, iY).computeVelocity(u);
                    T vThetaLb = (-y * u[0] + x * u[1]) / r;
                    plint b = (plint)r;
                    profSum[b] += vThetaLb * dx / dt; // m/s
                    profCnt[b] += (T)1;
                }
            }
        }
        if (iT % fieldEvery == 0) writeField(iT);
        if (iT % 4000 == 0) {
            pcout << "step " << iT << "/" << nSteps << " t=" << iT * dt << " s" << endl;
        }
    }
    writeField(nSteps - 1);

    // Profile from the steady window.
    std::ofstream profFile(outDir + "/profile.csv");
    profFile << "r_m,vtheta_m_s,vtheta_analytic_m_s\n";
    T tauFromProfile = (T)0;
    {
        // Shear rate at the inner wall from the two innermost fluid bins.
        plint b0 = (plint)RiLb + 1, b1 = b0 + 1;
        T r0 = ((T)b0 + 0.5) * dx, r1 = ((T)b1 + 0.5) * dx;
        T v0 = profSum[b0] / profCnt[b0], v1 = profSum[b1] / profCnt[b1];
        T omega0 = v0 / r0, omega1 = v1 / r1;
        T domegaDr = (omega1 - omega0) / (r1 - r0);
        // tau_rtheta = eta * r * d(omega)/dr is NEGATIVE (omega falls
        // with r); the torque magnitude uses -tau_rtheta.
        tauFromProfile = -etaPhys * r0 * domegaDr;
        pcout << "profile: vtheta(r0=" << r0 << ")=" << v0 << " m/s, r1=" << r1
              << " v1=" << v1 << ", tau_inner=" << tauFromProfile << " Pa" << endl;
    }
    for (plint b = 0; b < nBins; ++b) {
        if (profCnt[b] == 0) continue;
        T r = ((T)b + 0.5) * dx;
        profFile << r << "," << profSum[b] / profCnt[b] << "," << couetteVtheta(r)
                 << "\n";
    }

    T TqSteady = torqueSteadyCnt ? torqueSteadySum / torqueSteadyCnt : (T)0;
    T Tanalytic = 4 * PI * Lphys * etaPhys * Omega * Ri * Ri * Ro * Ro
                  / (Ro * Ro - Ri * Ri);
    pcout << "steady torque (momentum exchange): " << TqSteady * Lphys
          << " N*m (per-length " << TqSteady << " N*m/m)" << endl;
    pcout << "analytic torque: " << Tanalytic << " N*m, rel err "
          << std::abs(TqSteady * Lphys - Tanalytic) / Tanalytic * 100 << " %" << endl;
    pcout << "torque from profile: "
          << tauFromProfile * 2 * PI * Ri * Ri * Lphys << " N*m" << endl;

    return 0;
}
