// Couette rheometer — Palabos three-component (TC) case (2D annulus).
//
// Same geometry and sweep as the Newtonian example in
// examples/couette/couette.cpp; the fluid is the TC law instead:
//   tau = sigma_y + sigma_y*sqrt(gamma_dot/gamma_dot_c) + eta_bg*gamma_dot
// with sigma_y = 70 Pa, gamma_dot_c = 100 1/s, eta_bg = 0.05 Pa s
// (Marco, 2026-11-17). Torque on the inner wall is computed by momentum
// exchange at the inner-wall links, converted to N*m per unit length.
// couette_post.py validates the torque against the analytic TC torque
// (bisection + quadrature, ported from ~/workspace/tc-analytic.py)
// and enforces the 10% gate.
//
// Method — TC rheology as a LOCAL relaxation frequency (the documented
// Palabos pattern for non-Newtonian fluids; see the skill:
// "Non-Newtonian means making the relaxation frequency a function of
// the local shear rate — Carreau-Yasuda is the documented pattern").
// This file implements it as TCBGKdynamics, modelled line-for-line on
// Palabos' own BGKCarreauDynamics
// (src/complexDynamics/carreauDynamics.hh): in collide(),
//   1. PiNeq (the non-equilibrium stress) from the populations via
//      momentTemplates::compute_rhoBar_j_PiNeq;
//   2. local shear rate from PiNeq by Chapman-Enskog,
//      gamma_dot_lb = omega * sqrt(2*|PiNeq|^2) / (2*rho*cs^2)
//      (tensorNormSqr is the Frobenius norm squared; the factor
//      matches Palabos' own Carreau alpha = lambda^2*gamma_dot^2);
//   3. local apparent viscosity from the TC law,
//      eta = min(sigma_y/gamma_dot + sigma_y/sqrt(gamma_dot*gamma_dot_c)
//                + eta_bg, eta_max), eta_max = eta(1e-3 1/s) ~ 7.02e4 Pa s
//      (the cap regularizes the 1/gamma_dot divergence at the yield
//      surface; the plug then creeps instead of being exactly rigid);
//   4. local omega = 1/(3*nu_lb + 0.5), nu_lb = eta/rho * dt/dx^2.
//   Steps 2–4 are implicit (omega sets gamma_dot and vice versa), so
//   they are iterated to a fixed point (10 iterations max, 1e-10
//   relative tolerance), exactly as BGKCarreauDynamics iterates its
//   implicit omega. Then a standard BGK collision with that omega.
//   The dynamics object stores no per-cell state: omega is recomputed
//   from the cell's own populations at every collide, so it is correct
//   whether Palabos shares or clones the dynamics object per cell.
//   omega_min/omega_max seen are tracked in file-scope statics and
//   printed at the end as a sanity check on the viscosity range.
//
// Lattice (per sweep point — a single dt cannot serve the whole sweep,
// because eta at the label rate spans 722 Pa s at gamma=0.1 down to
// 1.45 Pa s at gamma=100; see case.yaml and lattice.py):
//   dx = 0.5 mm, omega_target = 1.0 AT THE LABEL VISCOSITY,
//   dt = nu_lb*dx^2/(eta_label/rho), nu_lb = 1/6,
//   n_steps = 24,000 per point (fixed step count, NOT fixed physical
//   time), ramp = first 20% of steps (quarter-sine), steady window =
//   last 20%. Initial condition: rest (equilibrium, u = 0) — the
//   Newtonian analytic profile is not the TC profile, so the run
//   spins up from rest under the ramp instead.
//
// Everything else (staircase annulus, VelocityBounceBack moving inner
// wall with setUwall() ramp, outer BounceBack, momentum-exchange
// torque, unit conversion rho*dx^4/dt^2 for the 2D per-length torque)
// is identical to examples/couette/couette.cpp — see its header and
// RUN_REPORT.md for the API references and the two bugs found there
// (dangling dynamics pointer; dx^5 vs dx^4 torque conversion).
//
// Build (serial, no MPI):
//   mkdir -p build && cd build && cmake .. -DPALABOS_ROOT=$HOME/palabos && make -j4
// Run:
//   ./couette --gamma 10 --out runs/gamma_10
//
// Defaults below are examples/couette-tc/case.yaml (do not change the
// case by editing numbers here — change case.yaml and keep this file
// in sync).

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "core/dynamicsIdentifiers.h"
#include "core/util.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
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
// TC fluid (Marco, 2026-11-17):
static T sigmaY = 70.0;     // Pa
static T gammaDotC = 100.0; // 1/s
static T etaBg = 0.05;      // Pa s (50 cP)
static T gammaDotReg = 1e-3; // 1/s, regularization rate; eta_max = eta(gammaDotReg)
static T etaMax = 70221.4;  // Pa s, = sigmaY/1e-3 + sigmaY/sqrt(1e-3*100) + etaBg
static T gammaTarget = 10.0; // 1/s, from --gamma
static int N_gap = 10;  // cells across the gap
static T omegaTarget = 1.0; // relaxation frequency AT THE LABEL VISCOSITY
static plint nSteps = 24000; // per point (fixed step count, see header)
static T rampFraction = 0.2; // quarter-sine ramp over the first 20% of steps
static string outDir = "runs/gamma_10";

// ---- derived lattice parameters (same formulas as lattice.py) ----
static T dx, dt, nuLb, Omega, OmegaLb, torqueConv, etaLabel;
static plint nx, ny, rampSteps;
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

static const T PI = std::acos((T)-1);

// ---- TC constitutive law (physical units) ----
// Apparent viscosity at shear rate gamma (1/s), capped at etaMax.
static T etaTC(T gamma)
{
    if (gamma <= (T)0) return etaMax;
    T eta = sigmaY / gamma + sigmaY / std::sqrt(gamma * gammaDotC) + etaBg;
    return std::min(eta, etaMax);
}

// ---- Analytic TC Couette reference (port of ~/workspace/tc-analytic.py) ----
// Invert the TC law: shear rate for a given shear stress (0 if unyielded).
static T gammaOfTau(T tau)
{
    if (tau <= sigmaY) return (T)0;
    T B = sigmaY / std::sqrt(gammaDotC);
    T disc = B * B - (T)4 * etaBg * (sigmaY - tau);
    T x = (-B + std::sqrt(disc)) / ((T)2 * etaBg);
    return x * x;
}

// Angular velocity for a given torque, by quadrature over the gap.
static T omegaOfT(T Tq, int n = 2000)
{
    if (Tq <= (T)0) return (T)0;
    T rYield = std::sqrt(Tq / ((T)2 * PI * Lphys * sigmaY));
    T rHi = std::min(Ro, rYield);
    if (rHi <= Ri) return (T)0;
    T lnLo = std::log(Ri), lnHi = std::log(rHi);
    T dln = (lnHi - lnLo) / n, total = (T)0;
    for (int i = 0; i < n; ++i) {
        T r = std::exp(lnLo + (i + (T)0.5) * dln);
        T tau = Tq / ((T)2 * PI * Lphys * r * r);
        total += gammaOfTau(tau) * dln; // gamma(r)/r * dr = gamma * d(ln r)
    }
    return total;
}

// Torque for a given angular velocity, by bisection on T.
static T torqueOfOmega(T Om)
{
    T Tlo = (T)2 * PI * Lphys * sigmaY * Ri * Ri; // tau(Ri) = sigma_y
    T Thi = Tlo * (T)2;
    while (omegaOfT(Thi) < Om) Thi *= (T)2;
    for (int i = 0; i < 200; ++i) {
        T Tmid = (T)0.5 * (Tlo + Thi);
        if (omegaOfT(Tmid) < Om) Tlo = Tmid;
        else Thi = Tmid;
        if ((Thi - Tlo) / Tmid < 1e-12) break;
    }
    return (T)0.5 * (Tlo + Thi);
}

static T yieldRadius(T Tq)
{
    return std::min(Ro, std::sqrt(Tq / ((T)2 * PI * Lphys * sigmaY)));
}

// Analytic azimuthal velocity (m/s) at radius r (m) for torque Tq:
// v_theta = r * int_r^{R_y} gamma_dot dr'/r'; zero in the plug (r >= R_y).
static T vthetaAnalytic(T r, T Tq)
{
    T rY = yieldRadius(Tq);
    if (r >= rY) return (T)0;
    int n = 400;
    T lnLo = std::log(r), lnHi = std::log(rY);
    T dln = (lnHi - lnLo) / n, total = (T)0;
    for (int i = 0; i < n; ++i) {
        T rp = std::exp(lnLo + (i + (T)0.5) * dln);
        total += gammaOfTau(Tq / ((T)2 * PI * Lphys * rp * rp)) * dln;
    }
    return r * total;
}

// ---- TC dynamics: BGK with a local relaxation frequency from the TC law ----
// Modelled on Palabos' own BGKCarreauDynamics
// (src/complexDynamics/carreauDynamics.hh): same collide structure,
// same fixed-point iteration on the implicit omega, same equilibrium.
static T omegaMinSeen = 1e30, omegaMaxSeen = -1e30;

template <typename T_, template <typename U> class Descriptor>
class TCBGKdynamics : public IsoThermalBulkDynamics<T_, Descriptor> {
public:
    TCBGKdynamics(T_ omega_) : IsoThermalBulkDynamics<T_, Descriptor>(omega_)
    { }
    TCBGKdynamics(HierarchicUnserializer &unserializer) :
        IsoThermalBulkDynamics<T_, Descriptor>(T_())
    {
        this->unserialize(unserializer);
    }
    virtual TCBGKdynamics<T_, Descriptor> *clone() const
    {
        return new TCBGKdynamics<T_, Descriptor>(*this);
    }
    virtual int getId() const
    {
        return id;
    }
    virtual void collide(Cell<T_, Descriptor> &cell, BlockStatistics &statistics)
    {
        T_ rhoBar;
        Array<T_, Descriptor<T_>::d> j;
        Array<T_, SymmetricTensor<T_, Descriptor>::n> PiNeq;
        momentTemplates<T_, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
        T_ rho = Descriptor<T_>::fullRho(rhoBar);
        T_ piNeqNormSqr = SymmetricTensor<T_, Descriptor>::tensorNormSqr(PiNeq);
        // Chapman-Enskog: gamma_dot_lb = omega*sqrt(2*|PiNeq|^2)/(2*rho*cs^2).
        // Implicit in omega (gamma_dot sets eta sets omega): fixed point.
        T_ omega = this->getOmega();
        for (int i = 0; i < 10; ++i) {
            T_ gammaLb = omega * std::sqrt((T_)2 * piNeqNormSqr)
                         / ((T_)2 * rho * Descriptor<T_>::cs2);
            T_ gammaPhys = gammaLb / dt;
            T_ eta = etaTC(gammaPhys);
            T_ nuLb = eta / rhoPhys * dt / (dx * dx);
            T_ omegaNew = (T_)1 / ((T_)3 * nuLb + (T_)0.5);
            if (std::fabs((omegaNew - omega) / omega) < 1e-10) {
                omega = omegaNew;
                break;
            }
            omega = omegaNew;
        }
        this->setOmega(omega);
        if (omega < omegaMinSeen) omegaMinSeen = omega;
        if (omega > omegaMaxSeen) omegaMaxSeen = omega;
        T_ uSqr = dynamicsTemplates<T_, Descriptor>::bgk_ma2_collision(
            cell, rhoBar, j, this->getOmega());
        if (cell.takesStatistics()) {
            gatherStatistics(statistics, rhoBar, uSqr);
        }
    }
    virtual T_ computeEquilibrium(
        plint iPop, T_ rhoBar, Array<T_, Descriptor<T_>::d> const &j, T_ jSqr,
        T_ = T_()) const
    {
        T_ invRho = Descriptor<T_>::invRho(rhoBar);
        return dynamicsTemplates<T_, Descriptor>::bgk_ma2_equilibrium(
            iPop, rhoBar, invRho, j, jSqr);
    }

private:
    static int id;
};

template <typename T_, template <typename U> class Descriptor>
int TCBGKdynamics<T_, Descriptor>::id =
    meta::registerGeneralDynamics<T_, Descriptor, TCBGKdynamics<T_, Descriptor> >(
        std::string("TCBGKdynamics"));

// Initial condition: rest (equilibrium, u = 0 everywhere).
struct RestIC {
    void operator()(plint, plint, T &rho, Array<T, 2> &u) const
    {
        rho = (T)1;
        u[0] = (T)0;
        u[1] = (T)0;
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

// Torque on the inner wall (N*m per unit length) from momentum exchange,
// sampled right after collideAndStream (populations just streamed in).
// Same code as examples/couette/couette.cpp (moving-wall generalisation
// of Palabos' MomentumExchangeBounceBack); the torque depends only on
// the populations and the wall velocity, not on the rheology.
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
    // comparable with the analytic TC torque and with couette_post.py.
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

    // eta_max from the regularization rate (same formula as lattice.py's
    // eta_max_check — do not hard-code it twice).
    etaMax = sigmaY / gammaDotReg + sigmaY / std::sqrt(gammaDotReg * gammaDotC) + etaBg;

    // Lattice parameters (identical to lattice.py).
    const T gap = Ro - Ri;
    dx = gap / N_gap;
    etaLabel = etaTC(gammaTarget); // TC viscosity at the label rate
    const T nuPhys = etaLabel / rhoPhys;
    nuLb = (1.0 / omegaTarget - 0.5) / 3.0;
    dt = nuLb * dx * dx / nuPhys;
    Omega = gammaTarget * (Ro * Ro - Ri * Ri) / (2.0 * Ro * Ro);
    OmegaLb = Omega * dt;
    rampSteps = (plint)std::ceil(rampFraction * (T)nSteps);
    torqueConv = rhoPhys * std::pow(dx, 4) / (dt * dt); // lattice torque -> N*m/m (2D: force-like)
    nx = (plint)std::ceil(2.0 * Ro / dx) + 2;
    ny = nx;
    cx = nx / 2.0;
    cy = ny / 2.0;
    RiLb = Ri / dx;
    RoLb = Ro / dx;

    T Tanalytic = torqueOfOmega(Omega);
    pcout << "Couette-TC (Palabos): gamma=" << gammaTarget << " 1/s, Omega=" << Omega
          << " rad/s, dx=" << dx << " m, dt=" << dt << " s, omega_target=" << omegaTarget
          << " (at eta_label=" << etaLabel << " Pa.s), steps=" << nSteps
          << ", run_time=" << nSteps * dt << " s, grid=" << nx << "x" << ny << endl;
    pcout << "analytic TC torque: " << Tanalytic << " N*m, yield radius "
          << yieldRadius(Tanalytic) * 1000 << " mm"
          << (yieldRadius(Tanalytic) >= Ro ? " (fully yielded)" : " (plug at outer wall)")
          << endl;

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
        nx, ny, new TCBGKdynamics<T, DESCRIPTOR>(omegaTarget));
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

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), RestIC());
    lattice.initialize();

    std::filesystem::create_directories(outDir);
    std::ofstream torqueFile(outDir + "/torque.csv");
    torqueFile << "time_s,torque_N_m_per_m\n";
    {
        // Per-run provenance (lattice params actually used by this binary).
        std::ofstream pj(outDir + "/params.json");
        pj << "{\n  \"gamma_1_s\": " << gammaTarget
           << ",\n  \"dx_m\": " << dx
           << ",\n  \"dt_s\": " << dt
           << ",\n  \"eta_label_Pa_s\": " << etaLabel
           << ",\n  \"omega_target\": " << omegaTarget
           << ",\n  \"n_steps\": " << nSteps
           << ",\n  \"run_time_s\": " << nSteps * dt
           << ",\n  \"analytic_torque_N_m\": " << Tanalytic
           << ",\n  \"yield_radius_m\": " << yieldRadius(Tanalytic)
           << "\n}\n";
    }

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
        // tau_rtheta = eta(gamma_inner) * r * d(omega)/dr is NEGATIVE
        // (omega falls with r); gamma_inner = -r*d(omega)/dr and the
        // torque magnitude uses eta at that local shear rate.
        T gammaInner = -r0 * domegaDr;
        tauFromProfile = etaTC(gammaInner) * gammaInner;
        pcout << "profile: vtheta(r0=" << r0 << ")=" << v0 << " m/s, r1=" << r1
              << " v1=" << v1 << ", gamma_inner=" << gammaInner << " 1/s, tau_inner="
              << tauFromProfile << " Pa" << endl;
    }
    for (plint b = 0; b < nBins; ++b) {
        if (profCnt[b] == 0) continue;
        T r = ((T)b + 0.5) * dx;
        profFile << r << "," << profSum[b] / profCnt[b] << "," << vthetaAnalytic(r, Tanalytic)
                 << "\n";
    }

    T TqSteady = torqueSteadyCnt ? torqueSteadySum / torqueSteadyCnt : (T)0;
    pcout << "steady torque (momentum exchange): " << TqSteady * Lphys
          << " N*m (per-length " << TqSteady << " N*m/m)" << endl;
    pcout << "analytic TC torque: " << Tanalytic << " N*m, rel err "
          << std::abs(TqSteady * Lphys - Tanalytic) / Tanalytic * 100 << " %" << endl;
    pcout << "torque from profile: "
          << tauFromProfile * 2 * PI * Ri * Ri * Lphys << " N*m" << endl;
    pcout << "omega seen in fluid cells: min=" << omegaMinSeen
          << " max=" << omegaMaxSeen << " (target 1.0 at label viscosity; "
          << "low end is the capped plug)" << endl;

    return 0;
}
