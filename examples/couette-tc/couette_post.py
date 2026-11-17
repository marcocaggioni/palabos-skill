#!/usr/bin/env python3
"""Couette post-processing for the TC case: torque.csv -> flow curve + TC validation.

Input per sweep point: runs/gamma_<value>/torque.csv with columns
    time_s, torque_N_m_per_m
(torque per unit length; the 2D case. Multiply by L for 3D torque.)

Conversions (same as the Newtonian Couette example):
    gamma_i = 2*Omega*Ro^2 / (Ro^2 - Ri^2)      Omega from gamma via case.yaml
    tau_i   = T / (2*pi*Ri^2*L),  T = torque_per_m * L
    eta_app = tau_i / gamma_i
Validation is TORQUE vs the analytic TC torque at the same Omega —
NOT eta_app vs the TC law at the label rate: the inner-wall rate
exceeds the label rate, and at gamma <= 10 part of the gap is an
unyielded plug attached to the outer wall. The analytic TC torque
comes from the reference implementation in ~/workspace/tc-analytic.py,
copied here so this repo is self-contained (bisection on T, with
Omega(T) by quadrature; the integral machinery was verified against
the known Newtonian analytic torque, 0.0000% diff):
    tau(r) = T/(2*pi*L*r^2); invert the TC law for gamma_dot(tau)
    (quadratic in sqrt(gamma_dot); gamma_dot = 0 where tau < sigma_y);
    Omega = integral_{Ri}^{R_y} gamma_dot dr/r,
    R_y = min(Ro, sqrt(T/(2*pi*L*sigma_y))).

Validation gate: the simulated torque must match the analytic TC
torque within 10% (default) at every point. Exits non-zero on
failure — do not proceed to a flow curve until it passes.

Usage: python3 couette_post.py [--runs runs] [--validation-tol 0.10]
Writes flow_curve.csv and flow_curve.png in examples/couette-tc/.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from lattice import load_case  # noqa: E402

HERE = Path(__file__).parent

# --- TC analytic reference (copy of ~/workspace/tc-analytic.py) ---
B_COEF = None  # set in main from case.yaml (sigma_y / sqrt(gamma_dot_c))


def gamma_of_tau(tau, sigma_y, eta_bg):
    """Invert TC law: shear rate for a given shear stress (0 if unyielded)."""
    if tau <= sigma_y:
        return 0.0
    disc = B_COEF * B_COEF - 4.0 * eta_bg * (sigma_y - tau)
    x = (-B_COEF + math.sqrt(disc)) / (2.0 * eta_bg)
    return x * x


def omega_of_T(T, Ri, Ro, L, sigma_y, eta_bg, n=2000):
    """Angular velocity for a given torque, by quadrature over the gap."""
    if T <= 0:
        return 0.0
    r_yield = math.sqrt(T / (2.0 * math.pi * L * sigma_y))
    r_hi = min(Ro, r_yield)
    if r_hi <= Ri:
        return 0.0
    total = 0.0
    ln_lo, ln_hi = math.log(Ri), math.log(r_hi)
    dln = (ln_hi - ln_lo) / n
    for i in range(n):
        ln_r = ln_lo + (i + 0.5) * dln
        r = math.exp(ln_r)
        tau = T / (2.0 * math.pi * L * r * r)
        total += gamma_of_tau(tau, sigma_y, eta_bg) * dln
    return total


def torque_of_omega(omega, Ri, Ro, L, sigma_y, eta_bg, tol=1e-12):
    """Torque for a given angular velocity, by bisection on T."""
    T_lo = 2.0 * math.pi * L * sigma_y * Ri * Ri
    T_hi = T_lo * 2.0
    while omega_of_T(T_hi, Ri, Ro, L, sigma_y, eta_bg) < omega:
        T_hi *= 2.0
    for _ in range(200):
        T_mid = 0.5 * (T_lo + T_hi)
        if omega_of_T(T_mid, Ri, Ro, L, sigma_y, eta_bg) < omega:
            T_lo = T_mid
        else:
            T_hi = T_mid
        if (T_hi - T_lo) / T_mid < tol:
            break
    return 0.5 * (T_lo + T_hi)


def tau_tc_label(gamma, sigma_y, gamma_dot_c, eta_bg):
    """TC law evaluated at the label rate (for reference only — the
    inner-wall rate exceeds the label rate, so tau_i != this)."""
    return (sigma_y + sigma_y * math.sqrt(gamma / gamma_dot_c)
            + eta_bg * gamma)


def read_torque(path):
    times, torques = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            times.append(float(row["time_s"]))
            torques.append(float(row["torque_N_m_per_m"]))
    return times, torques


def steady_mean(times, torques, frac):
    n = len(times)
    k = max(1, int(n * frac))
    window = torques[-k:]
    mean = sum(window) / len(window)
    sd = (sum((x - mean) ** 2 for x in window) / len(window)) ** 0.5
    return mean, sd


def main():
    global B_COEF
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--validation-tol", type=float, default=0.10)
    args = ap.parse_args()

    case = load_case()
    g, f, s = case["geometry"], case["fluid"], case["sweep"]
    Ri, Ro, L = g["R_i"], g["R_o"], g["L"]
    sigma_y = f["sigma_y_Pa"]
    gamma_dot_c = f["gamma_dot_c_1_s"]
    eta_bg = f["eta_bg_Pa_s"]
    B_COEF = sigma_y / math.sqrt(gamma_dot_c)
    frac = s["steady_window_fraction"]
    rows = []
    for gamma in s["shear_rate_1_s"]:
        path = Path(args.runs) / f"gamma_{gamma}" / "torque.csv"
        if not path.exists():
            print(f"missing {path} — run the sweep first")
            return 1
        times, torques = read_torque(path)
        T_per_m, sd = steady_mean(times, torques, frac)
        T = T_per_m * L
        Omega = gamma * (Ro * Ro - Ri * Ri) / (2.0 * Ro * Ro)
        tau = T / (2.0 * math.pi * Ri * Ri * L)
        T_an = torque_of_omega(Omega, Ri, Ro, L, sigma_y, eta_bg)
        tau_an = T_an / (2.0 * math.pi * Ri * Ri * L)
        err = abs(T - T_an) / T_an if T_an else float("nan")
        r_y = min(Ro, math.sqrt(T_an / (2.0 * math.pi * L * sigma_y)))
        plug = r_y < Ro - 1e-12
        rows.append({"gamma": gamma, "omega": Omega, "tau": tau,
                     "tau_analytic": tau_an, "eta_app": tau / gamma,
                     "eta_app_analytic": tau_an / gamma,
                     "T": T, "T_analytic": T_an, "rel_err": err,
                     "torque_sd": sd * L, "R_y_mm": r_y * 1000,
                     "plug": plug})
        print(f"gamma={gamma:>6} 1/s  tau={tau:9.4f} Pa  (analytic {tau_an:9.4f})"
              f"  T={T:.6e} N.m  analytic err={err*100:.2f}%"
              f"  R_y={r_y*1000:.2f} mm{' plug' if plug else ' fully yielded'}")

    max_err = max(r["rel_err"] for r in rows)
    passed = max_err <= args.validation_tol
    print(f"TC validation: max analytic-torque error {max_err*100:.2f}% "
          f"(tol {args.validation_tol*100:.0f}%) — {'PASSED' if passed else 'FAILED'}")

    with open(HERE / "flow_curve.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gamma_1_s", "omega_rad_s", "tau_Pa",
                                           "tau_analytic_Pa", "eta_app_Pa_s",
                                           "eta_app_analytic_Pa_s", "torque_N_m",
                                           "torque_analytic_N_m", "rel_err",
                                           "R_y_mm", "plug"])
        w.writeheader()
        for r in rows:
            w.writerow({"gamma_1_s": r["gamma"], "omega_rad_s": r["omega"],
                        "tau_Pa": r["tau"], "tau_analytic_Pa": r["tau_analytic"],
                        "eta_app_Pa_s": r["eta_app"],
                        "eta_app_analytic_Pa_s": r["eta_app_analytic"],
                        "torque_N_m": r["T"],
                        "torque_analytic_N_m": r["T_analytic"],
                        "rel_err": r["rel_err"], "R_y_mm": r["R_y_mm"],
                        "plug": r["plug"]})

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].loglog([r["gamma"] for r in rows], [r["tau_analytic"] for r in rows],
                 "-", label="analytic TC")
    ax[0].loglog([r["gamma"] for r in rows], [r["tau"] for r in rows],
                 "o", label="Palabos")
    ax[0].set_xlabel("shear rate (1/s)"); ax[0].set_ylabel("shear stress (Pa)")
    ax[0].set_title("Couette flow curve, TC fluid (Palabos)")
    ax[0].legend()
    ax[1].loglog([r["gamma"] for r in rows],
                 [r["eta_app_analytic"] for r in rows], "-", label="analytic TC")
    ax[1].loglog([r["gamma"] for r in rows], [r["eta_app"] for r in rows],
                 "o", label="Palabos")
    ax[1].set_xlabel("shear rate (1/s)"); ax[1].set_ylabel("apparent viscosity (Pa.s)")
    ax[1].set_title("Apparent viscosity")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(HERE / "flow_curve.png", dpi=150)
    print("wrote flow_curve.csv and flow_curve.png")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
