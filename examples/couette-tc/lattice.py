#!/usr/bin/env python3
"""Physical -> lattice parameters for the Palabos TC Couette case.

Reads examples/couette-tc/case.yaml. The TC apparent viscosity
    eta(gamma) = sigma_y/gamma + sigma_y/sqrt(gamma*gamma_dot_c) + eta_bg
spans ~1.5 Pa s (gamma=100) to ~722 Pa s (gamma=0.1), plus the cap
eta_max = eta(gamma_dot_reg) ~ 7.02e4 Pa s in the plug, so one dt
cannot serve the whole sweep: dt is solved PER POINT from the
viscosity at the label rate, at omega_target (default 1.0, so
nu_lb = 1/6), and the physical run time per point is n_steps*dt
(the step count, not the run time, is held fixed — 24,000 steps is
~40 viscous diffusion times at the label viscosity). The C++ case
recomputes the same per-point dt from gamma; the two must agree.

Formulas (D2Q9/D3Q19, cs^2 = 1/3):
    dx     = gap / N_gap
    nu_lb  = (1/omega - 0.5)/3
    dt     = nu_lb * dx^2 / (eta_label/rho)
    u_lb   = u_phys * dt / dx         ->  must stay < max_lattice_velocity
    torque conversion: lattice force -> N uses rho * dx^4 / dt^2
Local omega in the run: omega = 2/(1 + eta_local/eta_label), so
omega ~ 1 at the label rate, up to ~1.8 in the yielded zone near the
bob, down to omega(eta_max) in the capped plug (see per-point table).
"""
import json
import math
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # keep the file usable without PyYAML for a quick check
    yaml = None

CASE = Path(__file__).with_name("case.yaml")


def load_case():
    if yaml is None:
        raise SystemExit("PyYAML required: pip install pyyaml")
    return yaml.safe_load(CASE.read_text())


def eta_tc(gamma, f):
    """TC apparent viscosity at shear rate gamma (1/s), capped at eta_max."""
    if gamma <= 0:
        return f["eta_max_Pa_s"]
    eta = (f["sigma_y_Pa"] / gamma
           + f["sigma_y_Pa"] / math.sqrt(gamma * f["gamma_dot_c_1_s"])
           + f["eta_bg_Pa_s"])
    return min(eta, f["eta_max_Pa_s"])


def eta_max_check(f):
    """eta_max must equal eta(gamma_dot_reg); recompute, don't trust the file."""
    eta = (f["sigma_y_Pa"] / f["gamma_dot_reg_1_s"]
           + f["sigma_y_Pa"] / math.sqrt(f["gamma_dot_reg_1_s"]
                                         * f["gamma_dot_c_1_s"])
           + f["eta_bg_Pa_s"])
    return eta


def per_point_params(case):
    g = case["geometry"]
    f = case["fluid"]
    s = case["sweep"]
    lat = case["lattice"]
    Ri, Ro = g["R_i"], g["R_o"]
    gap = Ro - Ri
    N = lat["cells_across_gap"]
    dx = gap / N
    omega = lat["omega_target"]
    if not (0.5 < omega < 2.0):
        raise ValueError(f"omega_target={omega} outside (0.5, 2.0)")
    nu_lb = (1.0 / omega - 0.5) / 3.0
    n_steps = s["n_steps_per_point"]
    eta_max = eta_max_check(f)
    if abs(eta_max - f["eta_max_Pa_s"]) / eta_max > 1e-6:
        raise ValueError(f"case.yaml eta_max_Pa_s={f['eta_max_Pa_s']} != "
                         f"eta(gamma_dot_reg)={eta_max:.1f} — fix case.yaml")
    points = {}
    for gamma in s["shear_rate_1_s"]:
        eta_label = eta_tc(gamma, f)
        dt = nu_lb * dx * dx / (eta_label / f["rho_kg_m3"])
        Omega = gamma * (Ro * Ro - Ri * Ri) / (2.0 * Ro * Ro)
        u_wall = Omega * Ri
        u_lb = u_wall * dt / dx
        # omega in the capped plug: omega = 1/(3*nu_lb(eta_max)+0.5)
        nu_lb_max = (eta_max / f["rho_kg_m3"]) * dt / dx**2
        omega_plug = 1.0 / (3.0 * nu_lb_max + 0.5)
        points[str(gamma)] = {
            "gamma_1_s": gamma,
            "eta_label_Pa_s": eta_label,
            "dt_s": dt,
            "n_steps": n_steps,
            "run_time_s": n_steps * dt,
            "omega_target": omega,
            "omega_plug": omega_plug,
            "u_lb": u_lb,
            "mach_ok": u_lb < lat["max_lattice_velocity"],
        }
    return {
        "dx_m": dx,
        "cells_across_gap": N,
        "gap_m": gap,
        "eta_max_Pa_s": eta_max,
        "gamma_dot_reg_1_s": f["gamma_dot_reg_1_s"],
        "max_lattice_velocity": lat["max_lattice_velocity"],
        "force_conversion_N_per_lb_force_per_point": "rho*dx^4/dt^2 (dt varies)",
        "points": points,
    }


def main():
    case = load_case()
    p = per_point_params(case)
    out = CASE.with_name("lattice.json")
    out.write_text(json.dumps(p, indent=2) + "\n")
    print(f"dx = {p['dx_m']*1e3:.3f} mm   eta_max = {p['eta_max_Pa_s']:.1f} Pa.s "
          f"  steps/point = {list(p['points'].values())[0]['n_steps']:,}")
    print(f"{'gamma':>6} {'eta_label':>10} {'dt (s)':>11} {'run time':>10} "
          f"{'omega_plug':>10} {'u_lb':>8}")
    ok = True
    for pt in p["points"].values():
        print(f"{pt['gamma_1_s']:>6} {pt['eta_label_Pa_s']:>10.3f} "
              f"{pt['dt_s']:>11.3e} {pt['run_time_s']:>10.3e} "
              f"{pt['omega_plug']:>10.5f} {pt['u_lb']:>8.5f}")
        ok &= pt["mach_ok"]
    if not ok:
        print("FAIL: u_lb too high at some point — fix case.yaml first")
        return 1
    print("lattice check PASSED — wrote", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
