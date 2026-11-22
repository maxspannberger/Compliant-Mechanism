"""
analysis.py
===========
Glue module: takes the input displacement u (scalar or array), runs the
kinematic and static analysis and returns *all* quantities in one object.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass

import numpy as np

import kinematics as kin
import statics as sta
from parameters import PRBMParameters


@dataclass
class PRBMResult:
    p: PRBMParameters
    data: dict            # name -> 1-D numpy array (all the same length as u)

    def __getitem__(self, key):
        return self.data[key]

    def keys(self):
        return self.data.keys()


# Units used in the CSV header / console output
UNITS = {
    "u": "m", "R": "m", "beta": "rad", "alpha": "rad",
    "theta_mid": "rad", "theta_top": "rad", "theta_bot": "rad",
    "X": "m", "v": "m",
    "dtheta_mid_du": "rad/m", "dtheta_top_du": "rad/m",
    "dtheta_bot_du": "rad/m", "dv_du": "-",
    "Phi_top": "rad", "Phi_mid": "rad", "Phi_v_bot": "rad", "Phi_h_bot": "rad",
    "T_top": "Nm", "T_mid": "Nm", "T_v_bot": "Nm", "T_h_bot": "Nm",
    "F_top": "N", "F_mid": "N", "F_bot": "N", "F": "N", "F_two_legs": "N",
    "U_top": "J", "U_mid": "J", "U_v_bot": "J", "U_h_bot": "J",
    "U": "J", "U_integral": "J",
    "feasible": "bool",
}


def solve(u, p: PRBMParameters, Q_v: float = 0.0) -> PRBMResult:
    """Full PRBM solution for input displacement(s) u [m].

    Q_v : optional upward foot force [N] (eq. 75). Default 0 = free foot.
    """
    u = np.atleast_1d(np.asarray(u, dtype=float))

    # --- displacement analysis -----------------------------------------
    R, beta, alpha = kin.aux_triangle(u, p)
    th_mid = beta - alpha                                   # (18)/(19)
    th_top = kin.theta_top_of_theta_mid(th_mid, p)          # (21)
    X = kin.horizontal_separation(th_mid, p)                # (28)
    th_bot = kin.theta_bot_of_theta_mid(th_mid, p)          # (30)
    v = kin.foot_displacement(th_mid, p)                    # (32)

    # --- first-order coefficients ----------------------------------------
    with np.errstate(divide="ignore", invalid="ignore"):
        d_mid, d_top, d_bot = kin.kinematic_coefficients(th_mid, th_top, p)
        dv_du = kin.foot_velocity_coefficient(th_mid, d_mid, p)

    # --- statics -------------------------------------------------------
    phi = sta.spring_rotations(th_mid, th_top, th_bot, p)
    T = sta.spring_torques(phi, p)
    Fparts = sta.force_contributions(phi, d_top, d_mid, d_bot, p)
    F = sta.input_force(phi, d_top, d_mid, d_bot, p, Q_v=Q_v, dv_du=dv_du)
    Uparts = sta.spring_energy_parts(phi, p)
    U = sta.spring_energy(phi, p)
    # U_integral (eq. 73) needs a fine grid starting at u = 0: filled in sweep()
    U_int = np.full(len(u), np.nan)

    data = {
        "u": u, "R": R, "beta": beta, "alpha": alpha,
        "theta_mid": th_mid, "theta_top": th_top, "theta_bot": th_bot,
        "X": X, "v": v,
        "dtheta_mid_du": d_mid, "dtheta_top_du": d_top,
        "dtheta_bot_du": d_bot, "dv_du": dv_du,
        **phi, **T, **Fparts,
        "F": F, "F_two_legs": 2.0 * F,
        **Uparts, "U": U, "U_integral": U_int,
        "feasible": kin.is_feasible(u, p).astype(float),
    }
    return PRBMResult(p=p, data=data)


def sweep(p: PRBMParameters, u_max: float | None = None, n_points: int = 400,
          Q_v: float = 0.0, safety: float = 0.95) -> PRBMResult:
    """Solve on a uniform grid 0 ... u_max.

    If u_max is None, the largest admissible stroke is used. The grid is
    capped at `safety * stroke` because dtheta/du -> infinity at the
    lock-out position (upper links collinear).
    """
    stroke = kin.max_stroke(p)
    if u_max is None:
        u_max = safety * stroke
    elif u_max > stroke:
        print(f"[warning] requested u_max = {u_max*1e3:.2f} mm exceeds the "
              f"admissible stroke {stroke*1e3:.2f} mm; clipping.")
        u_max = safety * stroke
    res = solve(np.linspace(0.0, u_max, n_points), p, Q_v=Q_v)
    # eq. (73): U = integral of the free-foot force (independent check of eq. 74)
    F_free = res["F_top"] + res["F_mid"] + res["F_bot"]
    res.data["U_integral"] = sta.cumulative_trapezoid(F_free, res["u"])
    return res


def save_csv(result: PRBMResult, path: str) -> None:
    """Write all quantities (columns) vs u (rows) to a CSV file."""
    keys = list(result.keys())
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([f"{k} [{UNITS.get(k, '')}]" for k in keys])
        for i in range(len(result["u"])):
            w.writerow([f"{result[k][i]:.10g}" for k in keys])


def print_table(result: PRBMResult) -> None:
    """Pretty console output (used for single-u calls)."""
    n = len(result["u"])
    for i in range(n):
        print(f"\n--- u = {result['u'][i]*1e3:.4f} mm "
              f"({'feasible' if result['feasible'][i] else 'NOT FEASIBLE'}) ---")
        for k in result.keys():
            if k in ("u", "feasible"):
                continue
            val = result[k][i]
            if np.isnan(val):
                continue
            unit = UNITS.get(k, "")
            if unit == "rad":
                print(f"  {k:<14s} = {val:14.6f} rad  ({np.degrees(val):10.4f} deg)")
            elif unit == "m":
                print(f"  {k:<14s} = {val*1e3:14.6f} mm")
            else:
                print(f"  {k:<14s} = {val:14.6g} {unit}")
