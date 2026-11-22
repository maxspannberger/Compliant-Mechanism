"""
statics.py
==========
Force-displacement analysis by virtual work (Section 8 of the report).

All forces / torques are for ONE COMPLETE SYMMETRIC LEG (left + right side).
For two identical legs loaded together multiply F by 2.
"""
from __future__ import annotations

import numpy as np

from parameters import PRBMParameters


# ----------------------------------------------------------------------
# Step 2: spring rotations (Lagrange coordinates), eqs. 47-50
# ----------------------------------------------------------------------
def spring_rotations(theta_mid, theta_top, theta_bot, p: PRBMParameters) -> dict:
    d_bot = theta_bot - p.theta_bot0
    return {
        "Phi_top":   theta_top - np.pi / 2,     # (47)
        "Phi_mid":   theta_mid,                 # (48)
        "Phi_v_bot": d_bot,                     # (49)
        "Phi_h_bot": 2.0 * d_bot,               # (50)
    }


# ----------------------------------------------------------------------
# Step 3: spring torques, eqs. 51-55
# ----------------------------------------------------------------------
def spring_torques(phi: dict, p: PRBMParameters) -> dict:
    return {
        "T_top":   p.K_top * phi["Phi_top"],
        "T_mid":   p.K_mid * phi["Phi_mid"],
        "T_v_bot": p.K_v_bot * phi["Phi_v_bot"],
        "T_h_bot": p.K_h_bot * phi["Phi_h_bot"],   # = 2 K_h (theta_bot - theta_bot0)
    }


# ----------------------------------------------------------------------
# Spring energy, eq. 74, split per spring group (counting both sides)
# ----------------------------------------------------------------------
def spring_energy_parts(phi: dict, p: PRBMParameters) -> dict:
    return {
        "U_top":   4 * 0.5 * p.K_top * phi["Phi_top"] ** 2,            # 4 springs
        "U_mid":   4 * p.n_mid * 0.5 * p.K_mid * phi["Phi_mid"] ** 2,  # 4n springs
        "U_v_bot": 2 * 0.5 * p.K_v_bot * phi["Phi_v_bot"] ** 2,        # 2 springs
        "U_h_bot": 1 * 0.5 * p.K_h_bot * phi["Phi_h_bot"] ** 2,        # 1 spring
    }


def spring_energy(phi: dict, p: PRBMParameters):
    """Total stored energy, eq. (74)."""
    return sum(spring_energy_parts(phi, p).values())


# ----------------------------------------------------------------------
# Steps 5-7: virtual work -> F(u), eqs. 65 / 70
# ----------------------------------------------------------------------
def force_contributions(phi: dict, dtheta_top, dtheta_mid, dtheta_bot,
                        p: PRBMParameters) -> dict:
    """The three terms of eq. (65) (lower springs combined as in eq. 63)."""
    return {
        "F_top": 4 * p.K_top * phi["Phi_top"] * dtheta_top,
        "F_mid": 4 * p.n_mid * p.K_mid * phi["Phi_mid"] * dtheta_mid,
        "F_bot": 2 * (p.K_v_bot + 2 * p.K_h_bot) * phi["Phi_v_bot"] * dtheta_bot,
    }


def input_force(phi: dict, dtheta_top, dtheta_mid, dtheta_bot,
                p: PRBMParameters, Q_v: float = 0.0, dv_du=0.0):
    """Vertical input force F(u) on the upper bar, eq. (65).

    Q_v (optional): upward force on the foot, e.g. ground contact (eq. 75):
        F du + Q_v dv = sum(T dPhi)  ->  F = dU/du - Q_v dv/du
    """
    parts = force_contributions(phi, dtheta_top, dtheta_mid, dtheta_bot, p)
    return sum(parts.values()) - Q_v * dv_du


def cumulative_trapezoid(y, x):
    """Running integral of y dx (starts at 0); used for U(u) = int F du, eq. 73."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out
