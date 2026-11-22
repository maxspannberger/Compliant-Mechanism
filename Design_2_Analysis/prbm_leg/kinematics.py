"""
kinematics.py
=============
Displacement analysis (Sections 3-7 of the report).

Every function is vectorised: `u` may be a float or a numpy array [m].
"""
from __future__ import annotations

import numpy as np

from parameters import PRBMParameters


def _acos(x):
    """arccos that tolerates tiny round-off outside [-1, 1]."""
    return np.arccos(np.clip(x, -1.0, 1.0))


# ----------------------------------------------------------------------
# Vector loop 1: upper and middle PRBM links (Section 4, eqs. 12-21)
# ----------------------------------------------------------------------
def aux_triangle(u, p: PRBMParameters):
    """Auxiliary triangle O-P-C: returns R(u) (12), beta(u) (14), alpha(u) (17)."""
    u = np.asarray(u, dtype=float)
    R = np.sqrt(p.L_mid ** 2 + (p.L_top + u) ** 2)                      # (12)
    beta = np.arctan2(p.L_top + u, p.L_mid)                             # (14)
    cos_alpha = (p.L_mid ** 2 + R ** 2 - p.L_top ** 2) / (2 * p.L_mid * R)
    alpha = _acos(cos_alpha)                                            # (17)
    return R, beta, alpha


def theta_mid_of_u(u, p: PRBMParameters):
    """Middle link angle, eq. (19): theta_mid = beta - alpha."""
    _, beta, alpha = aux_triangle(u, p)
    return beta - alpha


def theta_top_of_theta_mid(theta_mid, p: PRBMParameters):
    """Upper link angle, eq. (21)."""
    return _acos(p.L_mid / p.L_top * (1.0 - np.cos(theta_mid)))


# ----------------------------------------------------------------------
# Vector loop 2: lower V-shaped linkage (Section 5, eqs. 28-32)
# ----------------------------------------------------------------------
def horizontal_separation(theta_mid, p: PRBMParameters):
    """X(u) = a + L_mid (cos(theta_mid) - 1), eq. (28)."""
    return p.a + p.L_mid * (np.cos(theta_mid) - 1.0)


def theta_bot_of_theta_mid(theta_mid, p: PRBMParameters):
    """Lower link angle, eq. (30) (third-quadrant branch)."""
    X = horizontal_separation(theta_mid, p)
    return np.pi + _acos(X / p.L_bot)


def foot_displacement(theta_mid, p: PRBMParameters):
    """Vertical foot displacement v(u), eq. (32)."""
    X = horizontal_separation(theta_mid, p)
    root = np.sqrt(np.maximum(p.L_bot ** 2 - X ** 2, 0.0))
    return p.b + p.L_mid * np.sin(theta_mid) - root


# ----------------------------------------------------------------------
# First-order kinematic coefficients (Section 7, eqs. 39-41)
# ----------------------------------------------------------------------
def kinematic_coefficients(theta_mid, theta_top, p: PRBMParameters):
    """Return dtheta_mid/du, dtheta_top/du, dtheta_bot/du  [rad/m].

    Eq. (39) is used in the equivalent form
        dtheta_mid/du = sin(theta_top) / ( L_mid sin(theta_top - theta_mid) ),
    since cos(tm) - sin(tm) cot(tt) = sin(tt - tm) / sin(tt).
    It is identical to (39) but has no 0/0 problem at theta_top = pi/2 and
    makes the singular configuration (theta_top = theta_mid, i.e. the two
    upper links stretched out in line) easy to see.
    """
    X = horizontal_separation(theta_mid, p)
    dmid = np.sin(theta_top) / (p.L_mid * np.sin(theta_top - theta_mid))   # (39)
    dtop = -(p.L_mid * np.sin(theta_mid)
             / (p.L_top * np.sin(theta_top))) * dmid                       # (40)
    root = np.sqrt(np.maximum(p.L_bot ** 2 - X ** 2, 1e-300))
    dbot = (p.L_mid * np.sin(theta_mid) / root) * dmid                     # (41)
    return dmid, dtop, dbot


def foot_velocity_coefficient(theta_mid, dtheta_mid_du, p: PRBMParameters):
    """dv/du, obtained by differentiating eq. (32) (needed only for eq. 75)."""
    X = horizontal_separation(theta_mid, p)
    root = np.sqrt(np.maximum(p.L_bot ** 2 - X ** 2, 1e-300))
    return p.L_mid * dtheta_mid_du * (np.cos(theta_mid)
                                      - X * np.sin(theta_mid) / root)


# ----------------------------------------------------------------------
# Feasibility (eq. 36 and the lower-link branch condition 0 < X <= L_bot)
# ----------------------------------------------------------------------
def is_feasible(u, p: PRBMParameters):
    """Boolean mask: is u inside the geometrically admissible motion range?"""
    u = np.asarray(u, dtype=float)
    R = np.sqrt(p.L_mid ** 2 + (p.L_top + u) ** 2)
    ok_upper = (R >= abs(p.L_mid - p.L_top)) & (R <= p.L_mid + p.L_top)   # (36)
    # lower-link condition is only meaningful where the upper loop closes
    with np.errstate(invalid="ignore"):
        th = theta_mid_of_u(u, p)
        X = horizontal_separation(th, p)
    ok_lower = (X > 0.0) & (X <= p.L_bot)
    return ok_upper & ok_lower


def max_stroke(p: PRBMParameters, n_scan: int = 20001) -> float:
    """Largest admissible input displacement u_max [m].

    Upper loop closes up to R = L_mid + L_top, i.e.
        u_geo = sqrt(L_top^2 + 2 L_mid L_top) - L_top.
    The lower-link condition X > 0 is checked numerically and may reduce it.
    At u_geo the two upper links are collinear: dtheta/du -> infinity
    (a toggle / lock-out position), so in practice stay below it.
    """
    u_geo = np.sqrt(p.L_top ** 2 + 2 * p.L_mid * p.L_top) - p.L_top
    u = np.linspace(0.0, u_geo, n_scan)
    ok = is_feasible(u, p)
    if ok.all():
        return float(u_geo)
    return float(u[np.argmax(~ok) - 1])
