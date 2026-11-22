"""
verify.py
=========
Independent self-checks of the implementation.  Run:  python verify.py
Every check prints PASS / FAIL.
"""
from __future__ import annotations

import numpy as np

import analysis
import kinematics as kin
import statics as sta
from parameters import default_parameters


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    return ok


def run(p=None):
    p = p or default_parameters()
    stroke = kin.max_stroke(p)
    u = np.linspace(0, 0.97 * stroke, 800)
    r = analysis.solve(u, p)
    ok = True

    # 1. initial configuration (eqs. 22, 33-35, 72)
    ok &= check("theta_mid(0) = 0", abs(r["theta_mid"][0]) < 1e-12)
    ok &= check("theta_top(0) = pi/2", abs(r["theta_top"][0] - np.pi / 2) < 1e-9)
    ok &= check("theta_bot(0) = theta_bot,0", abs(r["theta_bot"][0] - p.theta_bot0) < 1e-9)
    ok &= check("v(0) = 0", abs(r["v"][0]) < 1e-12)
    ok &= check("F(0) = 0 and U(0) = 0", abs(r["F"][0]) < 1e-12 and abs(r["U"][0]) < 1e-12)

    # 2. loop 1 closure (eqs. 7, 8)
    th_m, th_t = r["theta_mid"], r["theta_top"]
    res_x = p.L_mid * np.cos(th_m) + p.L_top * np.cos(th_t) - p.L_mid
    res_y = p.L_mid * np.sin(th_m) + p.L_top * np.sin(th_t) - (p.L_top + u)
    ok &= check("loop 1 closes", max(abs(res_x).max(), abs(res_y).max()) < 1e-12,
                f"(max residual {max(abs(res_x).max(), abs(res_y).max()):.1e} m)")

    # 3. loop 2 closure (eqs. 26, 27)
    th_b = r["theta_bot"]
    res_x2 = p.a + p.L_mid * (np.cos(th_m) - 1) + p.L_bot * np.cos(th_b)
    res_y2 = p.b + p.L_mid * np.sin(th_m) + p.L_bot * np.sin(th_b) - r["v"]
    ok &= check("loop 2 closes", max(abs(res_x2).max(), abs(res_y2).max()) < 1e-12,
                f"(max residual {max(abs(res_x2).max(), abs(res_y2).max()):.1e} m)")

    # 4. analytic kinematic coefficients vs. central finite differences (eqs. 39-41)
    for key, name in (("theta_mid", "dtheta_mid/du"), ("theta_top", "dtheta_top/du"),
                      ("theta_bot", "dtheta_bot/du"), ("v", "dv/du")):
        num = np.gradient(r[key], u)
        ana = r[{"theta_mid": "dtheta_mid_du", "theta_top": "dtheta_top_du",
                 "theta_bot": "dtheta_bot_du", "v": "dv_du"}[key]]
        err = np.max(abs(num[2:-2] - ana[2:-2])) / np.max(abs(ana[2:-2]))
        ok &= check(f"{name} analytic = numeric", err < 5e-3, f"(rel. err {err:.1e})")

    # 5. virtual-work force = derivative of the spring energy (eqs. 65 vs 74)
    F_num = np.gradient(r["U"], u)
    err = np.max(abs(F_num[2:-2] - r["F"][2:-2])) / np.max(abs(r["F"][2:-2]))
    ok &= check("F(u) = dU/du", err < 5e-3, f"(rel. err {err:.1e})")

    # 6. U(u) = integral of F (eq. 73)
    Ffree = r["F_top"] + r["F_mid"] + r["F_bot"]
    U_int = sta.cumulative_trapezoid(Ffree, u)
    err = abs(U_int[-1] - r["U"][-1]) / r["U"][-1]
    ok &= check("U = int F du", err < 1e-3, f"(rel. err {err:.1e})")

    # 7. foot force (eq. 75): F = dU/du - Qv dv/du
    Qv = 3.0
    rq = analysis.solve(u, p, Q_v=Qv)
    ok &= check("foot force term", np.allclose(rq["F"], r["F"] - Qv * r["dv_du"]))

    # 8. feasibility
    ok &= check("sweep inside feasible range", bool(r["feasible"].all()))

    print("\nALL CHECKS PASSED" if ok else "\nSOME CHECKS FAILED")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
