"""
plotting.py
===========
All graphs as functions of the input displacement u.  Each `draw_*` function
fills a given matplotlib Axes, so the same code is used for the individual
figures and for the overview sheet.
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np

from analysis import PRBMResult

MM = 1e3
DEG = np.degrees


def _style(ax, title, ylabel):
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("input displacement u [mm]")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)


# ----------------------------------------------------------------------
# draw_* functions (fill an Axes)
# ----------------------------------------------------------------------
def draw_angles(ax, r: PRBMResult):
    u = r["u"] * MM
    ax.plot(u, DEG(r["theta_mid"]), label=r"$\theta_{mid}$")
    ax.plot(u, DEG(r["theta_top"]), label=r"$\theta_{top}$")
    ax.plot(u, DEG(r["theta_bot"]), label=r"$\theta_{bot}$")
    ax.legend()
    _style(ax, "PRBM link angles", "angle [deg]")


def draw_rotations(ax, r: PRBMResult):
    u = r["u"] * MM
    ax.plot(u, DEG(r["Phi_top"]), label=r"$\Phi_{top}$")
    ax.plot(u, DEG(r["Phi_mid"]), label=r"$\Phi_{mid}$")
    ax.plot(u, DEG(r["Phi_v_bot"]), label=r"$\Phi_{v,bot}$")
    ax.plot(u, DEG(r["Phi_h_bot"]), "--", label=r"$\Phi_{h,bot}$")
    ax.legend()
    _style(ax, "Spring rotations (from stress-free state)", "rotation [deg]")


def draw_foot(ax, r: PRBMResult):
    u = r["u"] * MM
    ax.plot(u, r["v"] * MM, label="foot displacement $v$")
    ax.plot(u, (r["X"] - r.p.a) * MM, label=r"$X-a$ (half-gap change)")
    ax.legend()
    _style(ax, "Foot displacement and lower-link geometry", "[mm]")


def draw_coefficients(ax, r: PRBMResult):
    u = r["u"] * MM
    l1, = ax.plot(u, r["dtheta_mid_du"] / MM, label=r"$d\theta_{mid}/du$")
    l2, = ax.plot(u, r["dtheta_top_du"] / MM, label=r"$d\theta_{top}/du$")
    l3, = ax.plot(u, r["dtheta_bot_du"] / MM, label=r"$d\theta_{bot}/du$")
    ax2 = ax.twinx()                       # dv/du is dimensionless -> own axis
    l4, = ax2.plot(u, r["dv_du"], "k--", label=r"$dv/du$ (right axis)")
    ax2.set_ylabel("dv/du [-]")
    ax.legend(handles=[l1, l2, l3, l4], fontsize=8)
    _style(ax, "First-order kinematic coefficients", "dtheta/du [rad/mm]")


def draw_torques(ax, r: PRBMResult):
    u = r["u"] * MM
    ax.plot(u, r["T_top"] * MM, label=r"$T_{top}$")
    ax.plot(u, r["T_mid"] * MM, label=r"$T_{mid}$")
    ax.plot(u, r["T_v_bot"] * MM, label=r"$T_{v,bot}$")
    ax.plot(u, r["T_h_bot"] * MM, "--", label=r"$T_{h,bot}$")
    ax.legend()
    _style(ax, "Spring torques (one spring each)", "torque [N mm]")


def draw_force(ax, r: PRBMResult):
    u = r["u"] * MM
    ax.plot(u, r["F"], "k", lw=2, label="F(u), one leg")
    ax.plot(u, r["F_two_legs"], "k--", lw=1, label="2F(u), two legs")
    ax.plot(u, r["F_top"], alpha=0.7, label="top springs")
    ax.plot(u, r["F_mid"], alpha=0.7, label="middle springs")
    ax.plot(u, r["F_bot"], alpha=0.7, label="lower springs")
    ax.legend(fontsize=8)
    _style(ax, "Input force-displacement (virtual work)", "force [N]")


def draw_energy(ax, r: PRBMResult):
    u = r["u"] * MM
    ax.stackplot(u, r["U_top"], r["U_mid"], r["U_v_bot"], r["U_h_bot"],
                 labels=["top", "middle", "v,bot", "h,bot"], alpha=0.35)
    ax.plot(u, r["U"], "k", lw=2, label=r"$U=\sum \frac{1}{2}K\Phi^2$ (74)")
    ax.plot(u, r["U_integral"], "r--", lw=1.2, label=r"$\int F\,du$ (73)")
    ax.legend(fontsize=8)
    _style(ax, "Stored spring energy", "energy [J]")


# ----------------------------------------------------------------------
# Mechanism poses
# ----------------------------------------------------------------------
def draw_pose(axes, r: PRBMResult, n_poses: int = 4):
    """Draw the PRBM links of the right-hand side (and mirrored foot loop)."""
    p = r.p
    ax1, ax2 = axes
    idx = np.unique(np.linspace(0, len(r["u"]) - 1, n_poses).astype(int))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(idx)))
    for c, i in zip(colors, idx):
        u, tm = r["u"][i], r["theta_mid"][i]
        lab = f"u = {u*MM:.1f} mm"
        # loop 1 : O -> P -> C   (O = central attachment of mid link)
        O = np.array([0.0, 0.0])
        P = p.L_mid * np.array([np.cos(tm), np.sin(tm)])
        C = np.array([p.L_mid, p.L_top + u])
        ax1.plot(*zip(O, P, C), "-o", color=c, label=lab)
        # loop 2 : blue joints -> foot (origin = initial foot position)
        blue_r = np.array([p.a + p.L_mid * (np.cos(tm) - 1), p.b + p.L_mid * np.sin(tm)])
        blue_l = np.array([-blue_r[0], blue_r[1]])
        foot = np.array([0.0, r["v"][i]])
        ax2.plot(*zip(blue_l, foot, blue_r), "-o", color=c, label=lab)
    for ax, t in ((ax1, "Upper loop (right side): O - P - C"),
                  (ax2, "Lower V-linkage: blue joints - foot")):
        ax.set_aspect("equal")
        ax.set_title(t, fontsize=10)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)


# ----------------------------------------------------------------------
# Figure builders
# ----------------------------------------------------------------------
_SINGLE = {
    "angles": draw_angles,
    "spring_rotations": draw_rotations,
    "foot_displacement": draw_foot,
    "kinematic_coefficients": draw_coefficients,
    "spring_torques": draw_torques,
    "force_displacement": draw_force,
    "energy": draw_energy,
}


def plot_overview(r: PRBMResult):
    fig, axs = plt.subplots(2, 4, figsize=(20, 9))
    for ax, fn in zip(axs.flat, _SINGLE.values()):
        fn(ax, r)
    axs.flat[-1].axis("off")
    fig.suptitle("Compliant jumping leg - PRBM results vs. input displacement u", fontsize=13)
    fig.tight_layout()
    return fig


def plot_poses(r: PRBMResult, n_poses: int = 4):
    fig, axes = plt.subplots(1, 2, figsize=(11, 6))
    draw_pose(axes, r, n_poses)
    fig.tight_layout()
    return fig


def plot_all(r: PRBMResult, outdir: str, show: bool = False, dpi: int = 150):
    """Save every figure as PNG in `outdir`; return the list of file paths."""
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for name, fn in _SINGLE.items():
        fig, ax = plt.subplots(figsize=(7, 4.5))
        fn(ax, r)
        fig.tight_layout()
        path = os.path.join(outdir, f"{name}.png")
        fig.savefig(path, dpi=dpi)
        paths.append(path)
        if not show:
            plt.close(fig)
    for name, fig in (("poses", plot_poses(r)), ("overview", plot_overview(r))):
        path = os.path.join(outdir, f"{name}.png")
        fig.savefig(path, dpi=dpi)
        paths.append(path)
        if not show:
            plt.close(fig)
    if show:
        plt.show()
    return paths
