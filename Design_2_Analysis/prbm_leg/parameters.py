"""
parameters.py
=============
All inputs of the compliant jumping-leg PRBM live here, together with the
derived quantities of Sections 2 and 8 (Step 1) of the report.

Units: SI (metres, pascals, newtons, radians) everywhere inside the code.
Only the plots / console output convert to mm.

EDIT `default_parameters()` (or pass a JSON file with --params) and put in
YOUR dimensions and material data.  The numbers below are placeholders.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, fields

# --- constants of the PRBM for a fixed-guided beam (report eq. 2 and 42) ---
GAMMA = 0.85      # characteristic radius factor
K_THETA = 2.65    # stiffness coefficient


def rect_I(width: float, thickness: float) -> float:
    """Second moment of area of a rectangular flexure, I = w t^3 / 12.

    width     : out-of-plane depth (printed layer height direction) [m]
    thickness : in-plane thickness, i.e. the bending direction      [m]
    """
    return width * thickness ** 3 / 12.0


@dataclass(frozen=True)
class PRBMParameters:
    # ---- geometry (from the mechanism drawing) [m] ----
    l_mid: float      # physical length of a middle flexure
    l_top: float      # physical length of an upper flexure
    d_m_mid: float    # width of the central (hatched) body
    d_o_mid: float    # width of an outer side body
    l_v_bot: float    # length of the short vertical flexure at the V
    l_h_bot: float    # length of the short horizontal flexure at the V tip
    n_mid: int        # number of middle flexures per side

    # ---- material / cross-sections ----
    E: float          # Young's modulus [Pa]
    I_top: float      # second moment of area, upper flexure  [m^4]
    I_mid: float      # second moment of area, middle flexure [m^4]
    I_v_bot: float    # second moment of area, short vertical flexure [m^4]
    I_h_bot: float    # second moment of area, short horizontal flexure [m^4]

    # ------------------------------------------------------------------
    # PRBM lengths (Section 2)
    # ------------------------------------------------------------------
    @property
    def L_mid(self) -> float:
        return GAMMA * self.l_mid                                   # eq. (2)

    @property
    def L_top(self) -> float:
        return GAMMA * self.l_top                                   # eq. (2)

    @property
    def a(self) -> float:
        """Horizontal distance foot joint -> right blue joint, eq. (3)."""
        return (2 * self.l_mid + self.d_m_mid + self.d_o_mid) / 2.0

    @property
    def b(self) -> float:
        """Vertical distance foot joint -> right blue joint, eq. (4)."""
        return (self.l_mid + self.d_o_mid
                + (self.d_m_mid + self.l_v_bot - self.l_h_bot) / 2.0)

    @property
    def L_bot(self) -> float:
        return math.hypot(self.a, self.b)                           # eq. (5)

    @property
    def theta_bot0(self) -> float:
        """Initial angle of the lower link, eq. (34) (third quadrant)."""
        return math.pi + math.atan2(self.b, self.a)

    # ------------------------------------------------------------------
    # Spring constants (Section 8, Step 1)
    # ------------------------------------------------------------------
    @property
    def K_top(self) -> float:
        return 2 * GAMMA * K_THETA * self.E * self.I_top / self.l_top      # (43)

    @property
    def K_mid(self) -> float:
        return 2 * GAMMA * K_THETA * self.E * self.I_mid / self.l_mid      # (44)

    @property
    def K_v_bot(self) -> float:
        return self.E * self.I_v_bot / self.l_v_bot                        # (45)

    @property
    def K_h_bot(self) -> float:
        return self.E * self.I_h_bot / self.l_h_bot                        # (46)

    @property
    def n_springs(self) -> int:
        """Total number of torsional springs: 4 + 4n + 2 + 1 = 4n + 7."""
        return 4 * self.n_mid + 7

    # ------------------------------------------------------------------
    def summary(self) -> str:
        mm = 1e3
        return "\n".join([
            "PRBM parameters (derived)",
            f"  L_mid   = {self.L_mid*mm:9.3f} mm",
            f"  L_top   = {self.L_top*mm:9.3f} mm",
            f"  a       = {self.a*mm:9.3f} mm",
            f"  b       = {self.b*mm:9.3f} mm",
            f"  L_bot   = {self.L_bot*mm:9.3f} mm",
            f"  theta_bot,0 = {math.degrees(self.theta_bot0):8.3f} deg",
            f"  K_top   = {self.K_top:9.4f} N m/rad",
            f"  K_mid   = {self.K_mid:9.4f} N m/rad",
            f"  K_v,bot = {self.K_v_bot:9.4f} N m/rad",
            f"  K_h,bot = {self.K_h_bot:9.4f} N m/rad",
            f"  n_mid   = {self.n_mid}  ->  {self.n_springs} torsional springs",
        ])

    @classmethod
    def from_json(cls, path: str) -> "PRBMParameters":
        """Load parameters from a JSON file with the field names above (SI units)."""
        with open(path) as fh:
            data = json.load(fh)
        names = {f.name for f in fields(cls)}
        unknown = set(data) - names
        if unknown:
            raise KeyError(f"Unknown parameter(s) in {path}: {sorted(unknown)}")
        return cls(**data)


def default_parameters() -> PRBMParameters:
    """PLACEHOLDER VALUES - replace with your own design."""
    mm = 1e-3
    depth = 20 * mm                     # out-of-plane depth of the print
    return PRBMParameters(
        l_mid=30 * mm,
        l_top=25 * mm,
        d_m_mid=10 * mm,
        d_o_mid=10 * mm,
        l_v_bot=5 * mm,
        l_h_bot=5 * mm,
        n_mid=5,
        E=1.5e9,                        # PETG, order of magnitude
        I_top=rect_I(depth, 0.6 * mm),
        I_mid=rect_I(depth, 0.6 * mm),
        I_v_bot=rect_I(depth, 0.6 * mm),
        I_h_bot=rect_I(depth, 0.6 * mm),
    )
