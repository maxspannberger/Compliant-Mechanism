# Compliant jumping leg - PRBM analysis

Input: vertical displacement `u` of the upper bar. Output: every angle, the foot
displacement, kinematic coefficients, spring rotations/torques, the force F(u)
and the spring energy, as numbers (CSV / console) and graphs (PNG).

| file | content | report sections |
|---|---|---|
| `parameters.py` | inputs (geometry, E, I, n), L_mid, L_top, a, b, L_bot, spring constants | 2, 8 (step 1) |
| `kinematics.py` | vector loops, theta_mid/top/bot, X, v, dtheta/du, feasible range | 3-7 |
| `statics.py`    | spring rotations, torques, virtual-work force, energy | 8 |
| `analysis.py`   | `solve(u, p)` / `sweep(p)`: runs everything, CSV export, console table | - |
| `plotting.py`   | all graphs (individual PNGs, overview sheet, mechanism poses) | - |
| `main.py`       | command line | - |
| `verify.py`     | self-checks (loop closure, finite differences, F = dU/du, ...) | - |

## Use
    pip install -r requirements.txt
    python main.py                       # sweep 0 ... 95 % of max stroke -> ./results
    python main.py --u 5 10              # print every value at u = 5 mm and 10 mm
    python main.py --u-max 12 --n-points 600
    python main.py --params my_leg_example.json
    python main.py --Qv 2.0              # foot force, eq. (75)
    python verify.py

Or from Python:

    from parameters import default_parameters
    import analysis
    p = default_parameters()
    r = analysis.solve([0.005, 0.010], p)      # u in metres
    print(r["F"], r["theta_mid"], r["v"])

## Notes
* SI units inside the code (m, N, Pa, rad); plots are in mm / deg.
* The numbers in `default_parameters()` are PLACEHOLDERS - put in your own dimensions.
* F is for one complete leg (left + right); two legs need 2F (`F_two_legs`).
* At u_max the two upper links become collinear (theta_top = theta_mid):
  dtheta/du and F go to infinity. Stay below it.
