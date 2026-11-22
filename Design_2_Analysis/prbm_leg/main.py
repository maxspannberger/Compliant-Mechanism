"""
main.py
=======
Command-line entry point.

Examples
--------
    python main.py                       # full sweep 0 ... max stroke, saves CSV + PNGs in ./results
    python main.py --u-max 10            # sweep 0 ... 10 mm
    python main.py --u 8                 # single input u = 8 mm, prints every value
    python main.py --u 4 8 12            # several discrete inputs
    python main.py --params my_leg.json  # your own dimensions (SI units, see parameters.py)
    python main.py --Qv 2.0              # add an upward foot force of 2 N (eq. 75)
    python main.py --show                # also open the plot windows
"""
from __future__ import annotations

import argparse
import os

import matplotlib

import analysis
import kinematics as kin
from parameters import PRBMParameters, default_parameters


def parse_args():
    ap = argparse.ArgumentParser(description="PRBM analysis of the compliant jumping leg")
    ap.add_argument("--u", type=float, nargs="+", default=None,
                    help="one or more input displacements in mm (prints all values)")
    ap.add_argument("--u-max", type=float, default=None,
                    help="end of the sweep in mm (default: maximum admissible stroke)")
    ap.add_argument("--n-points", type=int, default=400, help="points in the sweep")
    ap.add_argument("--params", type=str, default=None,
                    help="JSON file with the parameters (SI units)")
    ap.add_argument("--Qv", type=float, default=0.0,
                    help="upward force on the foot in N (default 0 = free foot)")
    ap.add_argument("--outdir", type=str, default="results")
    ap.add_argument("--show", action="store_true", help="open plot windows")
    return ap.parse_args()


def main():
    args = parse_args()
    if not args.show:
        matplotlib.use("Agg")           # no GUI needed
    import plotting                      # after backend selection

    p = PRBMParameters.from_json(args.params) if args.params else default_parameters()
    print(p.summary())
    stroke = kin.max_stroke(p)
    print(f"\n  admissible stroke u_max = {stroke*1e3:.3f} mm "
          f"(lock-out: upper links collinear, force -> infinity)\n")

    os.makedirs(args.outdir, exist_ok=True)

    if args.u is not None:                                  # discrete inputs
        u = [x * 1e-3 for x in args.u]
        res = analysis.solve(u, p, Q_v=args.Qv)
        analysis.print_table(res)
        # sweep still used for the graphs, with the given points marked in the table only
        sweep_max = max(u)
    else:
        sweep_max = None if args.u_max is None else args.u_max * 1e-3

    res = analysis.sweep(p, u_max=sweep_max, n_points=args.n_points, Q_v=args.Qv)
    csv_path = os.path.join(args.outdir, "prbm_results.csv")
    analysis.save_csv(res, csv_path)
    paths = plotting.plot_all(res, args.outdir, show=args.show)

    print(f"\nSaved {csv_path}")
    for pth in paths:
        print(f"Saved {pth}")


if __name__ == "__main__":
    main()
