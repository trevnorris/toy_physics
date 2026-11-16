"""Signed radial and transverse held-mouth loads through one complete Cartesian cell."""

import argparse
import json
from pathlib import Path

import numpy as np

from cartesian_held_orientation import CartesianBath, measure


ARCHIVE = Path("docs/_measurements/cartesian_held_orientation.json")
GRIDS = ((4,16),(5,16),(4,24),(5,24))


def archived_rows():
    archive = json.loads(ARCHIVE.read_text())
    rows = []
    for label, shift in (("c4_base",0),("c4_half",.5)):
        run, = (r for r in archive["runs"] if r["name"]==label)
        if run["guard_exit"] or run["worker_exit"] or not run["limits_verified"]:
            raise ValueError("Unverified archived Cartesian control")
        row, = (r for r in map(json.loads,run["stdout"]) if r["d"]==8 and r["angle"]=="oblique")
        if row["cells"]!=4 or row["shift"]!=shift:
            raise ValueError("Wrong archived grid phase")
        rows.append(dict(radius=16,height=8,eighth=int(8*shift),**row))
    return rows


def compare(paths):
    rows = {}
    for r in archived_rows()+[json.loads(p.read_text()) for p in paths]:
        key = r["cells"], r["radius"], r["eighth"]
        if (key in rows or r["d"]!=8 or r["angle"]!="oblique" or
                r["height"]*2!=r["radius"] or r["shift"]*8!=r["eighth"]):
            raise ValueError(f"Mismatched or duplicate case {key}")
        rows[key] = r
    if set(rows)!={(c,R,i) for c,R in GRIDS for i in range(9)}:
        raise ValueError("Incomplete phase/domain/resolution comparison")
    for c,R in GRIDS:
        samples = [rows[c,R,i] for i in range(9)]
        u = np.array(samples[0]["direction"])
        axes = np.array([u,[-u[1],u[0]]])
        reflected_error = 0.
        for i in (1,2,3,5,6,7):
            a,b = samples[i],samples[8-i]
            # Inversion interchanges mouths and maps phase i/8 onto (8-i)/8.
            # Both geometries and both isolated baselines were solved independently.
            for src,tar in (("source_force",["parallel","perpendicular"]),
                            ("source_direct",["direct_parallel","direct_perpendicular"])):
                reflected_error = max(reflected_error, np.max(abs(axes@a[src]+np.array([b[t] for t in tar]))))
        if reflected_error>3e-12 or not np.isfinite(reflected_error):
            raise ArithmeticError(f"Complementary-phase inversion failed {reflected_error=}")
        for label,rad,side in (("mouth","parallel","perpendicular"),
                               ("direct","direct_parallel","direct_perpendicular")):
            force = np.array([[r[rad],r[side]] for r in samples])
            closure = max(abs(force[-1]-force[0]))
            if not np.isfinite(force).all() or closure>3e-12 or np.any(force[:,0]==0):
                raise ArithmeticError(f"Invalid full-cell readout {closure=}")
            # Use each of the eight unique phases once for the period mean.
            radial, transverse = force[:-1].T
            signed = 100*transverse/radial
            print(f"cells={c} R={R} H={R//2} readout={label} "
                  f"radial_min={radial.min():.12g} radial_mean={radial.mean():.12g} radial_max={radial.max():.12g} "
                  f"radial_span_percent={100*np.ptp(radial)/radial.mean():.9g} "
                  f"side_min={transverse.min():.12g} side_mean={transverse.mean():.12g} side_max={transverse.max():.12g} "
                  f"signed_min_percent={signed.min():.9g} signed_mean_percent={signed.mean():.9g} "
                  f"signed_max_percent={signed.max():.9g} end_minus_start_radial={force[-1,0]-force[0,0]:.9g} "
                  f"end_minus_start_side={force[-1,1]-force[0,1]:.9g} "
                  f"reflected_load_max_abs_diff={reflected_error:.9g}")
            print(f"cells={c} R={R} readout={label} phases_eighths=0,1,2,3,4,5,6,7,8 "
                  f"signed_percent_by_phase={','.join(f'{100*t/r:.9g}' for r,t in force)}")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action",required=True)
    run = sub.add_parser("measure")
    run.add_argument("--cells",type=int,choices=(4,5),required=True)
    run.add_argument("--radius",type=int,choices=(16,24),required=True)
    run.add_argument("--eighth",type=int,choices=range(9),required=True)
    comp = sub.add_parser("compare")
    comp.add_argument("readouts",nargs="+",type=Path)
    args = parser.parse_args()
    if args.action=="measure":
        shift = args.eighth/8
        # Evaluate 8/8 itself on the translated coordinates, without wrapping
        # phase or substituting the zero-phase result.
        slab = CartesianBath(args.cells,shift,radius=args.radius,height=args.radius//2)
        row = measure(slab,8,"oblique")
        print(json.dumps(dict(cells=args.cells,radius=args.radius,height=args.radius//2,
                              shift=shift,eighth=args.eighth,**row),allow_nan=False),flush=True)
    else:
        compare(args.readouts)
