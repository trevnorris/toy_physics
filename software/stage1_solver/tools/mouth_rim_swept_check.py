"""Compare the moving-cell energy/transpose with independent face traces."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_last_rings_fluid import GRIDS, read_source
from mouth_rim_disk_work_review import Disk
from mouth_rim_swept_cells import SweptCells
from mouth_shared_rim_shift import SharedRimChart
from mouth_supercore_joined_front import file_sha


def run(a):
    oracle = json.loads(a.faces.read_text())
    if file_sha(a.zero) != oracle["field"] or file_sha(a.wall) != oracle["wall"]:
        raise ValueError("Different held wall or wound states")
    expected = next(r["cases"] for r in oracle["readings"] if r["gauss"]==2)
    disk = Disk(a.state,a.wall,a.fields)
    model = disk.model
    x,d,_,_ = read_source(a.state,a.wall)
    front = Front(a.state,GRIDS["native"],1)
    native = Trial(front,a.state,x,(d.shape[1]-1)//3,d.shape[0]-1,False).model
    chart = SharedRimChart(model,native,a.axes)
    with np.load(a.zero) as f:
        driven,rest = f["wound"],f["rest"]
    rng = np.random.default_rng(3404)
    direction = rng.uniform(-1.,1.,driven.shape[1])
    records = []
    for name,n in (("fixed_wall_difference",0),("rim_four",4),("rim_eight",8)):
        material = chart.physical if n else chart.still
        velocity = (lambda b,r,z: chart.carrier(b,r,z,n)) if n else chart.stationary
        layer = SweptCells(model,front.c,chart.tents,material,velocity)
        for step, key in ((-2.e-5,"negative"),(2.e-5,"positive")):
            gradient = np.zeros_like(driven)
            energy = layer.evaluate(driven,step,gradient)
            correction = -(energy-layer.work(rest,step))/step
            errors = []
            for field in range(3):
                reference = gradient[field]@direction/step
                for delta in (1.e-7,5.e-8):
                    u = driven.copy()
                    u[field] += delta*direction
                    a1 = layer.work(u,step)
                    u[field] = driven[field]-delta*direction
                    a0 = layer.work(u,step)
                    numerical = (a1-a0)/(2*delta*step)
                    errors.append(abs(numerical-reference)/max(1.,abs(reference)))
            mismatch = correction-expected[name][key]
            if abs(mismatch)>2e-10 or max(errors)>1e-5:
                raise ArithmeticError(f"Unresolved face energy or derivative {name}: {mismatch}, {errors}")
            records.append(dict(case=name,step=step,points=len(layer.weight),correction=correction,
                                independent_mismatch=mismatch,normalized_gradient_errors=errors))
    return dict(field_sha256=oracle["field"],grid=disk.level,records=records)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("state","wall","fields","axes","zero","faces"):
        parser.add_argument(key,type=Path)
    print(json.dumps(run(parser.parse_args()),allow_nan=False))
