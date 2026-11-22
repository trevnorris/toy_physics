"""Bounded-order check where the smooth shallow radius-one solution meets its order floor."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_live_waist_active_wave import ActiveWave
from mouth_live_waist_bounded import bounded
from mouth_live_waist_order import Order
from mouth_live_waist_support import release
from mouth_live_waist_wave_location import share
from mouth_sealed_shape_relax import DELTA, Meridian
from mouth_shallow_puncture_fields import COMPONENTS, ROOT, emit, read_case
from mouth_shallow_puncture_order import loads


def at_shape(shape, n, mesh, chi, seed, target):
    order = Order(shape, n)
    material = ActiveWave.install(loads(shape, n, mesh, seed, target))
    if seed is None:
        weights = material.wave.project @ chi.T.ravel()
        out, _, _, _, _ = material.wave.pole(weights)
        if abs(out["frequency"] - target) > 1e-7:
            raise ArithmeticError("Failed to recover occupied shallow pole")
    chi, state = bounded(order, chi, material, 1.)
    return chi, state, material, order


def main(n):
    radius = 1
    prep, chi, _, _, _, _, source_sha = read_case(radius, n)
    mesh = prep["field_mesh"]
    source = ROOT/f"_scratch/s11c/shallow_puncture_fields01/r1_m{n}"
    guard = json.loads((source/"outcome.json").read_text())
    if not guard["limitsVerified"] or guard["exitCode"] or guard["childOutcome"]["exitCode"]:
        raise ValueError("Missing guarded transported-state result")
    raw = (source/"stdout").read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    rows = list(map(json.loads, raw.splitlines()))
    held, = [r for r in rows if r["kind"] == "center" and abs(r["depth"]-1) < 1e-10]
    x = np.array(prep["coordinates"])
    x[0] -= np.log(Meridian(radius, x).length)
    shape = Meridian(radius, x)
    emit(kind="preparation", mesh=n, radius=radius, source=source_sha, snapshot=rows[0]["snapshot_sha256"])
    if rows[0]["snapshot_sha256"] != source_sha:
        raise ValueError("Mismatched order preparation")
    center, state, mat, order = at_shape(shape,n,mesh,chi,None,held["wave"]["frequency"])
    freq, vector = mat.wave.seed
    weights = mat.wave.project @ center.T.ravel()
    split = min(shape.length, shape.length*order.t_waist + 2*DELTA)
    downstream = share(mat.wave,shape,vector,weights,mesh,split,shape.length)
    upstream = share(mat.wave,shape,vector,weights,mesh,0.,split)
    stored = state["wave_pocket_fraction"]
    if abs((upstream+downstream)/stored-1)>1e-5:
        raise ArithmeticError("Wave kinetic partition mismatch")
    corrections = release(order, mat, center, state["total_inside_and_drives"]) if state["fraction_at_zero"] else []
    emit(kind="bounded_stationary_check",state=state,reopening=corrections,
         downstream_share=downstream/stored)
    if any(c.get("undefined") or c["energy_slope"] < -.001 for c in corrections):
        raise ArithmeticError("Connected zero-order region opens downhill or without a defined pole")
    path = ROOT/f"_scratch/s11c/shallow_puncture_fields01/bounded_r1_m{n}/center.npz"
    np.savez_compressed(path, order=center, vector=vector, pole=freq, coordinates=x,
                        original_snapshot_sha256=source_sha, transported_readout_sha256=digest)
    emit(kind="center", radius=radius, mesh=n, saved_snapshot_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
         original_snapshot_sha256=source_sha,transported_readout_sha256=digest,
         decrease_from_transported=fsum(held["energy"].values())-
         (state["total_inside_and_drives"]+held["energy"]["order_outside"]),state=state,
         downstream_share=downstream/stored,max_order_change=float(np.max(abs(center-chi))))
    for axis, index in (("radial",2),("axial",0)):
        values = {}
        for step in (-.002,-.001,.001,.002):
            xx=x.copy()
            xx[index]+=step
            trial=Meridian(radius,xx)
            _, shifted, _, _=at_shape(trial,n,mesh,center,(freq,vector),freq.real)
            ext,=[r["energy"]["order_outside"] for r in rows if r["kind"]=="offset"
                  and (r["axis"],r["step"])==(axis,step)]
            energy=dict(order_inside=shifted["order_energy"],order_outside=ext,
                        carrier_two_sided=shifted["carrier_energy"],wave=shifted["wave_frequency"])
            coord=step if axis=="radial" else trial.length
            values[step]=(coord,energy)
            emit(kind="offset",axis=axis,step=step,state=shifted,coordinate=coord,energy=energy)
        for low,high in ((-.002,-.001),(.001,.002),(-.002,.002),(-.001,.001)):
            lo,a=values[low]
            hi,b=values[high]
            force={k:(a[k]-b[k])/(hi-lo) for k in COMPONENTS}
            emit(kind="secant",axis=axis,interval=(low,high),components=force,
                 force_before_static_crystal=fsum(force.values()))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh",type=int,choices=(48,72),required=True)
    main(parser.parse_args().mesh)
