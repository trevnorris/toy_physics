"""Release the internal order amplitude at each signed shallow shape, with both reservoirs present."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_exchange_shear import DEGREE
from mouth_live_waist_fields import CoupledLoads
from mouth_live_waist_order import Order
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import COMPONENTS, PORT, ROOT, emit, occupied, read_case
from mouth_unified_transverse import cell_difference


def loads(shape, n, mesh, seed, target):
    material = CoupledLoads(shape, n, 2*n, mesh, PORT, seed, target)
    # The stiffness minus the outgoing boundary block is unchanged by the
    # matched change of uniform lead operator and traction. Only the pocket
    # difference remains; it enters both the pole and its order derivative.
    material.wave.kc = material.wave.kc-cell_difference(material.wave, shape, mesh, DEGREE)
    return material


def equilibrate(shape, n, mesh, order_seed, pole_seed, ext_energy, target):
    order = Order(shape, n)
    material = loads(shape, n, mesh, pole_seed, target)
    if pole_seed is None:
        # The deep-to-shallow record supplies the target on this very shape.
        weights = material.wave.project@order_seed.T.ravel()
        out, _, _, _, _ = material.wave.pole(weights)
        if abs(out["frequency"]-target)>1e-7:
            raise ArithmeticError("Failed to recover transported pole before allowing order to react")
    chi, state = order.condition(order_seed, material, 1., BIAS, tolerance=4e-6)
    if state["min_order"] < -1e-10 or state["max_order"] > 1.+1e-9:
        raise ArithmeticError(f"Smooth stationary order left its allowed interval: {state}")
    energy = dict(order_inside=state["order_energy"], order_outside=ext_energy,
                  carrier_two_sided=state["carrier_energy"], wave=state["wave_frequency"])
    return chi, material.wave.seed, state, energy


def main(radius, n):
    prep, chi, _, _, _, _, source_sha = read_case(radius, n)
    mesh = prep["field_mesh"]
    source = ROOT/f"_scratch/s11c/shallow_puncture_fields01/r{radius}_m{n}"
    guard = json.loads((source/"outcome.json").read_text())
    if not guard["limitsVerified"] or guard["exitCode"] or guard["childOutcome"]["exitCode"]:
        raise ValueError("Missing guarded transported-state result")
    data = (source/"stdout").read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    rows = list(map(json.loads, data.splitlines()))
    info = rows[0]
    if (info["snapshot_sha256"], info["radius"], info["mesh"]) != (source_sha, radius, n):
        raise ValueError("Wrong transported source")
    frozen, = [r for r in rows if r["kind"] == "center" and abs(r["depth"]-1)<1e-10]
    x = np.array(prep["coordinates"])
    x[0] -= np.log(Meridian(radius, x).length)
    shape = Meridian(radius, x)
    center, seed, state, baseline = equilibrate(shape, n, mesh, chi, None, frozen["energy"]["order_outside"],
                                                frozen["wave"]["frequency"])
    freq, vector, wave = occupied(shape, n, mesh, center, *seed)
    if abs(state["wave_frequency"]-freq.real)>1e-8:
        raise ValueError("Reconstruction disagreed with the relaxed pole")
    path = ROOT/f"_scratch/s11c/shallow_puncture_fields01/relaxed_r{radius}_m{n}/center.npz"
    np.savez_compressed(path, order=center, vector=vector, pole=freq, coordinates=x,
                        original_snapshot_sha256=source_sha, transported_readout_sha256=digest)
    emit(kind="center", original_snapshot_sha256=source_sha, transported_readout_sha256=digest,
         saved_snapshot_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), radius=radius, mesh=n,
         energy=baseline, decrease_from_transported=fsum(frozen["energy"].values())-fsum(baseline.values()),
         max_order_change=float(np.max(abs(center-chi))), state=state, wave=wave)
    for axis, index, factor in (("radial", 2, 1/radius), ("axial", 0, 1.)):
        states = {}
        for step in (-.002, -.001, .001, .002):
            xx = x.copy()
            xx[index] += step*factor
            trial = Meridian(radius, xx)
            ext, = [r["energy"]["order_outside"] for r in rows if r["kind"]=="offset"
                    and (r["axis"], r["step"]) == (axis, step)]
            _, _, diagnostics, energies = equilibrate(trial, n, mesh, center, seed, ext, freq.real)
            coordinate = step if axis == "radial" else trial.length
            states[step] = (coordinate, energies)
            emit(kind="offset", axis=axis, step=step, state=diagnostics, coordinate=coordinate, energy=energies)
        for low, high in ((-.002, -.001), (.001, .002), (-.002, .002), (-.001, .001)):
            lo, a = states[low]
            hi, b = states[high]
            force = {k:(a[k]-b[k])/(hi-lo) for k in COMPONENTS}
            emit(kind="secant", axis=axis, interval=(low, high), components=force,
                 force_before_static_crystal=fsum(force.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=int, choices=(1, 2), required=True)
    parser.add_argument("--mesh", type=int, choices=(48, 72), required=True)
    args = parser.parse_args()
    main(args.radius, args.mesh)
