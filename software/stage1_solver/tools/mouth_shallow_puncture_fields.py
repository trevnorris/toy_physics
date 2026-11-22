"""Two-sided field work while shortening a broad occupied meridian to one slab thickness."""
import argparse
import hashlib
import json
from math import ceil, fsum, log
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from mouth_exchange_carrier import BIAS
from mouth_exchange_shear import DEGREE, solve
from mouth_live_waist_fields import Carrier
from mouth_live_waist_order import Order
from mouth_live_waist_wave import OccupiedWave
from mouth_live_waist_wave_location import share
from mouth_sealed_shape_relax import DELTA, Meridian
from mouth_unified_transverse import PORT, cell_difference
from mouth_waist_outer_order import exterior


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "research/pde_ledger_v3/_measurements"
COMPONENTS = ("order_inside", "order_outside", "carrier_two_sided", "wave")


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def read_case(radius, n):
    case = f"r{radius}_broad_{n}"
    record = json.loads((DATA/"mouth_live_waist_order.json").read_text())
    original = record["runs"][case]
    raw = original["literal_stdout"]
    if (not original["limits_verified"] or original["worker_exit"] or original["guard_exit"]
            or hashlib.sha256(raw.encode()).hexdigest() != original["stdout_sha256"]):
        raise ValueError("Invalid archived occupied state")
    rows = list(map(json.loads, raw.splitlines()))
    prep = rows[0]
    state, = [row for row in rows if row.get("kind") == "coupled" and row["action"] == 1.]
    path = Path("_scratch/s11c/live_waist_dynamics01/fields", case+".npz")
    digest = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    if (record["saved_field_sha256"][str(path)] != digest
            or (prep["radius"], prep["point"], prep["J"], prep["mesh"]) != (radius, "broad", 1., n)):
        raise ValueError("Snapshot/held geometry mismatch")
    symmetric = json.loads((DATA/"mouth_unified_transverse.json").read_text())["runs"]["main"]
    if not symmetric["limits_verified"] or symmetric["guard_exit"] or symmetric["child_exit"]:
        raise ValueError("Invalid archived symmetric pole")
    text = symmetric["literal_stdout"]
    if hashlib.sha256(text.encode()).hexdigest() != symmetric["transcript_sha256"]:
        raise ValueError("Archived symmetric-pole transcript hash mismatch")
    pole, = [item for item in map(json.loads, text.splitlines()) if item["case"] == case]
    if pole["sha256"] != digest:
        raise ValueError("Different order at the reference pole")
    archived = json.loads((DATA/"mouth_coupled_waist_balance.json").read_text())["readings"]
    base, = [r for r in archived if (r["radius"], r["order_mesh"], r["state"]) == (radius, n, "broad")]
    outer, = [v for v in base["energies"]["radial"] if v["displacement"] == 0.]
    ref = pole["endpoints"][-1]
    deep_values = dict(order_inside=state["order_energy"], order_outside=outer["exterior_order"],
                       carrier_two_sided=state["carrier_energy"], wave=ref["frequency"])
    with np.load(ROOT/path) as snapshot:
        return (prep, snapshot["order"], snapshot["vector"], ref, deep_values,
                base["exterior_mesh"], digest)


def occupied(shape, n, mesh, chi, prior, seed):
    op = OccupiedWave(shape, n, 2*n, mesh, PORT)
    weights = op.project @ chi.T.ravel()
    inside = op.matrix(weights[:, None, None]*op.mc)
    mass = op.cut_m + inside
    sym = op.kc-cell_difference(op, shape, mesh, DEGREE)
    stiff = op.cut_k+op.matrix(weights[:, None, None]*sym)
    out, vector, _ = solve(shape, mesh, PORT, 0., prior, incoming_vector=seed,
                           supplied_system=(op.uniform[0], stiff, mass, inside),
                           target_frequency=prior.real, shape_work=False)
    overlap = abs(np.vdot(seed, mass@vector))**2 / np.vdot(seed, mass@seed).real
    if overlap < .94 or out["pole_residual"] > 1e-8:
        raise ArithmeticError(f"Lost continuing occupied pole, overlap {overlap}, {out}")
    pole = complex(out["frequency"], -out["amplitude_loss"])
    waist = minimize_scalar(lambda t: shape.profile(t)[0], bounds=(0., 1.), method="bounded").x
    split = min(shape.length, waist*shape.length+2*DELTA)
    downstream = share(op, shape, vector, weights, mesh, split, shape.length)
    cavity = np.vdot(vector, inside@vector).real
    upstream = share(op, shape, vector, weights, mesh, 0., split)
    if abs((upstream+downstream)/cavity-1) > 1e-5:
        raise ArithmeticError("Wave kinetic partition mismatch")
    return pole, vector, dict(frequency=pole.real, loss=-pole.imag, residual=out["pole_residual"],
                              overlap=overlap, cavity_kinetic=cavity, downstream_share=downstream/cavity)


def fields(shape, n, mesh, exterior_mesh, chi, wave):
    inner, _ = Order(shape, n).energy_gradient(chi)
    outer = exterior(shape, exterior_mesh)
    kinetic, _, flow = Carrier(shape, n, 2*n, mesh).response(chi, BIAS)
    e = dict(order_inside=inner, order_outside=outer["energy"], carrier_two_sided=kinetic,
             wave=wave["frequency"])
    return e, dict(depth=shape.length, current=flow["current"], terminal=flow["terminal_current"],
                   receiver_order_residual=outer["residual"], wave=wave)


def run(radius, n):
    prep, chi, seed, ref, deep, exterior_mesh, digest = read_case(radius, n)
    mesh = prep["field_mesh"]
    x = np.asarray(prep["coordinates"])
    deep_shape = Meridian(radius, x)
    pole = complex(ref["frequency"], -ref["amplitude_loss"])
    pole, seed, start = occupied(deep_shape, n, mesh, chi, pole, seed)
    if abs(pole.real-ref["frequency"]) > 1e-8 or abs(start["downstream_share"]-ref["downstream_of_band"]) > 1e-5:
        raise ArithmeticError("Archived deep branch differs")
    emit(kind="preparation", radius=radius, mesh=n, reference_depth=deep_shape.length,
         reference_energy=deep, recovered_symmetric_pole=start, snapshot_sha256=digest)
    depth_change = -log(deep_shape.length)
    # Integrate along log depth on the same material mesh. An even count includes
    # the geometrical mean of initial and final depths as a separate energy reading.
    steps = 2*ceil(abs(depth_change)/.12)
    for j in range(1, steps+1):
        xx = x.copy()
        xx[0] += depth_change*j/steps
        shape = Meridian(radius, xx)
        pole, seed, branch = occupied(shape, n, mesh, chi, pole, seed)
        emit(kind="pole", depth=shape.length, **branch)
        if 2*j == steps or j == steps:
            e, metrics = fields(shape, n, mesh, exterior_mesh, chi, branch)
            emit(kind="center", **metrics, energy=e, energy_minus_deep={k:e[k]-deep[k] for k in COMPONENTS})
    if abs(shape.length-1.) > 1e-12:
        raise ValueError("Shallow comparison must use exactly the nominal slab thickness")
    baseline = xx.copy()
    for axis, index, factor in (("radial", 2, 1/radius), ("axial", 0, 1.)):
        states = {}
        for step in (-.002, -.001, .001, .002):
            xx = baseline.copy()
            xx[index] += step*factor
            trial = Meridian(radius, xx)
            _, _, wave = occupied(trial, n, mesh, chi, pole, seed)
            energies, diagnostics = fields(trial, n, mesh, exterior_mesh, chi, wave)
            physical = step if axis == "radial" else trial.length
            states[step] = (physical, energies)
            emit(kind="offset", axis=axis, step=step, **diagnostics, energy=energies)
        for low, high in ((-.002, -.001), (.001, .002), (-.002, .002), (-.001, .001)):
            lo, a = states[low]
            hi, b = states[high]
            force = {k:(a[k]-b[k])/(hi-lo) for k in COMPONENTS}
            emit(kind="secant", axis=axis, interval=(low, high), coordinate=(lo, hi),
                 components=force, force_before_static_crystal=fsum(force.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=int, choices=(1, 2), required=True)
    parser.add_argument("--mesh", type=int, choices=(48, 72), required=True)
    arguments = parser.parse_args()
    run(arguments.radius, arguments.mesh)
