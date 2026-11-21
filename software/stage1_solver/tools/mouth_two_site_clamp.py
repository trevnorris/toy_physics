"""Quasistatic two-pad force-gauge experiment at frozen local brane order.

The tagged rigid pads and elastic strands are trial substrate constituents.
Coordinates are in local pitches, areas in unit cells, and forces in model units.
The stage commands are feedback outputs. Only the equilibrated pad coordinates
and deflections of the gauges are available to the separate measurement reader.
"""
import hashlib
import json
from pathlib import Path

from mpmath import mp


mp.dps = 50
number = mp.mpf
ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "research/pde_ledger_v3/_measurements/mouth_offset_brane_readout.json"
SOURCE_HASH = "0e6500211a3cb7e255a3b0dcf02d8b758c6bcde638bef0b0f5593f0823298464"
REGISTRY = number("0.125")
GAP = number(1)
REST = number("0.95")
GAUGES = mp.diag([5, 5, 8, 8])


def packed(vector):
    return [mp.nstr(value, 46) for value in vector]


def preparation():
    data = SOURCE.read_bytes()
    if hashlib.sha256(data).hexdigest() != SOURCE_HASH:
        raise ValueError("Different archived half-order preparation")
    run = json.loads(data)["verified_runs"]["native_offset_brane48_01"]
    rows = [json.loads(line, parse_float=str) for line in run["literal_stdout"]]
    state = next(row for row in rows if row["test"] == "state")
    patch = next(row for row in rows if row["test"] == "contact_free_bulk")
    return state, patch


def material_stiffness(order):
    # Exact piecewise-linear minimizer between the grounded planes and pads.
    # The three gaps are 1, 1 and 1.5. mu0 and the patch area are both 1.
    stiffness = mp.zeros(4)
    for axis in (0, 1):
        stiffness[axis, axis] = 2 * order
        stiffness[axis+2, axis+2] = (1 + 1 / number("1.5")) * order
        stiffness[axis, axis+2] = stiffness[axis+2, axis] = -order
    return stiffness


def equilibrate(start, stage, bare_stiffness, reference, width, connected):
    pos = start.copy()
    for _ in range(30):
        force = bare_stiffness * (pos - reference) + GAUGES * (pos - stage)
        tangent = bare_stiffness + GAUGES
        lengths = []
        if connected:
            for span in (-width, width):
                v = [pos[2] - pos[0] + span, pos[3] - pos[1]]
                length = mp.sqrt(GAP**2 + sum(t*t for t in v))
                lengths.append(length)
                for i in (0, 1):
                    pull = (1 - REST/length) * v[i]
                    force[i] -= pull
                    force[i+2] += pull
                    for j in (0, 1):
                        entry = (i == j)*(1 - REST/length) + REST*v[i]*v[j]/length**3
                        tangent[i,j] += entry
                        tangent[i+2,j+2] += entry
                        tangent[i,j+2] -= entry
                        tangent[i+2,j] -= entry
        residual = max(abs(f) for f in force)
        if residual < number("1e-44"):
            return pos, residual, lengths
        pos -= mp.lu_solve(tangent, force)
    raise ArithmeticError("Pad equilibrium residual failed to converge")


def acquire(q, y, width, connected, stiffness, reference):
    target = mp.matrix([-q/2, -y/2, q/2, y/2])
    stage = target.copy()
    pos = target.copy()
    for iteration in range(1, 180):
        pos, residual, lengths = equilibrate(pos, stage, stiffness, reference, width, connected)
        error = target - pos
        if max(abs(e) for e in error) < number("1e-34"):
            pads = {}
            for name, first in (("wall", 0), ("ordered", 2)):
                pads[name] = {"xy": packed([pos[first], pos[first+1]]),
                              "stage_xy": packed([stage[first], stage[first+1]]),
                              "z": str(first//2), "gauge_k": str(int(GAUGES[first,first]))}
            return {"span": str(width), "fixture": "connected" if connected else "brane_only",
                    "q_target": str(q), "y_target": str(y), "pads": pads,
                    "strand_lengths": packed(lengths), "servo_iterations": iteration,
                    "equilibrium_residual": mp.nstr(residual, 12)}
        # Each calibrated stage sees only its own marker position.
        stage += error
    raise ArithmeticError("Measured pad positions failed to settle at target")


def main():
    state, patch = preparation()
    stiffness = material_stiffness(number(patch["local_order"]))
    ref_q = REGISTRY - number(patch["preload_strain"])
    reference = mp.matrix([-ref_q/2, 0, ref_q/2, 0])
    print(json.dumps({"source_sha256": SOURCE_HASH, "intake": state["intake"],
                      "position": patch["position"], "chi": patch["local_order"],
                      "normal_gap": str(GAP), "rest_length": str(REST), "strand_k": "1",
                      "zero_strain_registry": str(ref_q), "registry": str(REGISTRY),
                      "mode": "frozen order, fixed parallel pad orientations"}))
    full = number("1e-5")
    targets = [(REGISTRY, y) for y in (0, full, -full, full/2, -full/2)]
    targets += [(REGISTRY + dq, 0) for dq in (-full, full)]
    for width in (number("0.25"), number("1")):
        for q, y in targets:
            for connected in (False, True):
                print(json.dumps(acquire(q, number(y), width, connected, stiffness, reference)),
                      flush=True)


if __name__ == "__main__":
    main()
