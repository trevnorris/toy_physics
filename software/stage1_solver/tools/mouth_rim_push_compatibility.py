"""Check the finite-grid brane-only release condition on the saved resting drain."""

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from held_crossflow_budget import archived_packet
from mouth_feed_weak_flow import FeedPatch


def resting_trace(export_face):
    path = archived_packet("mouth_feed_weak_run.json", "n32_restoring_slope_u0.npz")
    patch = FeedPatch(32, 65, "restoring_slope")
    with np.load(path, allow_pickle=False) as packet:
        f, current = packet["field"], packet["current"]
        if (f.shape != patch.shape or float(packet["speed"]) != 0
                or float(packet["dx"]) != patch.dx or float(packet["dz"]) != patch.dz
                or not np.array_equal(packet["mouth"], patch.mouth)
                or not np.array_equal(current, patch.current)):
            raise ValueError("Mismatched zero-flow calibration")
    face = f[-1]
    # Raw checkpoint inputs only, in C order, for the separate Wolfram reader.
    with export_face.open("xb") as stream:
        np.asarray([face.real, face.imag, current], dtype="<f8").tofile(stream)
    n = abs(face)**2
    intake = current.sum() * patch.dx**3
    mouth = patch.mouth.copy()
    speed_shape = current/n
    loads = []
    for offset in (-1, 0, 1):
        # Translate the entire intake law by an exact grid cell, keeping the
        # lab carrier trace fixed. This controls sensitivity to X, not A.
        patch.mouth = np.roll(mouth, offset, axis=0)
        patch.restoring = np.roll(patch.weight, offset, axis=0)
        raw = np.roll(speed_shape, offset, axis=0) * n
        patch.current = intake * raw / (patch.dx**3 * raw.sum())
        reading = patch.read(f, 0)
        loads.append(float(reading["carried"] + reading["stress"]))
    print("inward_number_rate", intake)
    print("frozen_brane_loads", json.dumps(loads), flush=True)
    if (not np.isclose(intake, 0.25, rtol=0, atol=1e-13) or abs(loads[1]) > 1e-13
            or loads[0] == 0 or not np.isclose(loads[0], -loads[2], rtol=1e-10, atol=0)):
        raise ArithmeticError("Frozen-trace reflection or normalization check failed")


def release():
    t, tau, a, jerk, k, r0, r1, push = sp.symbols("t tau a jerk k r0 r1 F")
    # In the lab the regular finite-grid carrier has zdot = A z + b X.
    # The projected impulse kernel c.exp(A t).b begins r0 + r1 t.
    x = a*t**2/2 + jerk*t**3/6
    memory = sp.integrate((r0+r1*(t-tau))*x.subs(t, tau), (tau, 0, t))
    load = k*x + memory
    load_jet = sp.series(load, t, 0, 4)
    residual = (push + load).subs(t, 0)
    print("incremental carrier load through cubic order =", load_jet)
    print("brane-only release residual divided by signed nonzero push =", sp.cancel(residual/push))
    # The endpoint check tests what an endpoint displacement and traction
    # specify in the distributed dipole equation. Keep its coefficients live.
    w, length, amplitude, tension, gradient, inertia = sp.symbols("w L B T0 T1 mu0", nonzero=True)
    bump = amplitude*(w/length)**2*(1-w/length)**2
    endpoint = (bump.subs(w, 0), bump.diff(w).subs(w, 0),
                bump.subs(w, length), bump.diff(w).subs(w, length))
    acceleration = sp.simplify(sp.diff((tension+gradient*w)*bump.diff(w), w).subs(w, 0)/inertia)
    print("endpoint displacement/slope differences =", endpoint)
    print("endpoint acceleration difference =", acceleration, flush=True)
    if (residual != push or load.diff(t).subs(t, 0) != 0 or load.diff(t, 2).subs(t, 0) != a*k
            or any(endpoint) or sp.simplify(acceleration-2*amplitude*tension/(length**2*inertia)) != 0):
        raise ArithmeticError("Release expansion or endpoint variation failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export_face", type=Path)
    resting_trace(parser.parse_args().export_face)
    release()
