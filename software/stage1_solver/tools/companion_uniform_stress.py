"""Compare a held pair with a single held target in a matched uniform carrier stream."""

import argparse
import json
from pathlib import Path
import shlex

import numpy as np

from incident_flow_falloff import checked_grid, check_current
from lateral_reservoir_field import LateralSlab
from lateral_source_rate import verified_packet
from lateral_target_rate import force, profile
from mouth_feed_weak_check import dx4


def center_speed(slab, velocity):
    # Axis and center are between samples. Remove both quadratic errors by
    # even extrapolation in radius and symmetric interpolation in x.
    i = len(slab.x) // 2
    axis = (9 * velocity[:, 0] - velocity[:, 1]) / 8
    near = (axis[i-1] + axis[i]) / 2
    far = (axis[i-2] + axis[i+1]) / 2
    return (9 * near - far) / 8


def check_solution(s, f, normal, phase_step):
    # Direct discrete stationary equation. Across an axial link the physical
    # field acquires the uniform-stream phase, including at the two endcaps.
    end = np.concatenate((2-f[:, :1], f, 2-f[:, -1:]), axis=1)
    link = np.exp(1j*phase_step)
    lap = (link*end[:, 2:] - 2*f + link.conjugate()*end[:, :-2])/s.dx**2
    radii = np.arange(len(s.r)+1)*s.dx
    side = np.concatenate((f[:, :, :1], f, 2-f[:, :, -1:]), axis=2)
    lap += np.diff(radii*np.diff(side, axis=2), axis=2)/(s.r*s.dx**2)
    faces = np.concatenate((f[1:2], f, (f[-2]+2*s.dz*normal)[None]), axis=0)
    lap += np.diff(faces, n=2, axis=0)/s.dz**2
    bath_energy = 2*np.sin(phase_step/2)**2/s.dx**2
    residual = np.max(abs(-lap/2 + ((abs(f)**8-1)/4 - bath_energy)*f))
    left = (link * end[:, 0].conj() * f[:, 0]).imag
    right = (link.conjugate() * end[:, -1].conj() * f[:, -1]).imag
    ends = np.sum((left+right)*(2*np.pi*s.r*s.dx), axis=1)
    outer = 2*np.pi*radii[-1]*s.dx*2*f[:, :, -1].imag.sum(axis=1)
    inlet = s.dz*np.sum(s.sqrt_z**2*(ends+outer))/s.dx
    return residual, inlet


def uniform_solution(s, base, speed, weight, shape):
    # f is the envelope after removing the exact unit-density plane wave.
    # This drive changes the reservoir phase and leaves its density fixed.
    step_phase = -speed*s.dx
    plane = np.exp(-1j*speed*s.x)[:, None]
    def boundary(face):
        n = abs(face)**2
        raw = shape*n
        j = 0.25*raw/np.sum(s.dv*raw)
        return j, (-weight*(1-1/abs(face))+1j*j/n)*face

    f = base.copy()
    for step in range(1, 501):
        j, normal = boundary(f[-1])
        end = np.concatenate((2-f[:, :1], f, 2-f[:, -1:]), axis=1)
        correction = (-2*np.sin(step_phase/2)**2 * np.diff(end, n=2, axis=1)
                      + 1j*np.sin(step_phase)*(end[:, 2:]-end[:, :-2]))/(2*s.dx**2)
        rhs = correction - ((abs(f)**8-1)/4)*f
        re = s.elliptic(2*(f.real-1) + rhs.real, normal.real, 2)
        im = s.elliptic(rhs.imag.copy(), normal.imag, 0)
        updated = 1+re+1j*im
        error = np.max(abs(updated-f))
        f = updated
        if not np.isfinite(error):
            raise ValueError("Invalid uniform-stream envelope")
        if error < 2e-13:
            j, normal = boundary(f[-1])
            residual, inlet = check_solution(s, f, normal, step_phase)
            check_current(s, 0, j[None])
            print(f"uniform iterations={step} residual={residual:.8g} bath_inward={inlet:.14g}", flush=True)
            if not (residual < 5e-9 and abs(inlet-0.25) < 2e-8):
                raise ValueError("Uniform-stream balance failed")
            layers = plane*f[-4:]
            # Remove the target analytically: the empty bath is a plane wave.
            speed_read = center_speed(s, np.broadcast_to(
                -(plane.conj()*dx4(plane, s.dx)).imag, f[-1].shape))
            if abs(speed_read-speed) > 2e-14:
                raise ValueError("Empty-target incident calibration failed")
            return force(s, layers, j, 0.25)
    raise ArithmeticError("Uniform-stream field failed to converge")


def compare(distance):
    manifest = json.loads(Path("docs/_measurements/held_force_incident_flow.json").read_text())
    row, = [r for r in manifest["runs"] if r["name"] == f"L48C4_d{distance}"]
    args = shlex.split(row["command"])
    pins = {Path(k): v for k, v in (row["input_sha256"] | row["output_fields_sha256"]).items()}
    def packet(flag):
        return verified_packet(Path(args[args.index(flag)+1]), pins)

    s = LateralSlab(48, 4, 16)
    with packet("--baseline") as p:
        checked_grid(s, p)
        base = p["field"]
        current = check_current(s, 0, p["currents"])
    alone = force(s, base[-4:], current, 0.25)
    if np.max(abs(alone)) > 1e-11:
        raise ValueError("Isolated target fails reflection symmetry")
    weight, shape = profile(s, 0, abs(base[-1])**2)
    with packet("--source") as p:
        checked_grid(s, p)
        check_current(s, -distance, p["currents"])
        face = p["layers"][-1]
        velocity = -(face.conj()*dx4(face, s.dx)).imag/abs(face)**2
    center = center_speed(s, velocity)
    mean = np.sum(s.dv*current*velocity)/0.25
    with packet("--pair" if "--pair" in args else "--output") as p:
        checked_grid(s, p)
        check_current(s, -distance, p["currents"][:1])
        target = check_current(s, 0, p["currents"][1:])
        pair = force(s, p["layers"], target, 0.25) - alone
    archived = json.loads(row["stdout"][-1])
    np.testing.assert_allclose([mean, *pair], [archived[k] for k in
                               ("incident_speed", "carried", "stress", "direct")], rtol=2e-9, atol=1e-16)
    for match, speed in (("center", center), ("intake_mean", mean)):
        uniform = uniform_solution(s, base, speed, weight, shape) - alone
        diff = pair - uniform
        print(json.dumps(dict(d=distance, match=match, speed=float(speed),
                              pair=pair.tolist(), uniform=uniform.tolist(), difference=diff.tolist(),
                              combined_difference=float(sum(diff[:2])),
                              scale=0.25*speed)), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("distance", type=int, choices=(8, 16))
    compare(parser.parse_args().distance)
