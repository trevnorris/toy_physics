"""Cold torque of rotating a labeled interior against its mobile neighbors.

Each gate uses its archived positions as the zero angle. Labels are selected
once per law at entry. Both layers translate and deform freely. An exact local
gauge fixes only the whole patch's translation and Procrustes orientation.
The held angle is invariant to all three common rigid motions.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import root

from neutral_patch_resting_order_hold import LEVELS, digest, emit, hessian
from neutral_patch_sites import COUNT, LAWS, Patch


RECORD = Path("research/pde_ledger_v3/_measurements/neutral_patch_strain_hold.json")
TRACE_DIR = Path("_scratch/s11c/patch_resting_order_hold_v1/traces")
AMPLITUDES = (0.001, 0.0005)


def quarter_turn(xy):
    return np.column_stack((-xy[:, 1], xy[:, 0]))


def angular_terms(xy, inner, meter):
    """An individual angle is atan2(v.x,u.x) for centered reference offsets."""
    terms = []
    planes = (np.arange(COUNT), np.arange(COUNT, 2*COUNT))
    sets = (np.arange(2*COUNT),) if meter == "pooled" else planes
    for members in sets:
        for ids, weight in ((members[inner[members]], 1.), (members[~inner[members]], -1.)):
            offsets = xy[ids]-np.mean(xy[ids], axis=0)
            u, v = np.zeros_like(xy), np.zeros_like(xy)
            u[ids], v[ids] = offsets, quarter_turn(offsets)
            terms.append((weight/len(sets), u.ravel(), v.ravel()))
    return terms


def angle_response(flat, terms):
    angle = 0.
    grad = np.zeros_like(flat)
    curvature = np.zeros((len(flat), len(flat)))
    for weight, u, v in terms:
        a, b = u@flat, v@flat
        denom = a*a+b*b
        angle += weight*np.arctan2(b, a)
        grad += weight*(a*v-b*u)/denom
        curvature += weight*((2*a*b)*(np.outer(u, u)-np.outer(v, v)) +
                             (b*b-a*a)*(np.outer(u, v)+np.outer(v, u)))/denom**2
    return angle, grad, curvature


def free_basis(xy):
    x, y = np.zeros_like(xy), np.zeros_like(xy)
    x[:, 0], y[:, 1] = 1., 1.
    rot = quarter_turn(xy-xy.mean(axis=0))
    return null_space(np.stack([c.ravel() for c in (x, y, rot)])), rot.ravel()


def weak_holds(patch, ref, free, terms, stiffness, tangent):
    """Solve the nonlinear stationarity equations using the actual potential."""
    result = []
    for amplitude in (0., *[sign*a for a in AMPLITUDES for sign in (-1., 1.)]):
        def system(z, jacobian=False):
            flat = ref+free@z[:-1]
            angle, grad, curvature = angle_response(flat, terms)
            reduced = free.T@grad
            torque = z[-1]
            if jacobian:
                lagrangian = free.T@(hessian(patch, flat)-torque*curvature)@free
                return np.block([[lagrangian, -reduced[:, None]],
                                 [reduced[None, :], np.zeros((1, 1))]])
            force = patch.evaluate(flat.reshape(-1, 2))[1].ravel()
            return np.r_[free.T@force-torque*reduced, angle-amplitude]

        initial = np.r_[amplitude*tangent, amplitude*stiffness]
        solution = root(system, initial, jac=lambda z: system(z, True), tol=1.e-10)
        flat = ref+free@solution.x[:-1]
        torque = solution.x[-1]
        residual = max(abs(system(solution.x)))
        angle, angular_grad, angular_hess = angle_response(flat, terms)
        tangent_plane = free@null_space((angular_grad@free)[None, :])
        least_curve = np.linalg.eigvalsh(
            tangent_plane.T@(hessian(patch, flat)-torque*angular_hess)@tangent_plane)[0]
        energy, grad, _ = patch.evaluate(flat.reshape(-1, 2))
        # Full residual tests that the coordinate gauge supplied no real force.
        full_residual = np.max(abs(grad.ravel()-torque*angular_grad))
        result.append(dict(angle=angle, target=amplitude, energy=float(energy.sum()),
                           holding_torque=torque, stationarity_residual=residual,
                           full_residual=full_residual, least_constrained_curvature=least_curve))
    return result


def run_gate(law, level, xy, labels, step):
    patch = Patch(law)
    patch.depth[patch.cross] *= level/.95
    ref = xy.ravel()
    free, rot = free_basis(xy)
    full = hessian(patch, ref)
    def force(flat):
        return patch.evaluate(flat.reshape(-1, 2))[1].ravel()
    direct = np.stack([(force(ref+step*e)-force(ref-step*e))/(2*step)
                       for e in np.eye(len(ref))], axis=1)
    restricted = free.T@full@free
    min_mode = np.linalg.eigvalsh(restricted)[0]
    force_residual = max(abs(force(ref)))
    hessian_error = np.max(abs(full-direct))/np.max(abs(full))
    emit("gate_mechanics", law=law, chi=level, free_degrees=free.shape[1],
         force_residual=force_residual, min_eigenvalue=min_mode, hessian_fd_rel=hessian_error)
    if min_mode <= 0. or force_residual > 1.e-9 or hessian_error > 2.e-8:
        raise ArithmeticError("Unstable or unverified native equilibrium")
    for n, inner in labels.items():
        for meter in ("pooled", "within_plane"):
            terms = angular_terms(xy, inner, meter)
            _, derivative, _ = angle_response(ref, terms)
            reduced = free.T@derivative
            response = np.linalg.solve(restricted, reduced)
            stiffness = 1./(reduced@response)
            finite = weak_holds(patch, ref, free, terms, stiffness, stiffness*response)
            energy_zero = finite[0]["energy"]
            centered = []
            for amplitude, minus, plus in zip(AMPLITUDES, finite[1::2], finite[2::2], strict=True):
                torque_slope = (plus["holding_torque"]-minus["holding_torque"])/(2*amplitude)
                energy_slope = ((plus["energy"]-energy_zero)+(minus["energy"]-energy_zero))/amplitude**2
                centered.append(dict(amplitude=amplitude, torque_slope=torque_slope,
                                     energy_curvature=energy_slope,
                                     torque_rel_error=abs(torque_slope/stiffness-1),
                                     energy_rel_error=abs(energy_slope/stiffness-1)))
            translation_leak = max(abs(derivative.reshape(-1, 2).sum(axis=0)))
            rotation_leak = abs(derivative@rot)
            emit("neighbor_torque", law=law, chi=level, inner_per_plane=n, meter=meter,
                 stiffness=stiffness, common_translation_leak=translation_leak,
                 common_rotation_leak=rotation_leak, finite=finite, centered=centered)
            if (max(translation_leak, rotation_leak) > 1.e-12 or
                max(row['full_residual'] for row in finite) > 1.e-9 or
                min(row['least_constrained_curvature'] for row in finite) <= 0. or
                max(abs(row['angle']-row['target']) for row in finite) > 1.e-10 or
                max(centered[-1]['torque_rel_error'], centered[-1]['energy_rel_error']) > .01):
                raise ArithmeticError("Finite constraint, equilibrium or torque check failed")


def main(trace_dir):
    records = json.loads(RECORD.read_text())["guarded_measurements"]
    lines = records["_scratch/s11c/patch_resting_order_hold_v1/angle_guarded"]["literal_stdout"]
    hashes = {(row["law"], row["chi"]): row["trace_sha256"] for row in map(json.loads, lines)
              if row["test"] == "stationary_angle_control"}
    for law in LAWS:
        for level in LEVELS:
            path = trace_dir/f"{law}_{level:.2f}.npz"
            actual = digest(path)
            emit("input_trace", law=law, chi=level, path=str(path), sha256=actual)
            if hashes[law, level] != actual:
                raise ValueError("Cold preparation archive has changed")
            with np.load(path) as data:
                xy = data["rested_xy"]
            if level == LEVELS[0]:
                ranks = []
                for offset in (0, COUNT):
                    plane = xy[offset:offset+COUNT]
                    distances = np.linalg.norm(plane-plane.mean(axis=0), axis=1)
                    ranks.append(offset+np.argsort(distances, kind="stable"))
                labels = {}
                for n in (4, 6):
                    ids = np.sort(np.concatenate([rank[:n] for rank in ranks]))
                    labels[n] = np.isin(np.arange(2*COUNT), ids)
                    emit("entry_labels", law=law, inner_per_plane=n, inner_indices=ids.tolist())
            run_gate(law, level, xy, labels, 2.e-6)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-dir", type=Path, default=TRACE_DIR)
    main(parser.parse_args().trace_dir)
