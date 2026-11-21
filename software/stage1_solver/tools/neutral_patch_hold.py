"""Predict from neutral calibration first, then clamp a separately assembled patch."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator


def emit(**row):
    print(json.dumps(row, allow_nan=False), flush=True)


def predict(source, destination):
    # Only the reference patch metrology and axial tension readings enter.
    # The extrapolation from axial to tangential response is a hypothesis.
    raw = source.read_bytes()
    material = json.loads(raw)
    grid = material["normal_separations"]
    force = material["normal_tensions"]
    fine = PchipInterpolator(grid, force, extrapolate=False)
    coarse = PchipInterpolator(grid[::2], force[::2], extrapolate=False)
    gap = material["gap"]
    scale = material["force_unit_measured"]
    length = material["layered_crossover_length"]
    predictions = []
    for law, patch in material["patches"].items():
        pitch = patch["measured_pitch"]
        for fraction in (0., .0625, .125, -.125):
            q = fraction*pitch
            if law == "central":
                radius = np.hypot(gap, q)
                jacobian = q/radius
            elif law == "layered":
                radius = np.hypot(gap, q*q/length)
                jacobian = 2*q**3/(length**2*radius)
            else:
                raise ValueError("Unknown calibration")
            estimate = float(fine(radius)*jacobian)
            row = dict(law=law, fraction=fraction, pitch=pitch, q=q, gap=gap,
                       radial_read_at=float(radius), predicted=estimate,
                       normalized=estimate/scale, coarse_prediction=float(coarse(radius)*jacobian))
            emit(test="blind_prediction", **row)
            predictions.append(row)
    payload = dict(source_hash=hashlib.sha256(raw).hexdigest(), force_scale=scale, readings=predictions)
    destination.write_text(json.dumps(payload, allow_nan=False)+"\n")


def centered_force(patch, positions, ids, removed):
    # Virtual translation of both tagged sites at fixed midpoint. Spectators
    # and grips are at the identical configuration for the control and signal.
    left, right = ids
    diff = np.zeros_like(positions)
    diff[left,0], diff[right,0] = -.5, .5
    def gauge(step):
        energy = patch.evaluate(positions+step*diff, removed)[0]
        return float(energy.sum())
    h = 2e-5
    coarse = (gauge(h)-gauge(-h))/(2*h)
    fine = (gauge(h/2)-gauge(-h/2))/h
    return (4*fine-coarse)/3, coarse


def measure(source, seed, frames, selected_law=None):
    from neutral_patch_sites import LAWS, Patch, prepare

    raw = source.read_bytes()
    frozen = json.loads(raw)
    emit(test="prediction_file", sha256=hashlib.sha256(raw).hexdigest(), seed=seed)
    archive = {}
    for law in LAWS:
        if selected_law is not None and law != selected_law:
            continue
        patch = Patch(law)
        rested, tag, summary = prepare(patch, seed)
        ids = summary["tag"]
        archive[law] = {"rested_xy": rested.tolist(), "sites": ids, "held_states": []}
        frames.write_text(json.dumps(archive, allow_nan=False)+"\n")
        if tag is None:
            emit(test="clamp_unavailable", seed=seed, law=law, reason="no eligible initial tag")
            continue
        left, right = ids
        midpoint = (rested[left]+rested[right])/2
        calibration = [r for r in frozen["readings"] if r["law"] == law]
        for reading in calibration:
            q = reading["q"]
            start = rested.copy()
            start[ids] = midpoint+np.array([[-q/2, 0.], [q/2, 0.]])
            held, residual = patch.relaxed(start, ids)
            energy, gradient, distances = patch.evaluate(held)
            _, control_gradient, _ = patch.evaluate(held, removed=tag)
            intact = bool(patch.neighbors(energy, distances, .6)[tag])
            total, total_coarse = centered_force(patch, held, ids, None)
            control, control_coarse = centered_force(patch, held, ids, tag)
            measured = total-control
            site_force = gradient-control_gradient
            observed = held[right]-held[left]
            balance = float(max(abs(site_force[right]+site_force[left])))
            force_error = max(abs(measured-site_force[right,0]), abs(measured+site_force[left,0]))
            emit(test="held_sites", law=law, seed=seed, ids=ids, tagged_pair_persistent=summary["tag_survived"]["0.6"],
                 fraction=reading["fraction"], observed_xy=observed.tolist(),
                 check_patch_pitch=summary["measured_pitch"], reference_pitch=reading["pitch"],
                 total_holding=total, matched_other_material=control, single_contact=measured,
                 normalized=measured/frozen["force_scale"], predicted=reading["predicted"],
                 reading_minus_prediction=measured-reading["predicted"],
                 half_step_meter_change=measured-(total_coarse-control_coarse),
                 site_force_discrepancy=float(force_error), action_reaction_error=balance,
                 tag_still_registered=intact, max_spectator_gradient=residual)
            archive[law]["held_states"].append({"fraction": reading["fraction"], "xy": held.tolist()})
            if np.max(abs(observed-[q,0.])) > 1e-12 or force_error > 5e-8 or balance > 1e-10:
                raise ArithmeticError("Clamp alignment or virtual-work readout failed")
        frames.write_text(json.dumps(archive, allow_nan=False)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    train = sub.add_parser("predict")
    train.add_argument("calibration", type=Path)
    train.add_argument("output", type=Path)
    test = sub.add_parser("measure")
    test.add_argument("predictions", type=Path)
    test.add_argument("--seed", required=True, type=int)
    test.add_argument("--frames", required=True, type=Path)
    test.add_argument("--law", choices=("central", "layered"))
    args = parser.parse_args()
    if args.mode == "predict":
        predict(args.calibration, args.output)
    else:
        measure(args.predictions, args.seed, args.frames, args.law)
