"""Reduce signed, stationary R4 energies using equal physical supports on both grids."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_text())


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def validate_guard(folder):
    outcome = read(folder / "outcome.json")
    limits = read(folder / "limit-validation.json")
    duration = read(folder / "duration-validation.json")
    if not (
        outcome["exitCode"] == outcome["childOutcome"]["exitCode"] == 0
        and outcome["childOutcome"]["guardReason"] is None
        and outcome["limitsVerified"] and limits["verified"]
        and limits["actual"]["memory.swap.max"] == "0"
        and duration["verified"] and duration["wallDeadlineSeconds"] is None
        and duration["inactivityDeadlineSeconds"] is None
    ):
        raise ValueError(f"Unvalidated field solve at {folder}")


def paired(row, key):
    return math.fsum((row[0][key], -row[1][key]))


def slopes(records, magnitude, zero_values):
    """Axial reaction is minus the derivative of driven-minus-resting energy."""
    plus = records[magnitude]["pair"]
    minus = records[-magnitude]["pair"]
    outputs = {}
    for key, baseline in zero_values.items():
        pos = paired(plus, key)
        neg = paired(minus, key)
        outputs[key] = {
            "positive": -(pos - baseline) / magnitude,
            "negative": (neg - baseline) / magnitude,
            "centered": (neg - pos) / (2 * magnitude),
        }
    return outputs


def reduce_job(folder, prior, constraint_work):
    validate_guard(folder)
    result = read(folder / "result.json")
    meta = result["inputs"]
    mesh = meta["grid"]
    if read(folder.parent / "DONE.json") != {"double": 0, "local": 0} or (folder.parent / "ISSUE.json").exists():
        raise ValueError("Incomplete signed batch")
    if meta != read(folder / "inputs.json") or result["zero"] != read(folder / "zero.json"):
        raise ValueError(f"Signed batch summary disagrees with its individual records on {mesh}")
    command = read(folder / "invocation.json")["command"]
    for source, key in zip(command[4:9], ("state", "wall", "field", "load", "reference_axes"), strict=True):
        if digest(Path(source)) != meta[key + "_sha256"]:
            raise ValueError(f"Input file changed on {mesh}: {key}")
    old = prior["cases"][mesh + "_q2"]
    estimate = constraint_work["measurements"][mesh]
    for field in ("field", "wall", "load"):
        if old["source"][field] != meta[field + "_sha256"] or estimate["source"][field] != meta[field + "_sha256"]:
            raise ValueError(f"{mesh} changed the {field}")
    signed = {}
    zero = result["zero"]
    zero_energy = paired(zero, "rechecked_energy")
    baselines = dict(energy=zero_energy, frozen_energy=zero_energy,
                     order_pin_relative_to_carried=0., crystal_pin_relative_to_fixed_nodal_load=0.)
    histories = list(zero)
    with np.load(folder / "zero.npz", allow_pickle=False) as archive:
        reference_fields = [archive[k] for k in ("wound", "rest")]
    for row in result["signed"]:
        archive = folder / row["file"]
        if digest(archive) != row["sha256"] or row != read(archive.with_suffix(".json")):
            raise ValueError(f"Signed state hash mismatch at {archive}")
        with np.load(archive, allow_pickle=False) as saved:
            for key, reference, stats in zip(("wound", "rest"), reference_fields, row["pair"], strict=True):
                u = saved[key]
                if u.shape != reference.shape or not np.array_equal(
                    np.max(np.abs(u - reference), axis=1), stats["response_max"]
                ) or not np.isfinite(u).all() or stats["energy"] != stats["rechecked_energy"]:
                    raise ValueError(f"Saved {key} field does not reproduce its readout at {archive}")
        table = signed.setdefault(row["name"], {})
        if row["step"] in table or row["pin_info"]["minimum_axial_jacobian"] <= 0:
            raise ValueError(f"Unusable signed displacement in {archive}")
        table[row["step"]] = row
        histories.extend(row["pair"])
    step = meta["full_step"]
    expected = {"four", "eight", "four_minus_eight"}
    if signed.keys() != expected or any(
        table.keys() != {step, -step, step/2, -step/2} for table in signed.values()
    ):
        raise ValueError(f"Incomplete signed data in {folder}")
    max_stationarity = max(row["stationarity"] for row in histories)
    current_error = max(abs(row["feed"] - row["bulk"]) for row in histories)
    if max_stationarity >= meta["stationarity_tolerance"] or current_error >= 3.e-5:
        raise ValueError(f"Stationarity or reservoir current failed at {folder}")
    answers = {
        name: {label: slopes(table, h, baselines) for label, h in (("full", step), ("half", step/2))}
        for name, table in signed.items()
    }
    for item in estimate["banks"]:
        name = {(4, 8): "four_minus_eight", (4,): "four", (8,): "eight"}[tuple(item["widths"])]
        predicted = dict(zip(item["signed_steps"], item["predicted_pin_corrected_reaction"], strict=True))
        for label, h in (("full", step), ("half", step/2)):
            actual = answers[name][label]["energy"]
            answers[name][label]["settled_minus_fixed_state_prediction"] = {
                "positive": actual["positive"] - predicted[h],
                "negative": actual["negative"] - predicted[-h],
                "centered": actual["centered"] - (predicted[h] + predicted[-h])/2,
            }
    # A receiver table yields two absolute reactions. A zero-wall duct control
    # directly measures the derivative for their difference in the interior.
    gaps = {"receiver": {}, "duct_relative": {}}
    for label in ("full", "half"):
        gaps["receiver"][label] = {
            key: answers["four"][label]["energy"][key] - answers["eight"][label]["energy"][key]
            for key in ("centered", "positive", "negative")
        }
        gaps["duct_relative"][label] = answers["four_minus_eight"][label]["energy"]
    old_gaps = {}
    for row in old["probes"]:
        if row["strip"] == "near":
            work = {r["width"]: r["signed_energy"] for r in row["readings"]}
            old_gaps[row["bank"]] = work[4] - work[8]
    return dict(inputs=meta, outcomes=answers, four_minus_eight=gaps,
                previous_convected_gaps=old_gaps, max_stationarity=max_stationarity,
                max_current_imbalance=current_error, result_sha256=digest(folder / "result.json"),
                zero_state_sha256=digest(folder / "zero.npz"), zero=zero, signed=result["signed"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prior", type=Path)
    parser.add_argument("constraint_work", type=Path)
    parser.add_argument("jobs", type=Path, nargs=2)
    args = parser.parse_args()
    prior = read(args.prior)
    constraint_work = read(args.constraint_work)
    cases = {}
    for folder in args.jobs:
        case = reduce_job(folder, prior, constraint_work)
        meta = case["inputs"]
        cases[meta["grid"]] = case
    comparison = {}
    coarse, fine = (cases[mesh] for mesh in ("double", "local"))
    for key in ("wall_sha256", "state_sha256", "reference_axes_sha256", "terminal_knots",
                "receiver_bounds", "duct_bounds", "full_step"):
        if coarse["inputs"][key] != fine["inputs"][key]:
            raise ValueError(f"Mismatched physical inputs for {key}")
    for task in ("receiver", "duct_relative"):
        a, b = (case["four_minus_eight"][task]["half"]["centered"] for case in (coarse, fine))
        comparison[task] = dict(doubled=a, rim_refined=b, signed_change=b-a)
    print(json.dumps(dict(prior_sha256=digest(args.prior), pin_work_sha256=digest(args.constraint_work), comparison=comparison,
                          cases=[cases[k] for k in sorted(cases)]), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
