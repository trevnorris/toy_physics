"""Check paired archived fields and form signed force estimates independently."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_text())


def sha(path):
    with path.open("rb") as data:
        return hashlib.file_digest(data, "sha256").hexdigest()


def inspect(folder, estimates):
    inputs = read(folder / "inputs.json")
    outcomes = read(folder / "outcome.json")
    assert outcomes["limitsVerified"] and outcomes["stderrBytes"] == 0
    assert outcomes["exitCode"] == outcomes["childOutcome"]["exitCode"] == 0
    assert outcomes["childOutcome"]["guardReason"] is None
    limits = read(folder / "limit-validation.json")
    assert limits["verified"] and limits["actual"]["memory.swap.max"] == "0"
    duration = read(folder / "duration-validation.json")
    assert duration["verified"] and duration["wallDeadlineSeconds"] is None
    invocation = read(folder / "invocation.json")["command"]
    for index, key in enumerate(("state", "wall", "field", "load", "reference_axes"), 4):
        assert sha(Path(invocation[index])) == inputs[key + "_sha256"]
    fields = ("wound", "rest")
    with np.load(folder / "zero.npz") as archive:
        zero_fields = [archive[key] for key in fields]
    zeros = read(folder / "zero.json")
    result = read(folder / "result.json")
    assert result["zero"] == zeros and result["inputs"] == inputs
    summaries = {r["file"]: r for r in result["signed"]}
    assert len(summaries) == len(result["signed"]) == 12
    stationarity = [v["stationarity"] for v in zeros]
    mass_defects = [v["port_error"] for v in zeros]
    names = ("four", "eight", "four_minus_eight")
    readings = {}
    for name, prediction in zip(names, (estimates["banks"][1], estimates["banks"][2], estimates["banks"][0]), strict=True):
        assert tuple(prediction["widths"]) == {"four": (4,), "eight": (8,), "four_minus_eight": (4, 8)}[name]
        cells = {}
        for label, factor in (("plus", 1), ("plus_half", .5), ("minus", -1), ("minus_half", -.5)):
            row = read(folder / (name + "_" + label + ".json"))
            h = row["step"]
            assert h == factor * inputs["full_step"] and row == summaries[row["file"]]
            state = folder / row["file"]
            assert sha(state) == row["sha256"]
            assert row["pin_info"]["minimum_axial_jacobian"] > .99
            pair = row["pair"]
            with np.load(state) as archive:
                assert tuple(archive.files) == fields
                for j, key in enumerate(fields):
                    saved = archive[key]
                    assert saved.shape == zero_fields[j].shape and np.isfinite(saved).all()
                    measured = np.max(abs(saved - zero_fields[j]), axis=1)
                    assert np.array_equal(measured, pair[j]["response_max"])
                    assert pair[j]["energy"] == pair[j]["rechecked_energy"]
                    assert pair[j]["energy"] <= pair[j]["frozen_energy"] + 1e-12
                    assert (pair[j]["feed"] > 0) if j == 0 else (abs(pair[j]["feed"]) < 1e-12)
                    mass_defects.append(abs(pair[j]["feed"] - pair[j]["bulk"]))
                    stationarity.append(pair[j]["stationarity"])
            secants = {}
            for quantity in ("energy", "frozen_energy", "order_pin_relative_to_carried",
                             "crystal_pin_relative_to_fixed_nodal_load"):
                initial = [v["energy"] for v in zeros] if quantity in ("energy", "frozen_energy") else [0., 0.]
                # Subtract the zero level within each physical preparation before the
                # draining-minus-resting comparison. Force is the negative energy slope.
                secants[quantity] = (
                    (pair[1][quantity] - initial[1]) - (pair[0][quantity] - initial[0])
                ) / h
            entry = prediction["signed_steps"].index(h)
            secants["local_prediction"] = prediction["predicted_pin_corrected_reaction"][entry]
            cells[label] = secants
        for suffix in ("", "_half"):
            cells["centered" + suffix] = {
                key: (value + cells["minus" + suffix][key]) / 2
                for key, value in cells["plus" + suffix].items()
            }
        readings[name] = cells
    assert max(stationarity) < inputs["stationarity_tolerance"]
    assert max(mass_defects) < 3e-5
    return dict(result_sha256=sha(folder / "result.json"), maximum_residual=max(stationarity),
                maximum_current_difference=max(mass_defects), reactions=readings,
                receiving_four_minus_eight={key: v["energy"] - readings["eight"][key]["energy"]
                                            for key, v in readings["four"].items()})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("predictions", type=Path)
    args = parser.parse_args()
    assert read(args.run / "DONE.json") == {"double": 0, "local": 0}
    assert not (args.run / "ISSUE.json").exists()
    targets = read(args.predictions)["measurements"]
    print(json.dumps({level: inspect(args.run / level, targets[level]) for level in targets},
                     indent=2, allow_nan=False))
