"""Compare the validated signed R4 fields with the pre-displacement predictions."""
import argparse
import json
from math import fsum
from pathlib import Path

from mouth_rim_pinned_comparison import digest, read, validate_guard


def compare(batch, reduction, prediction):
    validate_guard(reduction)
    signed = read(reduction / "stdout")
    predictions = read(prediction)["cases"]
    summary = []
    for case, forecast in zip(signed["cases"], predictions, strict=True):
        grid = case["grid"]
        folder = batch / grid
        raw = read(folder / "result.json")
        if (forecast["grid"] != grid or forecast["zero_sha256"] != case["zero_sha256"]
                or forecast["surface_sha"] != signed["surface_sha256"]
                or digest(folder / "result.json") != case["result_sha256"]):
            raise ValueError(f"Unmatched fields or forecast on {grid}")
        command = read(folder / "invocation.json")["command"]
        for key, path in zip(("state_sha", "wall_sha", "fields", "load_sha", "axes_sha", "zero_sha"),
                             command[-8:-2], strict=True):
            if digest(Path(path)) != raw["input"][key]:
                raise ValueError(f"Different {key} on {grid}")
        rows = {(row["case"], row["step"]): row for row in raw["signed"]}
        work = {}
        four, eight = forecast["relative_centered"]
        for label, h in (("full", raw["input"]["step"]), ("half", raw["input"]["step"]/2)):
            widths = {}
            for i, name in enumerate(("four", "eight")):
                settled = case["steps"][name][label]["rechecked_energy"]
                predicted = dict(zip(("negative", "positive"), forecast["predicted_one_sided"][i], strict=True))
                predicted["centered"] = forecast["predicted_centered"][i]
                widths[name] = {k: settled[k] - predicted[k] for k in predicted}
            delta = case["four_minus_eight"][label]["centered"]
            common = [rows["four", s*h]["common_crystal_strain_energy_change"] for s in (1, -1)]
            other = [rows["eight", s*h]["common_crystal_strain_energy_change"] for s in (1, -1)]
            if common != other or any(rows["stationary", s*h]["common_crystal_strain_energy_change"] != 0
                                      for s in (1, -1)):
                raise ValueError(f"Different crystal boundary motion on {grid}")
            work[label] = dict(measured_minus_forecast=widths,
                               extension_gap_minus_forecast=fsum((delta, -four, eight)),
                               common_crystal_reference_energy_plus_minus=common,
                               common_reference_crystal_reaction=(common[1]-common[0])/(2*h))
        summary.append(dict(grid=grid, input=raw["input"], work=work,
                            signed_energies=[dict(case=r["case"], step=r["step"], sha256=r["sha256"],
                                rechecked_energy=[s["rechecked_energy"] for s in r["pair"]])
                                for r in raw["signed"]]))
    return dict(observable="Axial reaction is minus the derivative of draining-minus-resting energy. "
                           "All signed shifts carry the same bonded disk and rim except the fixed-wall control.",
                command=read(reduction / "invocation.json")["command"],
                prediction_sha256=digest(prediction), validated_readings=signed,
                comparison=summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("batch", "reduction", "prediction"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.batch, args.reduction, args.prediction), indent=2, allow_nan=False))
