"""Reduce force at fixed population displacement, retaining the separate spatial-cut reading."""
import argparse
import json
from pathlib import Path
from statistics import mean

import numpy as np

from mouth_bulk_temperature_readout import guarded, interval
from neutral_patch_sites import emit


ARMS = ("zero", "plus", "minus", "half_plus", "half_minus")


def read_history(root, rung, history):
    states, frames = {}, {}
    for arm in ARMS:
        folder = root / f"r{rung}h{history}_{arm}"
        setup, *observations = guarded(folder)
        states[arm], frames[arm] = setup, observations
        assert [x["frame"] for x in frames[arm]] == list(range(49))
        assert states[arm]["branch"] == arm and states[arm]["rung"] == rung and states[arm]["history"] == history
        with np.load(folder / "native/frames.npz", allow_pickle=False) as saved:
            assert saved["pos"].shape == (49, 128, 4) and saved["q"] == states[arm]["q"]
    base = states["zero"]
    for state in states.values():
        for key in ("source_sha256", "noise_seed", "fixed_box", "area_per_face",
                    "temperature", "dt", "steps_per_frame", "transit"):
            assert state[key] == base[key]
        np.testing.assert_allclose(state["baseline_coordinate"], base["baseline_coordinate"], rtol=0, atol=1.e-13)
    offsets = [states[arm]["q"] for arm in ARMS]
    np.testing.assert_allclose(offsets, [0., .0002, -.0002, .0001, -.0001], rtol=0, atol=0)
    for arm in ARMS:
        np.testing.assert_allclose([r["centroid"] for r in frames[arm]], states[arm]["target"], rtol=0, atol=1.e-11)
        np.testing.assert_allclose([r["interval"]["mean_thermal_traction"] for r in frames[arm][1:]],
                                   [r["interval"]["mean_thermal_traction"] for r in frames["zero"][1:]],
                                   rtol=0, atol=1.e-12)
    return base, frames, offsets[1]


def paired(responses, q):
    zero, plus, minus, half_plus, half_minus = responses
    full = (plus-minus)/(2*q)
    half = (half_plus-half_minus)/q
    return dict(full=full, half=half, full_minus_half=full-half,
                even_full=((plus-zero)+(minus-zero))/(2*q),
                even_half=((half_plus-zero)+(half_minus-zero))/q)


def reduce_history(root, rung, history):
    base, frames, q = read_history(root, rung, history)
    initial = paired([frames[a][0]["pair_tagged_traction"] for a in ARMS], q)
    emit(test="initial", rung=rung, history=history, temperature=base["temperature"], **initial)
    for quarter in range(4):
        start = 1+12*quarter
        ranges = [frames[a][start:start+12] for a in ARMS]
        # Every integration step has equal duration. Each interval mean covers the
        # whole interval, while spatial-plane stresses are sampled at frame ends.
        def interval_force(key):
            return paired([mean(r["interval"][key] for r in window) for window in ranges], q)
        response = interval_force("mean_conditional_traction")
        applied = interval_force("mean_actuator_traction")
        for key in response:
            np.testing.assert_allclose(response[key], applied[key], rtol=1.e-8, atol=1.e-7)
        plane = paired([mean(r["plane_mean_traction"] for r in window) for window in ranges], q)
        tagged = paired([mean(r["pair_tagged_traction"] for r in window) for window in ranges], q)
        emit(test="quarter", rung=rung, history=history, quarter=quarter+1,
             tagged_actuator_all_steps=response, tagged_pair_force_at_frames=tagged,
             spatial_pair_force_at_frames=plane)
    row = frames["zero"][-1]
    emit(test="zero_control", rung=rung, history=history, density=base["fixed_density"],
         changed_side=row["changed_side"], ever_changed_side_at_frames=row["ever_changed_side_at_frames"],
         original_neighbors=row["initial"], original_missing=row["original_missing"],
         minimum_cluster=min(r["largest_cluster"] for r in frames["zero"]),
         final_quarter_material_pressure=mean(r["material_pressure"] for r in frames["zero"][-12:]))
    return response


def run(root):
    cases = json.loads((root / "plan.json").read_text())["cases"]
    assert set(json.loads((root / "DONE.json").read_text())["verified"]) == {c[0] for c in cases}
    assert len(cases) == 20
    reduce_history(root, 0, 0)
    hot = [reduce_history(root, 4, history) for history in range(3)]
    for key in ("full", "half", "full_minus_half", "even_full", "even_half"):
        values = [h[key] for h in hot]
        emit(test="hot_late_summary", field=key, histories=values,
             mean=mean(values), approximate_95_across_three_heating_histories=interval(values))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch", type=Path)
    run(parser.parse_args().batch)
