"""Check signed hot shear over one transit and the separate flare-force scale."""
import hashlib
import json
from pathlib import Path
from statistics import mean

MEASUREMENTS = Path("research/pde_ledger_v3/_measurements")


def manifest(name):
    return json.loads((MEASUREMENTS / f"{name}.json").read_text())


def emit(kind, **values):
    print(json.dumps({"test": kind, **values}, allow_nan=False))


def frames(rec):
    path = Path(rec["stdout_path"])
    stream = path.read_bytes()
    assert hashlib.sha256(stream).hexdigest() == rec["stdout_sha256"]
    assert rec["guard_exit"] == rec["child_exit"] == 0 and rec["limits_verified"]
    setup, *rows = map(json.loads, stream.splitlines())
    assert [r["frame"] for r in rows] == list(range(49))
    return setup, rows


def signed_history(runs, rung, history):
    prefix = f"r{rung}h{history}_"
    base, undriven = frames(runs[prefix + "undriven"])
    assert base["shears"] == [0.] and base["history"] == history and base["rung"] == rung
    odd, even = {}, {}
    for phase, amplitude in (("shear", .004), ("half_shear", .002)):
        setup, driven = frames(runs[prefix + phase])
        for key in ("history", "rung", "source_sha256", "pressure_time", "transit", "dt", "axis"):
            assert setup[key] == base[key]
        assert setup["shears"] == [amplitude, -amplitude]
        assert setup["applied_pressures"] == [0., 0.]
        odd[phase], even[phase] = [], []
        for row, ref in zip(driven, undriven, strict=True):
            assert row["elapsed"] == ref["elapsed"]
            plus, minus = (s["shear"] for s in row["readings"])
            background = ref["readings"][0]["shear"]
            slope = (plus - minus) / (2 * amplitude)
            assert abs(slope - row["shear_slope"]) < 1.e-9
            odd[phase].append(slope)
            even[phase].append(((plus-background) + (minus-background)) / (2 * amplitude))
    for window, start, end in (("instant", 0, 1), ("q1", 1, 13), ("q2", 13, 25),
                               ("q3", 25, 37), ("q4", 37, 49)):
        full = mean(odd["shear"][start:end])
        half = mean(odd["half_shear"][start:end])
        emit("signed_window", temperature=base["temperature"], history=history, window=window,
             signed_full=full, signed_half=half, full_minus_half=full-half,
             even_full=mean(even["shear"][start:end]),
             even_half=mean(even["half_shear"][start:end]))


def saved_rows(record):
    assert record["guard_exit"] == record["child_exit"] == 0 and record["limits_verified"]
    return list(map(json.loads, record["literal_stdout"].splitlines()))


def flare_comparison():
    push = manifest("mouth_occupied_slow_push")["guarded_runs"]["occupied_push_comparison_01"]
    rows = saved_rows(push)
    slopes, = [r["slopes"] for r in rows if r["test"] == "fixed_entrance_and_middle_flare_probe"]
    gap = slopes[1] - slopes[0]
    emit("trial_flare_scale", slopes=slopes, midpoint=mean(slopes), half_gap=gap/2)
    slip = manifest("mouth_replacing_flare_slip")["guarded_runs"]
    rows = saved_rows(slip["_scratch/s11c/replacing_flare_34_2_256_01"])
    cases = [r for r in rows if r["test"] == "hold" and r["tau"] == 1/32
             and r["t_over_transit"] == 1.]
    assert {r["wall"] for r in cases} == {"smooth", "lattice"} and len(cases) == 2
    for r in cases:
        slopes = r["force_per_signed_displacement"]
        assert max(slopes)-min(slopes) < 1.e-8
        emit("matched_coating", wall=r["wall"], displacement=r["signed_displacements"],
             incremental_force_per_displacement=slopes[0])


if __name__ == "__main__":
    runs = manifest("mouth_bulk_temperature_steps")["runs"]
    for rung in (0, 4):
        for history in range(3):
            signed_history(runs, rung, history)
    flare_comparison()
