"""Replicate the two low-order layered hold forecasts on the same four parent frames."""
import argparse
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

import neutral_patch_order_hold as assay
from neutral_patch_strain_forecast import apply_hold, emit, kinematics, rms, sha


OUT = Path("_scratch/s11c/strain_hold_layered_ensemble")
ORIGINAL = Path("_scratch/s11c/strain_hold/predicted01.npz")
ORIGINAL_SHA = "3f25c98da8ac3cd330ad95e6743c1549b487c4e60f456b1e557b2a5398a41ebc"
KERNEL_SHA = "6224e323c8c86ef0192a4fd464620772eaaaa6b4bfd7200acea912d83a9372b3"
ASSAY_SHA = "20a51ff660a0853053790acad5a3371ab405511b6c35a5c78ebc6901276db9b5"
ORDERS = np.array([.1, .01])
PARENT_COUNT = 4
BATH_COUNT = 16
SEED_START = dict(hold=71000, pulse=81000)
BOOT_SEED = 92001
BOOT_DRAWS = 4000


def verify(path, expected):
    if not expected or sha(path) != expected:
        raise ValueError(f"Invalid digest for {path}")


def save(path, **data):
    with path.open("xb") as fp:
        np.savez_compressed(fp, **data)
    return sha(path)


def batch_file(phase, parent):
    return OUT / f"{phase}_parent_{parent}.npz"


def frozen_file():
    return OUT / "frozen64.npz"


def acquire(phase, parent, frozen_sha):
    if phase == "pulse":
        verify(frozen_file(), frozen_sha)
    verify(Path(assay.__file__), ASSAY_SHA)
    verify(Path(assay.__file__).with_suffix(".cpp"), KERNEL_SHA)
    path = batch_file(phase, parent)
    if path.exists():
        raise FileExistsError(path)
    prepared = assay.preparation("layered")
    gate_ids = [assay.ORDERS.tolist().index(order) for order in ORDERS]
    frames = prepared["gate_xy"][gate_ids, parent]
    scale = prepared["gate_order"][gate_ids] / prepared["gate_order"][0]
    seeds = SEED_START[phase] + BATH_COUNT*parent + np.arange(BATH_COUNT)
    time = np.linspace(0, prepared["transit"], assay.INTERVALS+1)
    stride = assay.STRIDES[-1]
    steps = stride * assay.INTERVALS
    grid = np.linspace(0, time[-1], steps+1)
    dt = grid[1]
    amplitudes = assay.AMPLITUDES[:3]
    drive = amplitudes[:, None]*assay.waveform(grid, time[-1], phase)[0]
    force = np.empty((len(ORDERS), BATH_COUNT, len(amplitudes), len(time)))
    reference = assay.Sites("layered", [[0, assay.COUNT]])
    final = np.empty((len(amplitudes), 2*assay.COUNT, 2))
    shift = np.repeat([-.5, .5], assay.COUNT)
    error = np.zeros(1)
    max_error = reference_error = 0.
    with tempfile.TemporaryDirectory() as tmp:
        kernel = assay.compile_kernel(tmp)
        for b, seed in enumerate(seeds):
            rng = np.random.default_rng(seed)
            noise = np.sqrt(2*assay.TEMPERATURE*dt)*rng.standard_normal((steps, 2*assay.COUNT, 2))
            for g, frame in enumerate(frames):
                start = np.repeat(frame[None], len(amplitudes), axis=0)
                start[:, :, 0] += drive[:, :1]*shift
                trace = force[g, b]
                for branch in range(len(amplitudes)):
                    status = kernel(assay.COUNT, 1, steps, stride, dt, scale[g], assay.GAP,
                                    assay.PITCH_UNIT, np.asarray(assay.WITHIN), np.asarray(assay.BETWEEN),
                                    frame, noise, drive[branch], trace[branch], final[branch], error)
                    if status:
                        raise ArithmeticError(f"Native kernel returned {status}")
                    max_error = max(max_error, error.item())
                reference_error = max(reference_error,
                                      np.max(np.abs(assay.raw_force(reference, start, scale[g])-trace[:, 0])),
                                      np.max(np.abs(assay.raw_force(reference, final, scale[g])-trace[:, -1])))
    emit(test="clamp_reference", phase=phase, parent=parent, max_plane_mean_error=max_error,
         max_endpoint_force_error=reference_error)
    if not np.isfinite(force).all() or max_error > 1e-12 or reference_error > 2e-11:
        raise ArithmeticError("Constraint or native-vs-reference check failed")
    # At fixed angular freedom, only the two x plane means and normal gap are held.
    odd = (force[:, :, 1]-force[:, :, 2])/(2*amplitudes[1])
    even = (force[:, :, 1]+force[:, :, 2]-2*force[:, :, 0])/(2*amplitudes[1])
    fingerprint = save(path, odd=odd, even=even, seeds=seeds, parent=parent, parent_seed=prepared["parent_seeds"][parent],
                       gate_xy=frames, orders=ORDERS, time=time, dt=dt, phase=phase, law="layered",
                       amplitudes=amplitudes, forecast_sha256=frozen_sha or "", archive_sha256=prepared["archive_sha256"])
    emit(test="layered_parent", phase=phase, parent=parent, bath_seeds=[int(seeds[0]), int(seeds[-1])],
         source_archive_sha256=prepared["archive_sha256"], pinned_assay=ASSAY_SHA, pinned_cpp=KERNEL_SHA,
         frozen_sha256=frozen_sha, output=str(path), sha256=fingerprint)


def collect(phase, expected=None):
    curves, pins = [], {}
    for parent in range(PARENT_COUNT):
        path = batch_file(phase, parent)
        pins[str(path)] = sha(path)
        with np.load(path) as data:
            if data["phase"].item() != phase or data["parent"].item() != parent:
                raise ValueError(f"Invalid phase or parent in {path}")
            if phase == "pulse" and data["forecast_sha256"].item() != expected:
                raise ValueError("Pulse run refers to a different prediction")
            start = SEED_START[phase]+BATH_COUNT*parent
            np.testing.assert_array_equal(data["seeds"], np.arange(start, start+BATH_COUNT))
            np.testing.assert_array_equal(data["orders"], ORDERS)
            if parent:
                np.testing.assert_array_equal(time, data["time"])
            time = data["time"]
            curves.append(data["odd"])
    return time, np.concatenate(curves, axis=1), pins


def freeze():
    verify(ORIGINAL, ORIGINAL_SHA)
    time, held, pins = collect("hold")
    shape = kinematics(time)[0]
    predicted = apply_hold(held, shape)
    baseline = held[:, :, 0, None]*shape
    digest = save(frozen_file(), predicted=predicted, held=held, baseline=baseline, time=time, orders=ORDERS,
                  parent_index=np.repeat(np.arange(PARENT_COUNT), BATH_COUNT),
                  seeds=SEED_START["hold"]+np.arange(PARENT_COUNT*BATH_COUNT),
                  old_forecast_sha256=ORIGINAL_SHA, hold_pins=json.dumps(pins))
    emit(test="layered_64_frozen", sha256=digest, hold_inputs=pins, old_forecast_sha256=ORIGINAL_SHA,
         full_record_tolerance=.1, late="after T/4, normalized to full measured RMS")


def metrics(predicted, actual, late):
    error = predicted-actual
    norm = rms(actual)
    return dict(full=rms(error)/norm, late=rms(error[..., late])/norm,
                late_absolute=rms(error[..., late]), late_signal=rms(actual[..., late])/norm)


def resample_scores(forecasts, pulse, late):
    """Conditional empirical percentile ranges with each of the four parents fixed."""
    rng = np.random.default_rng(BOOT_SEED)
    methods = ("prior8_fixed", "hold64_fixed", "hold64_new_baths", "no_relaxation")
    per_method = {name: {"full": [], "late": []} for name in methods}
    n = PARENT_COUNT*BATH_COUNT
    fixed = {name: curve.mean(axis=1) for name, curve in forecasts.items()}
    for first in range(0, BOOT_DRAWS, 50):
        batch = min(50, BOOT_DRAWS-first)
        # Equal parent weight, multinomial counts for the 16 independent future
        # baths within each parent. Both gates see the same selected histories.
        counts = rng.multinomial(BATH_COUNT, np.full(BATH_COUNT, 1/BATH_COUNT),
                                 size=(2, batch, PARENT_COUNT))
        w_hold, w_pulse = counts.reshape(2, batch, n)/n
        actual = w_pulse @ pulse
        trial = {"prior8_fixed": fixed["prior8"], "hold64_fixed": fixed["hold64"],
                 "hold64_new_baths": w_hold @ forecasts["hold64"], "no_relaxation": fixed["no_relaxation"]}
        for name, estimated in trial.items():
            values = metrics(estimated[:, None, :] if estimated.ndim == 2 else estimated, actual, late)
            for key in per_method[name]:
                per_method[name][key].append(values[key])
    return {name: {metric: np.quantile(np.concatenate(parts, axis=1), [.025, .975], axis=1).T
                   for metric, parts in results.items()} for name, results in per_method.items()}


def score(expected):
    verify(frozen_file(), expected)
    verify(ORIGINAL, ORIGINAL_SHA)
    time, actual, pulse_pins = collect("pulse", expected)
    with np.load(ORIGINAL) as prior, np.load(frozen_file()) as new:
        np.testing.assert_array_equal(prior["time"], time)
        np.testing.assert_array_equal(new["time"], time)
        for path, digest in json.loads(new["hold_pins"].item()).items():
            verify(Path(path), digest)
        if set(new["seeds"]) & set(SEED_START["pulse"]+np.arange(PARENT_COUNT*BATH_COUNT)):
            raise ValueError("Overlapping acquisition baths")
        locs = [prior["orders"].tolist().index(order) for order in ORDERS]
        forecasts = {"prior8": prior["layered_predicted"][-1, locs], "hold64": new["predicted"],
                     "no_relaxation": new["baseline"]}
        labels = {"prior8": prior["parent_index"], "hold64": new["parent_index"],
                  "no_relaxation": new["parent_index"]}
    late = time > .25*time[-1]
    pn = labels["hold64"]
    mean = actual.mean(axis=1)
    full_scores = {name: metrics(pred.mean(axis=1), mean, late) for name, pred in forecasts.items()}
    loo = {}
    for name, pred in forecasts.items():
        loo[name] = [metrics(pred[:, labels[name] != p].mean(axis=1), actual[:, pn != p].mean(axis=1), late)
                     for p in range(PARENT_COUNT)]
    spread = {}
    # Use the full ensemble's scale for split-half scores, as in the eight-bath acquisition.
    for name, curve in (("actual", actual), ("hold64", forecasts["hold64"])):
        split = curve[:, ::2].mean(axis=1)-curve[:, 1::2].mean(axis=1)
        spread[name] = {"full": rms(split)/rms(curve.mean(axis=1)),
                        "late": rms(split[..., late])/rms(curve.mean(axis=1))}
    bootstrap = resample_scores(forecasts, actual, late)
    paths = save(OUT / "pulse64.npz", time=time, orders=ORDERS, actual=actual,
                 parent_index=pn, seeds=SEED_START["pulse"]+np.arange(len(pn)), frozen_sha256=expected,
                 pulse_pins=json.dumps(pulse_pins))
    emit(test="layered_ensemble_validation", prior_forecast=ORIGINAL_SHA, enlarged_forecast=expected,
         direct_curves_sha256=paths, pulse_inputs=pulse_pins, bootstrap_seed=BOOT_SEED, bootstrap_draws=BOOT_DRAWS)
    for g, order in enumerate(ORDERS):
        row = {name: {k: float(v[g]) for k, v in result.items()} for name, result in full_scores.items()}
        for name in row:
            row[name]["loo_full"] = [float(v["full"][g]) for v in loo[name]]
            row[name]["loo_late"] = [float(v["late"][g]) for v in loo[name]]
        intervals = {name: {k: v[g].tolist() for k, v in values.items()} for name, values in bootstrap.items()}
        split = {name: {k: float(v[g]) for k, v in values.items()} for name, values in spread.items()}
        emit(test="layered_ensemble_score", order=order, methods=row, split_half=split,
             conditional_bath_percentiles=intervals)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("hold", "freeze", "pulse", "score"))
    parser.add_argument("--parent", type=int, choices=range(PARENT_COUNT))
    parser.add_argument("--forecast-sha")
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True, parents=True)
    emit(test="command", argv=sys.argv, source_sha256=sha(Path(__file__)),
         forecast_source_sha256=sha(Path(assay.__file__).with_name("neutral_patch_strain_forecast.py")))
    if args.phase in ("hold", "pulse"):
        if args.parent is None:
            parser.error("acquisition requires --parent")
        acquire(args.phase, args.parent, args.forecast_sha)
    elif args.phase == "freeze":
        freeze()
    else:
        score(args.forecast_sha)
