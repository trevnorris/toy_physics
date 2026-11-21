"""Constrain relative plane displacement at archived draining-patch gates.

Each invocation runs one phase and one law. The incremental actuator force is
obtained by differencing the raw internal forces between sister branches. The
deterministic direct bath reaction and the hold's initial impulse area are
recorded separately. The bath's random reaction on the mean is common to all
branches within a history and therefore disappears from their differences.
"""
import argparse
import ctypes
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from neutral_patch_conversion import Sites
from neutral_patch_sites import BETWEEN, COUNT, GAP, PITCH_UNIT, TEMPERATURE, WITHIN


ARCHIVES = Path("_scratch/s11c/conversion_pair_histories")
CONVERSION = Path("research/pde_ledger_v3/_measurements/neutral_patch_conversion.json")
TRACKING = Path("research/pde_ledger_v3/_measurements/mouth_contact_tracking.json")
OUTPUT = Path("_scratch/s11c/strain_hold")
ORDERS = np.array([.95, .5, .1, .01])
AMPLITUDES = np.array([0., .001, -.001, .0005, -.0005])
INTERVALS = 2048
# The grid has 6 (coarse) or 12 (fine) steps per output interval.
STRIDES = (6, 12)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preparation(law):
    manifest = json.loads(CONVERSION.read_bytes())
    archive = ARCHIVES / f"{law}.npz"
    if digest(archive) != manifest["trajectory_archives"][str(archive)]:
        raise ValueError(f"Archive changed: {archive}")
    run = manifest["verified_runs"][f"converting_{law}_p1_01"]
    setup = json.loads(run["literal_stdout"][0])
    if (str(archive) not in run["command"] or setup["law"] != law or setup["power"] != 1
            or setup["drive_order"][0] != 0 or not run["limits_verified"]):
        raise ValueError("Input must be a verified unforced linear-softening preparation")
    clock_runs = json.loads(TRACKING.read_bytes())["guarded_runs"]["contact_tracking_01"]
    clock, = (row for row in map(json.loads, clock_runs["literal_stdout_for_reported_cases"])
              if row["field_sha256"] == setup["fields"])
    marks = {r["order"]: r for r in clock["readings"]}
    with np.load(archive) as raw:
        chi, t = raw["chi"], raw["t"]
        gates = np.abs(chi[:, None]-ORDERS).argmin(axis=0)
        np.testing.assert_allclose(chi[gates], ORDERS, rtol=0., atol=1e-12)
        np.testing.assert_allclose(t[gates], [marks[c]["elapsed"] for c in ORDERS], atol=1e-10, rtol=0.)
        np.testing.assert_allclose(raw["w"][gates], [marks[c]["depth"] for c in ORDERS], atol=1e-10, rtol=0.)
        np.testing.assert_allclose(setup["intake"], .25, rtol=0., atol=1e-10)
        frames = raw["xy_zero"][gates]
        if frames.shape != (len(ORDERS), len(setup["seeds"]), 2*COUNT, 2):
            raise ValueError(f"Unexpected unforced frames: {frames.shape}")
        return dict(gate_xy=frames, gate_index=gates, gate_time=t[gates],
                    gate_w=raw["w"][gates], gate_order=chi[gates],
                    gate_dilation=np.array([marks[c]["angular_spacing_ratio"] for c in ORDERS]),
                    transit=t[-1]-t[0], intake=setup["intake"], parent_seeds=np.array(setup["seeds"]),
                    archive_path=str(archive), archive_sha256=manifest["trajectory_archives"][str(archive)])


def compile_kernel(directory):
    path = Path(directory) / "assay.so"
    source = Path(__file__).with_suffix(".cpp")
    subprocess.run(["c++", "-O3", "-fno-math-errno", "-std=c++17", "-shared", "-fPIC",
                    str(source), "-o", str(path)], check=True)
    kernel = ctypes.CDLL(str(path)).assay
    array = np.ctypeslib.ndpointer(dtype=np.float64, flags=("C_CONTIGUOUS", "ALIGNED"))
    kernel.argtypes = [ctypes.c_int]*4 + [ctypes.c_double]*4 + [array]*8
    kernel.restype = ctypes.c_int
    return kernel


def waveform(time, transit, phase):
    if phase == "hold":
        return np.ones_like(time), np.zeros_like(time)
    duration = transit / 4
    active = (time > 0) & (time < duration)
    angle = np.pi*time[active]/duration
    strain, rate = np.zeros_like(time), np.zeros_like(time)
    strain[active] = np.sin(angle)**2
    rate[active] = (np.pi/duration)*np.sin(2*angle)
    return strain, rate


def raw_force(reference, points, scale):
    grad = reference.read(points, scale)[0][:, :, 0]
    return (grad[:, COUNT:].sum(axis=1)-grad[:, :COUNT].sum(axis=1))/2


def simulate(phase, law, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    initial = preparation(law)
    frames = initial["gate_xy"]
    nparents = len(initial["parent_seeds"])
    parent_index = np.repeat(np.arange(nparents), 2)
    seeds = (51000 if phase == "hold" else 61000) + np.arange(len(parent_index))
    transit = initial["transit"]
    time = np.linspace(0, transit, INTERVALS+1)
    force = np.empty((len(STRIDES), len(ORDERS), len(seeds), len(AMPLITUDES), len(time)))
    strain_shape, rate_shape = waveform(time, transit, phase)
    strain = AMPLITUDES[:, None]*strain_shape
    bath = (COUNT/2)*AMPLITUDES[:, None]*rate_shape
    impulse = (COUNT/2)*AMPLITUDES if phase == "hold" else np.zeros_like(AMPLITUDES)
    grids = [np.linspace(0, transit, INTERVALS*stride+1) for stride in STRIDES]
    drives = [AMPLITUDES[:, None]*waveform(t, transit, phase)[0] for t in grids]
    dt = np.array([t[1]-t[0] for t in grids])
    last = np.empty(force.shape[:-1] + (2*COUNT, 2))
    clamp_error = np.zeros(1)
    max_error = reference_error = 0.
    within, between = np.asarray(WITHIN), np.asarray(BETWEEN)
    reference = Sites(law, [[0, COUNT]])
    shift = np.repeat([-.5, .5], COUNT)[None, :]
    with tempfile.TemporaryDirectory() as directory:
        kernel = compile_kernel(directory)
        for history, (parent, seed) in enumerate(zip(parent_index, seeds, strict=True)):
            rng = np.random.default_rng(seed)
            fine = np.sqrt(2*TEMPERATURE*dt[1])*rng.standard_normal((len(grids[1])-1, 2*COUNT, 2))
            increments = (fine.reshape(-1, 2, 2*COUNT, 2).sum(axis=1), fine)
            for g, order in enumerate(initial["gate_order"]):
                frame = frames[g, parent]
                scale = order/initial["gate_order"][0]
                start = np.repeat(frame[None], len(AMPLITUDES), axis=0)
                start[:, :, 0] += shift*strain[:, :1]
                expected = raw_force(reference, start, scale)
                for level, stride in enumerate(STRIDES):
                    trace = force[level, g, history]
                    final = last[level, g, history]
                    for branch in range(len(AMPLITUDES)):
                        status = kernel(COUNT, int(law == "layered"), len(grids[level])-1, stride,
                                        dt[level], scale, GAP, PITCH_UNIT, within, between, frame,
                                        increments[level], drives[level][branch], trace[branch],
                                        final[branch], clamp_error)
                        if status:
                            raise ArithmeticError(f"Nonfinite pair separation or force: {status}")
                        max_error = max(max_error, clamp_error.item())
                    reference_error = max(reference_error, np.max(abs(trace[:, 0]-expected)),
                                          np.max(abs(trace[:, -1]-raw_force(reference, final, scale))))
    checks = dict(max_plane_mean_error=float(max_error), max_endpoint_reference_error=float(reference_error))
    print(json.dumps(dict(test="endpoint_and_clamp_checks", **checks)), flush=True)
    if not np.isfinite(force).all() or max_error > 1e-12 or reference_error > 2e-11:
        raise ArithmeticError("Constraint or native-vs-Python force check failed")
    metadata = dict(phase=phase, law=law, orders=ORDERS, amplitudes=AMPLITUDES, parent_index=parent_index,
                    seeds=seeds, time=time, dt=dt, stride=np.array(STRIDES), pulse_duration=transit/4,
                    temperature=TEMPERATURE, normal_gap=GAP, direct_drag=COUNT/2, force=force,
                    strain=strain, bath=bath, bath_initial_impulse=impulse, final_xy=last,
                    cpp_sha256=digest(Path(__file__).with_suffix(".cpp")), python_sha256=digest(Path(__file__)),
                    **checks, **initial)
    # Publish only after the checks and the complete archive write.
    with tempfile.NamedTemporaryFile(dir=output.parent, suffix=".npz") as tmp:
        np.savez_compressed(tmp, **metadata)
        tmp.flush()
        output.hardlink_to(tmp.name)
    print(json.dumps(dict(test="acquisition", phase=phase, law=law, archive_sha256=initial["archive_sha256"],
                          output=str(output), output_sha256=digest(output), shape=list(force.shape),
                          parent_index=parent_index.tolist(), seeds=seeds.tolist(), dt=dt.tolist(),
                          transit=transit, temperature=TEMPERATURE)), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=("hold", "pulse"))
    parser.add_argument("law", choices=("central", "layered"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    simulate(args.phase, args.law, args.output or OUTPUT / f"{args.phase}_{args.law}.npz")
