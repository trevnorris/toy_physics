"""Direct torque-turn response of an unclamped thermal site patch.

Each gate starts from the archived intake-0.25 warm sites and stays at its
own order for one transit. The centroids and angle of each mobile layer are
free. The loads use the same interior labels as the cold torque assay.
"""
import argparse
import ctypes
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from neutral_patch_neighbor_torque import RECORD, TRACE_DIR, angle_response, angular_terms
from neutral_patch_order_hold import INTERVALS, ORDERS, STRIDES, digest, preparation
from neutral_patch_resting_order_hold import emit
from neutral_patch_sites import BETWEEN, COUNT, GAP, LAWS, PITCH_UNIT, TEMPERATURE, WITHIN, Patch


LOADS = np.array([0., .02, -.02, .01, -.01])
OUTPUT = Path("_scratch/s11c/patch_warm_neighbor_torque_v1")


def central_labels(law):
    prior = json.loads(RECORD.read_text())["guarded_measurements"]
    lines = prior["_scratch/s11c/patch_resting_order_hold_v1/angle_guarded"]["literal_stdout"]
    receipt, = (r for r in map(json.loads, lines) if r["law"] == law and r["chi"] == .95)
    path = TRACE_DIR/f"{law}_0.95.npz"
    fingerprint = digest(path)
    if fingerprint != receipt["trace_sha256"]:
        raise ValueError("Cold entry frame changed")
    with np.load(path) as raw:
        xy = raw["rested_xy"]
    group = np.ones(2*COUNT, dtype=np.int32)
    for offset in (0, COUNT):
        plane = xy[offset:offset+COUNT]
        distances = np.linalg.norm(plane-plane.mean(axis=0), axis=1)
        group[offset+np.argsort(distances, kind="stable")[:6]] = 0
    return group, fingerprint


def compile_kernel(directory):
    binary = Path(directory)/"warm_torque.so"
    subprocess.run(["c++", "-std=c++17", "-O3", "-fno-math-errno", "-fPIC", "-shared",
                    str(Path(__file__).with_suffix(".cpp")), "-o", str(binary)], check=True)
    kernel = ctypes.CDLL(str(binary)).warm_turn
    floats = np.ctypeslib.ndpointer(dtype=np.float64, flags=("C_CONTIGUOUS", "ALIGNED"))
    integers = np.ctypeslib.ndpointer(dtype=np.int32, flags=("C_CONTIGUOUS", "ALIGNED"))
    kernel.argtypes = [ctypes.c_int]*4 + [ctypes.c_double]*5 + [floats]*2 + [integers] + [floats]*7
    kernel.restype = ctypes.c_int
    return kernel


def read_results(law, trace, time):
    # Axes of the raw fit are grid, gate, history, torque, time and group.
    turns = trace[..., 0]-trace[..., 1]
    full = (turns[..., 1, :]-turns[..., 2, :])/(2*LOADS[1])
    half = (turns[..., 3, :]-turns[..., 4, :])/(2*LOADS[3])
    zero = turns[..., 0, :]
    def rel_rms(error, signal):
        return np.linalg.norm(error)/np.linalg.norm(signal)

    for g, chi in enumerate(ORDERS):
        fine = full[-1, g]
        mean = fine.mean(axis=0)
        positive = (turns[-1, g, :, 1]-zero[-1, g])/LOADS[1]
        negative = (turns[-1, g, :, 2]-zero[-1, g])/LOADS[2]
        ends = fine[:, -1]
        thermal_end = zero[-1, g, :, -1]
        integrated = np.trapz(fine, time, axis=1)/time[-1]
        emit("warm_neighbor_turn", law=law, chi=chi, sample_times=(time[[0,512,1024,1536,2048]]).tolist(),
             mean_signed_compliance=mean[[0,512,1024,1536,2048]].tolist(),
             endpoint_compliance_by_history=ends.tolist(), transit_average_compliance=float(integrated.mean()),
             transit_average_compliance_by_history=integrated.tolist(),
             endpoint_parent_pair_means=ends.reshape(-1,2).mean(axis=1).tolist(),
             endpoint_plus_minus_over_torque=[float(positive[:, -1].mean()),float(negative[:, -1].mean())],
             half_endpoint=float(half[-1, g, :, -1].mean()), coarse_endpoint=float(full[0,g,:,-1].mean()),
             mean_trace_half_vs_full_rel_rms=rel_rms(half[-1, g].mean(axis=0)-mean, mean),
             mean_trace_coarse_vs_fine_rel_rms=rel_rms(full[0,g].mean(axis=0)-mean, mean),
             max_history_trace_half_rel_rms=max(rel_rms(b-a, a) for a,b in zip(fine, half[-1,g], strict=True)),
             max_history_trace_grid_rel_rms=max(rel_rms(b-a, a) for a,b in zip(fine, full[0,g], strict=True)),
             no_load_end_turn_by_history=thermal_end.tolist(),
             no_load_max_abs_relative_turn=float(np.max(abs(zero[-1,g]))),
             driven_max_abs_increment=max(np.max(abs(positive)),np.max(abs(negative)))*LOADS[1])


def run(law, output):
    setup = preparation(law)
    group, cold_sha = central_labels(law)
    parents = np.repeat(np.arange(len(setup["parent_seeds"])), 2)
    seeds = 121000+np.arange(len(parents))
    time = np.linspace(0.,setup["transit"],INTERVALS+1)
    steps = np.array(STRIDES)*INTERVALS
    dt = time[-1]/steps
    trace = np.empty((len(STRIDES), len(ORDERS), len(seeds), len(LOADS), len(time), 2))
    final = np.empty(trace.shape[:4]+(2*COUNT,2))
    minimum_overlap, max_angle_step, max_force_defect, max_torque_defect = 1., 0., 0., 0.
    endpoint_error = wrapped_error = 0.
    emit("warm_preparation", law=law, intake=setup["intake"], orders=ORDERS.tolist(),
         gate_times=setup["gate_time"].tolist(), duration=time[-1], temperature=TEMPERATURE,
         torque=LOADS.tolist(), seed=seeds.tolist(), parent_index=parents.tolist(),
         parent_seeds=setup["parent_seeds"].tolist(), inner_labels=np.flatnonzero(group == 0).tolist(),
         cold_entry_sha256=cold_sha, warm_archive_path=setup["archive_path"], warm_archive_sha256=setup["archive_sha256"],
         dt=dt.tolist(), cpp_sha256=digest(Path(__file__).with_suffix(".cpp")), python_sha256=digest(Path(__file__)))
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    with tempfile.TemporaryDirectory() as directory:
        kernel = compile_kernel(directory)
        for history, (parent, seed) in enumerate(zip(parents, seeds, strict=True)):
            rng = np.random.default_rng(seed)
            fine = np.sqrt(2*TEMPERATURE*dt[-1])*rng.standard_normal((steps[-1], 2*COUNT, 2))
            noises = (fine.reshape(-1, 2, 2*COUNT, 2).sum(axis=1), fine)
            for g, chi in enumerate(ORDERS):
                frame = np.ascontiguousarray(setup["gate_xy"][g,parent])
                reference = frame.copy()
                for index in (0, 1):
                    reference[group == index] -= frame[group == index].mean(axis=0)
                terms = angular_terms(frame, group == 0, "pooled")
                scale = chi/ORDERS[0]
                patch = Patch(law)
                patch.depth[patch.cross] *= scale
                for grid, stride in enumerate(STRIDES):
                    for branch, load in enumerate(LOADS):
                        checks = np.zeros(4)
                        last_before = np.empty_like(frame)
                        angles = trace[grid,g,history,branch]
                        xy = final[grid,g,history,branch]
                        status = kernel(COUNT, int(law == "layered"), int(steps[grid]), stride, dt[grid],
                                        scale, GAP, PITCH_UNIT, load, np.array(WITHIN), np.array(BETWEEN),
                                        group, frame, reference, noises[grid], angles, last_before, xy, checks)
                        if status:
                            raise ArithmeticError(f"Invalid native update or fit ({status})")
                        minimum_overlap = min(minimum_overlap, checks[0])
                        max_angle_step = max(max_angle_step, checks[1])
                        max_force_defect = max(max_force_defect, checks[2])
                        max_torque_defect = max(max_torque_defect, checks[3])
                        _, angular_grad, _ = angle_response(last_before.ravel(), terms)
                        expected = last_before + noises[grid][-1] + dt[grid]*(
                            -patch.evaluate(last_before)[1] + load*angular_grad.reshape(-1,2))
                        endpoint_error = max(endpoint_error, np.max(abs(xy-expected)))
                        theta, _, _ = angle_response(xy.ravel(), terms)
                        difference = theta-(angles[-1,0]-angles[-1,1])
                        wrapped_error = max(wrapped_error, abs(np.arctan2(np.sin(difference),np.cos(difference))))
    emit("warm_native_checks", law=law, least_normalized_overlap=minimum_overlap,
         max_unwrap_step=max_angle_step, max_abs_group_angle=float(np.max(abs(trace))),
         unit_torque_group_force_defect=max_force_defect, unit_group_torque_defect=max_torque_defect,
         last_step_reference_error=endpoint_error, endpoint_angle_reference_error=wrapped_error)
    if (not np.isfinite(trace).all() or minimum_overlap <= 1.e-6 or max_angle_step > np.pi/2 or
        max(max_force_defect,max_torque_defect,endpoint_error,wrapped_error) > 1.e-10):
        raise ArithmeticError("Invalid torque balance, reference force, or angular overlap")
    read_results(law, trace, time)
    with tempfile.NamedTemporaryFile(dir=output.parent) as tmp:
        np.savez_compressed(tmp, group=group, angles=trace, final_xy=final, time=time, seeds=seeds,
                            parent_index=parents, orders=ORDERS, loads=LOADS, dt=dt, **setup)
        tmp.flush()
        output.hardlink_to(tmp.name)
    emit("warm_trace", law=law, path=str(output), sha256=digest(output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("law", choices=LAWS)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    run(args.law, args.output or OUTPUT/f"{args.law}.npz")
