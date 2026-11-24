"""Independently check both R4 clearance records against their saved signed measurements.

Run under s11c_guarded_run.py with the scratch and measurement directories as arguments.
This is a read-only reviewer. It does not import the evaluator or either production reducer.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(actual, expected, atol=1e-13):
    np.testing.assert_allclose(actual, expected, rtol=0, atol=atol)


def verify_guard(folder):
    outcome = read(folder / "outcome.json")
    assert outcome["limitsVerified"]
    assert outcome["exitCode"] == outcome["childOutcome"]["exitCode"] == 0
    assert int(read(folder / "effective-limits.json")["memory.swap.max"]) == 0


def review_shift(scratch, measurements):
    run = scratch / "r4_clearance_shift01"
    archive = read(measurements / "mouth_r4_clearance_shift.json")
    previous = scratch / "r4_bucket_unfolded01"
    old_force = read(previous / "depth/result.json")["rows"]
    sectors = []
    for index in range(6):
        verify_guard(run / f"guard_sector_{index}")
        sectors.append(read(run / f"sector_{index}/result.json")["rows"])
    for level in (1, 2, 4):
        verify_guard(run / f"guard_labels_{level}")
        for case in ("unshifted", "shifted"):
            path = run / f"labels_{level}/{case}.npz"
            report = read(path.parent / "result.json")["cases"][case]
            assert digest(path) == report["sha256"]
            with np.load(path) as saved:
                for tag in report["tests"]:
                    physical_order = saved[tag + "_order_t"]
                    expected = physical_order.copy()
                    if case == "shifted":
                        expected[1:-1] += np.diff(physical_order)[1:] / 4
                    np.testing.assert_array_equal(saved[tag + "_label_t"], expected)
                    if case == "unshifted":
                        old_path = (previous / f"{level}/labels.npz" if tag == "zero"
                                    else previous / f"depth/signed_{level}.npz")
                        with np.load(previous / f"{level}/labels.npz") as base, np.load(old_path) as old:
                            active = base["active_stencil"]
                            old_key = "labels" if tag == "zero" else tag + "_labels"
                            close(saved[tag + "_labels"][active], old[old_key][active],
                                  np.finfo(float).eps if tag == "zero" else 0)
    for degree in (6, 8):
        verify_guard(run / f"guard_static_{degree}")
        static = read(run / f"static_{degree}/result.json")
        integrated = np.array(static["raw_by_region"])
        for level in (1, 2, 4):
            for case in ("unshifted", "shifted"):
                key = f"{case}_{level}_{degree}"
                rows = [sector[key] for sector in sectors]
                pairs = np.array([r["energies"] for r in rows], dtype=np.longdouble).sum(axis=0)
                force = np.array([r["force"] for r in rows], dtype=np.longdouble).sum(axis=0)
                close(force, (pairs[0] - pairs[1]) / (2 * archive["h"]), 1e-8)
                local = np.array([r["local"] for r in rows], dtype=np.longdouble).sum(axis=(0, 1))
                close(force, local)
                bank = static["axes_runs"].index(f"{case}_{level}")
                energy, order = (integrated[:, bank] - integrated[:, 0]).sum(axis=0)
                depth = float(force[1] - force[0])
                published = archive["rows"][key]
                close([energy, order, depth, depth - order],
                      [published[k] for k in ("energy", "order", "depth", "remaining")])
                if case == "unshifted":
                    close(depth, old_force[f"{level}_{degree}"]["spatial_minus_native"])


def review_smooth(scratch, measurements):
    run = scratch / "r4_clearance_smooth01"
    reported = read(measurements / "mouth_r4_clearance_smooth.json")
    hold = read(measurements / "mouth_r4_two_metric_labels.json")["data"]["held_fields"]
    masks = read(measurements / "mouth_r4_clearance_regions.json")
    assert digest(scratch / "r4_clearance_regions01/masks.npz") == masks["masks_sha256"]
    verify_guard(run / "bath_guard")
    bath = read(run / "bath.json")
    for row in bath["rows"].values():
        close(sum(row["components"]), row["force"], 1e-15)
    for degree_name, results in reported["degrees"].items():
        degree = int(degree_name)
        shape = read(run / f"shape_{degree}/result.json")
        depth = read(run / f"depth/{degree}/result.json")
        verify_guard(run / f"shape_guard_{degree}")
        verify_guard(run / f"depth_guard_{degree}")
        assert shape["mask_sha256"] == depth["mask_sha256"] == masks["masks_sha256"]
        assert shape["checkpoint"] == depth["checkpoint_sha256"] == reported["checkpoint_sha256"]
        assert shape["common_quadrature"] == depth["partition_sha256"]
        for published in results["levels"]:
            level = str(published["level"])
            s, d = shape["results"][level], depth["rows"][level]
            banks = np.array(d["analytic_by_region"])
            regions = np.diff(banks, axis=1)[:, 0]
            fixed = np.array(s["synthetic_region_delta"])
            close(regions[:, [0, 2]], fixed[:, :2], 2e-12)
            close(regions[:, 1], regions[:, 2] + regions[:, 3])
            combined = np.c_[fixed[:, 0], regions[:, 1], fixed[:, 2:], regions[:, 2:]]
            close(combined, published["by_region"])
            close(combined.sum(axis=0), published["totals"])
            close(np.diff(s["material_spatial"], axis=0)[0], fixed.sum(axis=0))
            np.testing.assert_allclose(np.diff(d["full_force"]), [regions[:, 1].sum()], atol=1e-13, rtol=0)
            for i, raw in enumerate(d["fd"]):
                minus, plus = np.array(raw["signed_energies"])
                signed_total = ((minus - plus) / (2 * raw["h"])).sum(axis=0)
                signed_total += bath["rows"][str(raw["h"])]["force"]
                measured = d["full_force"] - signed_total
                close(measured, np.array(published["analytic_minus_fd_total_by_step"]["depth"])[:, i], 2e-12)
            h = np.array(shape["signed_steps"])
            pairs = np.array(s["signed_crystal_energies"])
            signed = (pairs[..., 0] - pairs[..., 1]) / (2 * h[None, :, None])
            separate_fields = [[m["central"][str(step)]["force"] for m in hold["modes"][1:]] for step in h]
            total = signed + separate_fields
            expected = np.array(s["analytic_total"])[:, None, :] - total
            close(expected, published["analytic_minus_fd_total_by_step"]["other_five"], 1e-11)
            assert np.max(np.abs(published["analytic_minus_smaller_central"])) < 1e-7
        totals = abs(np.array([row["totals"][:7] for row in results["levels"]]))
        assert np.all(np.diff(totals, axis=0) < 0)
        by_region = abs(np.array([row["by_region"] for row in results["levels"]]))
        observed_orders = np.log2(by_region[:-1] / by_region[1:])
        close(observed_orders, results["observed_orders_by_region"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scratch", type=Path)
    parser.add_argument("measurements", type=Path)
    args = parser.parse_args()
    review_shift(args.scratch, args.measurements)
    review_smooth(args.scratch, args.measurements)
    print("Verified both fixed-state R4 clearance controls.")
