"""Combine guarded interior and receiving-side shape energies on identical displacements."""
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUTER = ROOT / "_scratch/s11c/waist_outer01"
INSIDE = ROOT / "_scratch/s11c/live_waist_dynamics01"
STEPS = (-.002, -.001, 0., .001, .002)
AXES = {"radial": "waist", "axial": "depth"}


def checked(directory, sources):
    if (directory / "outcome.json").is_file():
        outcome = json.loads((directory / "outcome.json").read_text())
        normal = outcome["exitCode"] == 0 and outcome["limitsVerified"]
        status = "normal guard exit"
    elif directory == INSIDE / "r1_tight_48_support":
        adopted = json.loads((directory / "terminal_worker_verification.json").read_text())
        normal = (adopted["guard_exit"] is None and adopted["adopted_authorized_worker"]
                  and adopted["worker_exit"] == 0 and adopted["limits_verified"]
                  and all(hashlib.sha256((directory / n).read_bytes()).hexdigest() == digest
                          for n, digest in adopted["evidence_sha256"].items()))
        status = "verified terminal child; final guard receipt absent"
    else:
        raise ValueError(f"Missing terminal guard: {directory}")
    limits = json.loads((directory / "limit-validation.json").read_text())
    duration = json.loads((directory / "duration-validation.json").read_text())
    child = json.loads((directory / "child-outcome.json").read_text())
    if not normal or child["exitCode"] or child["guardReason"] or not (
        limits["verified"] and duration["verified"]
        and limits["actual"]["memory.swap.max"] == "0"
    ):
        raise ValueError(f"Incomplete guard: {directory}")
    raw = (directory / "stdout").read_bytes()
    invocation = json.loads((directory / "invocation.json").read_text())
    if invocation["wallDeadlineSeconds"] is not None or invocation["swapMax"] != 0:
        raise ValueError(f"Unexpected resource settings: {directory}")
    sources[directory.name] = {
        "path": str(directory.relative_to(ROOT)), "stdout_sha256": hashlib.sha256(raw).hexdigest(),
        "worker_exit": child["exitCode"], "status": status, "swap_bytes": 0,
        "memory_limit_bytes": int(limits["actual"]["memory.max"]),
        "command": invocation["command"],
    }
    return [json.loads(line) for line in raw.splitlines()]


def exactly_one(rows):
    if len(rows) != 1:
        raise ValueError(f"Expected one row, found {len(rows)}")
    return rows[0]


def outer_rows(radius, state, mesh, sources):
    name = f"r{radius}{state}{mesh}"
    path = exactly_one([d for d in (OUTER / name, OUTER / ("first_" + name))
                        if (d / "outcome.json").exists()])
    rows = checked(path, sources)
    principal = [r for r in rows if "coordinates" in r]
    if len(principal) != 10 or any((r["radius"], r["state"], r["mesh"]) != (radius, state, mesh)
                                    for r in principal):
        raise ValueError(f"Wrong exterior preparation at {path}")
    sources[path.name]["exterior_relaxation"] = principal
    return rows


def energy_table(center, signed, exterior, onset):
    table = {}
    for axis, coord in AXES.items():
        energies = []
        for h in STEPS:
            ext = exactly_one([r for r in exterior if r["axis"] == axis and r.get("step") == h])
            if h == 0:
                row = center
            elif onset:
                row = exactly_one([r for r in signed if r["axis"] == axis and r["offset"] == h])
            else:
                row = exactly_one([r for r in signed if r["coordinate"] == coord and r["offset"] == h])
            inside = row["coupled_energy" if onset else "total_inside_and_drives"]
            if not math.isclose(inside, sum(row[k] for k in ("order_energy", "carrier_energy", "wave_energy")),
                                abs_tol=1e-9, rel_tol=1e-11):
                raise ValueError(f"Component energies disagree: {axis} {h}")
            energies.append({
                "displacement": h, "inner_order": row["order_energy"],
                "carrier": row["carrier_energy"], "wave": row["wave_energy"],
                "exterior_order": ext["energy"],
                "total": inside + ext["energy"],
            })
        table[axis] = energies
    return table


def conjugate_forces(energies, axis, length):
    points = {r["displacement"]: r for r in energies}
    result = {}
    for key, (left, right) in {
        "tighter_or_shorter": (-.002, 0.), "center_half_step": (-.001, .001),
        "center_full_step": (-.002, .002), "wider_or_longer": (0., .002),
    }.items():
        distance = (right-left if axis == "radial" else
                    length * math.exp(left) * math.expm1(right-left))
        result[key] = {component: (points[left][component] - points[right][component]) / distance
                       for component in ("inner_order", "carrier", "wave", "exterior_order", "total")}
    return result


def measure(radius, state, order_mesh, outer_mesh, sources, onset=False):
    bounded = radius == 1 and state != "broad" and not onset
    name = f"r{radius}_{state}_{order_mesh}" + ("_support" if bounded else "")
    interior = checked(OUTER / f"onset_r1{state}" if onset else INSIDE / name, sources)
    exterior = outer_rows(radius, state, outer_mesh, sources)
    centers = [r for r in interior if r["kind"] == "center" or
               (r["kind"] == "coupled" and r["action"] == 1.)]
    center = exactly_one(centers)
    # The onset snapshot has the same material chart as the bounded static solves.
    if onset:
        reference = checked(INSIDE / f"r1_{state}_72_support", sources)
        geometry = exactly_one([r for r in reference if r["kind"] == "flow"])
    else:
        geometry = center
    table = energy_table(center, [r for r in interior if r["kind"] == "signed"], exterior, onset)
    return {
        "radius": radius, "state": state, "order_mesh": order_mesh, "exterior_mesh": outer_mesh,
        "preparation": "post-wave initial profile" if onset else
        ("bounded crystalline fraction, signed held branch" if bounded else "stationary smooth field"),
        "waist": geometry["waist_radius"], "depth": geometry["axial_depth"],
        "axis_order": center["chi_axis" if onset else "axis_order"],
        "current": center["current"], "terminal_current": center["terminal_current"],
        "zero_wave_cells": center.get("zero_cells", 0), "energies": table,
        "forces": {ax: conjugate_forces(table[ax], ax, geometry["axial_depth"]) for ax in AXES},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = {}
    readings = []
    for radius in (1, 2):
        for state in ("broad", "mid", "tight"):
            for pair in ((48, 40), (72, 70)):
                readings.append(measure(radius, state, *pair, sources))
    convergence = [measure(radius, state, 72, 108, sources)
                   for radius, state in ((1, "broad"), (2, "tight"))]
    initial = [measure(1, state, 72, 70, sources, onset=True) for state in ("mid", "tight")]
    controls = {name: checked(OUTER / name, sources) for name in ("range", "figure")}
    evidence = {
        "experiment": "fixed absolute J=1, sealed lateral side, reversible terminal, training radii only",
        "unit_and_sign": "Forces in model units. Positive q widens midprofile; positive L lengthens; F=-dE/dx.",
        "archived_geometry": "research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json",
        "readings": readings, "exterior_refinement": convergence, "initial_comparison": initial,
        "controls": controls, "sources": sources,
        "source_hashes": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [
            Path(__file__), *(Path(__file__).with_name(n) for n in
              ("mouth_waist_initial_force.py", "mouth_waist_outer_order.py",
               "mouth_waist_localization_figure.py", "mouth_waist_order_range.py",
               "mouth_sealed_shape_relax.py"))]},
    }
    args.output.write_text(json.dumps(evidence, indent=2, allow_nan=False) + "\n")
    for r in readings + convergence + initial:
        forces = [[round(x["total"], 6) for x in v.values()] for v in r["forces"].values()]
        print(r["radius"], r["state"], r["order_mesh"], r["exterior_mesh"], r["waist"], r["depth"],
              round(r["axis_order"],8), forces, r["preparation"])


if __name__ == "__main__":
    main()
