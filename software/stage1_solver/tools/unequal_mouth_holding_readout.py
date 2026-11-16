"""Read the two stationary brane loads, in one frame, from matched saved fields."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from lateral_paired_acceleration import target_load
from lateral_source_rate import verified_packet
from mouth_feed_weak_check import dx4


def archive(name):
    return json.loads((Path("docs/_measurements") / f"{name}.json").read_text())


def packet(record, basename):
    if record["guard_exit"] or record["worker_exit"] or not record["limits_verified"]:
        raise ValueError("Unverified input run")
    hashes = {Path(path): digest for path, digest in record["fields"].items()}
    path, = (path for path in hashes if path.name == basename)
    return verified_packet(path, hashes)


def reading(x, r, dz, layers, current, center, radius):
    dx = x[1] - x[0]
    dv = 2*np.pi*r*dx**2
    weight = np.maximum(1-((x[:, None]-center)/radius)**2-(r/radius)**2, 0)**3
    mask = weight > 0
    if not mask.any() or mask[:2].any() or mask[-2:].any() or mask[:, -1].any():
        raise ValueError("Mouth or derivative stencil crosses the edge of the saved patch")
    if np.max(np.abs(current[~mask])) > 1e-13:
        raise ValueError("Current extends outside the requested opening")
    face = layers[-1]
    density = abs(face)**2
    velocity_x = (face.conj() * dx4(face, dx)).imag / density
    carry = np.sum(dv * current * velocity_x)
    total = target_load(face, weight, current/density, SimpleNamespace(dx=dx, dv=dv))

    def normal(f):
        return (11*f[-1] - 18*f[-2] + 9*f[-3] - 2*f[-4]) / (6*dz)

    mixed = (dx4(face, dx).conj()*normal(layers)).real - dx4(normal(abs(layers)**2), dx)/4
    direct = np.sum((dv*mixed)[mask])
    return np.array([carry, total - carry, total, direct]), np.sum(dv*current)


def main():
    pairs = archive("lateral_paired_acceleration")["runs"]
    centered = archive("lateral_radius_momentum")["runs"]
    off_center = archive("lateral_target_rate")["runs"]
    target_archive = archive("lateral_mouth_radius")
    source_checks = archive("lateral_companion_load")["comparisons"]
    radius = 1.5
    for grid in ("e24_c4", "e24_c6", "e48_c4"):
        pair, = (r for r in pairs if r["grid"] == grid and r["radius"] == radius)
        right, = (r for r in centered if r["grid"] == grid and r["radius"] == radius)
        target_record, = (r for r in target_archive["runs"] if r["grid"] == grid and r["radius"] == radius)
        if pair["reference_fields"] != right["reference_fields"]:
            raise ValueError("The centered mouth has a different isolated calibration")
        reference, = pair["reference_fields"]
        left, = (r for r in off_center if r["reference_fields"].get(reference) == pair["reference_fields"][reference])
        with packet(right, "local_fields.npz") as data:
            xr, rr, dzr = data["x"], data["r"], float(data["dz"])
            if float(data["radius"]) != radius:
                raise ValueError("Incorrect isolated target radius")
            layers = data["field"][-4:]
            current = (layers[-1].conj()*data["normal"]).imag
            right_base, right_rate = reading(xr, rr, dzr, layers, current, 0, radius)
        for d in (4, 8):
            with packet(pair, f"d{d}.npz") as data:
                x, r, dz = data["x"], data["r"], float(data["dz"])
                if (int(data["d"]) != d or float(data["target_radius"]) != radius or dz != dzr
                        or not np.array_equal(x[abs(x)<10], xr) or not np.array_equal(r[r<10], rr)):
                    raise ValueError("Incorrect pair geometry or isolated target grid")
                layers, currents = data["field"][-4:], data["currents"]
                left_pair, left_rate = reading(x, r, dz, layers, currents[0], -d, 1)
                right_pair, paired_rate = reading(x, r, dz, layers, currents[1], 0, radius)
            with packet(left, f"source_d{d}.npz") as data:
                if (not np.array_equal(data["centers"], [-d]) or not np.array_equal(data["rates"], [0.25])
                        or not np.array_equal(data["x"], x) or not np.array_equal(data["r"], r)
                        or float(data["dz"]) != dz):
                    raise ValueError("Incorrect off-center isolated control")
                left_base, base_rate = reading(x, r, dz, data["layers"], data["currents"][0], -d, 1)

            first, second = left_pair - left_base, right_pair - right_base
            net = first[2] + second[2]
            mean = (abs(first[2]) + abs(second[2]))/2
            direct_net = first[3] + second[3]
            direct_mean = (abs(first[3]) + abs(second[3]))/2
            print(f"{grid} d={d} left_alone={left_base[2]:.12g} right_alone={right_base[2]:.12g} "
                  f"left_pair={left_pair[2]:.12g} right_pair={right_pair[2]:.12g}", flush=True)
            print(f"delta_left={first[2]:.12g} carry={first[0]:.12g} stress={first[1]:.12g} "
                  f"delta_right={second[2]:.12g} carry={second[0]:.12g} stress={second[1]:.12g}", flush=True)
            print(f"net={net:.12g} percent_of_mean={100*net/mean:.9g} direct_left={first[3]:.12g} "
                  f"direct_right={second[3]:.12g} direct_net={direct_net:.12g} "
                  f"direct_percent={100*direct_net/direct_mean:.9g}", flush=True)

            previous, = (row for row in target_archive["comparisons"]
                         if row["grid"] == grid and row["radius"] == radius and row["d"] == d)
            original_line, = (line for line in target_record["stdout"] if line.startswith(f"d={d} carry="))
            original_direct = float(dict(word.split("=") for word in original_line.split())["direct"])
            source, = (row for row in source_checks if row["grid"] == grid and row["radius"] == radius
                       and row["d"] == d and row["omega"] == 0)
            np.testing.assert_allclose([left_rate, base_rate, right_rate, paired_rate],
                                       [0.25, 0.25, previous["isolated_rate"], previous["paired_rate"]],
                                       rtol=1e-11, atol=1e-13)
            np.testing.assert_allclose([left_pair[2], -second[2], -second[3]],
                                       [source["static_companion_load"], previous["held"], original_direct],
                                       rtol=1e-8, atol=1e-13)


if __name__ == "__main__":
    main()
