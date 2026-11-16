"""Compare the two ends of the periodic feeding slab in signed sideways flow."""

import hashlib
import json
from io import BytesIO
from pathlib import Path

import numpy as np

from mouth_feed_weak_check import check_file


def archived_packet(manifest, basename):
    entries = json.loads((Path("docs/_measurements") / manifest).read_text())["runs"]
    (path, digest, outcome), = (
        (Path(name), digest, run["outcome"])
        for run in entries for name, digest in run.get("fields", {}).items()
        if Path(name).name == basename and "depth01" in name
    )
    guard = json.loads((path.parent.parent / "outcome.json").read_text())
    if (outcome["guard_exit"] or outcome["worker_exit"] or not guard["limitsVerified"]
            or guard["exitCode"] or guard["childOutcome"]["exitCode"]):
        raise ValueError(f"Unverified run for {path}")
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"Checksum mismatch for {path}")
    return path


def signed_readings(path, zero_density, speed):
    law = "restoring_slope"
    with np.load(path, allow_pickle=False) as data:
        field, current, mouth = data["field"], data["current"], data["mouth"]
        dx, dz = float(data["dx"]), float(data["dz"])
        if field.shape != (65, 32, 32, 32) or float(data["speed"]) != speed:
            raise ValueError("Unexpected field or drive")
    coord = (np.arange(current.shape[0]) + 0.5 - current.shape[0] / 2) * dx
    x, y, z = np.meshgrid(coord, coord, coord, sparse=True, indexing="ij")
    weight = np.maximum(1 - x*x - y*y - z*z, 0)**3
    readings = []
    for sign in (1, -1):
        if sign == -1:
            # The envelope has the feed plane wave factored out. Reflection
            # changes its x argument and the feed speed, leaving drainage inward.
            field = field[:, ::-1]
            current = current[::-1]
        raw = weight * abs(field[-1])**2 / zero_density
        np.testing.assert_allclose(current, 0.25 * raw / (dx**3 * raw.sum()), rtol=1e-11, atol=1e-14)
        if not np.all(current[mouth] > 0) or np.any(current[~mouth] != 0):
            raise ValueError("Intake left its support or reversed")
        if sign == 1:
            row = check_file(path, law)
        else:
            with BytesIO() as reflected:
                np.savez(reflected, field=field, current=current, mouth=mouth, dx=dx, dz=dz, speed=-speed)
                reflected.seek(0)
                row = check_file(reflected, law)
        residual, _, incident, carry, stress, direct_carry, direct_stress, _ = row
        # The normal current on the first depth link uses the independently
        # solved field next to the reservoir. The holder does not enter here.
        rate = dx**3 * np.imag(field[0].conj() * field[1]).sum() / dz
        feed = incident * rate
        opening = carry + stress
        direct = direct_carry + direct_stress
        print(f"U={incident:+.5f} feed_rate={rate:.12g} feed_momentum={feed:+.12g} "
              f"carried={carry:+.12g} stress={stress:+.12g} aperture={opening:+.12g} "
              f"aperture_layers={direct:+.12g} residual={residual:.4g}", flush=True)
        readings.append(np.array([feed, opening, direct]))
    pos, neg = readings
    odd = (pos - neg) / (2 * speed)
    print(f"abs_U={speed:g} per_speed_feed={odd[0]:.12g} per_speed_aperture={odd[1]:.12g} "
          f"per_speed_layers={odd[2]:.12g} aperture_minus_feed_per_speed={odd[1]-odd[0]:+.12g}",
          flush=True)
    np.testing.assert_allclose(neg, -pos, rtol=1e-11, atol=1e-14)


def main():
    stem = "n32_restoring_slope_u"
    base = archived_packet("mouth_feed_weak_run.json", stem + "0.npz")
    with np.load(base, allow_pickle=False) as data:
        zero_density = abs(data["field"][-1])**2
        if float(data["speed"]) != 0:
            raise ValueError("Calibration is driven")
    np.testing.assert_allclose(zero_density, zero_density[::-1], rtol=1e-12, atol=1e-14)
    for speed in (0.01, 0.02):
        path = archived_packet("mouth_intake_redistribution.json", stem + f"{speed:g}.npz")
        signed_readings(path, zero_density, speed)


if __name__ == "__main__":
    main()
