"""Integrate the signed feeding-carrier load on the intact faces near one drain."""

import json
from pathlib import Path

import numpy as np

from drain_width_traction import normal_stress
from lateral_source_rate import verified_packet


def read_faces(tag, hashes):
    path = Path(f"_scratch/s11c/lateral_bath_{tag}_01/fields/d0.npz")
    with verified_packet(path, hashes) as packet:
        field, currents = packet["field"], packet["currents"]
        x, r, dz = packet["x"], packet["r"], float(packet["dz"])
        if int(packet["d"]) != 0 or len(currents) != 1:
            raise ValueError("An isolated mouth is required")

    dx = x[1] - x[0]
    radius_squared = x[:, None]**2 + r**2
    aperture = radius_squared < 1
    measure = 2*np.pi*r*dx**2
    rate = np.sum(currents[0]*measure)
    if (abs(rate - 0.25) > 1e-12 or np.any((currents[0] > 0) != aperture)
            or not np.isclose((len(field)-1)*dz, 1)):
        raise ValueError("Unexpected aperture, depth or intake")

    # normal_stress discards the outermost lateral cells. Every ball must fit.
    if min(x[-2], -x[1], r[-2]) <= 16:
        raise ValueError("Integration reaches the reservoir stencil")
    upper = normal_stress(field[-1], field[-2], field[-3], r, dx, dz)
    lower = normal_stress(field[0], field[1], field[2], r, dx, dz)
    aperture = aperture[1:-1, :-1]
    radius_squared = radius_squared[1:-1, :-1]
    print(f"grid={tag} intake={rate:.12g}")
    for ball in (2, 4, 8, 16):
        inside = radius_squared <= ball**2
        loads = []
        for top, bottom in zip(upper, lower):
            # The lower face is intact even directly under the opening.
            # The upper aperture carries flux into the receiver and is a separate cut.
            net = np.where(aperture, 0.0, top) - bottom
            loads.append(np.sum((net*measure[:-1])[inside]))
        if not np.all(np.isfinite(loads)):
            raise ValueError("Nonfinite stress on intact faces")
        print(f"  radius={ball} intact_direct={loads[0]:+.12g} intact_pde={loads[1]:+.12g}")


if __name__ == "__main__":
    runs = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())["runs"]
    hashes = {Path(k): v for run in runs for k, v in run.get("fields", {}).items()}
    for grid in ("e24_c4", "e24_c6", "e48_c4"):
        read_faces(grid, hashes)
