"""Compare fixed-feed phase and ideal receiver phase for identical archived currents."""

import json
from pathlib import Path

import numpy as np
from scipy.special import xlogy

from lateral_source_rate import verified_packet


def receiver_phase(x0, x, r, dx, current):
    # On the symmetry axis the Neumann kernel 1/distance**2 has an exact
    # volume integral over every cylindrical ring cell, including at the probe.
    def primitive(d, s):
        return xlogy(d, d*d + s*s) - 2*d + 2*s*np.arctan2(d, s)

    left, right = x[:, None] - x0 - dx/2, x[:, None] - x0 + dx/2
    inner, outer = r-dx/2, r+dx/2
    rings = np.pi*((primitive(right, outer)-primitive(left, outer))
                   - (primitive(right, inner)-primitive(left, inner)))
    return -np.sum(rings*current)/(2*np.pi**2)


def read(tag, hashes):
    path = Path(f"_scratch/s11c/lateral_bath_{tag}_01/fields/d0.npz")
    with verified_packet(path, hashes) as p:
        if int(p["d"]) != 0 or len(p["currents"]) != 1:
            raise ValueError("Expected an isolated mouth")
        x, r, face, current = p["x"], p["r"], p["field"][-1], p["currents"][0]
    dx = 2*r[0]
    rate = np.sum(current*(2*np.pi*r*dx**2))
    if abs(rate-0.25) > 1e-12 or np.any(current < 0):
        raise ValueError("Unexpected number current")
    # Density is positive and the trace stays on a single phase branch.
    if np.min(face[:, :2].real) <= 0:
        raise ValueError("Feed phase needs unwrapping")
    phase = np.angle(face[:, :2])
    axis = (9*phase[:, 0] - phase[:, 1])/8  # even quadratic through two rings
    def contrast(values):
        return values[1] - (values[0]+values[2])/2
    probes = [-0.5, 0, 0.5]
    feeding = contrast(np.interp(probes, x, axis))
    first_ring = contrast(np.interp(probes, x, phase[:, 0]))
    bulk = contrast([receiver_phase(q, x, r, dx, current) for q in probes])
    print(f"{tag} rate={rate:.12g} feed_axis={feeding:+.12g} "
          f"feed_first_ring={first_ring:+.12g} receiver_axis={bulk:+.12g} "
          f"receiver_minus_feed={bulk-feeding:+.12g}")


if __name__ == "__main__":
    archive = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    hashes = {Path(k): v for run in archive["runs"] for k, v in run.get("fields", {}).items()}
    for tag in ("e24_c4", "e24_c6", "e48_c4"):
        read(tag, hashes)
