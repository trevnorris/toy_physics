"""Independent archive readout using cellwise constants on exact ball intersections."""
import hashlib
import json
from math import pi
from pathlib import Path

import numpy as np


def cut_moment_weights(x, r, h):
    # On positive x, transform x dx = d(x²)/2.  The annular area is
    # pi times the difference of the two positive-part linear functions.
    left, right = (x[:, None]-h/2)**2, (x[:, None]+h/2)**2
    def integral(radius):
        interior = 1-radius**2
        return np.maximum(interior-left, 0)**2-np.maximum(interior-right, 0)**2
    return pi/4*(integral(r-h/2)-integral(r+h/2))


def trace(path, digest):
    with open(path, "rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == digest
    with np.load(path) as z:
        x, r = z["x"], z["r"]
        h = 2*r[0]
        ix = np.abs(x) < 1 + h/4
        jr = r < 1 + h/4
        face = z["layers" if "layers" in z else "field"][-1][ix][:, jr]
        density = face.real**2 + face.imag**2
        p = (density**5-1)/5
        speed = z["currents"][-1][ix][:, jr]/density
        speed_asym = np.max(np.abs(speed-speed[::-1]))
        return x[ix], r[jr], h, p, speed_asym


for line in Path("_scratch/s11c/mouth_front_pressure_dipole02/stdout").read_text().splitlines():
    row = json.loads(line)
    base, pair = [trace(path, pin) for path, pin in row["inputs"].items()]
    x, r, h, p, speed_asym = pair
    assert np.array_equal(x, base[0]) and np.array_equal(r, base[1])
    weights = cut_moment_weights(x[x>0], r, h)
    assert abs(weights.sum()-pi/4)<2e-14
    q = p-base[3]
    contrast = q[::-1][x>0]-q[x>0]
    result = np.sum(weights*contrast)
    mask = weights>2e-15
    print(f"grid={row['grid']} d={row['d']} cell_pressure={result:.13g}"
          f" ratio_to_bilinear={result/row['moment']['pressure']:.10g}"
          f" cells_with_negative_contrast={np.count_nonzero(contrast[mask]<0)}/{mask.sum()}"
          f" paired_v_odd_max={speed_asym:.6g}")
