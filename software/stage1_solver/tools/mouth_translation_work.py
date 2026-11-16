"""Virtual translation of the density boundary and a reciprocity check for the open intake."""

import json
from hashlib import sha256
from pathlib import Path

import numpy as np


def read_face(path, digest):
    if sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"Unverified field {path}")
    with np.load(path, allow_pickle=False) as data:
        field = data["layers"] if int(data["d"]) else data["field"]
        return data["x"], data["r"], abs(field[-1]), data["currents"][-1]


def conservative_work(x, r, amplitude):
    # X>0 moves the target toward the source, located on the negative x axis.
    dx = x[1]-x[0]
    dv = 2*np.pi*r*dx**2
    def energy(offset):
        base = np.maximum(1-(x[:, None]+offset)**2-r**2, 0)
        return np.sum(dv*base**3*(amplitude-1)**2)/2
    base = np.maximum(1-x[:, None]**2-r**2, 0)
    analytic = np.sum(dv*3*x[:, None]*base**2*(amplitude-1)**2)
    h = 1e-4*dx
    differenced = (energy(-h)-energy(h))/(2*h)
    if not np.isfinite([analytic,differenced]).all() or abs(analytic-differenced)>1e-12:
        raise ArithmeticError("Boundary energy displacement derivative disagrees")
    return np.array([analytic, differenced])


def open_intake_check(x, r, amplitude, current):
    # For the density-responsive continuation from the isolated calibration,
    # v_normal is fixed in shape up to one common factor at fixed total rate.
    dx = x[1]-x[0]
    dv = 2*np.pi*r*dx**2
    weight = np.maximum(1-x[:, None]**2-r**2, 0)**3
    dipole = x[:, None]*weight
    rate = np.sum(dv*current)
    speed_shape = current/amplitude**2
    def intake(eps):
        trial_amplitude = amplitude*(1+eps*dipole)
        raw = speed_shape*trial_amplitude**2
        return 0.25*raw/np.sum(dv*raw)
    eps = 1e-4
    dj = (intake(eps)-intake(-eps))/(2*eps)
    mean = np.sum(dv*current*dipole)/rate
    expected = 2*np.sum(dv*current*(dipole-mean)**2)
    mixed = np.sum(dv*dipole*dj)
    if (abs(rate-0.25)>2e-13 or abs(np.sum(dv*dj))>2e-11
            or not 0 < expected or not np.isclose(mixed,expected,rtol=1e-8,atol=1e-14)):
        raise ArithmeticError("Density variation failed the fixed-rate or mixed-variation check")
    print(f"  rate={rate:.12g} integral_dj={np.sum(dv*dj):.4g} "
          f"phase_test_times_dj={mixed:.12g} twice_weighted_variance={expected:.12g}")


def main():
    archive = json.loads(Path("docs/_measurements/lateral_reservoir_pair.json").read_text())
    hashes = {Path(name): digest for run in archive["runs"]
              for name, digest in run.get("fields", {}).items()}
    for grid in ("e24_c6", "e48_c4"):
        folder = Path(f"_scratch/s11c/lateral_bath_{grid}_01/fields")
        path = folder/"d0.npz"
        x, r, amplitude, current = read_face(path, hashes[path])
        print(grid)
        open_intake_check(x, r, amplitude, current)
        baseline = conservative_work(x, r, amplitude)
        for d in (3, 4, 6, 8):
            path = folder/f"d{d}.npz"
            pair_x, pair_r, pair_amplitude, _ = read_face(path, hashes[path])
            if not np.array_equal(pair_x, x) or not np.array_equal(pair_r, r):
                raise ValueError("Mismatched paired and isolated geometry")
            load, check = conservative_work(x, r, pair_amplitude)-baseline
            print(f"  d={d} restoring_work_force={load:.12g} energy_difference_force={check:.12g}")


if __name__ == "__main__":
    main()
