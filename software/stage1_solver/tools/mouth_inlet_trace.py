"""Compare the full feeding-mouth trace with a uniform receiver inlet."""

import json
from pathlib import Path

import numpy as np

from mouth_carrier_momentum import read_state
from mouth_feed_weak_check import dx4
from mouth_release_input_check import manifest_fields


def aperture_state(run, grid, speed, record):
    folder = Path("_scratch/s11c") / run
    outcome = json.loads((folder / "outcome.json").read_text())
    if not outcome["limitsVerified"] or outcome["exitCode"] or outcome["childOutcome"]["exitCode"]:
        raise ValueError(f"Unverified archived run {run}")
    hashes = manifest_fields(Path("docs/_measurements") / record)
    field, current, dx, _ = read_state(folder / "fields", grid, "restoring_slope", speed, hashes)
    if dx != 4/grid or np.min(field[-1].real) <= 0:
        raise ValueError("Unexpected grid or phase branch")
    return field[-1], current, dx


def phase_fit(phase, basis, weights):
    root = np.sqrt(weights)
    coeff, _, rank, _ = np.linalg.lstsq(root[:, None]*basis, root*phase, rcond=None)
    if rank != basis.shape[1]:
        raise ValueError("Unresolved uniform phase or drift")
    error = phase-basis@coeff
    return coeff, np.sqrt(np.average(error**2, weights=weights))


def compare(grid, rest_run, driven_run):
    base, base_j, dx = aperture_state(rest_run, grid, 0, "mouth_feed_weak_run.json")
    coord = (np.arange(grid)+0.5-grid/2)*dx
    xyz = np.stack(np.meshgrid(coord, coord, coord, indexing="ij"), axis=-1)
    aperture = np.sum(xyz**2, axis=-1) < 1
    weights = base_j[aperture]
    basis = np.column_stack((np.ones(len(weights)), xyz[aperture]))
    np.testing.assert_allclose(dx**3*weights.sum(), 0.25, rtol=0, atol=1e-12)
    if np.any(weights <= 0) or np.any(base_j[~aperture]):
        raise ValueError("Archived geometry or current support differs")
    _, base_rms = phase_fit(np.angle(base[aperture]), basis, weights)
    print(f"grid={grid} rest phase RMS after best affine fit={base_rms:.12g} rad")
    x = xyz[..., 0]
    for magnitude in (0.02, 0.01):
        face, current, _ = aperture_state(driven_run, grid, magnitude,
                                          "mouth_intake_redistribution.json")
        # Negative drive is the documented reflected partner about x=0.
        phases, carried, load = [], [], []
        for speed, f, j in ((magnitude, face, current), (-magnitude, face[::-1], current[::-1])):
            np.testing.assert_allclose(dx**3*j.sum(), 0.25, rtol=0, atol=1e-12)
            if np.any(j[~aperture]) or np.any(j[aperture] <= 0):
                raise ValueError("Drive changed the aperture or flow direction")
            phase = np.angle(f[aperture])+speed*x[aperture]
            phases.append(phase)
            density = abs(f)**2
            vx = speed + np.imag(f.conj()*dx4(f, dx))/density
            carried.append(dx**3*np.sum(j*vx))
            # Full-mouth integral of mixed stress by parts; the restoring
            # density slope has compact support and vanishes at the rim.
            weight = np.maximum(1-np.sum(xyz**2, axis=-1), 0)**3
            slope = -2*weight*(1-1/abs(f))
            mixed = dx**3*np.sum(slope*dx4(density, dx))/4
            load.append(carried[-1]+mixed)
        odd = (phases[0]-phases[1])/2
        fit, rms = phase_fit(odd, basis, weights)
        _, const_rms = phase_fit(odd, basis[:, :1], weights)
        incoming = (load[0]-load[1])/2
        transport = (carried[0]-carried[1])/2
        # For the prescribed homogeneous unit-density inlet: zero in-plane
        # carriage; integrated mixed stress = -(pi/3)*s1 on the unit ball.
        required_slope = -3*incoming/np.pi
        print(f"grid={grid} |U|={magnitude:g} odd phase constant_rms={const_rms:.12g} "
              f"affine_x={fit[1]:.12g} affine_rms={rms:.12g} rad")
        print(f"grid={grid} |U|={magnitude:g} incoming={incoming:.12g} "
              f"carried={transport:.12g} carried_fraction={transport/incoming:.12g} "
              f"required_rim_slope={required_slope:.12g}")
        np.testing.assert_allclose(carried[0], -carried[1], rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(load[0], -load[1], rtol=1e-12, atol=1e-14)
        # Check the least-squares residual independently via orthogonality.
        np.testing.assert_allclose(basis.T@(weights*(odd-basis@fit)), 0, atol=1e-11)


if __name__ == "__main__":
    compare(24, "mouth_feed_weak01", "mouth_redistribution_coarse01")
    compare(32, "mouth_feed_depth01", "mouth_redistribution_depth01")
