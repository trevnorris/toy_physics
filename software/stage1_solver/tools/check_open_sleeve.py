"""Independently read stationary currents and calibrate the lateral probe on saved grids."""

import argparse
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from stage1_solver.mouth_probe import aperture_center_x, spherical_rim_probe


def check_field(path):
    with np.load(path, allow_pickle=False) as data:
        r, w = data["r_index"], data["w_index"]
        nr, nw = int(r.max()) + 1, int(w.max()) + 1
        dr, dw = float(data["dr"]), float(data["dw"])
        radial_faces = np.arange(nr + 1) * dr
        areas = 4 * np.pi / 3 * np.diff(radial_faces**3)
        psi = np.zeros((nr, nw), dtype=complex)
        psi[r, w] = data["psi"]
        bath = np.sqrt(data["bath_density"])
        occupied = np.zeros_like(psi, dtype=bool)
        occupied[r, w] = True
        if np.any(abs(psi[occupied]) == 0):
            raise ArithmeticError("An occupied cell has zero density")
        scale = data["hbar"] / data["constituent_mass"]
        rate = float(data["rate"])
        inlet = bath * np.exp(1j * data["feed_phase"])

        def current(left, right, spacing):
            # Face-average density amplitude times the face-normal gradient.
            # Links to a reflecting wall have zero current.
            return scale * np.imag(np.conj((left + right) / 2) * (right - left) / spacing)

        axial = np.zeros((nr, nw + 1))
        axial[:, 0] = current(inlet, psi[:, 0], dw / 2) * occupied[:, 0]
        axial[:, -1] = current(psi[:, -1], bath, dw / 2) * occupied[:, -1]
        axial[:, 1:-1] = current(psi[:, :-1], psi[:, 1:], dw) * (occupied[:, :-1] & occupied[:, 1:])
        rates = areas @ axial
        crossing_error = max(abs(rates - rate))
        print(f"{path.name} rate={rate:g} min_plane_rate={min(rates):.12g} "
              f"max_plane_rate={max(rates):.12g} max_abs_plane_error={crossing_error:.3e}")

        saved_sleeve = int(data["sleeve_cells"])
        mouth_radius = int(data["mouth_cells"]) * dr
        if not np.all(occupied[:, :saved_sleeve] == (radial_faces[:-1] < mouth_radius)[:, None]):
            raise ValueError("Saved sleeve does not match the declared aperture")
        if not np.isfinite(crossing_error) or crossing_error > 1e-10:
            raise ArithmeticError("An integrated plane current disagrees with the prescribed rate")
        return mouth_radius


def check_probe(radius):
    mu, weights = leggauss(32)
    for force in (0.0, 0.015, 0.0075):
        traction = spherical_rim_probe(radius, mu, force)
        measured = 2 * np.pi * radius**2 * np.dot(weights, mu * traction)
        print(f"rim_load={force:.12g} integrated_force_x={measured:.12g}")
        if abs(measured - force) > 1e-13:
            raise ArithmeticError("Rim force calibration failed")
    for offset in (0, 0.1*radius, -0.1*radius):
        rim = offset*mu + np.sqrt(radius**2 - offset**2 * (1 - mu**2))
        measured = aperture_center_x(rim, mu, weights)
        print(f"test_geometry_x={offset:.12g} centroid_x={measured:.12g}")
        if abs(offset - measured) > 1e-12 * radius:
            raise ArithmeticError("Geometric readout failed")
    distorted = radius * (1 + 0.05 * (3 * mu**2 - 1) / 2)
    measured = aperture_center_x(distorted, mu, weights)
    print(f"centered_quadrupole_centroid_x={measured:.3e}")
    if abs(measured) > 1e-12 * radius:
        raise ArithmeticError("A centered quadrupole produced a false translation")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fields", nargs="+", type=Path)
    radii = [check_field(path) for path in parser.parse_args().fields]
    for radius in sorted(set(radii)):
        check_probe(radius)
