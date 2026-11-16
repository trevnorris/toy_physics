"""Translate each restoring profile on the stationary paired and isolated traces."""

import numpy as np

from unequal_mouth_holding_readout import archive, packet


def aperture(x, r, center, radius):
    dx = x[1] - x[0]
    # Include the support of both displaced profiles in the energy difference.
    ix = abs(x - center) < radius + dx
    ir = r < radius + dx
    return x[ix] - center, r[ir], np.ix_(ix, ir)


def probe(x, r, radius, pair, alone, current, isolated_current, fixed_rate):
    dx = x[1] - x[0]
    dv = 2 * np.pi * r * dx**2
    axial = x[:, None]

    def base(offset):
        return np.maximum(1 - ((axial - offset)**2 + r**2) / radius**2, 0)

    difference = ((pair - 1)**2 - (alone - 1)**2) / 2
    force = -np.sum(dv * difference * 6 * axial * base(0)**2 / radius**2)

    ring_volume = np.pi * dx * ((r + dx/2)**2 - (r - dx/2)**2)
    def energy(offset, amplitude):
        return np.sum(ring_volume * (amplitude - 1)**2 * base(offset)**3) / 2

    step = dx * 1e-4
    def displaced_force(amplitude):
        return (energy(-step, amplitude) - energy(step, amplitude)) / (2 * step)
    difference_force = displaced_force(pair) - displaced_force(alone)

    # Vary the amplitude by u*eta at fixed profile and test the induced
    # phase flux with eta. The amplitude boundary law is phase-independent.
    eta = axial / radius * base(0)**3
    rate = np.sum(dv * current)
    mean = np.sum(dv * current * eta) / rate if fixed_rate else 0
    mixed = 2 * np.sum(dv * current * (eta - mean)**2)

    def prepared_current(epsilon, variation):
        result = isolated_current * (pair * (1 + epsilon * variation) / alone)**2
        return result * (0.25 / np.sum(ring_volume * result)) if fixed_rate else result

    eps = 1e-5
    np.testing.assert_allclose(current, prepared_current(0, eta), rtol=1e-11, atol=1e-14)
    derivative = (prepared_current(eps, eta) - prepared_current(-eps, eta)) / (2 * eps)
    mixed_fd = np.sum(dv * eta * derivative)
    # A uniform amplitude change has nonzero current-weighted mean, unlike
    # the nearly centered dipole. This distinguishes the two boundary laws.
    uniform_fd = np.sum(ring_volume * (prepared_current(eps, 1) - prepared_current(-eps, 1))) / (2*eps)
    return force, difference_force, mixed, mixed_fd, uniform_fd


def main():
    pairs = archive("lateral_paired_acceleration")["runs"]
    centered = archive("lateral_radius_momentum")["runs"]
    off_center = archive("lateral_target_rate")["runs"]
    for grid in ("e24_c4", "e24_c6", "e48_c4"):
        pair, = (run for run in pairs if run["grid"] == grid and run["radius"] == 1.5)
        right, = (run for run in centered if run["grid"] == grid and run["radius"] == 1.5)
        reference, = pair["reference_fields"]
        left, = (run for run in off_center
                 if run["reference_fields"].get(reference) == pair["reference_fields"][reference])
        if right["reference_fields"] != pair["reference_fields"]:
            raise ValueError(f"Isolated calibration differs for {reference}")
        for distance in (4, 8):
            with packet(pair, f"d{distance}.npz") as data:
                if int(data["d"]) != distance or float(data["target_radius"]) != 1.5:
                    raise ValueError("Unexpected paired geometry")
                currents, face = data["currents"], abs(data["field"][-1])
                for i, (center, radius) in enumerate(((-distance, 1), (0, 1.5))):
                    x, r, indices = aperture(data["x"], data["r"], center, radius)
                    if i == 0:
                        with packet(left, f"source_d{distance}.npz") as control:
                            if control["centers"].tolist() != [center] or control["rates"].tolist() != [0.25]:
                                raise ValueError("Isolated left mouth differs")
                            xc, rc, local = aperture(control["x"], control["r"], center, radius)
                            alone = abs(control["layers"][-1][local])
                            isolated_current = control["currents"][0][local]
                    else:
                        with packet(right, "local_fields.npz") as control:
                            if float(control["radius"]) != radius:
                                raise ValueError("Isolated right mouth differs")
                            xc, rc, local = aperture(control["x"], control["r"], center, radius)
                            isolated_face = control["field"][-1][local]
                            alone = abs(isolated_face)
                            isolated_current = (isolated_face.conj()*control["normal"][local]).imag
                    np.testing.assert_array_equal(x, xc)
                    np.testing.assert_array_equal(r, rc)
                    current = currents[i][indices]
                    f, fd, mixed, mixed_fd, uniform = probe(x, r, radius, face[indices], alone,
                                                           current, isolated_current, fixed_rate=(i == 0))
                    rate = np.sum(2*np.pi*r*(x[1]-x[0])**2 * current)
                    uniform_expected = 0 if radius == 1 else 2*rate
                    print(f"{grid} d={distance} a={radius:g} restoring={f:.12g} energy_fd={fd:.12g} "
                          f"mixed={mixed:.12g} current_fd={mixed_fd:.12g} "
                          f"uniform_rate_fd={uniform:.12g} uniform_expected={uniform_expected:.12g}", flush=True)
                    np.testing.assert_allclose(uniform, uniform_expected, rtol=1e-9, atol=1e-11)
                    np.testing.assert_allclose(f, fd, rtol=3e-6, atol=1e-15)
                    np.testing.assert_allclose(mixed, mixed_fd, rtol=3e-8, atol=1e-13)
                    if not mixed > 0:
                        raise ArithmeticError("Degenerate current variation")


if __name__ == "__main__":
    main()
