"""Mechanical energy flow of a cubic lattice of marked angular elements."""

import numpy as np


PLANES = ((0, 1), (1, 2), (2, 0))
EDGE_INERTIA, CORNER_INERTIA, SPRING = 2.0, 1.0, 8.0
MASS = np.array([EDGE_INERTIA] * 3 + [CORNER_INERTIA] * 12)


def spring_rows(k):
    rows = np.zeros((12, 15), dtype=complex)
    crossings = []
    for plane, (i, j) in enumerate(PLANES):
        sides = ((i, 1, None), (j, 1, i), (i, -1, j), (j, -1, None))
        for corner, (axis, sign, shift) in enumerate(sides):
            row = 4 * plane + corner
            rows[row, 3 + row] = -1
            rows[row, 3 + 4 * plane + (corner + 1) % 4] = 1
            phase = 1 if shift is None else np.exp(1j * k[shift])
            rows[row, axis] = -sign * phase
            if shift is not None:
                crossings.append((row, shift, axis, sign, phase))
    return rows, crossings


def spectrum(k):
    rows, _ = spring_rows(k)
    scaled = rows / np.sqrt(MASS)
    values, vectors = np.linalg.eigh(SPRING * scaled.conj().T @ scaled)
    return values, vectors / np.sqrt(MASS[:, None])


def measure(k, index, direction):
    values, vectors = spectrum(k)
    freq, coords = np.sqrt(values[index]), vectors[:, index]
    rows, crossings = spring_rows(k)
    strains = rows @ coords
    velocity = -1j * freq * coords
    kinetic = MASS * abs(velocity) ** 2 / 4
    elastic = SPRING * np.vdot(strains, strains).real / 4
    energy = sum(kinetic) + elastic
    current = np.zeros(3)
    for row, shift, axis, sign, phase in crossings:
        # Torque exerted by the spring in this cell on the edge in the next.
        torque = SPRING * sign * strains[row]
        current[shift] += (torque * np.conj(phase * velocity[axis])).real / 2
    h = 2e-5
    freq_minus = np.sqrt(spectrum(k - h * direction)[0][index])
    freq_plus = np.sqrt(spectrum(k + h * direction)[0][index])
    slope = (freq_plus - freq_minus) / (2 * h)
    return freq, current @ direction / energy, slope, sum(kinetic[:3]) / energy


def main():
    speed = np.sqrt(SPRING / (4 * (EDGE_INERTIA + 2 * CORNER_INERTIA)))
    print(f"I={EDGE_INERTIA:g} J={CORNER_INERTIA:g} K={SPRING:g} predicted_speed={speed:.12f}")
    worst_difference = 0.
    for ray in ((1, 0, 0), (1, 1, 1), (1, 2, 3)):
        direction = np.array(ray) / np.linalg.norm(ray)
        for magnitude in (.32, .16, .08, .04):
            k = magnitude * direction
            rows, _ = spring_rows(k)
            open_rows = np.delete(rows, [3, 7, 11], axis=0)
            nullities = tuple(15 - np.linalg.matrix_rank(r) for r in (rows, open_rows))
            observations = [measure(k, index, direction) for index in (4, 5)]
            print(f"ray={ray} k={magnitude:g} null_intact_open={nullities} "
                  f"phase_speeds={[round(o[0]/magnitude, 9) for o in observations]} "
                  f"torque_power_over_energy={[round(o[1], 9) for o in observations]} "
                  f"edge_kinetic_fractions={[round(o[3], 9) for o in observations]}")
            difference = max(abs(o[1] - o[2]) for o in observations)
            worst_difference = max(worst_difference, difference)
            assert nullities == (4, 6)
            if magnitude == .04:
                assert all(abs(o[1] - speed) < .0005 and abs(o[3] - .25) < .0002
                           for o in observations)
    print(f"max_abs_speed_difference_from_frequency_derivative={worst_difference:.10g}")
    assert worst_difference < 1e-6


if __name__ == "__main__":
    main()
