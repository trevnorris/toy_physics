"""Primitive angular springs outside a mouth with a prescribed radial texture.

R is the half-width of the inner cubic sampling cut in lattice units.  All
exterior links are varied in the energy comparison. The radial links below
only provide one admissible upper bound for that comparison.
"""

import numpy as np


def spring_energy(flux):
    # Four identical periodic contacts in series, K=1 per contact per face.
    return 8 * np.sin(flux / 8) ** 2


def trial_face_flux(outer, sign):
    axis_points = np.arange(-outer, outer + 1)
    vertices = np.stack(np.meshgrid(*([axis_points] * 3), indexing="ij"), axis=-1)
    radius = np.linalg.norm(vertices, axis=-1, keepdims=True)
    normal = sign * vertices / np.maximum(radius, 1)
    normal[outer, outer, outer] = [0, 0, 1]  # unused interior vertex
    x, y, z = np.moveaxis(normal, -1, 0)
    azimuth = np.exp(1j * np.arctan2(y, x))
    state = np.stack((np.sqrt((1 + z) / 2), azimuth * np.sqrt((1 - z) / 2)), axis=-1)
    flux = []
    for axis in range(3):
        i, j = (axis + 1) % 3, (axis + 2) % 3
        corners = []
        for di, dj in ((0, 0), (1, 0), (1, 1), (0, 1)):
            slices = [slice(None)] * 3
            slices[i] = slice(di, di + 2 * outer)
            slices[j] = slice(dj, dj + 2 * outer)
            corners.append(state[tuple(slices)])
        loop = np.ones(corners[0].shape[:-1], dtype=complex)
        for a, b in zip(corners, corners[1:] + corners[:1]):
            loop *= np.sum(a.conj() * b, axis=-1)
        flux.append(np.angle(loop))
    return flux


def shell_count(flux, outer, cut):
    total = 0.
    for axis, field in enumerate(flux):
        lower = [slice(outer-cut, outer+cut)] * 3
        upper = lower.copy()
        lower[axis], upper[axis] = outer-cut, outer+cut
        total += field[tuple(upper)].sum() - field[tuple(lower)].sum()
    return total / (2 * np.pi)


def trial(radius, sign):
    outer = 2 * radius
    fields = trial_face_flux(outer, sign)
    energy = 0.
    max_ratio = 0.
    for axis, field in enumerate(fields):
        coords = np.indices(field.shape, dtype=float)
        for i in range(3):
            coords[i] += -outer if i == axis else 0.5 - outer
        shell = np.floor(np.abs(coords).max(axis=0))
        use = shell >= radius
        energy += np.sum(spring_energy(field[use]))
        max_ratio = max(max_ratio, np.max(abs(field[use]) * 2 * shell[use] ** 2))
    net = sum(np.diff(field, axis=i) for i, field in enumerate(fields))
    cells = (abs(np.indices(net.shape) + .5 - outer).max(axis=0) >= radius)
    max_defect = abs(net[cells]).max() / (2 * np.pi)
    counts = [shell_count(fields, outer, r) for r in (radius, radius+1, outer)]
    return energy, max_ratio, max_defect, counts


def main():
    attachment = 4.0
    # All following thresholds include both faces with exact counterrotation.
    cube_floor = 12 * spring_energy(np.pi / 3)
    branch_floor = 2 * spring_energy(np.pi)
    print(f"K=1 G={attachment:g} unit_cube_floor={cube_floor:.12f} "
          f"branch_tie_floor={branch_floor:.12f}")
    for radius in (5, 10):
        upper = 2 * (72 + 48 / radius + 12 / radius**2) / (32 * (radius - 1))
        fmax = 1 / (2 * radius**2)
        # The unnormalized nearest-site overlap multiplies the link energy.
        rim_floor = (2 * attachment * np.cos(1 / (2 * radius))
                     * (1 - np.cos((np.pi - fmax) / 4)))
        print(f"R={radius} exterior_trial_upper={upper:.12f} "
              f"wrong_rim_count_floor={rim_floor:.12f}")
        energies = []
        for sign in (+1, -1):
            energy, ratio, defect, counts = trial(radius, sign)
            energies.append(energy)
            print(f"R={radius} face_sign={sign:+d} trial_one_face={energy:.12f} "
                  f"counts_inner_next_outer={[round(q, 12) for q in counts]} "
                  f"max_local_flux_over_bound={ratio:.9f} max_cube_count_abs={defect:.3g}")
            assert np.allclose(counts, sign, rtol=0, atol=1e-12)
            assert ratio <= 1 + 1e-12 and defect < 1e-12
        assert sum(energies) < upper < min(rim_floor, cube_floor, branch_floor)
    print("G=0: flat trial, both compact energies and both mouth counts are zero.")


if __name__ == "__main__":
    main()
