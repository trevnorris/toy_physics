"""Static compact torsion on one square and on a closed oriented cubical mouth."""

from itertools import product
from math import asin, cos, pi, sin, sqrt

import numpy as np


def principal(angle):
    return np.angle(np.exp(1j * angle))


def ring_energy(flux):
    # Each of four independent rotatable attachments has energy 1 - cos(strain).
    return 4 * (1 - np.cos(principal(flux) / 4))


def stationary_square(flux):
    # Equal transmitted torque implies each bond has angle a or pi-a.
    # k counts the latter and a lies in [-pi/2, pi/2].
    candidates = []
    for k in range(5):
        multiplicity = 4 - 2 * k
        if multiplicity == 0:
            if abs(principal(flux)) < 1e-12:
                candidates.append((k, None, 4.0))
            continue
        for winding in range(-5, 6):
            angle = (2 * pi * winding - flux - k * pi) / multiplicity
            if abs(angle) <= pi / 2 + 1e-13:
                candidates.append((k, angle, 4 - multiplicity * cos(angle)))
    return candidates


def check_square():
    samples = [0, pi / 3, pi, 1.9 * pi, 2 * pi]
    for flux in samples:
        stationary = stationary_square(flux)
        minimum = min(energy for _, _, energy in stationary)
        print(f"square flux/pi={flux/pi:.6f} stationary_min={minimum:.9f}"
              f" compact_energy={ring_energy(flux):.9f}"
              f" marked_harmonic_energy={flux**2/8:.9f}")
    sweep = np.linspace(-2*pi, 2*pi, 513)
    largest = max(abs(min(e for _, _, e in stationary_square(f)) - ring_energy(f))
                  for f in sweep)
    print(f"all-stationary sweep discrepancy={largest:.3g}"
          f" equilibrium_cut_torque={sin(pi/4):.9f}")
    assert largest < 1e-12


def cube():
    vertices = np.array(list(product((-1., 1.), repeat=3)))
    faces = []
    for axis, sign in product(range(3), (-1, 1)):
        j, k = (axis+1) % 3, (axis+2) % 3
        face = [i for i, v in enumerate(vertices) if v[axis] == sign]
        face.sort(key=lambda i: np.arctan2(sign * vertices[i, k], vertices[i, j]))
        faces.append(face)
    edges = sorted({tuple(sorted((i, j))) for face in faces
                    for i, j in zip(face, face[1:] + face[:1])})
    incidence = np.zeros((len(faces), len(edges)))
    for p, face in enumerate(faces):
        for i, j in zip(face, face[1:] + face[:1]):
            incidence[p, edges.index(tuple(sorted((i, j))))] = 1 if i < j else -1
    assert np.all(incidence.sum(axis=0) == 0)
    return vertices / sqrt(3), edges, incidence


def overlaps(normals, edges):
    z = []
    for x, y, n in normals:
        # North chart is regular at all eight vertices.
        z.append(np.array([1+n, x+1j*y]) / sqrt(2*(1+n)))
    return np.array([np.vdot(z[i], z[j]) for i, j in edges])


def oriented_counts(incidence, top, bottom):
    curvatures = [principal(incidence @ links) for links in (top, bottom, top-bottom)]
    counts = [sum(curvature)/(2*pi) for curvature in curvatures]
    assert np.max(np.abs(np.array(counts) - np.round(counts))) < 1e-11
    return [round(value) for value in counts], curvatures


def mouth_sector():
    normals, edges, incidence = cube()
    for sign in (1, -1, 0):
        texture = sign * normals if sign else np.tile([1., 0., 0.], (8, 1))
        top = np.angle(overlaps(texture, edges))
        bottom = np.angle(overlaps(-texture, edges))
        counts, curvature = oriented_counts(incidence, top, bottom)
        print(f"texture={sign:+d} counts_top_bottom_difference={counts}"
              f" first_plaquette/pi={np.round([p[0]/pi for p in curvature], 8)}"
              f" total_torsion={sum(ring_energy(curvature[2])):.9f}")
        assert counts == [sign, -sign, 2*sign]
    # Keeping the mouth order while introducing untied, new mechanical rotors.
    unconstrained = np.zeros(len(edges))
    _, curvature = oriented_counts(incidence, unconstrained, unconstrained)
    print(f"independent_flat_rotors total_torsion={sum(ring_energy(curvature[2])):.9f}")

    # The ideal configuration above is a variational trial. The surface contains
    # only the local site, link and differential square terms, all nonnegative.
    tilt, link_error = 1/30, 1/12
    site_gain, link_gain = 10000, 2000
    ideal_overlap = min(abs(overlaps(normals, edges)))
    # In the site-wise shortest spinor gauge, ||z'-z|| <= 2 sin(tilt/4).
    complex_shift = 4*sin(tilt/4)
    link_shift = asin(complex_shift / ideal_overlap) + link_error
    worst_diff_phase = 2*pi/3 + 8 * link_shift
    site_floor = site_gain * (1-cos(tilt))
    link_floor = link_gain * (ideal_overlap-complex_shift) * (1-cos(link_error))
    trial = 6*ring_energy(2*pi/3)
    print(f"finite attachment J={site_gain} K={link_gain} tangent_spin_tilt={tilt:.9f}"
          f" attachment_deviation={link_error:.9f} trial={trial:.9f}"
          f" site_floor={site_floor:.9f} bond_floor={link_floor:.9f}"
          f" face_link_error_bound={link_shift:.9f} worst_diff_face/pi={worst_diff_phase/pi:.9f}")
    assert trial < min(site_floor, link_floor) and worst_diff_phase < pi

    for delta in (0.1, 0.01, 0.001):
        # Any six cube phases with sum 2pi have link representatives. Subtract
        # a full turn before solving the exact, zero-sum real incidence problem.
        primitive = np.array([pi-delta, pi-delta] + [delta/2]*4)
        lifted = primitive.copy()
        lifted[0] -= 2*pi
        top = np.linalg.lstsq(incidence, lifted, rcond=None)[0]
        counts, curvature = oriented_counts(incidence, top, -top)
        mismatch = np.max(abs(curvature[0] - primitive))
        print(f"branch-control delta={delta:g} counts={counts}"
              f" torsion={sum(ring_energy(curvature[2])):.12g}"
              f" primitive_gauge_action={sum(1-np.cos(curvature[0])):.9f}"
              f" reconstruction={mismatch:.2g}")
        assert counts == [1, -1, 0] and mismatch < 1e-12


def torque_frames():
    # Signed angle of two marked directions in the same material bearing plane.
    axis, reference, mark = np.eye(3)
    mark = cos(0.7)*reference + sin(0.7)*mark
    direction = np.array([2., 1., -3.]) / sqrt(14)
    generator = np.array([np.cross(direction, vector) for vector in np.eye(3)]).T
    rotation = np.eye(3) + sin(1.3)*generator + (1-cos(1.3))*(generator @ generator)

    def material_angle(n, u, v):
        return np.arctan2(n @ np.cross(u, v), u @ v)

    original = material_angle(axis, reference, mark)
    rotated = material_angle(rotation @ axis, rotation @ reference, rotation @ mark)
    # Equal scalar torques on perpendicular shafts need a vector reaction.
    torque_first, torque_second = np.eye(3)[:2]
    bearing = -(torque_first-torque_second)
    print(f"material_mark_angle={original:.9f} rigidly_rotated={rotated:.9f}"
          f" unit_right_angle_bearing_torque_norm={np.linalg.norm(bearing):.9f}")
    assert abs(original-rotated) < 1e-12


if __name__ == "__main__":
    check_square()
    mouth_sector()
    torque_frames()
