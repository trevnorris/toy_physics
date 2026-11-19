"""Directed material connections on a closed conversion section, with free site and bond states."""

from itertools import product

import numpy as np
import sympy as sp


def unit(vector):
    return vector / np.linalg.norm(vector)


def section(level):
    """Subdivide an outward-oriented octahedron on a unit section."""
    points = list(np.concatenate((np.eye(3), -np.eye(3))))
    faces = []
    for signs in product((0, 3), repeat=3):
        i, j, k = [axis + shift for axis, shift in enumerate(signs)]
        if np.dot(points[i], np.cross(points[j], points[k])) < 0:
            j, k = k, j
        faces.append((i, j, k))
    for _ in range(level):
        midpoints = {}

        def midpoint(i, j):
            edge = tuple(sorted((i, j)))
            if edge not in midpoints:
                midpoints[edge] = len(points)
                points.append(unit(points[i] + points[j]))
            return midpoints[edge]

        children = []
        for i, j, k in faces:
            a, b, c = midpoint(i, j), midpoint(j, k), midpoint(k, i)
            children.extend(((i, a, c), (a, j, b), (c, b, k), (a, b, c)))
        faces = children
    return np.array(points), np.array(faces)


def state(direction):
    x, y, z = direction
    operator = np.array([[z, x - 1j*y], [x + 1j*y, -z]])
    return np.linalg.eigh(operator)[1][:, -1]


def preferred_bonds(states, faces):
    edges = {}
    for i, j, k in faces:
        for a, b in ((i, j), (j, k), (k, i)):
            edge = tuple(sorted((a, b)))
            if edge not in edges:
                overlap = np.vdot(states[edge[0]], states[edge[1]])
                if abs(overlap) < 1e-12:
                    raise ValueError("Neighbor order states are orthogonal")
                edges[edge] = overlap / abs(overlap)
    return edges


def fluxes(bonds, faces):
    def oriented(i, j):
        return bonds[i, j] if i < j else np.conj(bonds[j, i])
    return np.angle([oriented(i, j)*oriented(j, k)*oriented(k, i) for i, j, k in faces])


def count(angles):
    total = sum(angles) / (2*np.pi)
    assert np.isclose(total, round(total), rtol=0, atol=1e-11)
    return round(total)


def spherical_areas(points, faces):
    a, b, c = points[faces].transpose(1, 0, 2)
    triple = np.einsum("ij,ij->i", a, np.cross(b, c))
    pairs = (a*b + b*c + c*a).sum(axis=1)
    return 2*np.arctan2(triple, 1 + pairs)


def mouth(radius, level):
    directions, faces = section(level)
    areas = spherical_areas(directions, faces)
    assert min(areas) > 0 and np.isclose(sum(areas), 4*np.pi)
    rng = np.random.default_rng(741)
    rephasings = np.exp(1j*rng.uniform(-np.pi, np.pi, len(directions)))
    positions = radius*directions
    length = np.linalg.norm(positions, axis=1)
    currents = 0.25 * areas / (4*np.pi)
    records = {}
    for side in (1, -1):
        # w_s(r)=s R^2/r. Use the positive-w coorientation for this first comparison.
        gradient = -side*radius**2*positions / length[:, None]**3
        normal = np.array([unit(-v) for v in gradient])
        states = np.array([state(v) for v in normal])
        bonds = preferred_bonds(states, faces)
        angles = fluxes(bonds, faces)
        winding = count(angles)
        physical_flux_density = angles / (radius**2 * areas)
        changed = preferred_bonds(states*rephasings[:, None], faces)
        assert np.max(np.abs(angles - fluxes(changed, faces))) < 1e-13
        assert np.max(np.abs(angles - side*areas/2)) < 1e-13
        print(f"R={radius} level={level} faces={len(faces)} side={side:+d}"
              f" inward={sum(currents):.9f} count={winding:+d}"
              f" mean-density={physical_flux_density.mean():+.9f}"
              f" max-loop={max(abs(angles)):.9f}")
        # The full link space admits other counts even at the same site states.
        sectors = []
        for power in (-2, -1, 0, 1, 2):
            varied = {edge: phase**power for edge, phase in bonds.items()}
            sectors.append(count(fluxes(varied, faces)))
        assert sectors == [power*side for power in (-2, -1, 0, 1, 2)]
        # Turn the slope past the edge of the wrapping domain at the same scalar intake.
        tilted_states = [state(unit(v + 2*np.array([0., 0., 1.]))) for v in normal]
        tilted = count(fluxes(preferred_bonds(tilted_states, faces), faces))
        uniform = count(fluxes(preferred_bonds([state((0,0,1))]*len(normal), faces), faces))
        conjugate = count(fluxes(preferred_bonds(states.conj(), faces), faces))
        # Each outward material normal is evaluated on its own face of w_s +/- thickness/2.
        outward_counts = {}
        for bank in (1, -1):
            outward_states = [state(bank*v) for v in normal]
            outward_counts[bank] = count(fluxes(preferred_bonds(outward_states, faces), faces))
        assert tilted == uniform == 0 and conjugate == -winding
        print(f"  accessible={sectors} uniform={uniform} tilted={tilted} conjugate={conjugate:+d}"
              f" upper/lower={outward_counts[1]:+d},{outward_counts[-1]:+d}"
              f" protruding={outward_counts[side]:+d} both-sum={sum(outward_counts.values())}")
        records[side] = (winding, outward_counts[side])
    assert records[1] == (1,1) and records[-1] == (-1,1)


def finite_locking():
    # E contains J(1-n.p) per site, K(|v|-Re[exp(-ia) v]) per bond with
    # v=<z_i|z_j>, and nonnegative exterior energy; both local terms can vary.
    cos_alpha = sp.Rational(99,100)
    sin_double = 2*cos_alpha*sp.sqrt(1-cos_alpha**2)
    min_modulus = sp.sqrt((1-sin_double)/2)
    assert 1-3*sin_double > 0 and min_modulus > sp.Rational(3,5)
    # Positive real trace(P_i P_j P_k) throughout the short-geodesic deformation
    # implies |preferred plaquette phase| < pi/2. Three link errors < pi/6
    # keep the physical plaquette below the principal branch cut under interpolation.
    site, bond = 150*(1-cos_alpha), 20*sp.Rational(3,5)*(1-sp.cos(sp.pi/6))
    floor = min(site, bond)
    print("finite locking: sin(2alpha), minimum overlap, site floor, bond floor, minimum:",
          sin_double, min_modulus, site, bond, floor)
    assert floor == sp.Rational(3,2)
    # C=1 denotes the exterior energy of one uniform unit at R=1.
    for radius in (1, 2):
        print("radius, preferred trial upper, different-sector floor, margin:",
              radius, sp.Rational(1,radius), floor, floor-sp.Rational(1,radius))
        assert floor > sp.Rational(1,radius)


def exterior_recipient():
    # Auxiliary model: a thin spherical layer at L may carry integer compensating flux.
    # Its given cost per unit includes exchange with a distant charge reservoir.
    for radius in (1, 2):
        cut = 4*radius
        for cost in (sp.Rational(1,16), sp.Rational(3,16),
                     sp.Rational(1,cut), sp.Rational(1,2)):
            face_signs = []
            for inner in (1,-1):
                def energy(outer):
                    return (sp.Rational(1,radius) - sp.Rational(1,cut)
                            + sp.Rational(outer**2,cut) + cost*abs(outer-inner))
                # For integer outer other than 0 or inner, both terms are
                # bounded below by those at inner, with at least one strict.
                candidates = [0, inner]
                winners = [v for v in candidates if energy(v) == min(map(energy,candidates))]
                face_signs.append(winners)
            print(f"screen R={radius} L={cut} cost={cost} far for +,-: {face_signs}")
            assert face_signs[1] == [-x for x in face_signs[0]]


if __name__ == "__main__":
    for size, resolution in ((1, 1), (2, 1), (2, 2)):
        mouth(size, resolution)
    finite_locking()
    exterior_recipient()
