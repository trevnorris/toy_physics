"""Interval energy comparisons for seven shells of prepared orientations.

The gauge uses the shared meridian of a given stream line, so radial overlap
phases are zero. K=1, the outer vacuum coefficient is C=1, its radius is 16,
and the inner angular attachment has strength 64. All decimal inputs are exact. mpmath.iv encloses the angle,
energy and attachment calculations, including mesh construction and rotation.
"""
from itertools import product
from math import floor

from mpmath import iv


iv.dps = 27
R = iv.mpf
ZERO = R(0)
ONE = R(1)
PI = iv.pi
TAU = 2 * PI


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return tuple(a[(i + 1) % 3] * b[(i + 2) % 3]
                 - a[(i + 2) % 3] * b[(i + 1) % 3] for i in range(3))


def normalized(v):
    size = iv.sqrt(sum(x ** 2 for x in v))
    return tuple(x / size for x in v)


def mesh(level):
    vertices = [tuple(R(s * int(i == j)) for j in range(3))
                for s in (1, -1) for i in range(3)]
    faces = []
    for signs in product((1, -1), repeat=3):
        i, j, k = [axis + 3 * int(s < 0) for axis, s in enumerate(signs)]
        if signs[0] * signs[1] * signs[2] < 0:
            k, j = j, k
        faces.append((i, j, k))
    for _ in range(level):
        cache = {}

        def midpoint(i, j):
            e = tuple(sorted((i, j)))
            if e not in cache:
                cache[e] = len(vertices)
                vertices.append(normalized([x + y for x, y in zip(vertices[i], vertices[j])]))
            return cache[e]

        children = []
        for i, j, k in faces:
            a, b, c = midpoint(i, j), midpoint(j, k), midpoint(k, i)
            children.extend(((i, a, c), (a, j, b), (c, b, k), (a, b, c)))
        faces = children
    return vertices, faces


def rotated(vertices, rotation):
    rotation = tuple(map(R, rotation))
    angle = iv.sqrt(dot(rotation, rotation))
    axis = normalized(rotation)
    cosine, sine = iv.cos(angle), iv.sin(angle)
    directions = []
    for v in vertices:
        along = dot(v, axis) * (1 - cosine)
        across = cross(axis, v)
        directions.append(tuple(cosine * x + sine * z + along * y
                                for x, y, z in zip(v, axis, across)))
    return directions


def ray(direction, retained, k, eta):
    if k == 0:
        return (ONE, ZERO, ZERO), ONE
    nx, ny, nz = (eta * x for x in direction)
    t = retained ** k
    a = iv.sqrt(t)
    b = (1 - a) * (1 + a - a * nz)
    vector = (b * nx, b * ny, a + b * nz)
    length = iv.sqrt(sum(c ** 2 for c in vector))
    x, y, z = (c / length for c in vector)
    # The mesh is rotated, so its vertices avoid both poles.
    h = iv.sqrt((1 + z) / 2)
    return (h, x / (2 * h), y / (2 * h)), length


def bond(u, v):
    h, x, y = u
    H, X, Y = v
    re, im = h * H + x * X + y * Y, x * Y - y * X
    modulus = iv.sqrt(re ** 2 + im ** 2)
    assert modulus > 0
    return iv.atan2(im, re), modulus


def principal(angle):
    n = (angle + PI) / TAU
    low, high = floor(float(n.a)), floor(float(n.b))
    assert low == high
    result = angle - low * TAU
    assert result > -PI and result < PI
    return result


def spring(flux, contacts):
    return 2 * contacts * iv.sin(flux / (2 * contacts)) ** 2


def min_bound(xs):
    # Choose an endpoint whose value is provably at or below the true minimum.
    return min(x.a for x in xs)


def rounded_bound(x, digits, upper=False):
    # int converts a point interval to its truncated integer. Directed outputs
    # therefore do not depend on binary-to-decimal rounding of a float.
    scale = 10 ** digits
    y = (x * scale).b if upper else (x * scale).a
    return (int(y) + int(upper)) / scale


def check(level, rotation, eta):
    vertices, faces = mesh(level)
    normals = rotated(vertices, rotation)
    edges = sorted({tuple(sorted((i, j))) for f in faces for i, j in zip(f, f[1:] + f[:1])})
    retained = iv.cos(iv.sqrt(R(2)) * R(".33")) ** 2
    shells = []
    side_energy = ZERO
    all_bonds = []
    for k in range(7):
        prepared = [ray(n, retained, k, eta) for n in normals]
        bonds = {}
        for i, j in edges:
            u, p = prepared[i]
            v, q = prepared[j]
            phase, overlap = bond(u, v)
            bonds[i, j] = (phase, iv.sqrt(p * q) * overlap)
        angles = []
        slips = []
        for f in faces:
            loop = [(tuple(sorted((i, j))), 1 if i < j else -1)
                    for i, j in zip(f, f[1:] + f[:1])]
            angle = principal(sum(s * bonds[e][0] for e, s in loop))
            angles.append(angle)
            # A change of the cut integer implies a branch crossing on some
            # triangle along the shortest-edge interpolation. At a crossing
            # sum |edge slips| >= pi - |preferred face|. Restrict their sum
            # to this gap (<=pi) by monotonicity. Pairwise equalizing improves
            # the bound since cos(x)+cos(y)=2cos((x+y)/2)cos((x-y)/2).
            gap = PI - abs(angle)
            slips.append(6 * min_bound(bonds[e][1] for e, _ in loop) * iv.sin(gap / 6) ** 2)
        winding = sum(angles) / TAU
        degree = round(float(winding.mid))
        assert winding.a <= degree <= winding.b and winding.delta < R("1e-16")
        if k:
            side_energy += sum(spring(principal(bonds[e][0] - previous[e][0]), 4)
                               for e in edges)
        all_bonds.append(bonds)
        shells.append((degree, sum(spring(a, 4) for a in angles), min_bound(slips)))
        previous = bonds
    trial = sum(s[1] for s in shells) + side_energy
    # Carry the innermost angular connection through every cut. In this
    # meridional gauge preferred radial phases and all trial radial phases
    # vanish. The radial vacuum test field F_i/r^2 is interpreted
    # using the solid angle of its triangle. Its flux per solid angle is
    # held constant with radius, and its divergence away from the cut is zero.
    areas = [2 * iv.atan2(dot(normals[i], cross(normals[j], normals[k])),
                          1 + dot(normals[i], normals[j]) + dot(normals[i], normals[k])
                          + dot(normals[j], normals[k])) for i, j, k in faces]
    tail = sum(f ** 2 / area for f, area in zip(angles, areas)) / 32
    mismatch = sum(2 * weight * iv.sin((all_bonds[-1][e][0] - phase) / 2) ** 2
                   for previous in all_bonds[:-1] for e, (phase, weight) in previous.items())
    carried_base = 7 * shells[-1][1] + tail
    through = carried_base + R(".0001") * mismatch
    charged = 5 * spring(TAU / 5, 4)
    tie = spring(PI, 4)
    # Convex W on [-pi,pi] gives sum_face W >= 5W(2pi/5) on a
    # five-face prism with nonzero index. A face at a branch tie costs W(pi).
    # In both regimes the inner angular links are unconstrained at strength 64.
    forbidden_carried = min_bound([charged, tie, 64 * shells[-1][2]])
    weak_upper = (forbidden_carried - carried_base) / mismatch
    strong_lower = trial / shells[0][2]
    degrees = [s[0] for s in shells]
    floors = [rounded_bound(s[2], 6) for s in shells]
    print(f"level={level} rotation={','.join(rotation)} eta={eta:+d} "
          f"degrees={degrees} preferred_trial_upper={rounded_bound(trial, 6, True):.6f} "
          f"wrong_cut_lower_per_G={floors} carried_trial_upper={rounded_bound(through, 6, True):.6f} "
          f"rim_field_per_cut_upper={rounded_bound(shells[-1][1], 8, True):.8f} "
          f"radial_vacuum_upper={rounded_bound(tail, 8, True):.8f} "
          f"carried_attachment_per_feed_upper={rounded_bound(mismatch, 7, True):.7f} "
          f"charged_prism_lower={rounded_bound(charged, 6):.6f} "
          f"tie_lower={rounded_bound(tie, 6):.6f} "
          f"transmission_for_g_below={rounded_bound(weak_upper, 6):.6f} "
          f"screening_for_g_above={rounded_bound(strong_lower, 6, True):.6f}", flush=True)
    assert trial < 8 * shells[0][2] and trial < 64 * shells[-1][2]
    assert through < forbidden_carried and weak_upper > R(".015")
    assert strong_lower < R("1.76")


if __name__ == "__main__":
    for level in (2, 3):
        for rotation in ((".26", "-.19", ".17"), ("-.13", ".24", ".31")):
            for eta in (1, -1):
                check(level, rotation, eta)
