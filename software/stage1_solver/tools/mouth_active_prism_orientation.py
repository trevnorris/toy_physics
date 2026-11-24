"""Certify the registered R4 prism chart where the spatial absolute Q1 order is positive.

The determinant of P1(triangle) x P1(angle) is affine in barycentric position
and quadratic in angular fraction. Positive Bernstein coefficients certify the
whole wedge. On the remaining wedges the binary input numbers are interpreted
as exact rationals. A clipped triangular enclosure covers every possible
nonpositive determinant at every angle in a subprism. Nonnegative Q1 basis
functions identify the spatial nodes whose open support intersects its box.
Unseparated cases are subdivided in material coordinates.
"""
from bisect import bisect_left, bisect_right
from fractions import Fraction
from itertools import product

import numpy as np


def _interval(value):
    return value, value


def _down(value):
    return np.nextafter(value, -np.inf)


def _up(value):
    return np.nextafter(value, np.inf)


def _add(a, b):
    return _down(a[0] + b[0]), _up(a[1] + b[1])


def _sub(a, b):
    return _down(a[0] - b[1]), _up(a[1] - b[0])


def _mul(a, b):
    products = [x*y for x in a for y in b]
    return _down(np.minimum.reduce(products)), _up(np.maximum.reduce(products))


def _component(v, k):
    return v[0][..., k], v[1][..., k]


def _cross(a, b):
    return [_sub(_mul(_component(a, i), _component(b, j)),
                 _mul(_component(a, j), _component(b, i)))
            for i, j in ((1, 2), (2, 0), (0, 1))]


def _dot_vertices(cross, delta):
    terms = [_mul((lo[..., None], hi[..., None]), _component(delta, k))
             for k, (lo, hi) in enumerate(cross)]
    return _add(_add(terms[0], terms[1]), terms[2])


def _easy_wedges(bottom, top):
    """Outward-rounded Bernstein lower bounds on all nine coefficients."""
    edges = [[_sub(tuple(v[..., i, :] for v in plane),
                   tuple(v[..., 0, :] for v in plane)) for i in (1, 2)]
             for plane in (bottom, top)]
    (x0, t0), (x1, t1) = edges
    da = _sub(top, bottom)
    lo = _dot_vertices(_cross(x0, t0), da)[0]
    hi = _dot_vertices(_cross(x1, t1), da)[0]
    mixed = [_add(a, b) for a, b in zip(_cross(x0, t1), _cross(x1, t0))]
    mid_twice = _dot_vertices(mixed, da)[0]
    return np.all((lo > 0) & (hi > 0) & (mid_twice > 0), axis=-1)


def _radial_bound(body):
    """Enclose the native R4 piecewise cubic plus quartic meridian on each knot interval."""
    c = body.buffers['spline']
    e, m, a, extra = body.buffers['shape'][:4]
    n = len(c)
    knots = np.arange(n + 1) / n
    t = (np.maximum(0., _down(knots[:-1])), np.minimum(1., _up(knots[1:])))
    h = (0., _up(1./n))
    base = _interval(c[:, 0])
    for column in c[:, 1:].T:
        base = _add(_mul(base, h), _interval(column))
    prod = _mul(t, _sub(_interval(1.), t))
    odd = _sub(_mul(_interval(2.), t), _interval(1.))
    mid = _mul(_mul(_interval(4.), prod), _add(_interval(m), _mul(_interval(a), odd)))
    quart = _mul(_mul(_interval(16.), _interval(extra)), _mul(prod, prod))
    profile = _mul(_interval(4.), _add(_add(base, _mul(_interval(e), t)), _add(mid, quart)))
    # 0 <= bump <= 1, 0 <= x^2 <= 1, -1/2 <= P2(cos(theta)) <= 1.
    polar = _mul(_mul(_interval(3.), _interval(body.epsilon)), (-.5, 1.))
    return float(_add(profile, (min(0., polar[0]), max(0., polar[1])))[0].min())


def _rational(x):
    return Fraction.from_float(float(x))


def _exact_array(x):
    return np.vectorize(_rational, otypes=[object])(x)


def _det(a, b, z):
    return (a[0]*(b[1]*z[2] - b[2]*z[1]) + a[1]*(b[2]*z[0] - b[0]*z[2])
            + a[2]*(b[0]*z[1] - b[1]*z[0]))


def _coefficients(points):
    (x, t), (v, w) = [plane[1:] - plane[0] for plane in points]
    return np.array([[_det(x, t, d), (_det(x, w, d) + _det(v, t, d))/2, _det(v, w, d)]
                     for d in points[1] - points[0]], dtype=object)


def _quadratic_min(row):
    lo, middle, hi = row
    a, b = lo - 2*middle + hi, 2*(middle-lo)
    if a > 0 and 0 < -b < 2*a:
        u = -b/(2*a)
        return lo + u*(b+a*u), u
    return (lo, 0) if lo <= hi else (hi, 1)


UNIT = np.array([[Fraction(int(i == j)) for j in range(3)] for i in range(3)], dtype=object)
HALF = Fraction(1, 2)


def _possible_fold(minima):
    """Vertices of the triangle clipped by sum_i lambda_i min_angle D_i <= 0."""
    bary = [UNIT[i] for i in range(3) if minima[i] <= 0]
    for i, j in ((0, 1), (1, 2), (2, 0)):
        if minima[i]*minima[j] < 0:
            α = minima[i]/(minima[i]-minima[j])
            bary.append(UNIT[i]+α*(UNIT[j]-UNIT[i]))
    return np.array(bary, dtype=object)


class _Support:
    def __init__(self, body, order):
        self.axes = [tuple(map(_rational, body.buffers['field_'+k])) for k in ('x', 't', 'a')]
        self.ids = np.asarray(body.ids).transpose(0, 2, 1)
        self.values = np.asarray(order)[self.ids]
        if self.ids.shape != tuple(map(len, self.axes)) or any(
                len(v) < 2 or any(a >= b for a, b in zip(v, v[1:])) for v in self.axes):
            raise ValueError('Malformed native order axes')
        if not np.all(np.isfinite(self.values) & (0 <= self.values) & (self.values <= 1)):
            raise ValueError('Absolute order outside [0,1]')
        self.positive = self.values > 0

    def box(self, points):
        bounds = zip(points.min(axis=(0, 1)), points.max(axis=(0, 1)), self.axes)
        return tuple(slice(bisect_right(axis, lo)-1, bisect_left(axis, hi)+1)
                     for lo, hi, axis in bounds)

    def at(self, point):
        basis = []
        for axis, x in zip(self.axes, point):
            i = min(len(axis)-2, bisect_right(axis, x)-1)
            f = (x-axis[i])/(axis[i+1]-axis[i])
            basis.append(((i, 1-f), (i+1, f)))
        value = Fraction(0)
        for (i, wi), (j, wj), (k, wk) in product(*basis):
            if wi and wj and wk:
                value += wi*wj*wk*_rational(self.values[i, j, k])
        return value


def _split(points, coeff, reference, support, box):
    # Refine in the direction of spatial order variation. Where the support
    # enclosure has only positive coefficients, determinant variation is sufficient.
    changing = [k for k in range(3) if np.any(np.diff(support.positive[box], axis=k))]
    def score(d):
        return max((float(np.max(abs(d[..., k])) /
                          (support.axes[k][box[k].stop-1]-support.axes[k][box[k].start]))
                    for k in changing), default=0.)
    edges = ((0, 1), (1, 2), (2, 0))
    scores = [score(points[:, i]-points[:, j]) for i, j in edges] + [score(points[1]-points[0])]
    # Break ties using the longest reference edge in units of the original prism.
    sizes = [np.max(abs(reference[:, i]-reference[:, j])) for i, j in edges]
    sizes.append(np.max(abs(reference[1]-reference[0])))
    axis = max(range(4), key=lambda i: (scores[i], sizes[i]))
    if axis == 3:
        mid, refmid = (points[0]+points[1])/2, (reference[0]+reference[1])/2
        ab, bc = (coeff[:, 0]+coeff[:, 1])/2, (coeff[:, 1]+coeff[:, 2])/2
        center = (ab+bc)/2
        return [(np.stack((points[0], mid)), np.stack((coeff[:, 0], ab, center), axis=1),
                 np.stack((reference[0], refmid))),
                (np.stack((mid, points[1])), np.stack((center, bc, coeff[:, 2]), axis=1),
                 np.stack((refmid, reference[1])))]
    i, j = edges[axis]
    children = []
    for k in (i, j):
        p, c, r = points.copy(), coeff.copy(), reference.copy()
        p[:, k], c[k], r[:, k] = (points[:, i]+points[:, j])/2, (coeff[i]+coeff[j])/2, (reference[:, i]+reference[:, j])/2
        children.append((p, c, r))
    return children


def _certify_wedge(points, scale, support, activation, counts, max_depth):
    coeff = _coefficients(points)
    ref = np.array([[list(b)+[Fraction(y)] for b in UNIT] for y in (0, 1)], dtype=object)
    stack = [(points, coeff, ref, 0)]
    uncertain = None
    while stack:
        p, c, r, depth = stack.pop()
        counts['exact_subprisms'] += 1
        counts['max_depth_used'] = max(depth, counts['max_depth_used'])
        vertices = [_quadratic_min(row) for row in c]
        mins = [v for v, _ in vertices]
        if min(mins) > 0:
            continue
        possible = _possible_fold(mins)
        box = support.box(possible @ p)
        if not support.positive[box].any():
            counts['inactive_fold_enclosures'] += 1
            if activation is not None:
                activation.update(support.ids[box].ravel().tolist())
            continue

        # Test rational witnesses, including the exact angular minima, and
        # points inside the potentially folded polygon for an interior overlap.
        candidates = [(UNIT[i], u) for i, (v, u) in enumerate(vertices) if v <= 0]
        center = possible.sum(axis=0)/len(possible)
        for b in [center, *possible]:
            row = b @ c
            v, angle = _quadratic_min(row)
            if v <= 0:
                candidates.append((b, angle))
                for y in (0, HALF, 1):
                    if (1-y)**2*row[0]+2*y*(1-y)*row[1]+y*y*row[2] <= 0:
                        candidates.append((b, y))
        for b, y in candidates:
            position = b @ ((1-y)*p[0]+y*p[1])
            chi = support.at(position)
            if chi > 0:
                at = b @ ((1-y)*r[0]+y*r[1])
                det = (b @ c) @ np.array([(1-y)**2, 2*y*(1-y), y*y])/scale
                return 'active_fold', dict(barycentric=list(map(float, at[:3])), angular_fraction=float(at[3]),
                    placement_x_t_angle=list(map(float, position)), chart_determinant=float(det),
                    exact_determinant_sign=int(det < 0)*-1, order=float(chi), exact_order_positive=True)
        if depth == max_depth:
            if uncertain is None:
                uncertain = dict(reference_vertices=np.asarray(r, dtype=float).tolist(),
                                 determinant_lower_bound=float(min(mins)/scale),
                                 possible_spatial_node_ranges=[[s.start, s.stop-1] for s in box])
            continue
        for child in _split(p, c, r, support, box):
            stack.append((*child, depth+1))
    return ('unresolved', uncertain) if uncertain is not None else ('certified', None)


def active_orientation(body, solid, order, *, max_depth=12, first_failure=True, collect_activation=True):
    """Check one unpacked placement and the corresponding full carrier order vector.

    The carrier IDs in potential_activation_nodes cover the Q1 basis functions
    that could turn on in nonpositive-J material. A certified result with
    collect_activation=True gives a complete (possibly conservative) list.
    The signed chart determinant has the sign of physical J on the registered
    interior with the separately checked radial multiplier and positive depth.
    Symmetry faces use one-sided prism determinants. Bounds use the native
    normalized spatial chart, including every angular interval.
    """
    if max_depth < 0 or np.ndim(order) != 1 or np.any(body.buffers['base'] != 0):
        raise ValueError('Expected absolute order, a flat carrier vector and nonnegative subdivision depth')
    axes = [np.asarray(v, dtype=float) for v in (body.x, body.t, body.a)]
    if any(len(v) < 2 or not np.all(np.isfinite(v)) or np.any(np.diff(v) <= 0) for v in axes):
        raise ValueError('Malformed material axes')
    grid = np.stack(np.broadcast_arrays(axes[0][None,:,None], axes[1][:,None,None], axes[2][None,None,:]), axis=-1)
    solid = np.asarray(solid, dtype=float)
    if solid.shape != grid.shape or not np.all(np.isfinite(solid)):
        raise ValueError('Malformed unpacked material displacement')
    support = _Support(body, order)
    if any((v[0], v[-1]) != (f[0], f[-1]) for v, f in zip(axes, support.axes)):
        raise ValueError('Mismatched material and spatial chart domain')
    if [v[0] for v in axes] != [0., 0., 0.] or [v[-1] for v in axes] != [1., 1., np.pi/2]:
        raise ValueError('Expected the registered R4 chart')
    coords = grid + solid
    counts = dict(total_prisms=2*int(np.prod(np.array(grid.shape[:3])-1)), positive_bernstein_prisms=0,
                  exact_prisms=0, exact_subprisms=0, inactive_fold_enclosures=0, max_depth_used=0,
                  active_fold_prisms=0, unresolved_prisms=0)
    activation = set() if collect_activation else None
    radial = _radial_bound(body)
    def result(status, witness):
        return dict(status=status, all_active_cells_oriented=status == 'certified', witness=witness,
                    physical_radial_chart_derivative_lower_bound=radial, counts=counts,
                    potential_activation_nodes=sorted(activation) if activation is not None else None)
    if not np.isfinite(radial) or not np.isfinite(body.shape.length) or min(radial, body.shape.length) <= 0:
        return result('unresolved', dict(reason='Positive physical chart multiplier unverified'))
    for k, axis in enumerate(axes):
        x = coords[..., k]
        outside = (x < axis[0]) | (x > axis[-1])
        # Exact boundary ties in the rounded sum need the lost low bits.
        for idx in np.argwhere(((x == axis[0]) | (x == axis[-1])) & (solid[..., k] != 0)):
            key = (*idx, k)
            value = _rational(grid[key]) + _rational(solid[key])
            if not support.axes[k][0] <= value <= support.axes[k][-1]:
                outside[tuple(idx)] = True
        if outside.any():
            idx = tuple(map(int, np.argwhere(outside)[0]))
            return result('outside_chart', dict(node_t_x_angle=list(idx), component=k, value=float(x[idx])))
    failure, witness = 'certified', None
    for triangle, offsets in enumerate((((0, 0), (0, 1), (1, 0)), ((1, 1), (1, 0), (0, 1)))):
        # offsets here are (t,x); the two triangle derivatives both change
        # sign on the upper triangle, leaving the determinant sign unchanged.
        def corners(values):
            return np.stack([values[j:j+len(body.t)-1, i:i+len(body.x)-1] for j, i in offsets], axis=-2)
        v = corners(coords)
        lo, hi = _down(v), _up(v)
        easy = _easy_wedges((lo[:,:,:-1], hi[:,:,:-1]), (lo[:,:,1:], hi[:,:,1:]))
        counts['positive_bernstein_prisms'] += int(easy.sum())
        for j, i, k in np.argwhere(~easy):
            def planes(values):
                return np.stack([[values[j+b, i+d, k+y] for b, d in offsets] for y in (0, 1)])
            p = _exact_array(planes(grid)) + _exact_array(planes(solid))
            scale = Fraction(1)
            for idx, axis in zip((i, j, k), axes):
                scale *= _rational(axis[idx+1])-_rational(axis[idx])
            counts['exact_prisms'] += 1
            status, item = _certify_wedge(p, scale, support, activation, counts, max_depth)
            if status != 'certified':
                counts[status+'_prisms'] += 1
                if failure == 'certified' or (status == 'active_fold' and failure != status):
                    failure = status
                    witness = dict(cell_t_x_angle=list(map(int, (j, i, k))), triangle=triangle, **item)
                if first_failure:
                    return result(failure, witness)
    return result(failure, witness)
