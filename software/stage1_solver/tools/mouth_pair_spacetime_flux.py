"""Compact geometric flux of translated, prescribed radial spinor preparations.

Axes are t,x,y,z. The pair occupies a tensor product of two trial states, so
its normalized overlap is the product of their normalized overlaps. No energy
minimization, torque, propagation law or electric normalization enters here.
"""

from itertools import combinations

import numpy as np


TAU = 2 * np.pi
TOL = 1e-11
COORDS = [np.arange(-18, 19), np.arange(-8, 9), np.arange(-8, 9)]
GRID = np.stack(np.meshgrid(*COORDS, indexing="ij"), axis=-1)
START_X = np.array([-7.69, 8.31])
STEPS = np.arange(5) / 4


def end(field, axis, high):
    index = [slice(None)] * field.ndim
    index[axis] = slice(1, None) if high else slice(None, -1)
    return field[tuple(index)]


def preparation(x_path, sign, core_shift=0.):
    center = np.stack((x_path, np.full_like(x_path, .37), np.full_like(x_path, .43)), axis=-1)
    relative = GRID[None, ...] - center[:, None, None, None, :]
    radius = np.max(abs(relative), axis=-1)
    # This bounded deformation of the interior is identically zero for r_inf >= 4.
    relative[..., 0] -= core_shift * np.clip((4 - radius) / 2, 0, 1)
    normal = sign * relative / np.linalg.norm(relative, axis=-1, keepdims=True)
    x, y, z = np.moveaxis(normal, -1, 0)
    spinor = np.stack((np.sqrt((1 + z) / 2),
                       np.exp(1j * np.arctan2(y, x)) * np.sqrt((1 - z) / 2)), axis=-1)
    return spinor, radius >= 4


def overlap_links(spinor):
    links = []
    minimum = 1.
    for axis in range(4):
        v = np.sum(end(spinor, axis, False).conj() * end(spinor, axis, True), axis=-1)
        minimum = min(minimum, np.min(abs(v)))
        if np.any(abs(v) < 1e-12):
            raise ValueError(f"undefined overlap phase on axis {axis}, minimum {minimum}")
        links.append(v / abs(v))
    return links, minimum


def plaquettes(links):
    field = {}
    for a, b in combinations(range(4), 2):
        loop = (end(links[a], b, False) * end(links[b], a, True)
                * end(links[a], b, True).conj() * end(links[b], a, False).conj())
        field[a, b] = np.angle(loop)
    return field


def charges_current(field):
    def derivative(a, b, c):
        return (np.diff(field[b, c], axis=a) - np.diff(field[a, c], axis=b)
                + np.diff(field[a, b], axis=c)) / TAU

    # In these orientations partial_t q + div J = 0 on each 4-cell.
    raw = [derivative(1, 2, 3), -derivative(0, 2, 3),
           derivative(0, 1, 3), -derivative(0, 1, 2)]
    integer = [np.rint(component).astype(int) for component in raw]
    residual = max(np.max(abs(x - n)) for x, n in zip(raw, integer))
    closure = sum(np.diff(component, axis=i) for i, component in enumerate(integer))
    return integer[0], integer[1:], residual, np.max(abs(closure))


def cut_charge(field, center_x, half_width):
    # Read the six exterior faces directly, independently of the interior fill.
    origin = np.array([center_x, 0, 0]) - np.array([axis[0] for axis in COORDS])
    flux = 0.
    for i, face in enumerate((field[2, 3], -field[1, 3], field[1, 2])):
        low = [slice(None)] + [slice(v-half_width, v+half_width) for v in origin]
        high = low.copy()
        low[i+1], high[i+1] = origin[i] - half_width, origin[i] + half_width
        flux += (face[tuple(high)] - face[tuple(low)]).sum(axis=(1, 2))
    return flux / TAU


def exterior_faces(sites, a, b):
    adjacent = end(sites, a, False) & end(sites, a, True)
    return end(adjacent, b, False) & end(adjacent, b, True)


def trial(x_paths, core_shift=0.):
    states, good = zip(*(preparation(x, sign, core_shift)
                         for x, sign in zip(x_paths, (+1, -1))))
    (plus, mplus), (minus, mminus) = map(overlap_links, states)
    field = plaquettes([a * b for a, b in zip(plus, minus)])
    # Components agree modulo 2 pi. A principal-value wrap can change local
    # counts relative to a sum of the individual compact currents.
    fa, fb = plaquettes(plus), plaquettes(minus)
    wraps = sum(np.count_nonzero(abs(fa[k] + fb[k]) > np.pi) for k in field)
    margin = min(np.pi - np.max(abs(face)) for face in field.values())
    return field, good[0] & good[1], min(mplus, mminus), wraps, margin


def transport(name, directions, steps=STEPS, core_shift=0.):
    paths = START_X[:, None] + np.asarray(directions)[:, None] * steps
    field, exterior, overlap, wraps, margin = trial(paths, core_shift)
    q, current, residual, closure = charges_current(field)
    count_error = max(np.max(abs(cut_charge(field, center, r) - sign))
                      for center, sign in ((-8, 1), (8, -1)) for r in (6, 7))
    exterior_field = [np.max(abs(face[exterior_faces(exterior, *key)]))
                      for key, face in field.items() if key[0] == 0]
    coarse_x = COORDS[0][:-1] + .5
    dipole = np.sum(q * coarse_x[None, :, None, None], axis=(1, 2, 3))
    flow = np.stack([j.sum(axis=(1, 2, 3)) for j in current], axis=-1)
    # Integer current integrated on two oriented planes swept by these paths.
    planes = [int(current[0][:, i - COORDS[0][0]].sum()) for i in (-7, 8, 9)]
    boundary = max(np.max(abs(np.take(j, (0, -1), axis=i+1))) for i, j in enumerate(current))
    x_indices = [[int(index) + COORDS[0][0] for index in np.where(qq.sum(axis=(1, 2)))[0]]
                 for qq in q[[0, -1]]]
    current_steps = flow[:, 0].tolist()
    print(f"{name} end_cells_x={x_indices} dPx={dipole[-1]-dipole[0]:g} "
          f"sumJx_by_step={current_steps} absJx={abs(current[0]).sum()} "
          f"Jx_planes[-7,8,9]={planes} other_absJ={sum(abs(j).sum() for j in current[1:])}")
    print(f"  min_overlap={overlap:.8f} pair_wraps={wraps} min_branch_margin={margin:.8f} "
          f"integer_error={residual:.3g} continuity={closure} boundary_current={boundary} "
          f"cut_error_R6_R7={count_error:.3g} outer_max_abs_Fti={max(exterior_field):.9g}")
    assert np.max(abs(q.sum(axis=(1, 2, 3)))) == 0
    assert np.all(abs(q).sum(axis=(1, 2, 3)) == 2)
    assert np.array_equal(np.diff(dipole), flow[:, 0]) and boundary == closure == 0
    assert wraps == 0 and margin > 1e-5 and max(residual, count_error) < TOL
    assert all(np.count_nonzero(j) == 0 for j in current[1:])
    # Independently compare to geometric crossings of the displaced core centers.
    expected = np.floor(paths + core_shift).astype(int)
    assert np.array_equal(np.diff(expected[0] - expected[1]), flow[:, 0])
    for i, plane in enumerate((-7, 8, 9)):
        right = (expected >= plane).astype(int)
        net_right = right[0] - right[1]
        assert planes[i] == net_right[-1] - net_right[0]
    return field, exterior, current


def gauge_control():
    state, _ = preparation(np.full_like(STEPS, START_X[0]), +1)
    links, _ = overlap_links(state)
    base = plaquettes(links)
    rng = np.random.default_rng(410)
    for label, phase in (("common_time", .7 * np.arange(len(STEPS))[:, None, None, None]),
                         ("local_spacetime", rng.uniform(-np.pi, np.pi, state.shape[:-1]))):
        transformed, _ = overlap_links(state * np.exp(1j * phase)[..., None])
        correct = plaquettes(transformed)
        temporal_change = np.max(abs(np.angle(transformed[0] * links[0].conj())))
        error = max(np.max(abs(correct[k] - f)) for k, f in base.items())
        # A changed spatial phase convention with the old temporal link retained
        # defines different plaquettes. On this static reference this can give
        # apparent currents around closed paths, with no spin orientation motion.
        wrong = plaquettes([links[0], *transformed[1:]])
        _, j, _, _ = charges_current(wrong)
        false_current = sum(abs(part).sum() for part in j)
        print(f"{label} max_temporal_link_change={temporal_change:.8f} "
              f"covariant_F_error={error:.3g} old_time_link_absJ={false_current}")
        assert error < TOL and temporal_change > .1
        assert (false_current > 0) == (label == "local_spacetime")


def main():
    print("Pair = normalized product-state overlap, lattice spacing 1, current in turns per step.")
    relative, exterior, standard_j = transport("relative", (+1, -1))
    transport("common", (+1, +1))
    transport("fixed", (0, 0))
    shifted, _, alternate_j = transport("relative_core_shift_0.4", (+1, -1), core_shift=.4)
    diff = max(np.max(abs((shifted[key] - face)[exterior_faces(exterior, *key)]))
               for key, face in relative.items())
    change = sum(abs(j - alt).sum() for j, alt in zip(standard_j, alternate_j))
    print(f"core_fill_change max_exterior_plaquette_difference={diff:.3g} abs_J_difference={change}")
    assert diff < TOL and change > 0
    subcell, good, currents = transport("subcell", (+1, -1), np.array([0., .11]))
    exterior_dt = max(np.max(abs(np.diff(face, axis=0))[exterior_faces(good, *key).all(axis=0)[None]])
                      for key, face in subcell.items() if key[0] > 0)
    print(f"subcell exterior_max_delta_spatial_F={exterior_dt:.9g}")
    assert exterior_dt > 1e-6 and all(np.count_nonzero(j) == 0 for j in currents)
    gauge_control()
    # At sublattice separation the compact product need not resolve both units.
    paths = [np.full(2, x) for x in (.15, -.15)]
    close, _, _, wraps, margin = trial(paths)
    q, _, _, _ = charges_current(close)
    single_counts = []
    for path, sign in zip(paths, (+1, -1)):
        spin, _ = preparation(path, sign)
        single, _, _, _ = charges_current(plaquettes(overlap_links(spin)[0]))
        single_counts.append(int(single[0].sum()))
    print(f"close_pair separation=0.30 single_Q={single_counts} pair_abs_Q={abs(q[0]).sum()} "
          f"pair_wraps={wraps} min_branch_margin={margin:.8f}")
    assert single_counts == [1, -1] and wraps > 0 and np.count_nonzero(q) == 0


if __name__ == "__main__":
    main()
