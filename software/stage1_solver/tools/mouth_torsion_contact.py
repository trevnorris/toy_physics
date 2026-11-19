"""Incremental torque on registered links from one tilted material patch."""

from itertools import product

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve


def cut_graph(radius):
    sites = [v for v in product(range(-radius, radius + 1), repeat=3)
             if max(map(abs, v)) == radius]
    ids = {v: i for i, v in enumerate(sites)}
    edges = []
    for v in sites:
        for axis in range(3):
            w = list(v)
            w[axis] += 1
            if tuple(w) in ids:
                edges.append((ids[v], ids[tuple(w)]))
    ends = np.array(edges)
    incidence = coo_matrix((np.tile([-1., 1.], len(ends)),
                            (np.repeat(np.arange(len(ends)), 2), ends.ravel())),
                           shape=(len(ends), len(sites))).tocsr()
    return np.array(sites), ends, incidence, ids[(radius, 0, 0)]


def mouth_response(radius, sign, flat=False):
    sites, ends, grad, patch = cut_graph(radius)
    normals = sign * sites / np.linalg.norm(sites, axis=1)[:, None]
    if flat:
        normals[:] = (sign, 0, 0)
    tilt = np.zeros_like(normals)
    tilt[patch, 1] = sign
    i, j = ends.T
    dot = np.einsum("ij,ij->i", normals[i], normals[j])
    overlap = np.sqrt((1 + dot) / 2)
    # Derivative of arg(z_i^dag z_j) when each direction follows the
    # shortest material rotation, with no extra rotation about the normal.
    beta = -np.sum((tilt[i]+tilt[j]) * np.cross(normals[i], normals[j]), axis=1) / (2*(1+dot))
    weight = 4 * overlap
    rigid = weight * beta
    reduced = grad[:, 1:]  # common roll costs no energy
    hess = reduced.T @ diags(weight) @ reduced
    roll = spsolve(hess.tocsc(), -reduced.T @ rigid)
    free = weight * (beta + reduced @ roll)

    # Independent derivative from two actual spinor preparations.
    x, y, z = normals.T
    states = np.stack((np.sqrt((1+z)/2), np.exp(1j*np.arctan2(y, x))*np.sqrt((1-z)/2)), axis=1)
    axis = np.cross(normals[patch], tilt[patch])
    ax, ay, az = axis
    pauli = np.array([[az, ax-1j*ay], [ax+1j*ay, -az]])
    eps = 1e-4
    links = []
    for angle in (-eps, eps):
        moved = states.copy()
        moved[patch] = (np.cos(angle/2)*np.eye(2) - 1j*np.sin(angle/2)*pauli) @ states[patch]
        links.append(np.sum(moved[i].conj() * moved[j], axis=1))
    numeric_beta = np.angle(links[0].conj()*links[1]) / (2*eps)
    # For zero bulk wave number a uniform edge angle is the gradient of a
    # Cartesian coordinate on this closed finite surface.
    uniform_coupling = np.linalg.norm((sites[j]-sites[i]).T @ free)
    print(f"R={radius} texture={0 if flat else sign:+d} "
          f"rigid_torque_norm={np.linalg.norm(rigid):.9f} "
          f"free_torque_norm={np.linalg.norm(free):.9f} "
          f"rigid_div_norm={np.linalg.norm(grad.T @ rigid):.9f} "
          f"free_div_norm={np.linalg.norm(grad.T @ free):.3g} "
          f"free_uniform_coupling={uniform_coupling:.3g} "
          f"spinor_derivative_error={max(abs(numeric_beta-beta)):.3g}")
    assert max(abs(numeric_beta-beta)) < 1e-8 and uniform_coupling < 1e-10
    assert np.linalg.norm(grad.T @ free) < 1e-10
    return free, rigid


if __name__ == "__main__":
    for radius in (5, 10):
        a = mouth_response(radius, 1)
        b = mouth_response(radius, -1)
        control = mouth_response(radius, 1, flat=True)
        print(f"R={radius} mirrored_torque_sum_norm={np.linalg.norm(a[0]+b[0]):.3g}")
        assert np.linalg.norm(a[0]) > 0 and np.allclose(a, -np.array(b))
        assert not np.any(control)
