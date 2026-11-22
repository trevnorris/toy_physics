"""Resolve the neutral equilibrium order inside a saved, axisymmetric, held throat."""
import argparse
import json

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

from mouth_sealed_shape_relax import DELTA, Meridian


def order_grid(mesh, shape):
    # r=x b(w), with the S² integration done exactly. Q1 gradient operator
    # and positive lumped integration of the local double-well energy.
    x = np.linspace(0, 1, mesh+1)
    w = np.linspace(0, shape.length, 2*mesh+1)
    rr = shape.radius(w[:, None])[0] * x
    ids = np.arange(rr.size).reshape(rr.shape)
    elements = np.stack([ids[:-1, :-1], ids[:-1, 1:], ids[1:, :-1], ids[1:, 1:]], -1)
    dx, dw = x[1], w[1]
    lump = np.zeros(rr.size)
    local = np.zeros((*elements.shape[:2], 4, 4))
    nodes, weights = np.polynomial.legendre.leggauss(3)
    for i, a in enumerate((nodes+1)/2):
        u = x[:-1] + a*dx
        for j, b in enumerate((nodes+1)/2):
            radius, slope = shape.radius(w[:-1, None]+b*dw)
            dr = np.array([b-1, 1-b, -b, b]) / (dx*radius[:, :, None])
            dz = (np.array([a-1, -a, 1-a, a])/dw
                  - (slope*u)[:, :, None]*dr)
            basis = np.array([(1-a)*(1-b), a*(1-b), (1-a)*b, a*b])
            volume = np.pi*(u*radius)**2*radius*dx*dw*weights[i]*weights[j]
            local += (DELTA**2*volume[:, :, None, None])*(
                dr[..., :, None]*dr[..., None, :] + dz[..., :, None]*dz[..., None, :])
            np.add.at(lump, elements, volume[:, :, None]*basis)
    rows = np.broadcast_to(elements[..., :, None], local.shape)
    cols = np.broadcast_to(elements[..., None, :], local.shape)
    stiff = coo_matrix((local.ravel(), (rows.ravel(), cols.ravel())),
                       shape=(rr.size, rr.size)).tocsr()
    return w, rr, ids, lump, stiff


def material(mesh, shape, axial_radius):
    w, rr, ids, lump, stiff = order_grid(mesh, shape)
    free = np.ones_like(rr, dtype=bool)
    free[[0, -1], :] = False
    free[:, -1] = False
    free = free.ravel()
    distance = np.minimum(rr[:, -1:] - rr, shape.length-w[:, None])
    u = (.5*np.tanh(distance/(2*DELTA))).ravel()  # chi-1/2
    # Reproducible and independent of the starting guess in the interior.
    u[ids[0]] = .5*np.tanh((shape.neck-rr[0])/(2*DELTA))
    u[ids[-1]] = 0
    u[ids[:, -1]] = 0
    reduced = 2*stiff[free][:, free]

    def energy(y):
        return y@(stiff@y) + lump@((y*y-.25)**2)

    err = None
    for iteration in range(60):
        g = (2*(stiff@u) + lump*(4*u**3-u))[free]
        hess = reduced + diags(lump[free]*(12*u[free]**2-1))
        du = spsolve(hess, -g)
        err = max(abs(du))
        if err < 1e-10:
            break
        if g@du >= 0:
            raise ArithmeticError("Order correction does not lower the patch energy")
        old = energy(u)
        alpha = 1
        while alpha > 1e-8:
            trial = u.copy()
            trial[free] += alpha*du
            if energy(trial) <= old + 1e-12 and min(trial) > -1e-7:
                u = trial
                break
            alpha /= 2
        else:
            raise ArithmeticError("Order line search failed")
    else:
        raise ArithmeticError("Order profile not stationary")
    axis = .5 + u[ids[:, 0]]
    if np.min(axis) < .5-1e-7 or np.max(axis) > 1+1e-6:
        raise ArithmeticError("Order left the prepared branch")
    cut_w = axial_radius*shape.length
    return dict(mesh=mesh, waist_depth=cut_w, waist_radius=float(shape.radius(cut_w)[0]),
                axis_order_at_waist=float(np.interp(cut_w, w, axis)),
                ordered_side_energy=energy(u), iterations=iteration, last_correction=err)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=int, required=True, choices=(1, 2))
    parser.add_argument("--shape", nargs=4, type=float, required=True)
    args = parser.parse_args()
    shape = Meridian(args.radius, args.shape)
    ts = np.linspace(0, 1, 4097)
    i = np.argmin(shape.profile(ts)[0])
    if i in (0, len(ts)-1):
        raise ValueError("No internal waist in sample")
    rows = [material(mesh, shape, ts[i]) for mesh in (64, 112)]
    print(json.dumps(dict(radius=args.radius, coordinates=args.shape, rows=rows), allow_nan=False),
          flush=True)
    if abs(rows[0]["axis_order_at_waist"]-rows[1]["axis_order_at_waist"]) > .015:
        raise ArithmeticError("Centerline reading is mesh sensitive")


if __name__ == "__main__":
    main()
