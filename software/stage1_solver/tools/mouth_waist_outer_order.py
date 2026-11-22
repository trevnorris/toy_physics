"""Neutral receiving-side half of the registered order front on held meridians."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.special import expit

from mouth_sealed_shape_relax import DELTA, Meridian


ROOT = Path(__file__).resolve().parents[3]
GEOMETRY = ROOT / "research/pde_ledger_v3/_measurements/mouth_sealed_shape_relax.json"


def archived_shape(radius, state):
    rows = json.loads(GEOMETRY.read_text())["training"][str(radius)]
    if state == "broad":
        sample = rows["states"][-2]
    elif state == "mid":
        sample = rows["extension_states"][1]
    else:
        sample = rows["extension_states"][-1]
    return np.array(sample["coordinates"])


def exterior(shape, mesh, cut=4.6, *, return_side_field=False, side_order_work=None):
    # Three material blocks meet at the disk rim. The receiver wraps the
    # outside of the duct and occupies the whole domain above the terminal disk.
    t = np.linspace(0, 1, mesh+1)
    xout = .10*t + .90*t**2
    xin = 1 - (.16*(1-t)+.84*(1-t)**2)
    wa = shape.length*(1-np.cos(np.linspace(0, np.pi, 2*mesh+1)))/2
    wb = shape.length+16*DELTA*(.36*t+.64*t**2)
    index = {}
    coordinates = []
    weights = []
    blocks = []
    boundary = {}
    q, qw = np.polynomial.legendre.leggauss(3)
    for xt, zt, inner, above in ((xout, wa, False, False),
                                 (xout, wb, False, True), (xin, wb, True, True)):
        b, _ = shape.radius(zt[:, None])
        r = b*xt if inner else b + (cut-b)*xt
        ids = np.empty_like(r, dtype=int)
        for j, i in np.ndindex(r.shape):
            key = tuple(np.round([zt[j], r[j, i]], 12))
            node = index.get(key)
            if node is None:
                node = index[key] = len(index)
                coordinates.append((r[j, i], zt[j]))
            ids[j, i] = node
            if (not above and (i == 0 or j == 0)) or (inner and j == 0):
                boundary[node] = .5
            elif above and j == len(zt)-1:
                boundary[node] = 0.
            elif i == len(xt)-1 and not inner:
                boundary[node] = expit(-zt[j]/DELTA)
        if (return_side_field or side_order_work is not None) and not above:
            side_ids = ids
        e = np.stack([ids[:-1, :-1], ids[:-1, 1:], ids[1:, :-1], ids[1:, 1:]], -1)
        dx = np.diff(xt)[None, :]
        dz = np.diff(zt)[:, None]
        local = np.zeros((*e.shape, 4))
        lumps = np.zeros(e.shape)
        for a, ga in zip((q+1)/2, qw/2):
            x = xt[None, :-1]+a*dx
            for bq, gb in zip((q+1)/2, qw/2):
                z = zt[:-1, None]+bq*dz
                radius, slope = shape.radius(z)
                if inner:
                    width, alpha, radial = radius, x, radius*x
                else:
                    width, alpha = cut-radius, 1-x
                    radial = radius+x*width
                nr = np.array([bq-1, 1-bq, -bq, bq])/(width*dx)[..., None]
                nw = np.array([a-1, -a, 1-a, a])/dz[..., None] - (alpha*slope)[..., None]*nr
                measure = 4*np.pi*radial**2*width*dx*dz*ga*gb
                n = np.array([(1-a)*(1-bq), a*(1-bq), (1-a)*bq, a*bq])
                local += (DELTA**2*measure)[..., None, None]*(
                    nr[..., :, None]*nr[..., None, :] + nw[..., :, None]*nw[..., None, :])
                lumps += measure[..., None]*n
        blocks.append((e.reshape(-1, 4), local.reshape(-1, 4, 4)))
        weights.append((e, lumps))
    ids = np.concatenate([b[0] for b in blocks])
    matrices = np.concatenate([b[1] for b in blocks])
    stiff = sparse.coo_matrix((matrices.ravel(), (
        np.broadcast_to(ids[:, :, None], matrices.shape).ravel(),
        np.broadcast_to(ids[:, None, :], matrices.shape).ravel())),
        shape=(len(index), len(index))).tocsr()
    vol = np.zeros(len(index))
    for e, load in weights:
        np.add.at(vol, e, load)
    r, w = np.asarray(coordinates).T
    b, _ = shape.radius(w)
    # Initial distance to disk edge, plane or side; the order then relaxes.
    side = np.hypot(np.maximum(r-b, 0), np.maximum(w-shape.length, 0))
    near = np.minimum(w, side)
    chi = expit(-near/DELTA)
    pins = np.array(list(boundary))
    chi[pins] = list(boundary.values())
    free = np.ones(len(index), dtype=bool)
    free[pins] = False
    reduced = 2*stiff[free][:, free]

    def energy(c):
        return c@(stiff@c) + vol@((c*(1-c))**2)

    if side_order_work is not None:
        from scipy.optimize import minimize
        # The supplied coefficient multiplies the physical nodal order once.
        source = np.zeros_like(chi)
        source[side_ids] = side_order_work(wa, xout)
        scale = np.sqrt((2*stiff.diagonal()+2*vol)[free])
        well = 3/DELTA

        def objective(x):
            c = chi.copy()
            c[free] = x/scale
            g = 2*(stiff@c) + vol*2*c*(1-c)*(1-2*c) + source/well
            return energy(c) + source@c/well, g[free]/scale

        sol = minimize(objective, chi[free]*scale, jac=True, method="L-BFGS-B",
                       bounds=list(zip(np.zeros_like(scale), scale)),
                       options=dict(gtol=2e-8, ftol=3e-16, maxiter=2000, maxls=40))
        chi[free] = sol.x/scale
        value, gradient = objective(sol.x)
        projected = sol.x-np.clip(sol.x-gradient, 0, scale)
        residual = np.max(abs(projected))
        # A line search can hit objective roundoff at an already stationary
        # solution. Verify the projected gradient and an actual feasible step.
        plateau = (sol.status == 2 and residual < 2e-7 and
                   objective(sol.x-projected)[0] >= value-16*np.spacing(abs(value)))
        if (not np.isfinite(value + residual) or residual > 7e-6
                or not (sol.success or plateau)):
            raise ArithmeticError(f"Receiving order search failed at residual {residual}: {sol.message}")
        result = dict(energy=well*energy(chi), static_work=source@chi,
                      residual=residual, iterations=sol.nit,
                      minimum_order=min(chi), maximum_order=max(chi))
        if return_side_field:
            return result, (wa, xout, chi[side_ids])
        return result

    for iteration in range(60):
        residual = (2*(stiff@chi) + vol*(2*chi*(1-chi)*(1-2*chi)))[free]
        hess = reduced+sparse.diags(vol[free]*(2-12*chi[free]+12*chi[free]**2))
        step = spsolve(hess, -residual)
        if max(abs(step)) < 1e-10:
            break
        old = energy(chi)
        alpha = 1.
        while alpha > 1e-8:
            trial = chi.copy()
            trial[free] += alpha*step
            if energy(trial) < old+1e-12 and max(trial) <= .50001:
                chi = trial
                break
            alpha /= 2
        else:
            raise ArithmeticError("Exterior order failed its energy descent check")
    else:
        raise ArithmeticError("Exterior order did not settle")
    a_b = 3/DELTA
    result = dict(energy=a_b*energy(chi), minimum_order=min(chi), maximum_order=max(chi),
                  residual=max(abs(residual)), last_step=max(abs(step)), iterations=iteration)
    if return_side_field:
        return result, (wa, xout, chi[side_ids])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", required=True, type=int, choices=(1, 2))
    parser.add_argument("--state", required=True, choices=("broad", "mid", "tight"))
    parser.add_argument("--mesh", required=True, type=int)
    a = parser.parse_args()
    x = archived_shape(a.radius, a.state)
    for axis, i, scale in (("radial", 2, 1/a.radius), ("axial", 0, 1.)):
        energies = {}
        for step in (0., -.002, -.001, .001, .002):
            dx = x.copy()
            dx[i] += scale*step
            row = exterior(Meridian(a.radius, dx), a.mesh)
            print(json.dumps(dict(radius=a.radius, state=a.state, axis=axis, step=step,
                                  coordinates=dx.tolist(), mesh=a.mesh, **row)), flush=True)
            energies[step] = row["energy"]
        for h in (.001, .002):
            print(json.dumps(dict(axis=axis, h=h, slope=(energies[h]-energies[-h])/(2*h),
                                  stiffness=(energies[h]+energies[-h]-2*energies[0.])/h**2)),
                  flush=True)


if __name__ == "__main__":
    main()
