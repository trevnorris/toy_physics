"""Receiving-order equilibrium and traces for a joint shape update.

Use the same three native blocks as mouth_supercore_terminal_order. Resolve
near-stationary constrained forces with exact quartic energy increments.
"""
import numpy as np
from scipy import sparse
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
from scipy.sparse.linalg import spsolve
from scipy.special import expit

from mouth_sealed_shape_relax import DELTA


def receiver(shape, mesh, cut, side_order_work, *, polar_shift=0., material_path=None):
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
        if not above:
            side_ids = ids
        if above and inner:
            terminal_ids=ids
        if above and not inner:
            outer_ids=ids
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
                rr = ww = volume = np.ones_like(measure)
                rw = np.zeros_like(measure)
                if polar_shift and not above:
                    from mouth_polar_bonded import outer_order_metric
                    rr, ww, rw, volume = outer_order_metric(shape, radial, z, polar_shift)
                local += (DELTA**2*measure)[..., None, None]*(
                    rr[...,None,None]*nr[..., :, None]*nr[..., None, :]
                    + ww[...,None,None]*nw[..., :, None]*nw[..., None, :]
                    + rw[...,None,None]*(nw[..., :, None]*nr[..., None, :]
                                        + nr[..., :, None]*nw[..., None, :]))
                lumps += (measure*volume)[...,None]*n
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
    well = 3 / DELTA
    source = np.zeros_like(chi)
    source[side_ids] = side_order_work(wa, xout)
    if polar_shift:
        from mouth_internal_collar_check import read
        from mouth_polar_bonded import outer_crystal_source_delta
        pair, labels, _, _ = read(material_path)
        pair = type(pair)(shape, pair.co, pair.na, pair.nt)
        source[side_ids] += outer_crystal_source_delta(pair, labels, wa, xout, polar_shift)
    scale = np.sqrt((2 * stiff.diagonal() + 2 * vol)[free])

    def gradient(c):
        return 2 * (stiff @ c) + 2 * vol * c * (1 - c) * (1 - 2*c) + source / well

    def objective(x):
        c = chi.copy()
        c[free] = x / scale
        return c @ (stiff @ c) + vol @ (c * (1-c))**2 + source @ c / well, gradient(c)[free] / scale

    result = minimize(objective, chi[free]*scale, jac=True, method="L-BFGS-B",
                      bounds=list(zip(np.zeros_like(scale), scale)),
                      options=dict(gtol=2e-8, ftol=3e-16, maxiter=2000, maxls=40))
    chi[free] = result.x / scale
    for iteration in range(35):
        g = gradient(chi)
        x = chi[free] * scale
        pg = x - np.clip(x-g[free]/scale, 0, scale)
        residual = np.max(abs(pg))
        if iteration == 0:
            initial = residual
        if residual <= 1e-8:
            break
        active = free & ~(((chi <= 0) & (g >= 0)) | ((chi >= 1) & (g <= 0)))
        h = 2*stiff[active][:, active] + sparse.diags((2*vol*(1-6*chi+6*chi**2))[active])
        step = np.zeros_like(chi)
        step[active] = spsolve(h.tocsc(), -g[active])
        descent = np.zeros_like(chi)
        descent[free] = -pg/scale
        accepted = False
        for direction in (step, descent):
            for halving in range(45):
                proposed = np.clip(chi + 2.**(-halving)*direction, 0, 1)
                d = proposed-chi
                slope = g @ d
                change = slope + d @ (stiff @ d) + (vol*d**2) @ (
                    1-6*chi+6*chi**2 + (4*chi-2)*d + d**2)
                if slope < 0 and np.isfinite(change) and change <= 1e-4*slope:
                    chi = proposed
                    accepted = True
                    break
            if accepted:
                break
        if not accepted:
            raise ArithmeticError(f"Receiver has no feasible energy descent at {residual}")
    else:
        raise ArithmeticError(f"Receiving-order stationarity unresolved at {residual}")
    if not np.isfinite(residual) or min(chi) < 0 or max(chi) > 1:
        raise ArithmeticError("Invalid bounded receiving order")
    disk = RegularGridInterpolator((wb-shape.length, xin), chi[terminal_ids])
    side = RegularGridInterpolator((wa, xout), chi[side_ids])
    outer = RegularGridInterpolator((wb, xout), chi[outer_ids])

    def profile(x, z):
        return disk((np.clip(z, 0, wb[-1]-shape.length), x))

    def outside(r, z, above):
        b, _ = shape.radius(z)
        x = (r-b)/(cut-b)
        clipped = np.clip(x, 0, 1)
        local = (outer((np.clip(z, wb[0], wb[-1]), clipped)) if above else
                 side((np.clip(z, 0, shape.length), clipped)))
        return np.where(x > 1., expit(-z/DELTA), local)

    profile.outside = outside
    profile.side = (wa, xout, chi[side_ids])
    info = dict(energy=float(well*(chi@(stiff@chi) + vol@((chi*(1-chi))**2))),
                static_work=float(source@chi), residual=float(residual),
                before_newton=float(initial), newton_steps=iteration,
                minimum_order=float(min(chi)), maximum_order=float(max(chi)))
    return profile, info
