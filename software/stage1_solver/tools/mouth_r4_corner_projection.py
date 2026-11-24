"""Project a packed active-label direction onto one-sided determinant tangent planes."""
from itertools import product

import numpy as np
from scipy.optimize import nnls

from mouth_cp10_label_condition import physical_frame, reference_frame

BITS = np.array(list(product(range(2), repeat=3)))


def packed_lift(opt, direction, harmonic=True):
    """Full nodal increment with held groups fixed and optional harmonic dead lift."""
    result = np.zeros_like(opt.base)
    start = 0
    for m, cm in enumerate(opt._components()):
        stop = start + int(cm.free.sum())
        block = direction[start:stop]
        if harmonic:
            increments = cm.lift(block)
        else:
            values = np.zeros(len(cm.free))
            values[cm.free] = block
            increments = values[cm.groups]
        result[..., m] = increments.reshape(opt.stencil.shape)
        start = stop
    if start != len(direction):
        raise ValueError('Wrong packed dimension')
    return result


def point_normal(opt, labels, cell, local, reference=None):
    """Packed active-label normal at fixed chart position; use metric Jinv with reference.

    The chart normal is regular on the polar axis. Physical normals require a
    regular physical frame and a one-sided differentiable material embedding.
    The corner gradient sums over active registry groups with held/dead DOFs fixed.
    """
    cell = np.asarray(cell, dtype=int)
    t = np.asarray(local, dtype=float)
    ids = BITS + cell
    nodes = labels[tuple(ids.T)]
    widths = np.array([ax[i+1]-ax[i] for ax, i in zip(opt.axes, cell)])
    point = np.array([ax[i] for ax, i in zip(opt.axes, cell)]) + widths*t
    factors = np.where(BITS, t, 1-t)
    basis = factors.prod(axis=1)
    basis_gradient = np.column_stack([(2*BITS[:,k]-1)*factors[:,(k+1)%3]*factors[:,(k+2)%3]/widths[k]
                                      for k in range(3)])
    xi = basis@nodes
    # Subtract common offset before differentiation to avoid cancellation at shallow cells.
    D = (nodes-nodes[0]).T @ basis_gradient
    cofactor = np.column_stack([np.cross(D[:,(k+1)%3], D[:,(k+2)%3]) for k in range(3)])
    chart_det = float(D[:,0] @ cofactor[:,0])
    gradient = basis_gradient @ cofactor.T
    scale = 1.
    if reference is not None:
        B, h = physical_frame(opt.body, point[None])
        C, h0 = reference_frame(opt.body, reference, xi[None])
        scale = float(abs(np.linalg.det(C[0]))*h0[0]/(np.linalg.det(B[0])*h[0]))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError('Physical metric ratio is singular here; use the chart determinant')
        x, s, a = xi
        if x < opt.body.cut:
            radius, slope = reference.radius(s*reference.length)
            log_grad = np.array([2/x, 3*reference.length*slope/radius, 1/np.tan(a)])
        else:
            r = C[0,2,2]
            log_grad = np.array([-1/x+2*C[0,0,0]/r, 2*C[0,0,1]/r, 1/np.tan(a)])
        gradient = scale*(gradient + chart_det*basis[:,None]*log_grad)
    flat = np.ravel_multi_index(tuple(ids.T), opt.stencil.shape)
    normal = np.concatenate([np.bincount(cm.groups[flat], weights=gradient[:,k],
                                        minlength=len(cm.free))[cm.free]
                             for k, cm in enumerate(opt._components())])
    return normal, dict(cell=cell.tolist(), local=t.tolist(), point=point.tolist(), label=xi.tolist(),
                        chart_det=chart_det, metric_factor=scale, determinant=scale*chart_det,
                        nodes=ids.tolist(), node_gradient=gradient.tolist(), basis=basis.tolist())


def project_halfspaces(direction, normals, metric=None):
    """Project in inverse-metric norm onto n·d>=0; unit steps still need positivity checks."""
    N = np.atleast_2d(normals).astype(float)
    if not N.size:
        return direction.copy()
    sizes = np.linalg.norm(N, axis=1)
    if np.any(sizes == 0):
        raise ValueError('Empty active constraint normal')
    N = N/sizes[:,None]
    response = np.array([metric(n) for n in N]) if metric else N
    values = N@direction
    if np.min(values) >= 0:
        return direction.copy()
    grams = N@response.T
    if len(N) == 1:
        result = direction + response[0]*(-values[0]/grams[0,0])
    else:
        eig, vec = np.linalg.eigh((grams+grams.T)/2)
        keep = eig > eig.max()*max(N.shape[0], 1)*np.finfo(float).eps
        sq = np.sqrt(eig[keep])
        matrix = sq[:,None]*vec[:,keep].T
        rhs = -(vec[:,keep].T@values)/sq
        multipliers = nnls(matrix, rhs)[0]
        result = direction + multipliers@response
    if np.min(N@result) < -1e-10*max(1., np.linalg.norm(direction)):
        raise ArithmeticError('Unresolved corner projection')
    return result
