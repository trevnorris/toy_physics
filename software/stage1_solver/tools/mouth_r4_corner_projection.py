"""Project a packed active-label direction onto one-sided determinant tangent planes."""
from itertools import combinations, product

import mpmath as mp

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


class ProjectionError(ArithmeticError):
    """An unresolved projection, with the rows needed to reproduce it."""
    def __init__(self, message, diagnostics):
        super().__init__(message)
        self.diagnostics = diagnostics


class NonDescendingProjection(ProjectionError):
    """The verified projection of this particular input is non-descending."""


def _extended_dual(rows, response, direction, rhs, seed, diagonal):
    """Solve the small dual on independent binding-row sets using original row products."""
    used = np.any(rows != 0, axis=0)
    with mp.workdps(90):
        A = mp.matrix(rows[:,used].tolist())
        if diagonal is None:
            B = mp.matrix(response[:,used].tolist())
        else:
            weights = diagonal[used]
            B = mp.matrix([[A[i,k]*float(w) for k,w in enumerate(weights)] for i in range(A.rows)])
        G = A * B.T
        b = A * mp.matrix(direction[used].tolist()) - mp.matrix(rhs.tolist())
        # A redundant row need never be chosen as an independent multiplier.
        # Inspect the NNLS support first; the exhaustive fallback is for small m.
        def supports():
            yield seed
            for size in range(A.rows+1):
                for bank in combinations(range(A.rows),size):
                    if bank != seed:
                        yield bank
        for bank in supports():
            lam = mp.zeros(A.rows,1)
            if bank:
                block = mp.matrix([[G[i,j] for j in bank] for i in bank])
                try:
                    active = mp.lu_solve(block,mp.matrix([-b[i] for i in bank]))
                except (ValueError,ZeroDivisionError):
                    continue
                if min(active)<0:
                    continue
                for i,x in zip(bank,active):
                    lam[i] = x
            slack = b+G*lam
            if min(slack) >= -mp.mpf('1e-65')*max(1,max(map(abs,b))):
                # Form the projected direction before rounding, including any
                # motion outside constraint support for a nonlocal legacy metric.
                out = direction.copy()
                moved = np.flatnonzero(np.any(response != 0,axis=0))
                for j in moved:
                    if diagonal is None:
                        correction = mp.fdot(lam, response[:,j].tolist())
                    else:
                        correction = float(diagonal[j])*mp.fdot(lam,rows[:,j].tolist())
                    out[j] = float(mp.mpf(float(direction[j]))+correction)
                return out,np.array([float(x) for x in lam])
    return None


def project_halfspaces(direction, normals, metric=None, *, gradient=None, return_info=False,
                       rhs=None, diagonal=None):
    """Project p onto A d >= b in the inverse-metric norm; directly certify the answer."""
    p = np.asarray(direction,dtype=float)
    normals = np.asarray(normals,dtype=float).reshape((-1,len(p)))
    lengths = np.linalg.norm(normals,axis=1)
    if not np.all(np.isfinite(lengths) & (lengths>0)):
        raise ValueError('Nonfinite or empty constraint row')
    A = normals/lengths[:,None]
    b = np.zeros(len(A)) if rhs is None else np.asarray(rhs,dtype=float)/lengths
    if b.shape != (len(A),) or not np.all(np.isfinite(b)):
        raise ValueError('Invalid corner constraint right-hand side')
    if diagonal is not None:
        diagonal = np.asarray(diagonal,dtype=float)
        if metric is not None or diagonal.shape!=p.shape or not np.all(np.isfinite(diagonal) & (diagonal>0)):
            raise ValueError('Expected one positive diagonal inverse metric')
        response = A*diagonal
    else:
        response = np.array([metric(a) for a in A]).reshape(A.shape) if metric else A
    relevant = np.any((A!=0) | (response!=0),axis=0)
    # Verify against the returned binary64 direction, in extended arithmetic.
    # Outside this support the direction is copied bit for bit.
    L = np.longdouble
    al = A[:,relevant].astype(L)
    pl = p[relevant].astype(L)
    rl = response[:,relevant].astype(L)

    def check(out, lam, method):
        dl = out[relevant].astype(L)
        ll = lam.astype(L)
        slack = al@dl-b
        change = dl-pl
        if diagonal is not None:
            station = change/diagonal[relevant] - ll@al
        else:
            station = change - ll@rl
        primal = np.minimum(slack,0)
        dual = np.minimum(ll,0)
        comp = ll*slack
        norm = float(np.linalg.norm(np.concatenate((station,primal,dual,comp))))
        slope = None if gradient is None else float(np.dot(gradient,out))
        info = dict(method=method,multipliers=lam.tolist(),slacks=slack.astype(float).tolist(),
                    primal_violation=float(-np.min(primal,initial=0)),dual_violation=float(-np.min(dual,initial=0)),
                    stationarity=float(np.linalg.norm(station)),complementarity=float(np.max(abs(comp),initial=0)),
                    kkt_norm=norm,slope=slope)
        ok = (np.isfinite(norm) and norm<=1e-10 and info['dual_violation']<=1e-12
              and (slope is None or np.isfinite(slope) and slope<0))
        return info,ok

    def diagnostics(info):
        return dict(check=info,columns=np.flatnonzero(relevant).tolist(),normalized_rows=A[:,relevant].tolist(),
                    normalized_rhs=b.tolist(),input_on_columns=p[relevant].tolist())

    lam = np.zeros(len(A))
    values = A@p-b
    try:
        if np.any(values<0):
            gram = A@response.T
            if len(A)==1:
                lam[0] = -values[0]/gram[0,0]
            else:
                eig,vec = np.linalg.eigh((gram+gram.T)/2)
                keep = eig > eig.max()*len(A)*np.finfo(float).eps
                sq = np.sqrt(eig[keep])
                mat = sq[:,None]*vec[:,keep].T
                lam = nnls(mat,-(vec[:,keep].T@values)/sq)[0]
        out = p + lam@response
        info,passed = check(out,lam,'float64')
    except (RuntimeError,ArithmeticError,ValueError,np.linalg.LinAlgError) as err:
        info,passed = dict(method='float64',error=str(err)),False
    if not passed:
        initial = info
        if len(A):
            seed = tuple(int(i) for i in np.flatnonzero(lam>0))
            polished = _extended_dual(A,response,p,b,seed,diagonal)
            if polished is None:
                raise ProjectionError('No verified extended-precision KKT solution',diagnostics(initial))
            out,lam = polished
        else:
            out = p.copy()
        info,passed = check(out,lam,'extended')
        info['float_check'] = initial
        if not passed:
            if info['kkt_norm']<=1e-10 and info['dual_violation']<=1e-12 and gradient is not None and info['slope']>=0:
                raise NonDescendingProjection('Projected input has no strict descent',diagnostics(info))
            raise ProjectionError('Returned extended projection failed direct KKT verification',diagnostics(info))
    return (out,info) if return_info else out
