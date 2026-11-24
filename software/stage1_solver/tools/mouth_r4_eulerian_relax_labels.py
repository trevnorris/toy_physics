"""Constrained minimization of the welded R4 spatial labels on a fixed physical grid.

The two physical metrics and their moving weld belong to solid.evaluate. Here
only the native registry, polar regularity and the zero-order harmonic gauge
are imposed. Strictly melted, unconstrained degrees of freedom are eliminated
by their unweighted harmonic equation. Gauge strength never enters that solve.
"""
from itertools import product

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_cp10_label_grid import BITS, matrix_for
from mouth_r4_eulerian_solid import InvalidSpatialCrystal


FORCE_TOL = 3.e-5
HARMONIC_TOL = 3.e-10
UROUND = np.finfo(float).eps


class LabelStall(RuntimeError):
    def __init__(self, reason, labels, diagnostics):
        self.labels = labels.copy()
        self.diagnostics = dict(diagnostics, stationary=False, reason=reason)
        super().__init__(reason)


class InfeasibleLabels(ValueError):
    def __init__(self, reason, **witness):
        self.witness = dict(reason=reason, **witness)
        super().__init__(reason)


def absmax(x):
    return float(np.max(np.abs(x), initial=0.))


def cell_order(order):
    windows = np.lib.stride_tricks.sliding_window_view(order, (2, 2, 2))
    return np.any(windows > 0, axis=(-3, -2, -1))


def slices(bits, dims):
    return tuple(slice(k, k+n) for k,n in zip(bits, dims))


def bisect_quadratic(coefficients, axis):
    b = np.moveaxis(coefficients, axis, -1)
    half_lo = (b[...,0] + b[...,1])*.5
    half_hi = (b[...,1] + b[...,2])*.5
    mid = (half_lo + half_hi)*.5
    return [np.moveaxis(q, -1, axis) for q in
            (np.stack((b[...,0], half_lo, mid), axis=-1),
             np.stack((mid, half_hi, b[...,2]), axis=-1))]


def certify_subcell(bern, margin, cell):
    """Use Bernstein subdivision only if the parent enclosure straddles zero."""
    stack = [(bern, np.zeros(3), np.ones(3), 0)]
    bound = np.inf
    while stack:
        b, lo, hi, depth = stack.pop()
        lower = np.min(b)-margin
        if lower > 0:
            bound = min(bound, lower)
            continue
        corners = b[::2,::2,::2]
        k = np.unravel_index(corners.argmin(), corners.shape)
        if corners[k] <= margin or depth >= 12:
            where = lo + np.asarray(k)*(hi-lo)
            raise InfeasibleLabels('Nonpositive or roundoff-unresolved chart determinant',
                cell=cell, chart_fractions=where.tolist(), chart_det=float(corners[k]),
                chart_lower=float(lower), roundoff_margin=float(margin), subdivisions=depth)
        middle = (lo+hi)*.5
        for i, u in enumerate(bisect_quadratic(b, 0)):
            for j, v in enumerate(bisect_quadratic(u, 1)):
                for k, w in enumerate(bisect_quadratic(v, 2)):
                    upper = np.array([i,j,k], bool)
                    stack.append((w, np.where(upper,middle,lo), np.where(upper,hi,middle), depth+1))
    return bound


def chart_orientation(axes, labels, active):
    """Whole-cell lower bounds on det(dxi/dp), not only values at Gauss sites."""
    if not np.isfinite(labels).all():
        raise InfeasibleLabels('Nonfinite label coordinates')
    for d, axis in enumerate(axes):
        coord = labels[...,d]
        # Q1 shape functions are positive and form a partition of unity.
        if coord.min() < axis[0]-2.e-13 or coord.max() > axis[-1]+2.e-13:
            raise InfeasibleLabels('Label leaves the saved material domain',
                component=d, lowest=float(coord.min()), highest=float(coord.max()))
    dims = active.shape
    vertices = {b:labels[slices(b,dims)] for b in BITS}
    cols = []
    norm_product = np.ones(dims)
    for d, axis in enumerate(axes):
        tangents = [k for k in range(3) if k != d]
        width = np.diff(axis).reshape(tuple(len(axis)-1 if k==d else 1 for k in range(3))+(1,))
        corners = []
        for l,m in product(range(2),repeat=2):
            b = [0,0,0]; b[tangents[0]]=l; b[tangents[1]]=m
            start = vertices[tuple(b)]; b[d] = 1
            corners.append((vertices[tuple(b)]-start)/width)
        rows = []
        for u,v in product((0.,.5,1.),repeat=2):
            rows.append((1-u)*((1-v)*corners[0]+v*corners[1])
                        + u*((1-v)*corners[2]+v*corners[3]))
        cols.append(rows)
        norm_product *= np.maximum.reduce([np.sqrt(np.sum(c*c,axis=-1)) for c in corners])
    bern = np.empty(dims+(3,3,3))
    for i,j,k in product(range(3),repeat=3):
        x,y,z = cols[0][3*j+k], cols[1][3*i+k], cols[2][3*i+j]
        bern[...,i,j,k] = np.einsum('...i,...i->...', x, np.cross(y,z))
    for d in (3,4,5):
        b = np.moveaxis(bern,d,-1)
        b[...,1] = 2*b[...,1]-.5*(b[...,0]+b[...,2])
    margin = 128*UROUND*np.maximum(norm_product,np.finfo(float).tiny)
    lower = bern.min(axis=(-1,-2,-3))-margin
    for row in np.argwhere(lower<=0):
        c = tuple(row)
        lower[c] = certify_subcell(bern[c],margin[c],row.tolist())
    idx = np.unravel_index(lower.argmin(), dims)
    return dict(chart_det_global_lower=float(lower[idx]),chart_lowest_cell=list(map(int,idx)),
                chart_det_active_lower=float(lower[active].min()) if active.any() else None)


class Component:
    """Physical registry and a sparse harmonic extension for one label component."""
    def __init__(self, ids, idx, held, stencil, K):
        _, self.groups = np.unique(idx,return_inverse=True)
        self.groups = self.groups.ravel()
        count = self.groups.max()+1
        self.C = sparse.csr_matrix((np.ones(ids.size),(ids.ravel(),self.groups)),shape=(ids.size,count))
        self.K = (self.C.T @ K @ self.C).tocsr()
        self.active = np.bincount(self.groups,weights=stencil.ravel(),minlength=count)>0
        self.held = np.bincount(self.groups,weights=held.ravel(),minlength=count)>0
        self.free = self.active & ~self.held
        self.dead = ~self.active & ~self.held
        self.representative = np.full(count,-1,dtype=np.intp)
        self.representative[self.groups] = ids.ravel()
        self.source = np.full(count,-1,dtype=np.intp)
        self.source[self.groups[stencil.ravel()]] = ids[stencil]
        self.extension = splu(self.K[self.dead][:,self.dead].tocsc()) if self.dead.any() else None
        self.cross = self.K[self.dead][:,self.free].tocsr()

    def pack(self, delta, stencil):
        v = np.zeros(len(self.held))
        flat = delta.ravel()
        v[self.free] = flat[self.source[self.free]]
        if absmax((v[self.groups]-flat)[stencil.ravel()])>2.e-9:
            raise ValueError('Conflicting regularity or materially displaced registered label')
        return v[self.free]

    def lift(self, x):
        v = np.zeros(len(self.held))
        v[self.free] = x
        if self.extension is not None:
            v[self.dead] = self.extension.solve(-(self.cross@x))
        return v[self.groups]

    def sum(self, values):
        return np.bincount(self.groups,weights=values.ravel(),minlength=len(self.held))

    def preconditioner(self, Kactive, scale):
        # This damped elliptic inverse is an optimization metric only. It is
        # never included in the energy, force, orientation or stopping tests.
        a = (self.C.T @ Kactive @ self.C).tocsr()[self.free][:,self.free]
        if a.shape[0]==0:
            return lambda rhs: rhs
        w = scale.ravel()[self.source[self.free]]
        stiffness = a + sparse.diags(np.maximum(a.diagonal(),1.e-12)*1.e-6)
        inverse = splu(stiffness.tocsc())
        return lambda rhs: inverse.solve(rhs/w)/w


class LabelRelaxation:
    def __init__(self, body, labels, order, epsilon):
        self.body = body
        self.axes = tuple(body.buffers['field_'+d] for d in ('x','t','a'))
        self.base = np.stack(np.meshgrid(*self.axes,indexing='ij'),axis=-1)
        self.epsilon = float(epsilon)
        # The two larger strengths are numerical continuation initializers.
        if self.epsilon not in (0.,1.e-10,1.e-8,1.e-6,1.e-4):
            raise ValueError('Unregistered harmonic gauge strength')
        if np.shape(labels) != self.base.shape:
            raise ValueError('Unexpected spatial-label shape')
        self.active = self.K = self.components = None
        self._update(order)

    def _update(self, order):
        if order.shape != self.base.shape[:-1] or not np.isfinite(order).all() or order.min()<0 or order.max()>1:
            raise ValueError('Order must be the bounded native x,t,a nodal field')
        active = cell_order(order)
        if self.active is not None and np.array_equal(active,self.active):
            return
        self.active = active
        self.stencil = np.zeros_like(order,dtype=bool)
        for bits in BITS:
            self.stencil[slices(bits,active.shape)] |= active
        self.K = matrix_for(self.body,self.axes,~active)
        self.components = None
        self.Kactive = None

    def _components(self):
        if self.components is not None:
            return self.components
        ids = np.arange(self.stencil.size).reshape(self.stencil.shape)
        comps = []
        for c in range(3):
            idx = ids.copy()
            held = np.zeros(self.stencil.shape,bool)
            held[:,[0,-1]] = True
            held[-1] = True
            if c==0:
                held[0] = True
            if c==1:
                idx[0] = ids[0,:,:1]
            if c==2:
                held[:,:,[0,-1]] = True
                mask = ~self.stencil[0] & ~held[0]
                idx[0][mask] = ids[1][mask]
            comps.append(Component(ids,idx,held,self.stencil,self.K))
        self.components = comps
        return comps

    def gauge_energy(self, labels):
        d = (labels-self.base).reshape(-1,3)
        return self.epsilon * .5 * float(np.sum(d*(self.K@d)))

    def _point(self, solid, labels, order):
        try:
            energy, elastic_grad, _, physical = solid.evaluate(labels,order,gradient=True)
        except InvalidSpatialCrystal as error:
            raise InfeasibleLabels(str(error), **error.diagnostics) from error
        # The native kernel bounds the orientation polynomial over every cell,
        # including faces and exact-zero-order cells, prior to integration.
        if physical['min_chart_Jinv_bound'] <= 0 or physical['min_Jinv_global'] <= 0:
            raise InfeasibleLabels('Invalid signed physical volume', **physical)
        # At determinant scales approaching the square root of roundoff,
        # independently test whole-cell orientation including cancellation
        # margins before permitting a step to cross a numerical fold.
        if physical['min_chart_Jinv_bound'] < np.sqrt(UROUND):
            chart_orientation(self.axes,labels,self.active)
        if not np.isfinite(energy) or not np.isfinite(elastic_grad).all():
            raise ArithmeticError('Nonfinite elastic energy or derivative at oriented labels')
        delta = (labels-self.base).reshape(-1,3)
        gg = self.K@delta
        gauge = self.epsilon * .5 * float(np.sum(delta*gg))
        force = []; gauge_force = []; dead = []; diag_elastic = []; projected = []
        for c, cm in enumerate(self._components()):
            eg = cm.sum(elastic_grad[...,c]); unit = cm.sum(gg[:,c]); g = eg + self.epsilon*unit
            force.append(absmax(g[cm.free])); gauge_force.append(absmax(unit[cm.dead]))
            dead.append(absmax(eg[cm.dead])); diag_elastic.append(absmax(eg[cm.free]))
            projected.append(g[cm.free])
        grad = np.concatenate(projected)
        residual = max(force); raw_elastic = max(diag_elastic); harmonic = max(gauge_force)
        stationary = residual<=FORCE_TOL and raw_elastic<=FORCE_TOL and harmonic<=HARMONIC_TOL and max(dead)<=1.e-12
        info = dict(physical,stationary=bool(stationary),elastic_energy=float(energy),gauge_energy=gauge,
                epsilon=self.epsilon,gauge_to_elastic=gauge/energy if energy else None,
                max_free_label_force=residual,max_free_elastic_force=raw_elastic,
                max_unit_harmonic_residual=harmonic,max_dead_elastic_force=max(dead),
                elastic_force_by_component=diag_elastic,free_force_by_component=force,
                harmonic_free_by_component=gauge_force)
        return float(energy)+gauge, grad, info

    def measure(self, solid, labels, order):
        self._update(order)
        delta = (labels-self.base).reshape(-1,3)
        reg = 0.
        for c,cm in enumerate(self._components()):
            d = delta[:,c]
            reg = max(reg, absmax(d[cm.held[cm.groups]]),
                      absmax(d-d[cm.representative][cm.groups]))
        # Following an order update an axis node can become strictly melted.
        # Its regular harmonic continuation may then require another sweep.
        # A conflict on the active trace is a distinct inadmissible state.
        self._pack(labels)
        info = self._point(solid,labels,order)[2]
        if reg>2.e-9:
            info['stationary'] = False
        return dict(info,max_registry_error=reg,independent_chart=chart_orientation(self.axes,labels,self.active))

    def _pack(self, labels):
        delta = labels-self.base
        return np.concatenate([cm.pack(delta[...,c],self.stencil) for c,cm in enumerate(self._components())])

    def _lift(self, x):
        labels = self.base.copy(); k = 0
        for c,cm in enumerate(self.components):
            n = np.count_nonzero(cm.free)
            labels[...,c] += cm.lift(x[k:k+n]).reshape(self.stencil.shape)
            k += n
        return labels

    def _precondition(self, order):
        if self.Kactive is None:
            self.Kactive = matrix_for(self.body,self.axes,self.active)
        chi_mean = np.lib.stride_tricks.sliding_window_view(order,(2,2,2)).mean(axis=(-3,-2,-1))
        scales = np.zeros_like(order); n = np.zeros_like(order)
        for bits in BITS:
            s = slices(bits,self.active.shape)
            scales[s] += chi_mean; n[s] += self.active
        scales = np.sqrt(np.maximum(scales/np.maximum(1.,n),1.e-4))
        x,t,_ = np.moveaxis(self.base,-1,0)
        R = self.body.sh.radius(t*self.body.sh.length)[0]
        length = (R, self.body.sh.length, np.maximum(x*R,.5*np.diff(self.axes[0])[0]*R))
        pre = [cm.preconditioner(self.Kactive,scales*length[c]) for c,cm in enumerate(self.components)]
        sections = np.cumsum([np.count_nonzero(cm.free) for cm in self.components])[:-1]
        return lambda v: np.concatenate([fn(part) for fn,part in zip(pre,np.split(v,sections))])

    def solve(self, solid, labels, order, log):
        self._update(order)
        x = self._pack(labels)
        initial = labels
        labels = self._lift(x)
        try:
            energy, grad, info = self._point(solid,labels,order)
        except InfeasibleLabels as error:
            raise LabelStall('Harmonic extension of the new active trace is infeasible',initial,error.witness) from error
        difference = absmax((labels-initial)[self.stencil])
        if difference>2.e-9:
            raise LabelStall('Harmonic extension altered the active trace', initial,dict(info,max_active_change=difference))
        log(dict(event='label_start',max_active_change=difference,**info))
        precondition = None; history = []; alphas = []; iteration = 0
        while not info['stationary']:
            if info['max_unit_harmonic_residual']>HARMONIC_TOL or info['max_dead_elastic_force']>1.e-12:
                raise LabelStall('Eliminated melted DOF has a resolved residual',labels,info)
            if precondition is None:
                precondition = self._precondition(order)
            v = grad.copy(); weights = []
            for s,y,inv in reversed(history):
                a = inv*np.dot(s,v); v -= a*y; weights.append(a)
            v = precondition(v)
            for (s,y,inv),a in zip(history,reversed(weights)):
                v += (a-inv*np.dot(y,v))*s
            direction = -v
            derivative = float(np.dot(grad,direction))
            if not np.isfinite(derivative) or derivative>=0:
                history.clear(); direction = -precondition(grad); derivative = float(np.dot(grad,direction))
            if not derivative<0:
                raise LabelStall('No descending direction at unresolved label residual',labels,info)
            alpha = 1.; rejected = None
            while True:
                trial = x+alpha*direction
                if np.array_equal(trial,x):
                    raise LabelStall('No representable descending oriented label step',labels,
                        dict(info,accepted_steps=len(alphas),trial_alpha=alpha,last_rejection=rejected))
                trial_labels = self._lift(trial)
                try:
                    new_energy,new_grad,new_info = self._point(solid,trial_labels,order)
                    if new_energy<energy and new_energy <= energy+1.e-4*alpha*derivative:
                        break
                    rejected = dict(reason='strict Armijo decrease failed',energy_delta=new_energy-energy)
                except InfeasibleLabels as error:
                    rejected = error.witness
                alpha *= .5
                if alpha<1.e-4 and -alpha*derivative<=16*UROUND*max(1.,abs(energy)):
                    raise LabelStall('Tiny oriented search interval with energetically unresolved descent',labels,
                        dict(info,accepted_steps=len(alphas),trial_alpha=alpha,last_rejection=rejected))
            iteration += 1
            diff = new_energy-energy
            s = trial-x; y = new_grad-grad; sy = float(np.dot(s,y))
            if sy>1.e-12*np.linalg.norm(s)*np.linalg.norm(y):
                history.append((s,y,1./sy)); history=history[-9:]
            if alpha<1.e-4 and -diff<=16*UROUND*max(1.,abs(energy)) and alphas and alphas[-1]<1.e-4:
                raise LabelStall('Repeated tiny label steps without resolved energy decrease',trial_labels,
                    dict(new_info,accepted_steps=iteration,last_alpha=alpha,energy_delta=diff,last_rejection=rejected))
            alphas.append(alpha)
            x,labels,energy,grad,info = trial,trial_labels,new_energy,new_grad,new_info
            log(dict(event='label_accepted',iteration=iteration,alpha=alpha,energy_delta=diff,
                     last_rejection=rejected,**info))
        return labels,dict(info,independent_chart=chart_orientation(self.axes,labels,self.active),accepted_steps=iteration,min_accepted_alpha=min(alphas,default=None),
                           last_accepted_alpha=alphas[-1] if alphas else None)
