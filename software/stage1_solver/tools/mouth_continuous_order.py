"""Absolute crystalline order on the whole held R4 field chart.

The carrier keeps its six-block slit topology. Order has a single spatial
node wherever the banks coincide. The only prescribed order is the distant
flat material trace and the bulk roof. At the lower slab face the variational
natural condition applies. No internal interface work or new solid is added.
"""
import numpy as np
from scipy import sparse
from scipy.interpolate import RegularGridInterpolator
from scipy.sparse.linalg import cg, LinearOperator, spsolve
from scipy.special import expit

from mouth_breathing_mobile_patch import KAPPA, ORDER_INERTIA as I, WELL, WIDTH
from mouth_exchange_carrier import BIAS
from mouth_free_shape_solid import Solid
from mouth_internal_collar import CollarPair
from mouth_last_rings_fluid import RimCarrier
from mouth_polar_bonded import exterior_crystal_sites
from mouth_sealed_shape_relax import Meridian


def norm(a):
    return float(np.max(np.abs(a), initial=0.))


def mat(A, weight):
    return A.T @ A.multiply(weight[:, None])


def topology(m):
    """Identify coinciding physical nodes, retaining separate carrier banks."""
    pairs = np.zeros(m.size, dtype=np.int64)-1
    zall = np.r_[m.z[0], m.z[1][1:], m.z[2][1:]]
    for (slab, side), ids in m.ids.items():
        ir, angle, iz = np.indices(ids.shape)
        ir += side*(len(m.rr[0])-1)
        angle[ir == 0] = 0
        iz += sum(len(z)-1 for z in m.z[:slab])
        keys = (ir*len(m.a)+angle)*len(zall)+iz
        if np.any((pairs[ids] != -1) & (pairs[ids] != keys)):
            raise ArithmeticError('Inconsistent physical charts')
        pairs[ids] = keys
    distinct, ix = np.unique(pairs, return_inverse=True)
    if distinct[0] < 0:
        raise ArithmeticError('Orphan carrier node')
    return ix, zall[distinct % len(zall)]


def exterior_load(base, shape, shift, m, angular_degree=4):
    """Deposit the full W of the existing exterior label triangles on order Q1.

    The material quadrature is the CollarPair triangle rule. Integrate each FE
    angular interval separately. Its polar pullback uses the complete 4D F.
    """
    pair = CollarPair(shape, base.co, base.na, base.nt)
    wt, jac, K, hoop, d, s, v = exterior_crystal_sites(pair, base.labels)
    ext = ~pair.inside
    z, r = pair.z[ext]*shape.length, pair.r[ext]
    b = shape.radius(z)[0]
    x = (r-b)/(28.-b)
    ids = m.ids[1, 1]
    def locate(nodes, values):
        i = np.clip(np.searchsorted(nodes, values)-1, 0, len(nodes)-2)
        f = (values-nodes[i])/(nodes[i+1]-nodes[i])
        if f.min() < -1.e-12 or f.max() > 1.+1.e-12:
            raise ValueError('Exterior crystal left its field chart')
        return i, f
    ix, fx = locate(m.rr[1], x)
    iz, fz = locate(m.z[1], z)
    q, qw = np.polynomial.legendre.leggauss(angular_degree)
    out = np.zeros(m.size)
    for k, (angle, da) in enumerate(zip(m.a, np.diff(m.a))):
        for u, weight in zip((q+1)/2, qw/2):
            a = angle+u*da
            P = (3*np.cos(a)**2-1)/2
            Ptheta = -3*np.sin(a)*np.cos(a)
            ar, arad, az = 1+shift*d*P, 1+shift*s*P, shift*v*P
            J = jac * ar*arad**2
            tr = (ar**2*K[:, 0, 0]+2*ar*az*K[:, 0, 1]+(1+az**2)*K[:, 1, 1]
                  + hoop**2*((shift*s*Ptheta)**2+2*arad**2))
            W = (tr-4)/2-np.log(J)+(J-1)**2/4
            deposit = wt*W*(da*weight*np.sin(a))
            for di, fi in ((0, 1-fx), (1, fx)):
                for dj, fj in ((0, 1-fz), (1, fz)):
                    for dk, fk in ((0, 1-u), (1, u)):
                        out += np.bincount(ids[ix+di, k+dk, iz+dj], deposit*fi*fj*fk,
                                           minlength=m.size)
    return out


def seed_order(m, base_order, profile, old):
    """Project the old charts to one initial trace; this trace is then free."""
    number, _ = topology(m)
    votes = np.zeros(number.max()+1)
    total = np.zeros_like(votes)
    interpolation = RegularGridInterpolator(tuple(np.linspace(0, 1, n) for n in base_order.shape), base_order)
    for (slab, side), ids in m.ids.items():
        x, a, z = np.meshgrid(m.rr[side], m.a, m.z[slab], indexing='ij')
        if slab == 0:
            # Extend the prepared flat half-front under the intact face. Away
            # from the entrance in the ordered chamber the seed approaches one.
            b = m.shape.radius(z)[0]
            r = b+x*(28-b) if side else x*b
            c = expit(np.hypot(z, np.maximum(4-r, 0.))/WIDTH)
        elif side == 1:
            b = m.shape.radius(z)[0]
            c = profile.outside(b+x*(28.-b), z, slab == 2)
        elif slab == 1:
            c = interpolation((z/m.shape.length, x))
        else:
            c = profile(x, z-m.shape.length)
        ix = number[ids].ravel()
        total += np.bincount(ix, (c+old[2, ids]).ravel(), minlength=len(total))
        votes += np.bincount(ix, minlength=len(total))
    new = old.copy()
    new[2] = (total/votes)[number]
    return new


class ContinuousTrial:
    def __init__(self, base, state, x, n, k, angles, degree, subdivisions, *, precondition=True, axes=None):
        self.x = np.array(x)
        self.shape = Meridian(4., base.sh.params+self.x[:4]*[1., .25, .25, .25],
                              extra=base.sh.extra+.25*self.x[4])
        self.polar_shift = float(self.x[5])
        self.model = m = RimCarrier(self.shape, base.c, 1, base.grid, polar_shift=self.polar_shift,
                                    precondition=precondition, axes=axes)
        # In this chart u[2] is the absolute Q1 order; no prepared gradient or
        # well remains anywhere in the kinetic, gradient or potential term.
        m.ref = np.zeros_like(m.ref)
        self.b = Solid(state, self.shape, self.polar_shift, n, k, angles, m,
                       degree=degree, subdivisions=subdivisions)
        self.b.buffers['base'].fill(0.)
        self.number, self.depth = topology(m)
        count = self.number.max()+1
        self.pull = sparse.csr_matrix((np.ones(m.size), (np.arange(m.size), self.number)),
                                      shape=(m.size, count))
        self.V = m.V @ self.pull
        self.G = [g @ self.pull for g in m.G]
        self.fixed_order = np.unique(self.number[np.r_[m.feed, m.bulk]])
        self.port_value = expit(-self.depth[self.fixed_order]/WIDTH)
        roof = np.unique(self.number[np.r_[m.ids[2,0][:,:,-1].ravel(), m.ids[2,1][:,:,-1].ravel()]])
        self.port_value[np.isin(self.fixed_order, roof)] = 0.
        self.free = np.ones(count, bool)
        self.free[self.fixed_order] = False
        self.ext = exterior_load(base, self.shape, self.polar_shift, m)
        self.first = np.unique(self.number, return_index=True)[1]
        self.K = None
        if precondition:
            self.K = sum(mat(g, KAPPA*m.hq) for g in self.G).tocsr()
            metric = self.K + mat(self.V, 2*WELL*m.hq)
            self.diag = metric.diagonal()
            m.solvers[2] = None
            self.Kv = (m.B+sum(mat(g, m.hq*m.D[:,4]**2) for g in m.G)).tocsr()
        self.f = self

    def transported(self, u):
        u = u.copy()
        c = u[2, self.first].copy()
        if norm(c[self.number]-u[2]) > 2.e-13:
            raise ValueError('Split absolute material trace')
        c[self.fixed_order] = self.port_value
        u[2] = c[self.number]
        return u

    def fixed(self, q, u):
        u = self.transported(u)
        m = self.model
        e, _, _, _ = self.b.evaluate(q, u[2], load=False, gradient=False)
        m.material_load = self.ext
        return float(e+m.evaluate(u)[0])

    def projected(self, c, g):
        # Diagonally scaled feasible gradient mapping; a zero at either bound
        # is counted only for the correct sign of the discrete obstacle force.
        return c[self.free]-np.clip(c[self.free]-g[self.free]/self.diag[self.free], 0., 1.)

    def order_terms(self, c, n, gp, load, rest):
        m = self.model
        C = self.V @ c
        w = m.hq
        gc = [g @ c for g in self.G]
        ec = (.5*c@(self.K@c) + w@(WELL*(C*(1-C))**2) + load@c)
        grad = self.K@c+self.V.T@(2*WELL*w*C*(1-C)*(1-2*C))+load
        if rest:
            return ec, grad, C
        den = n+I*sum(d*d for d in gc)
        dot = sum(a*b for a, b in zip(gc,gp))
        pi = I*n*dot/den
        # Constant 1/2 n |grad p|² is independent of order.
        ec -= .5*(w*pi)@dot
        for g, a, b in zip(self.G, gc, gp):
            grad -= g.T@(w*pi*(b-I*dot/den*a))
        return ec, grad, C

    def solve_order(self, u, load, rest, tol, log=None):
        m = self.model
        c = u[2,self.first].copy()
        n = (m.D[:,4]*(m.V @ u[0]))**2
        gp = [] if rest else [g@u[1] for g in m.G]
        e, g, C = self.order_terms(c,n,gp,load,rest)
        iteration = 0
        while True:
            projected = self.projected(c,g)
            resid = norm(projected)
            if resid < tol:
                break
            # Reduced-space Newton direction. In spinodal regions its Hessian
            # may have a negative eigenvalue. Shift the well curvature by its
            # full negative lower bound if the projected proposal cannot descend.
            free = self.free & ~(((c==0)&(g>=0)) | ((c==1)&(g<=0)))
            local = self.V[:, free]
            K = self.K[free][:, free]
            accepted = False
            for damping in (0., 1.10, 3.0, None):
                if damping is None:
                    direction = -g[free]/self.diag[free]
                else:
                    h = K + mat(local, WELL*m.hq*(2-12*C*(1-C)+damping))
                    direction = spsolve(h.tocsc(), -g[free], permc_spec='MMD_AT_PLUS_A')
                s = min(1., .35/max(1.e-16, norm(direction)))
                for _ in range(42):
                    new = c.copy()
                    new[free] = np.clip(c[free]+s*direction,0.,1.)
                    dx = new-c
                    slope = g@dx
                    if slope < 0:
                        en, gn, Cn = self.order_terms(new,n,gp,load,rest)
                        if abs(slope) < 1.e-8*max(1.,abs(e)):
                            _, gm, _ = self.order_terms((c+new)/2,n,gp,load,rest)
                            delta = ((g+gn+4*gm)@dx)/6
                            ok = delta < 1.e-4*slope and en<=e+32*np.spacing(abs(e))
                        else:
                            ok = en-e < 1.e-4*slope
                        if ok:
                            accepted = True
                            break
                    s *= .5
                if accepted:
                    break
            if not accepted:
                raise ArithmeticError(f'Bounded order has no energy-lowering step, residual {resid}')
            c,e,g,C = new,en,gn,Cn
            iteration += 1
            if log:
                log(iteration, e, resid, c, s, damping)
        u[2] = c[self.number]
        return iteration, resid

    def carrier(self, u, rest, tol):
        m = self.model
        F, w = m.D[:,4], m.hq
        gc = [g@u[2] for g in m.G]
        gc2 = sum(d*d for d in gc)
        free = m.masks[0]
        for iteration in range(100):
            t = m.V@u[0]
            n = (F*t)**2
            den = n+I*gc2
            if rest:
                u[1].fill(0.)
                p_error = 0.
            else:
                def phase(p):
                    gs = [g@p for g in m.G]
                    d = sum(a*b for a,b in zip(gs,gc))
                    return sum(g.T@(w*n*(v-I*d/den*c)) for g,v,c in zip(m.G,gs,gc))
                rhs = -phase(u[1])[free]
                p_error = norm(m.solvers[1](rhs))
                if p_error>=tol:
                    def mult(step):
                        p = np.zeros(m.size)
                        p[free]=step
                        return phase(p)[free]
                    size = (free.sum(),)*2
                    dx, status = cg(LinearOperator(size, mult, dtype=float), rhs, rtol=2e-9, atol=0.,
                                    M=LinearOperator(size,m.solvers[1], dtype=float))
                    if status:
                        raise ArithmeticError(f'Phase linear solve {status}')
                    u[1,free] += dx
            gp = [g@u[1] for g in m.G]
            dot = sum(a*b for a,b in zip(gc,gp))
            speed = sum(p*p for p in gp)-I**2*gc2*dot**2/den**2
            g = self.Kv@u[0]+m.V.T@(w*(.5*F**10*t*(t**8-1)+F**2*t*speed))
            a_error = norm(m.solvers[0](g[free]))
            if max(a_error,p_error) < tol:
                return iteration, max(a_error,p_error)
            if a_error >= tol:
                curvature = (.5*F**10*(9*t**8-1)+F**2*speed
                             +4*F**4*t*t*I**2*gc2*dot**2/den**3)
                H = self.Kv[free][:,free]+mat(m.V[:,free],w*curvature)
                step, status = cg(H,-g[free],rtol=1e-9,atol=0.,
                                  M=LinearOperator(H.shape, m.solvers[0], dtype=float))
                if status:
                    raise ArithmeticError(f'Amplitude linear solve {status}')
                scale = min(1., .1/max(1.e-15,norm(step)))
                u[0,free] += scale*step
        raise ArithmeticError('Fixed-order carrier failed to settle')

    def match(self, q, u, rest=False, log=None):
        if self.K is None:
            raise ValueError('Field match requires preconditioning')
        u = self.transported(u)
        m = self.model
        u[0,~m.masks[0]]=1.
        u[1,m.feed] = 0.
        u[1,m.bulk] = 0. if rest else BIAS
        if rest:
            u[1].fill(0.)
        _, _, inner, _ = self.b.evaluate(q, gradient=False)
        load = inner+self.ext
        joined = self.pull.T @ load
        order_steps, carrier_steps = [], []
        for sweep in range(80):
            count, _ = self.solve_order(u, joined, rest, 2.e-9, log)
            order_steps.append(count)
            count, cr = self.carrier(u, rest, 8.e-8)
            carrier_steps.append(count)
            m.material_load = load
            energy, g = m.evaluate(u)
            gres = self.pull.T@g[2]
            pr = norm(self.projected(u[2,self.first],gres))
            if max(cr,pr)<2.e-7:
                break
        else:
            raise ArithmeticError('Coupled order and carrier residual unresolved')
        solid_energy, gq, _, _ = self.b.evaluate(q,u[2],load=False)
        work_error = solid_energy-inner@u[2]
        if abs(work_error) > 2e-8*max(1.,abs(solid_energy)):
            raise ArithmeticError(f'Material deposition mismatch {work_error}')
        meta = dict(field=dict(stationarity=max(cr,pr), carrier_metric_residual=cr,
                               order_projected_residual=pr, order_steps=order_steps, carrier_steps=carrier_steps),
                    solid_energy=solid_energy, exchange_error=work_error, energy=float(energy))
        return u, float(energy), gq, load, meta
