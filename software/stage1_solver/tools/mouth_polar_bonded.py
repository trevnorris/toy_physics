"""Finite quadrupolar wall variation for the registered radius-four material.

The same spatial map follows the internal meridional displacement, transports
the Eulerian order chart, and extends across the wall to the fixed r=28 bath.
The body releases two meridional displacements per free node. The receiver's
lattice labels and both reference order fields are carried along this one
polar coordinate. They have no additional polar displacement/order variables.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_bonded_material_motion import Bonded
from mouth_internal_collar import interpolation
from mouth_internal_collar_check import read
from mouth_joined_angular_work import bump
from mouth_poly_shallow_release import CUT
from mouth_sealed_shape_relax import DELTA


BATH = 28.


def angles(n):
    if n < 4:
        raise ValueError('At least four angular quadrature points required')
    nodes, weights = leggauss(n)
    # Equatorial reflection and azimuthal symmetry reduce the solid-angle
    # mean to the unit interval in cos(theta). The parent measures include 4π.
    mu = (1 + nodes) / 2
    return (3 * mu**2 - 1) / 2, -3 * mu * np.sqrt(1 - mu**2), weights / 2


def geometry(shape, r, z, exterior):
    """Coefficients d,s,v of R_r=1+dP, R/r=1+sP, R_z=vP per unit amplitude."""
    b, slope = shape.radius(z)
    h, hz = bump(z, shape.length)
    if exterior:
        width = BATH - b
        k = (BATH-r) / width
        return -h / width, k*h/r, k*(hz + h*slope/width)
    k = (r/b)**3
    s = k*h/r
    return 3*s, s, k*(hz - 3*h*slope/b)


def require_map(d, s):
    # Check the whole angular range, including both endpoints.
    if any(np.any((a <= -1) | (a >= 2)) for a in (d, s)):
        raise ValueError('Folded quadrupolar map')


def crystal_delta(j, K, hoop, d, s, v, quad, derivatives=False):
    """Angular mean of W(Fpolar F)-W(F) for the full four-dimensional solid.

    Fpolar sends e_theta to s P_theta e_r + (1+sP)e_theta and e_phi to
    (1+sP)e_phi. The meridional row is (1+dP,vP). Using <P²>=1/5 and
    <P_theta²>=6/5 integrates the entire trace change exactly.
    """
    require_map(d, s)
    rr, rz, zz = K[:, 0, 0], K[:, 0, 1], K[:, 1, 1]
    delta = (d*d*rr + 2*d*v*rz + v*v*zz + 8*s*s*hoop**2) / 10
    if derivatives:
        # The interior chart has d=3s; derivatives below vary both together.
        vs = (3*d*rr + 3*v*rz + 8*s*hoop**2) / 5
        vv = (d*rz + v*zz) / 5
        vj = np.zeros_like(j)
    P, _, weights = quad
    for p, w in zip(P, weights):
        a, g = 1+d*p, 1+s*p
        # Product-minus-one without subtracting nearly equal Jacobians.
        df = (d+2*s)*p + (s*s+2*d*s)*p**2 + d*s*s*p**3
        delta += w*(-np.log1p(d*p) - 2*np.log1p(s*p)
                    + .25*j*df*(2*(j-1)+j*df))
        if derivatives:
            vs += w*p*(-3/a - 2/g + .5*j*(j*(1+df)-1)*(3*g*g+2*a*g))
            vj += w*.5*df*(j*(2+df)-1)
    return (delta, vs, vv, vj) if derivatives else delta


def order_at(body, r, z):
    shape = body.chart.shape
    b, slope = shape.radius(z)
    nt, nx = np.subtract(body.grid.shape, 1)
    t, x = z/body.L*nt, r/b*nx
    j, i = t.astype(int), x.astype(int)
    v, a = t-j, x-i
    c = body.grid
    lo, hi = c[j, i], c[j+1, i]
    dl, dh = c[j, i+1]-lo, c[j+1, i+1]-hi
    chi = (1-v)*(lo+a*dl) + v*(hi+a*dh)
    cr = nx*((1-v)*dl + v*dh)/b
    cw = nt*(hi-lo+a*(dh-dl))/body.L - r/b*slope*cr
    return chi, cr, cw


class PolarBonded(Bonded):
    def __init__(self, path, n, k, joined=None, *, shape=None, full_depth=True,
                 epsilon=0., angular_order=8):
        super().__init__(path, n, k, joined, shape=shape, full_depth=full_depth)
        self.epsilon = epsilon
        self.quad = angles(angular_order)

    def get_material_load(self, q, model):
        """Exact on the mobile terminal band, where the polar map is identity.

        The historical assembler returns entries on held nodes as well. Only
        its restriction to the free order coordinates is a polar adjoint.
        """
        from mouth_bonded_front_adjoint import get_material_load
        if self.epsilon != 0 and .62*self.L > self.L-model.h:
            raise ValueError('Polar support overlaps the mobile terminal order')
        return get_material_load(self, q, model)

    def energy(self, q, diagnose=False):
        base = super().energy(q, diagnose)
        if self.epsilon == 0:
            return base
        if self.front is not None:
            tz, _, dc = self.front
            overlap = (tz[:-1] < .62) & (tz[1:] > .06)
            if np.any(dc[:-1][overlap]) or np.any(dc[1:][overlap]):
                raise ValueError('Polar shape overlaps a nonzero terminal order correction')
        u = self.free.copy()
        u[self.mask] = q/self.scale
        nodes = u.reshape(-1, 2)[self.id]
        val = np.einsum('ij,ijk->ik', self.basis, nodes)
        r, z = self.r + val[:, 0], self.w + val[:, 1]
        active = (z > .06*self.L) & (z < .62*self.L)
        B = np.stack([np.einsum('ij,ijk->ik', d, nodes) for d in (self.dr, self.dw)],
                     axis=2) + np.eye(2)
        BK = B @ self.K
        K = BK @ B.transpose(0, 2, 1)
        hoop = self.h*r/self.r
        j = self.j*np.linalg.det(B)*(r/self.r)**2
        sh = self.chart.shape
        d, s, v = [self.epsilon*a for a in geometry(sh, r[active], z[active], False)]
        delta, vs, vv, vj = crystal_delta(j[active], K[active], hoop[active], d, s, v,
                                          self.quad, derivatives=True)
        chi, cr, cw = order_at(self, r[active], z[active])
        weights = self.measure[active]
        dE = weights @ (chi*delta)
        if diagnose:
            # The meridional solution has its own full-depth extrema. Rebuild
            # the polar eigensystem on quadrature points plus pole and equator.
            limits = [np.inf, -np.inf, np.inf, -np.inf]
            if np.any(~active):
                ev = np.linalg.eigvalsh(K[~active])
                limits = [j[~active].min(), j[~active].max(),
                          min(hoop[~active].min(), np.sqrt(ev[:, 0].min())),
                          max(hoop[~active].max(), np.sqrt(ev[:, 1].max()))]
            P, Pa, _ = self.quad
            for p, pa in zip(np.r_[P, -.5, 1.], np.r_[Pa, 0., 0.]):
                A = np.zeros((sum(active), 4, 4))
                A[:, 0, 0], A[:, 0, 1], A[:, 0, 2] = 1+d*p, v*p, s*pa
                A[:, 1, 1], A[:, 2, 2], A[:, 3, 3] = 1., 1+s*p, 1+s*p
                C = np.zeros_like(A)
                C[:, :2, :2] = K[active]
                C[:, 2, 2] = C[:, 3, 3] = hoop[active]**2
                ev = np.linalg.eigvalsh(A @ C @ A.transpose(0, 2, 1))
                jm = j[active]*(1+d*p)*(1+s*p)**2
                limits = [min(limits[0], jm.min()), max(limits[1], jm.max()),
                          min(limits[2], np.sqrt(ev[:, 0].min())),
                          max(limits[3], np.sqrt(ev[:, -1].max()))]
            inc = weights*chi*delta
            core = self.xx[active] < self.x[self.cut_node]
            base.update(energy=base['energy']+dE, min_J=limits[0], max_J=limits[1],
                        min_stretch=limits[2], max_stretch=limits[3],
                        core_polar_change=inc[core].sum(), inner_collar_polar_change=inc[~core].sum())
            return base

        r, z, j, hoop, B, BK = (a[active] for a in (r, z, j, hoop, B, BK))
        b, bw, bww = sh.profile(z/self.L)
        h, hw = bump(z, self.L)
        y = (z/self.L-.06)/.56
        hww = 32*(1 - 6*y + 6*y*y)/(.56*self.L)**2
        k = (r/b)**3
        vw = self.epsilon*k*(hww - 6*hw*bw/b + h*(12*(bw/b)**2 - 3*bww/b))
        # r,z dependence of Dpolar at the displaced material point. The order
        # chart is polar-transported, so χ,χ_r,χ_w here are its meridional pullback.
        wr = (2*s*vs + 3*v*vv + 2*j*vj + 8*s*s*hoop**2/5)/r
        wz = vs*v/r + vv*vw
        pos = np.column_stack((delta*cr + chi*wr, delta*cw + chi*wz))
        vec = np.column_stack((d, v))
        gradB = (vec[:, :, None]*vec[:, None, :] @ BK)/5
        gradB += (j*vj)[:, None, None]*np.linalg.inv(B).transpose(0, 2, 1)
        ids = self.id[active].ravel()
        for k in (0, 1):
            row = weights[:, None]*(pos[:, k, None]*self.basis[active] + chi[:, None]*(
                gradB[:, k, 0, None]*self.dr[active] + gradB[:, k, 1, None]*self.dw[active]))
            self.native[..., k] += np.bincount(ids, weights=row.ravel(),
                                              minlength=self.npts).reshape(self.mask.shape[:2])
        return base[0]+dE, self.native[self.mask]/self.scale


def order_sites(sh, grid, exterior):
    t, x, c = grid
    well = (c*(1-c))**2
    parts = []
    tail = []
    q, w = leggauss(3)
    dt, dx = np.diff(t)[:, None], np.diff(x)[None, :]
    L = sh.length
    for y, wy in zip((q+1)/2, w/2):
        z = (t[:-1, None]+dt*y)*L
        b, slope = sh.radius(z)
        lo = (1-y)*c[:-1]+y*c[1:]
        well_line = (1-y)*well[:-1]+y*well[1:]
        cd = np.diff(c, axis=0)/(dt*L)
        if exterior:
            # Continue the FEM's flat boundary column to the bath radius.
            # Its radial derivative vanishes, so only boundary volume moves.
            _, shift, _ = geometry(sh, np.full_like(z, CUT), z, True)
            flat_density = (3/DELTA)*(DELTA**2*cd[:, -1]**2+well_line[:, -1])
            tail.append((shift.ravel()*CUT, (4*np.pi*L*dt.ravel()*wy)*flat_density))
        for a, wa in zip((q+1)/2, w/2):
            xx = x[None, :-1]+a*dx
            width = CUT-b if exterior else b
            r = b+xx*width if exterior else xx*b
            jac = 4*np.pi*r*r*width*dx*L*dt*wa*wy
            gr = np.diff(lo, axis=1)/(dx*width)
            gz = (1-a)*cd[:, :-1]+a*cd[:, 1:] - slope*(1-xx if exterior else xx)*gr
            potential = (1-a)*well_line[:, :-1]+a*well_line[:, 1:]
            d, s, v = geometry(sh, r, z, exterior)
            parts.append(tuple(u.ravel() for u in
                               np.broadcast_arrays(jac*3/DELTA, gr, gz, potential, d, s, v)))
    data = tuple(np.concatenate(u) for u in zip(*parts))
    return data, tuple(np.concatenate(u) for u in zip(*tail))


def order_delta(data, epsilon, quad):
    wt, gr, gz, potential, d, s, v = data
    d, s, v = epsilon*d, epsilon*s, epsilon*v
    require_map(d, s)
    P, Pa, weights = quad
    total = np.zeros_like(gr)
    for p, pa, w in zip(P, Pa, weights):
        a, g = 1+d*p, 1+s*p
        f = a*g*g
        df = (d+2*s)*p+(s*s+2*d*s)*p*p+d*s*s*p**3
        # grad chi moves by Dpolar^{-T}, with polar shear even for initial
        # meridional chi. Difference form avoids cancelling total field energies.
        radial = gr/a
        az = -v*p*radial
        angular = -s*pa*radial/g
        dr2 = -d*p*(2+d*p)*radial**2
        dg2 = dr2 + angular**2 + az*(2*gz+az)
        total += w*(DELTA**2*f*dg2+df*(DELTA**2*(gr*gr+gz*gz)+potential))
    return wt@total


def outer_order_metric(shape, r, z, epsilon, angular_order=8):
    """Return <J D^{-1}D^{-T}>_{rr,ww,rw} and <J> for the exterior side.

    These multiply physical meridional Q1 gradients on the undeformed chart.
    The third term multiplies nr⊗nw+nw⊗nr. The last scales the mass lump.
    """
    d, s, v = (epsilon*c for c in geometry(shape, r, z, True))
    require_map(d, s)
    rr = np.zeros(np.broadcast_shapes(np.shape(r), np.shape(z)))
    rw = np.zeros_like(rr)
    volume = np.zeros_like(rr)
    for p, pa, weight in zip(*angles(angular_order)):
        a, g = 1+d*p, 1+s*p
        f = a*g*g
        rr += weight*f*(1+(s*pa/g)**2+(v*p)**2)/a**2
        rw -= weight*f*v*p/a
        volume += weight*f
    return rr, volume, rw, volume


def exterior_crystal_sites(pair, labels):
    sh = pair.shape
    ext = ~pair.inside
    nodes = labels.reshape(-1, 2)[pair.id[ext]]
    S = sh.neck+np.sum(pair.B[ext]*nodes[:, :, 0], axis=1)
    inv = np.stack([np.einsum('qi,qik->qk', g[ext], nodes) for g in (pair.E, pair.T)], axis=-1)
    D = np.linalg.det(inv)
    if np.min(D) <= 0 or np.min(S) <= 0:
        raise ValueError('Folded exterior material coordinates')
    spatial = np.zeros_like(inv)
    spatial[:, 0, 0], spatial[:, 0, 1], spatial[:, 1, 1] = pair.re[ext], pair.rt[ext], sh.length
    F = spatial @ np.linalg.inv(inv)
    r, z = pair.r[ext], pair.z[ext]*sh.length
    hoop = r/S
    j = -pair.re[ext]*sh.length*hoop**2/D
    wt = pair.qw[ext]*S**2*D
    return wt, j, F@F.transpose(0, 2, 1), hoop, *geometry(sh, r, z, True)


def outer_crystal_source_delta(pair, labels, wa, xout, epsilon, angular_order=8):
    """Deposit ∂(E_polar−E_spherical)/∂χ on native receiver nodes (wa physical)."""
    zero = np.zeros((len(wa), len(xout)))
    if epsilon == 0:
        return zero
    wt, j, K, hoop, d, s, v = exterior_crystal_sites(pair, labels)
    change = wt*crystal_delta(j, K, hoop, epsilon*d, epsilon*s, epsilon*v, angles(angular_order))
    ext = ~pair.inside
    _, ids, basis = interpolation(pair.z[ext], pair.coordinate[ext],
                                  (wa/pair.shape.length, xout, zero))
    return np.bincount(ids.ravel(), weights=(basis*change[:, None]).ravel(),
                       minlength=zero.size).reshape(zero.shape)


class StaticPolar:
    """Changes relative to the spherical static energy, excluding PolarBonded.

    The duct and receiving order energies use the existing FEM gradient and
    lumped well on their own meshes. The feed and both post-terminal blocks
    have identically zero shape increments. The physical receiver is continued
    from CUT to BATH with its fixed flat trace, using its nodal FEM interpolant
    so the continuation joins continuously. Its constant base energy drops out.
    """
    def __init__(self, path, *, shape=None, outer_order=None, angular_order=8):
        pair, labels, inner, outer = read(path)
        if shape is not None:
            pair = type(pair)(shape, pair.co, pair.na, pair.nt)
        sh = pair.shape
        if outer_order is not None:
            if isinstance(outer_order, tuple):
                wa, x, c = outer_order
                outer = (wa/sh.length, x, c)
            else:
                outer = (*outer[:2], outer_order)
        for t, x, c in (inner, outer):
            if c.shape != (len(t), len(x)) or not np.allclose(
                    [t[0], t[-1], x[0], x[-1]], [0, 1, 0, 1]):
                raise ValueError('Invalid reference order chart')
        self.quad = angles(angular_order)
        self.inner, _ = order_sites(sh, inner, False)
        self.outer, self.tail = order_sites(sh, outer, True)
        ext = ~pair.inside
        wt, *fields = exterior_crystal_sites(pair, labels)
        chi = interpolation(pair.z[ext], pair.coordinate[ext], outer)[0]
        self.elastic = (wt*chi, *fields)

    def inner_energy(self, epsilon):
        return float(order_delta(self.inner, epsilon, self.quad))

    def tail_energy(self, epsilon):
        shift, density_weight = self.tail
        # The changing flat-shell volume is -<((CUT+e*shift*P)^3-CUT^3)>/3.
        # <P>=0, <P²>=1/5, <P³>=2/35.
        move = epsilon*shift
        return float(-density_weight @ (CUT*move**2/5+2*move**3/105))

    def terms(self, epsilon):
        weight, j, K, hoop, d, s, v = self.elastic
        outer_W = weight @ crystal_delta(j, K, hoop, epsilon*d, epsilon*s, epsilon*v, self.quad)
        outside = order_delta(self.outer, epsilon, self.quad)
        return dict(inner_order=self.inner_energy(epsilon), outer_order=float(outside+self.tail_energy(epsilon)),
                    outer_crystal=float(outer_W))
