"""Finite quadrupolar displacement of the six-domain radius-four carrier and front.

The polar coordinate and all node labels follow the radial pullback. The
physical azimuthal vortex factor is evaluated at the moved cylindrical radius.
The mobile order band and both bath boundaries keep their original coordinates.
Use a material quadrature deposition for material_load. The native static order
and material energies enter separately in the whole-mouth action.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from mouth_bonded_front_adjoint import Coupled
from mouth_exchange_carrier import BIAS
from mouth_joined_angular_work import bump
from mouth_joined_front_stationarity import solve as solve_stationary
from mouth_wound_carrier import straight_vortex


class PolarCoupled(Coupled):
    def __init__(self, shape, order, winding, grid, *, polar_shift=0.,
                 terminal_profile=None, precondition=True):
        super().__init__(shape, order, winding, grid, terminal_profile=terminal_profile,
                         precondition=precondition)
        self.warm_path = None
        self.material_load = np.zeros(self.size)
        self._undeformed = (self.D, self.G, self.ref, self.B)
        self._vortex = straight_vortex(winding, 28.)
        _, r, angle, z, _, _, _, _, live, _, _ = self.D.T
        h, dh = bump(z, shape.length)
        self._active = h > 0
        if np.any(self._active & (live != 0)) or .62 * shape.length >= shape.length - self.h:
            raise ValueError('The side mode requires a separate, unperturbed terminal band')
        r, angle, z, h, dh = (a[self._active] for a in (r, angle, z, h, dh))
        b, slope = shape.radius(z)
        inside = r < b
        width = np.where(inside, b, 28. - b)
        x = np.where(inside, r / b, (r - b) / width)
        k = np.where(inside, x**3, 1. - x)
        kx = np.where(inside, 3. * x**2, -1.)
        xw = -slope * np.where(inside, x, 1. - x) / width
        c = np.cos(angle)
        P = (3.*c*c - 1.) / 2.
        self._mode = (k*h*P, kx*h*P/width, -3.*k*h*c*np.sin(angle),
                      (kx*xw*h + k*dh)*P)
        self.polar_shift = 0.
        self.set_shift(polar_shift)

    def _side_matrix(self, shift):
        # Integrated oriented area vector on r=R(theta,z) is proportional to
        # (1, -R_theta/R, -R_z). The straight vortex has zero axial gradient.
        # Sealed traces have opposite normals and independent amplitudes.
        quad = (1. + np.array([-1., 1.])/np.sqrt(3.))/2.
        dz = np.diff(self.z[1])[None, :]
        da = np.diff(self.a)[:, None]
        matrix = sparse.csr_matrix((self.size, self.size))
        for a in quad:
            angle = self.a[:-1, None] + a*da
            sn, co = np.sin(angle), np.cos(angle)
            for t in quad:
                z = self.z[1][None, :-1] + t*dz
                b, _ = self.shape.radius(z)
                h, _ = bump(z, self.shape.length)
                radius = b + shift*h*(3.*co*co-1.)/2.
                rtheta = -3.*shift*h*co*sn
                f, df = self._vortex(radius*sn)
                weight = np.pi*da*dz*sn*radius*(radius*sn-rtheta*co)*f*df
                v = np.kron([1.-a, a], [1.-t, t])
                local = weight[..., None, None]*np.outer(v, v)
                for side in (0, 1):
                    face = self.ids[1, side][0 if side else -1]
                    nodes = np.stack((face[:-1, :-1], face[:-1, 1:],
                                      face[1:, :-1], face[1:, 1:]), axis=-1)
                    rows = np.broadcast_to(nodes[..., :, None], local.shape).ravel()
                    cols = np.broadcast_to(nodes[..., None, :], local.shape).ravel()
                    matrix += sparse.coo_matrix(((1.-2.*side)*local.ravel(), (rows, cols)),
                                                shape=matrix.shape).tocsr()
        return matrix

    def set_shift(self, shift):
        """Set an absolute displacement on this meridian. Keep its SPD solver metric."""
        shift = float(shift)
        if not np.isfinite(shift):
            raise ValueError('Nonfinite polar displacement')
        if shift == self.polar_shift:
            return
        base, gradient, reference, _ = self._undeformed
        if shift == 0.:
            self.D, self.G, self.ref, self.B = self._undeformed
        else:
            take = self._active
            r, theta = base[take, 1], base[take, 2]
            mode, dr, da, dw = self._mode
            radius = r + shift*mode
            radial_metric = 1. + shift*dr
            if np.min(radial_metric) <= 0 or np.min(radius) <= 0:
                raise ValueError('Nonpositive quadrature Jacobian on the polar chart')
            # At each sampled depth the extrema of dR/dr over x and angle
            # occur at the side or axis and at a pole or the equator.
            depths = np.unique(base[take, 3])
            b, _ = self.shape.radius(depths)
            h, _ = bump(depths, self.shape.length)
            for wall_shift in (shift*h, -.5*shift*h):
                if np.any(1.+3.*wall_shift/b <= 0.) or np.any(1.-wall_shift/(28.-b) <= 0.):
                    raise ValueError('Folded radial chart at the wall')
            scale_r, scale_a = 1./radial_metric, r/radius
            shear_a, shear_w = -shift*da*scale_r/radius, -shift*dw*scale_r
            det = (radius/r)**2 * radial_metric
            data = base.copy()
            data[take, 0] *= det
            data[take, 9] *= det
            data[take, 1] = radius
            data[take, 6] = radius*np.sin(theta)
            data[take, 4], data[take, 5] = self._vortex(data[take, 6])
            ref = reference.copy()
            ref[take, 1] *= scale_r
            ref[take, 2] = scale_a*reference[take, 2] + shear_a*reference[take, 1]
            ref[take, 3] += shear_w*reference[take, 1]
            def full(active, rest):
                factors = np.full(len(base), rest)
                factors[take] = active
                return factors[:, None]
            gr, ga, gw = gradient
            G = [gr.multiply(full(scale_r, 1.)).tocsr(),
                 (ga.multiply(full(scale_a, 1.))+gr.multiply(full(shear_a, 0.))).tocsr(),
                 (gw+gr.multiply(full(shear_w, 0.))).tocsr()]
            self.D, self.G, self.ref, self.B = data, G, ref, self._side_matrix(shift)
        self.hq = self.D[:, 0]
        self.polar_shift = shift

    def solve(self, rest=False, initial=None):
        """Rematch both reservoirs and the front with the held shared material load."""
        if not self.precondition:
            raise ValueError('A preconditioned instance is required for field relaxation')
        if initial is not None:
            if initial.shape != self.seed.shape:
                raise ValueError('Incompatible warm field topology')
            self.seed = initial.copy()
        self.seed[0, ~self.masks[0]] = 1.
        self.seed[1, self.feed] = 0.
        self.seed[1, self.bulk] = 0. if rest else BIAS
        return solve_stationary(self, rest)


def check(path, joint):
    """Compare the new finite chart to an independent signed-work calculation."""
    from mouth_bonded_front_adjoint import get_material_load
    from mouth_bonded_material_motion import Bonded
    from mouth_joined_angular_work import deformed_work
    from mouth_joined_front_audit import Front
    from mouth_poly_shallow_release import CUT
    from mouth_supercore_joined_front import check_guard, emit, file_sha
    from mouth_supercore_terminal_order import terminal

    check_guard(path)
    check_guard(joint)
    rows = list(map(json.loads, (joint.parent/'stdout').read_text().splitlines()))
    prep, outcome = rows[0], rows[-1]
    if (outcome['kind'] != 'held_equilibrium' or outcome['sha256'] != file_sha(joint)
            or prep['source_sha256'] != file_sha(path) or prep['winding'] != 1):
        raise ValueError('Unverified radius-four wound preparation')
    base = Front(path, tuple(prep['grid']), 1)
    with np.load(joint) as fields:
        u, d = fields['field'], fields['displacement']
    body = Bonded(path, (d.shape[1]-1)//3, d.shape[0]-1, full_depth=True)
    n = base.order_mesh
    t, x = np.linspace(0., 1., 2*n+1), np.linspace(0., 1., n+1)
    source = body.chart.evaluate(base.labels, (t, x, base.c),
                                 ((1.-np.cos(np.pi*t))/2., .1*x+.9*x*x,
                                  base.outer))['outer_source']
    profile = terminal(base.sh, n, CUT, lambda w, xx: source)
    args = base.sh, base.c, base.N, base.grid
    old = Coupled(*args, terminal_profile=profile, precondition=False)
    model = PolarCoupled(*args, terminal_profile=profile)
    model.material_load = old.material_load = get_material_load(body, d[body.mask]*body.scale, old)
    e0, g0 = old.evaluate(u)
    en, gr = model.evaluate(u)
    gradient_error = float(abs(gr-g0).max())
    emit(kind='zero_baseline', state_sha256=file_sha(path), joint_sha256=file_sha(joint),
         grid=base.grid, energy=en, difference=en-e0, gradient_error=gradient_error,
         stationarity=float(abs(model.norm(model.vars(gr))).max()))
    if abs(en-e0) > 2.e-12 or gradient_error > 2.e-12:
        raise ArithmeticError('Zero-displacement model mismatch')
    for shift in (-.01, -.0002, .0002, .01):
        model.set_shift(shift)
        old_delta = deformed_work(old, u, shift)
        en, gr = model.evaluate(u)
        delta_side = .5*((u[0]-1.)@((model.B-old.B)@(u[0]+1.)))
        error = [en-e0-old_delta.sum(), delta_side-old_delta[1]]
        emit(kind='signed_mapping', shift=shift, field_energy_change=en-e0,
             wall_energy_change=delta_side, independent=old_delta.tolist(),
             differences=error)
        if np.max(np.abs(error)) > 1.e-10:
            raise ArithmeticError('Moved energies disagree with the direct signed pull')
    # Weak force checks on all three fields include the shared material work.
    rng = np.random.default_rng(72402)
    for i, mask in enumerate(model.masks):
        direction = np.zeros_like(u)
        direction[i, mask] = rng.normal(size=mask.sum())
        direction /= np.linalg.norm(direction)
        derivative = (gr*direction).sum()
        h = 1.e-5
        signed = (model.evaluate(u+h*direction)[0]-model.evaluate(u-h*direction)[0])/(2*h)
        emit(kind='finite_field_derivative', field=i, shift=model.polar_shift,
             step=h, analytic=derivative, signed=signed)
        if abs(signed-derivative) > 3.e-7:
            raise ArithmeticError('Finite field derivative check failed')
    rematched, info = model.solve(initial=u)
    read = model.measure(rematched)
    emit(kind='rematched', shift=model.polar_shift, **info, **read)
    if read['port_error'] > 3.e-5 or abs(read['weak_disk_current']-read['feed_current']) > 3.e-5:
        raise ArithmeticError('Rematched currents disagree')
    model.set_shift(0.)
    eback = model.evaluate(u)[0]
    emit(kind='return_to_zero', energy=eback, difference=eback-e0)
    if eback != e0:
        raise ArithmeticError('Returning to zero changed the field objective')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('state', type=Path)
    parser.add_argument('joint', type=Path)
    options = parser.parse_args()
    check(options.state, options.joint)
