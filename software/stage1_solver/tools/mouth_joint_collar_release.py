"""Relax the two material labels and five meridian coordinates of the ordered R=4 mouth.

The steady carrier here has the original sharp transmitting end. The joined
reservoir/layer calculation provides a separate test of the end work. A
stationary point of this meridional subspace still needs angular virtual work.
"""
import hashlib
import json

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_internal_collar import CollarPair, transfer
from mouth_poly_bonded_core import core_source, core_stretches
from mouth_poly_shallow_release import CUT, Order, exterior
from mouth_polycrystal_material import stable_energy
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_shape import emit, preparation
from mouth_wound_carrier import Wound


def joined_collar_density(shape, co, ref, pair, labels):
    """Stored elastic energy per current four-volume, per unit of local order."""
    def density(r, w, theta):
        r, w, _ = np.broadcast_arrays(r, w, theta)
        out = np.zeros(r.shape)
        within = (w >= 0) & (w <= shape.length)
        if np.any(within):
            t = w[within] / shape.length
            rr = r[within]
            b = shape.profile(t)[0]
            cut = b * np.exp(-.5 * np.exp(np.polynomial.legendre.legval(2*t-1, co)))
            inner = rr < cut
            val = np.zeros_like(rr)
            if np.any(inner):
                q, shear, k, jac = core_stretches(shape, co, ref, co, t[inner], rr[inner])
                val[inner] = stable_energy(3*q*q + shear*shear + k*k, jac)/jac
            if np.any(~inner):
                val[~inner] = pair.energy_density(labels, rr[~inner], w[within][~inner])
            out[within] = val
        return out
    return density


def resample(field, t, x):
    old = tuple(np.linspace(0., 1., n) for n in field.shape)
    return RegularGridInterpolator(old, field)(np.stack(np.meshgrid(t, x, indexing='ij'), -1))


class Release:
    def __init__(self, args):
        self.a = args
        self.ref, _, self.co = preparation(args.order)
        self.t = np.linspace(0., 1., 2*args.order + 1)
        self.x = np.linspace(0., 1., args.order + 1)
        self.ot = (1-np.cos(np.pi*self.t))/2
        self.ox = .10*self.x + .90*self.x**2
        self.scale = np.array([1., .25, .25, .25])
        digest = hashlib.sha256(args.state.read_bytes()).hexdigest()
        validation = json.loads((args.state.parent/'outcome.json').read_bytes())
        limits = json.loads((args.state.parent/'limit-validation.json').read_bytes())
        if (validation['exitCode'] or validation['childOutcome']['exitCode'] or
                not limits['verified'] or int(limits['actual']['memory.swap.max']) or
                digest not in (args.state.parent/'stdout').read_text()):
            raise ValueError('Unverified input state')
        with np.load(args.state) as z:
            if int(z['radius']) != 4 or not np.array_equal(z['registration'], self.co):
                raise ValueError('Preparation mismatch')
            self.start = np.r_[(z['shape_params']-self.ref.params)/self.scale, 4*z['shape_extra']]
            self.initial_chi = resample(z['order'], self.t, self.x)
            self.initial_labels = (transfer(z['labels'], args.na, args.nt)
                                   if 'labels' in z else None)
        emit(kind='preparation', sha256=digest, order_mesh=args.order, carrier_mesh=args.flow,
             polar_mesh=args.polar, half_collar_radial_cells=args.na, collar_depth_cells=args.nt,
             registration=self.co.tolist(), entrance=4, winding=1, disk='sharp reversible')
        self.count = 0

    def model(self, v):
        return Meridian(4., self.ref.params + v[:4]*self.scale, extra=v[4]/4)

    def pair(self, shape):
        return CollarPair(shape, self.co, self.a.na, self.a.nt)

    def outside(self, shape, source, field=False):
        return exterior(shape, self.a.order, CUT, return_side_field=field,
                        side_order_work=lambda w, x: source)

    def inner(self, c):
        return self.t, self.x, c

    def outer(self, c):
        return self.ot, self.ox, c

    def wound(self, sh, c):
        return Wound(sh, c, 1, self.a.flow, self.a.polar, 28., 22.)

    def fixed_fields(self, v, c, y, amp, p):
        sh = self.model(v)
        pair = self.pair(sh)
        # dE_el/dchi is exactly independent of chi at fixed inverse labels.
        work = pair.evaluate(y, self.inner(c), self.outer(np.zeros_like(c)))
        outside = self.outside(sh, work['outer_source'])
        core = core_source(sh, self.co, self.ref, self.co, self.t, self.x)
        return (Order(sh, self.a.order).energy_gradient(c)[0] + outside['energy'] +
                outside['static_work'] + ((work['inner_source']+core)*c).sum() +
                self.wound(sh, c).amplitude(amp, p)[0])
