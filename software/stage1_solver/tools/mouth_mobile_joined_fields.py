"""Two connected reservoirs, one mobile conversion band and a live solid load at R4.

Geometry coordinates transport the physical carrier and order elements together.
The internal solid samples this order chart at its current placement and returns
a nodal load from the same interpolation.  E_field includes load @ dc.  Add the
solid energy at reference order and JointFields.static exactly once for the wound
shape objective.  Rest and zero-winding fields use the same held lattice.
"""
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_exchange_carrier import BIAS, mapped_nodes
from mouth_joined_front_stationarity import solve as solve_stationary
from mouth_last_rings_fluid import GRIDS, RimCarrier, local_axes
from mouth_polar_bonded import StaticPolar
from mouth_poly_shallow_release import CUT, Order
from mouth_receiver_order_control import side_source
from mouth_receiver_order_equilibrium import receiver
from mouth_sealed_shape_relax import Meridian


def field_axes(shape, grid, halfwidth):
    """Carry the archived local refinements with the current exit and band edges."""
    if grid != GRIDS['local']:
        return None
    nr, _, nz = GRIDS['double']
    length = shape.length
    rr, parent_z = mapped_nodes(nr, length, 28., 22.)
    low, high = length-halfwidth, length+halfwidth
    below = np.r_[np.linspace(0., low, nr+1)[:-1], np.linspace(low, length, nz+1)]
    fraction = (np.sinh(2*np.linspace(0., 1., nr+1))/np.sinh(2.))**2
    above = np.r_[np.linspace(length, high, nz+1)[:-1], high+(22.-high)*fraction]
    return local_axes(SimpleNamespace(shape=shape, h=halfwidth, rr=rr,
                                      z=(parent_z[0], below, above)))


class JointFields:
    def __init__(self, base, path, x, *, winding=1, grid=None, precondition=True):
        self.x = np.asarray(x, dtype=float).copy()
        if self.x.shape != (6,) or not np.isfinite(self.x).all():
            raise ValueError('Expected six finite radius-four shape coordinates')
        self.grid = tuple(base.grid if grid is None else GRIDS[grid] if isinstance(grid, str) else grid)
        self.order = base.c
        self.shape = Meridian(4., base.sh.params+self.x[:4]*[1., .25, .25, .25],
                              extra=base.sh.extra+.25*self.x[4])
        self.polar_shift = float(self.x[5])
        # The adopted receiver is a constrained native equilibrium at this physical
        # wall. Its local angularly dependent correction then relaxes with the
        # carrier in the two-sided terminal band. The interior base order is
        # supplied by the radius-four preparation.
        source = side_source(base, self.shape)
        self.profile, self.receiver = receiver(
            self.shape, base.order_mesh, CUT, lambda _w, _r: source,
            polar_shift=self.polar_shift, material_path=path)
        polar = StaticPolar(path, shape=self.shape, outer_order=self.profile.side)
        self.static_terms = dict(
            inner_order=float(Order(self.shape, base.order_mesh).energy_gradient(self.order)[0]),
            inner_polar=polar.inner_energy(self.polar_shift), tail=polar.tail_energy(self.polar_shift),
            receiver_order=self.receiver['energy'], receiver_crystal=self.receiver['static_work'])
        self.static = sum(self.static_terms.values())
        # A single set of basis functions defines all terms on either side of
        # the disk. There is no secondary carrier extension in this functional.
        self.axes = field_axes(self.shape, self.grid, .12)
        self.model = self._model(winding, precondition)
        for side in (0, 1):
            if not np.array_equal(self.model.ids[1,side][:,:,-1], self.model.ids[2,side][:,:,0]):
                raise ArithmeticError('Terminal traces do not share field nodes')

    def _model(self, winding, precondition):
        return RimCarrier(self.shape, self.order, winding, self.grid, axes=self.axes,
                          polar_shift=self.polar_shift, terminal_profile=self.profile,
                          precondition=precondition)

    @staticmethod
    def _install(model, load):
        load = np.asarray(load)
        if load.shape != (model.size,) or not np.isfinite(load).all():
            raise ValueError('Invalid current material load on the field chart')
        model.material_load = load

    @staticmethod
    def _check_field(model, u):
        if u.shape != model.seed.shape or not np.isfinite(u).all():
            raise ValueError('Field topology or coefficients invalid on the current chart')
        if np.any(u[2, ~model.masks[2]]):
            raise ValueError('Order correction must vanish on prescribed material traces')

    def energy(self, u, load):
        """Excess coherent/order energy including the current material coupling."""
        self._check_field(self.model, u)
        self._install(self.model, load)
        return self.model.evaluate(u)[0]

    def gradient(self, u, load):
        """Unscaled gradient on all nodes, including the prescribed port reactions."""
        self._check_field(self.model, u)
        self._install(self.model, load)
        return self.model.evaluate(u)[1]

    def residual(self, u, load):
        """Largest metric-scaled residual on the native free field coordinates."""
        if not self.model.precondition:
            raise ValueError('Checking the metric residual requires a solver metric')
        free = self.model.vars(self.gradient(u, load))
        return float(np.max(abs(self.model.norm(free))))

    @classmethod
    def _solve(cls, model, u, load, rest, tolerance):
        if not model.precondition:
            raise ValueError('Field relaxation requires a solver metric')
        cls._install(model, load)
        # A warm start supplies free coefficients only. In particular, an
        # incoming drive may have a different prescribed phase from a rest run.
        seed = model.seed.copy() if u is None else np.array(u, dtype=float, copy=True)
        cls._check_field(model, seed)
        seed[0, ~model.masks[0]] = 1.
        seed[1, model.feed] = 0.
        seed[1, model.bulk] = 0. if rest else BIAS
        model.seed = seed
        return solve_stationary(model, rest, tolerance=tolerance)

    def solve(self, u, load, rest=False, *, tolerance=2e-7):
        return self._solve(self.model, u, load, rest, tolerance)

    def _measure(self, model, u, load):
        self._check_field(model, u)
        self._install(model, load)
        info = model.measure(u)
        # Q1 order is multi-affine on each of the common physical elements.
        # Nodal extrema include the band pins and edges missed at Gauss sites.
        inner = RegularGridInterpolator((np.linspace(0., 1., len(self.order)),
                                         np.linspace(0., 1., self.order.shape[1])), self.order)
        extremes = []
        for slab in (1, 2):
            z = model.z[slab]
            band = abs(z-model.shape.length) <= model.h+1e-13
            radius, depth = np.meshgrid(model.rr[0], z[band], indexing='ij')
            c0 = (inner((depth/model.shape.length, radius)) if slab == 1
                  else self.profile(radius, depth-model.shape.length))
            chi = c0[:, None, :]+u[2, model.ids[slab,0][:,:,band]]
            extremes.append((float(chi.min()), float(chi.max())))
        info['nodal_order_range_in_band'] = [min(e[0] for e in extremes), max(e[1] for e in extremes)]
        info['disk_minus_feed_current'] = float(info['weak_disk_current']-info['feed_current'])
        return info

    def measure(self, u, load):
        return self._measure(self.model, u, load)

    def controls(self, wound, load, rest_initial=None, unwound_initial=None, *,
                 tolerance=2e-7, return_model=False):
        """Rematch both controls on this shape, at the same physical solid displacement."""
        if self.model.winding != 1:
            raise ValueError('The controls require a wound driven endpoint')
        rest, rs = self.solve(wound if rest_initial is None else rest_initial,
                              load, rest=True, tolerance=tolerance)
        zero = self._model(0, True)
        if any(not np.array_equal(ids, zero.ids[key]) for key, ids in self.model.ids.items()):
            raise ArithmeticError('Different material node numbering on unwound chart')
        plain, zs = self._solve(zero, wound if unwound_initial is None else unwound_initial,
                                load, False, tolerance)
        driven = self.measure(wound, load)
        resting = self.measure(rest, load)
        unwound = self._measure(zero, plain, load)
        denominator = unwound['feed_current']
        info = dict(wound=driven, rest=resting, unwound=unwound,
                    retained_fraction=driven['feed_current']/denominator if denominator > 0 else None,
                    rest_solver=rs, unwound_solver=zs)
        return (rest, plain, info, zero) if return_model else (rest, plain, info)
