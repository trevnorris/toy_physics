"""Reduced shape work for the R4 melt-only determinant-barrier gauge."""
from itertools import product
from pathlib import Path
from tempfile import mkdtemp

import numpy as np

from mouth_cp10_field_shape import Bath
from mouth_cp10_label_gauge import VERTICES, depth_grid, profile, work
from mouth_free_shape_solid import Solid
from mouth_r4_eulerian_relax_labels import InfeasibleLabels, LabelRelaxation, chart_orientation
from mouth_r4_eulerian_release import block_axes, order_grid
from mouth_r4_eulerian_solid import InvalidSpatialCrystal, SpatialCrystal
from mouth_r4_fold_order import axes_model
from mouth_r4_harmonic_shape import HarmonicShapeDriver
from mouth_r4_melt_barrier import MeltBarrier


def unit_work(shape, polar, reference, axes, labels, melt, dq):
    """Analytic continued quadrature at fixed nodal offsets, including the changing corner determinants."""
    x,ts,a = axes
    length = shape.length*np.exp(dq[0])
    t = depth_grid(length,reference)
    moved = labels.astype(np.result_type(labels,dq),copy=True)
    moved[...,1] += (t-ts)[None,:,None]
    harmonic = work(shape,polar,reference,axes,moved,melt,dq)
    i,j,k = np.nonzero(melt)
    widths = np.column_stack((np.diff(x)[i],np.diff(t)[j],np.diff(a)[k]))
    weight = np.zeros(len(i),dtype=moved.dtype)
    for q in product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)),repeat=3):
        xx,tt,aa = (np.stack((x[i],t[j],a[k]),axis=1)+widths*q).T
        rm,_ = profile(shape,tt,dq)
        y = (tt-.06)/.56
        band = np.where((y.real>0)&(y.real<1),16*y**2*(1-y)**2,0.)
        hp = (polar+dq[5])*band*(3*np.cos(aa)**2-1)/2
        r = xx*rm+xx**3*hp
        dr = rm+3*xx**2*hp
        weight += 4*np.pi/64 * widths.prod(axis=1)*length*dr*r**2*np.sin(aa)
    value = harmonic
    for vertex in VERTICES:
        first = moved[i+vertex[0],j+vertex[1],k+vertex[2]]
        columns = []
        for d in range(3):
            other = vertex.copy()
            other[d] = 1-other[d]
            second = moved[i+other[0],j+other[1],k+other[2]]
            columns.append((second-first)*((1-2*vertex[d])/widths[:,d,None]))
        det = np.linalg.det(np.stack(columns,axis=-1))
        if np.min(det.real)<=0:
            raise ValueError('Analytic barrier evaluation must start with positive corner determinants')
        value += np.sum(weight*(det-1-np.log(det)))
    return value


class BarrierShapeWork:
    def __init__(self, driver, trial, solid):
        self.driver, self.trial, self.solid = driver,trial,solid
        self.initial_shape = driver.x.copy()
        self.labels = driver.labels.copy()
        self.order = order_grid(trial,driver.fields)
        self.axes = block_axes(trial.b)
        self.bath = Bath(trial,driver.base,driver.fields,driver.reference_depth)
        self.last_dx = self.last_labels = self.last_info = None

    def __call__(self, dx):
        d, trial = self.driver,self.trial
        if np.any(dx):
            shape = d.shape(self.initial_shape+dx)
            body = Solid(d.prep,shape,trial.polar_shift+dx[5],8,24,6,
                         axes_model(shape,d.reference_depth),degree=3,subdivisions=2)
            solid = SpatialCrystal(body)
        else:
            body,solid = trial.b,self.solid
        warm = self.labels.copy()
        warm[...,1] += (block_axes(body)[1]-self.axes[1])[None,:,None]
        opt = LabelRelaxation(body,warm,self.order,d.epsilon)
        folder = Path(mkdtemp(prefix='barrier_shape_',dir=d.folder))
        barrier = MeltBarrier(body,opt,folder)
        try:
            chart_orientation(opt.axes,warm,opt.active)
            labels,inner = barrier.solve(warm,'inner',folder)
        except InfeasibleLabels as exc:
            raise InvalidSpatialCrystal(str(exc),exc.witness) from exc
        if not inner['stationary']:
            raise ArithmeticError(f'Barrier shape work has a nonstationary inner solution: {inner}')
        if not np.array_equal(warm[opt.stencil],labels[opt.stencil]):
            raise ArithmeticError('A shape probe changed its active trace during inner equilibration')
        chart_orientation(opt.axes,labels,opt.active)
        E,_,_,info = solid.evaluate(labels,self.order,False)
        gauge = barrier.sample(labels)[0]
        self.last_dx,self.last_labels,self.last_info = dx.copy(),labels,dict(info,inner=inner,
                gauge_unit=float(gauge),melted_label_change=float(np.max(abs(labels-warm))))
        return np.r_[self.bath.energy(dx),E,d.epsilon*gauge].astype(np.longdouble)


class BarrierShapeDriver(HarmonicShapeDriver):
    """Use with a barrier-aware coupled field solver; step installs the evaluated endpoint labels."""
    def parts_at(self,t,spatial):
        return BarrierShapeWork(self,t,spatial)
