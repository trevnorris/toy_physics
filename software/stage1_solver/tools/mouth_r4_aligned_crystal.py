"""R4 Newton functional with a moving-point determinant barrier on all 12 edges.

The physical crystal, order interpolation, harmonic term and whole-cell
certificate are inherited from the held solver. Only the determinant mean
changes. The same local C++ energy propagates nodal, shape, tangent and finite
work increments; geometric transport and the harmonic term are applied outside
the cell kernel.
"""
import copy
import ctypes as ct
from functools import lru_cache
from hashlib import sha256
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from mouth_r4_eulerian_relax_labels import InfeasibleLabels
from mouth_r4_newton_crystal import JointCrystal, WeightedGauge
from mouth_r4_smooth_static import BITS, Gauge, LD


@lru_cache(maxsize=1)
def kernel():
    source = Path(__file__).with_name('mouth_r4_aligned_stencil.cpp')
    cache = source.parents[3]/'_scratch/s11c/r4_aligned_build'
    cache.mkdir(parents=True, exist_ok=True)
    path = cache/(sha256(source.read_bytes()).hexdigest()+'.so')
    if not path.exists():
        fd, name = tempfile.mkstemp(dir=cache,suffix='.so')
        os.close(fd)
        try:
            subprocess.run(['g++','-std=c++17','-O3','-fPIC','-shared',str(source),'-o',name],check=True)
            os.replace(name,path)
        finally:
            Path(name).unlink(missing_ok=True)
    lib = ct.CDLL(str(path))
    for name in ('real_values','complex_values','hessian_action','signed_work'):
        fn = getattr(lib,name)
        fn.argtypes = [ct.c_int]+[ct.c_void_p]*5
        fn.restype = ct.c_int
    return lib


def cell_call(mode, points, widths, change=None, gradient=False):
    """Return per-cell unit barrier means and nodal derivatives or increments."""
    dtype = np.clongdouble if np.iscomplexobj(points) or np.iscomplexobj(widths) else LD
    points = np.ascontiguousarray(points,dtype)
    widths = np.ascontiguousarray(widths,dtype)
    n = len(points)
    if points.shape != (n,8,3) or widths.shape != (n,3) or not np.all(widths.real>0):
        raise ValueError('Invalid eight-node cell or chart widths')
    grad = None
    if mode == 'values':
        fn = kernel().complex_values if dtype == np.clongdouble else kernel().real_values
        value = np.zeros(n,dtype)
        grad = np.zeros_like(points) if gradient else None
        args = [points,widths,value,grad]
    elif mode in ('hessian_action','signed_work'):
        if dtype != LD or change is None or np.shape(change) != points.shape:
            raise ValueError('A real, cellwise label displacement is required')
        fn = getattr(kernel(),mode)
        change = np.ascontiguousarray(change,LD)
        value = np.zeros_like(points) if mode=='hessian_action' else np.zeros(n,LD)
        args = [points,change,widths,value]
    else:
        raise ValueError(f'Unknown barrier operation {mode}')
    bad = ct.c_int(-1)
    rc = fn(n,*[a.ctypes.data if a is not None else None for a in args],ct.byref(bad))
    if rc:
        raise InfeasibleLabels('Nonpositive or unresolved aligned determinant',
                               selected_cell_index=bad.value,stencil_index=rc-1)
    return value, grad


class AlignedGauge(Gauge):
    def cell_barrier(self, offsets, widths, gradient):
        # Q1 translation cancels. offsets are the changes in xi-background
        # relative to vertex 000, transported with the moving geometry.
        return cell_call('values',offsets+BITS*widths[:,None,:],widths,gradient=gradient)


class AlignedWeightedGauge(WeightedGauge):
    def tangent(self, labels, step=None):
        if step is None:
            # Used only at machine epsilon in inverse_diagonal's positive floor.
            # The Jacobi diagonal itself is read from the actual joint Hessian.
            return np.broadcast_to(self.K.diagonal().reshape(labels.shape[:-1]+(1,)),labels.shape)
        moved = step.reshape(-1,3)[self.ids]
        touched = np.any(moved,axis=(1,2))
        ids = self.ids[touched]
        action,_ = cell_call('hessian_action',labels.reshape(-1,3)[ids],self.width[touched],moved[touched])
        weighted = action*(8*self.weight[touched])[:,None,None]
        out = np.zeros(labels.shape,LD)
        np.add.at(out.reshape(-1,3),ids.ravel(),weighted.reshape(-1,3))
        harmonic = (self.K@step.reshape(-1,3)).reshape(labels.shape)
        return np.asarray(out+harmonic,float)

    def work(self, old, new):
        a = old.reshape(-1,3).astype(LD)
        dv = new.reshape(-1,3).astype(LD)-a
        moved = dv[self.ids]
        touched = np.any(moved,axis=(1,2))
        signed,_ = cell_call('signed_work',a[self.ids[touched]],self.width[touched],moved[touched])
        midpoint = (a-self.base.reshape(-1,3))+dv/2
        return np.sum(signed*(8*self.weight[touched]),dtype=LD)+np.sum(dv*(self.K@midpoint),dtype=LD)


class AlignedCrystal(JointCrystal):
    gauge_class = AlignedGauge
    weighted_class = AlignedWeightedGauge


def held_forces(driver, trial, solid, labels, fields, workdir):
    """Held six physical-plus-gauge shape forces, including moving-knot work."""
    from mouth_r4_eulerian_release import order_grid
    from mouth_r4_smooth_endpoint import SHAPE_WORK_TOL, shape_work
    anchor = copy.copy(driver)
    anchor.folder,anchor.labels,anchor.fields = Path(workdir),labels,fields
    order = order_grid(trial,fields)
    energy,gradient,_,_ = solid.evaluate(labels,order)
    gauge = AlignedGauge(trial,driver.reference_depth,order)
    rows = shape_work(anchor,trial,solid,order,energy,gradient,gauge=gauge)
    return dict(six_total=[r['analytic'] for r in rows],rows=rows,
                shape_work_tolerance=SHAPE_WORK_TOL,
                shape_work_passed=(len(rows)==6 and all(abs(bank['minus_analytic']) <= SHAPE_WORK_TOL
                    for row in rows for bank in row['central'])))
