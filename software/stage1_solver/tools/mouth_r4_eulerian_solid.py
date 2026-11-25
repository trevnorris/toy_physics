"""Production crystal variation in a continuous spatial chart with two natural metrics."""
import ctypes as ct
from functools import lru_cache
from hashlib import sha256
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from mouth_free_shape_solid import Data, Grid, Pointer, pointer


class InvalidSpatialCrystal(ValueError):
    def __init__(self, message, diagnostics):
        self.diagnostics=diagnostics
        super().__init__(message)


@lru_cache(maxsize=1)
def _library():
    source=Path(__file__).with_suffix('.cpp')
    files=[source,source.with_name('mouth_cp10_spatial_crystal.cpp'),source.with_name('mouth_free_shape_solid.cpp')]
    key=sha256(b''.join(p.read_bytes() for p in files)).hexdigest()
    folder=source.parents[3]/'_scratch'/'s11c'/'r4_eulerian_crystal_build'
    folder.mkdir(parents=True,exist_ok=True)
    so=folder/(key+'.so')
    if not so.is_file():
        fd,name=tempfile.mkstemp(dir=folder,suffix='.so')
        os.close(fd)
        try:
            subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC',str(source),'-o',name],check=True)
            os.replace(name,so)
        finally:
            Path(name).unlink(missing_ok=True)
    lib=ct.CDLL(str(so))
    chart=lib.r4_chart_orientation
    chart.argtypes=[Grid,Grid,Grid,Pointer,ct.c_void_p,Pointer]
    chart.restype=ct.c_int
    f=lib.r4_eulerian_material
    f.argtypes=[ct.POINTER(Data),Pointer,Pointer,Grid]+[Pointer]*6
    f.restype=ct.c_int
    return lib


def _kernel():
    return _library().r4_eulerian_material


def chart_bounds(axes, labels, active):
    buffers=tuple(np.ascontiguousarray(a,dtype=float) for a in axes)
    labels=np.ascontiguousarray(labels,dtype=float)
    cells=np.ascontiguousarray(active,dtype=np.uint8)
    if (len(buffers)!=3 or any(a.ndim!=1 or a.size<2 or not np.isfinite(a).all()
                               or not np.all(np.diff(a)>0) for a in buffers)):
        raise ValueError('The label axes must be finite and strictly increasing')
    dims=tuple(len(a) for a in buffers)
    if labels.shape!=dims+(3,) or cells.shape!=tuple(n-1 for n in dims):
        raise ValueError('Label or active cell shape differs from chart axes')
    grids=[Grid(len(a),pointer(a)) for a in buffers]
    out=np.zeros(11)
    rc=_library().r4_chart_orientation(*grids,pointer(labels),cells.ctypes.data,pointer(out))
    cell=out[2:5].astype(int).tolist()
    if rc:
        raise InvalidSpatialCrystal('Nonpositive or roundoff-unresolved chart determinant',
            dict(code=rc,cell=cell,chart_fractions=out[5:8].tolist(),chart_det=float(out[8]),
                 chart_lower=float(out[0]),roundoff_margin=float(out[9]),subdivisions=int(out[10])))
    return dict(chart_det_global_lower=float(out[0]),chart_lowest_cell=cell,
                chart_det_active_lower=float(out[1]) if cells.any() else None)


class SpatialCrystal:
    """The body supplies the geometry and both registered reference embeddings."""
    def __init__(self,body,degree=3):
        self.body=body
        self.axes=tuple(body.buffers['field_'+k] for k in 'xta')
        self.shape=tuple(len(a) for a in self.axes)
        if not isinstance(degree,(int,np.integer)) or degree<2:
            raise ValueError('Crystal quadrature must have degree at least two')
        gp,gw=np.polynomial.legendre.leggauss(degree)
        self.points=np.ascontiguousarray((gp+1)/2)
        self.weights=np.ascontiguousarray(gw/2)
        self.rule=Grid(len(gp),pointer(self.points))
        self.function=_kernel()

    def evaluate(self,labels,order,gradient=True):
        labels=np.ascontiguousarray(labels,dtype=np.float64)
        order=np.ascontiguousarray(order,dtype=np.float64)
        if labels.shape!=self.shape+(3,) or order.shape!=self.shape or not np.isfinite(labels).all():
            raise ValueError('Expected registered chart labels and full spatial order')
        if not np.isfinite(order).all() or order.min()<0 or order.max()>1:
            raise ValueError('Order outside the unit interval')
        if any(labels[...,k].min()<ax[0]-2.e-13 or labels[...,k].max()>ax[-1]+2.e-13 for k,ax in enumerate(self.axes)):
            raise InvalidSpatialCrystal('Label leaves the registered chart',{'code':'label_domain'})
        grad=np.zeros_like(labels) if gradient else None
        load=np.zeros_like(order);energy=np.zeros(1);diagnostics=np.zeros(23)
        intact=np.full((self.shape[1]-1,4),-1.)
        rc=self.function(ct.byref(self.body.data),pointer(labels),pointer(order),self.rule,pointer(self.weights),
                         pointer(energy),pointer(grad) if gradient else Pointer(),pointer(load),pointer(diagnostics),pointer(intact))
        if rc:
            reason={2:'nonpositive sampled chart determinant',3:'unresolved chart orientation bound',
                    4:'invalid physical reference or signed volume',5:'nonfinite crystal energy',6:'indeterminate material face motion'}[rc]
            item=dict(code=rc,chart_bound=float(diagnostics[2]),cell=diagnostics[5:8].astype(int).tolist(),
                      local=diagnostics[8:11].tolist(),chi=float(diagnostics[11]))
            raise InvalidSpatialCrystal(reason,item)
        def minimum(index,start):
            if not np.isfinite(diagnostics[index]):return None
            cell=diagnostics[start:start+3].astype(int)
            local=diagnostics[start+3:start+6]
            chart=[float(ax[c]+s*(ax[c+1]-ax[c])) for ax,c,s in zip(self.axes,cell,local)]
            return dict(cell=cell.tolist(),local=local.tolist(),chart=chart,chi=float(diagnostics[start+6]),
                        Jinv=float(diagnostics[index]),grad_xi_norm=float(diagnostics[start+7]),W=float(diagnostics[start+8]))
        info=dict(min_Jinv_global=float(diagnostics[0]),min_Jinv_active=float(diagnostics[1]),
                  min_chart_Jinv_bound=float(diagnostics[2]),min_chart_Jinv_active_bound=float(diagnostics[3]),
                  max_Jinv_sampled=float(1./diagnostics[4]),worst_global=minimum(0,5),worst_active=minimum(1,14),
                  intact_by_depth=[None if v[0]<0 else dict(peak_W=v[0],physical_depth=v[1],chi=v[2],Jinv=v[3]) for v in intact])
        if not (np.isfinite(energy).all() and np.isfinite(load).all() and (grad is None or np.isfinite(grad).all())):
            raise InvalidSpatialCrystal('Nonfinite crystal variation',info)
        return float(energy[0]),grad,load,info
