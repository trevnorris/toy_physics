"""Six frozen total shape forces on nested physical grids, at CP10's material motion.

The true depth motion moves the original order knots and linearly subdivides those
same intervals. The Q1 order values of the refined grid are therefore fixed.
Transport of the exact inverse labels has a one-sided tangent at element faces.
Central energy differences transport those labels before evaluating the original
two-metric crystal energy with registered triangle cuts in the physical volume.
"""
import argparse
import ctypes as ct
import json
from pathlib import Path
import subprocess

import numpy as np

from mouth_cp10_spatial_crystal import MODES, load_context, signed_geometry
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest
from mouth_free_shape_solid import Data, Grid, Pointer, pointer


def subdivide(axis, level):
    """Subdivide without moving the parent knots, including roundoff."""
    frac = np.arange(level)/level
    return np.r_[(axis[:-1,None] + np.diff(axis)[:,None]*frac).ravel(), axis[-1]]


def kernel(area):
    source = Path(__file__).with_name('mouth_cp10_spatial_crystal.cpp')
    lib = area/'kernel.so'
    subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC','-DCP10_SPLIT_REGISTRY_DIAGONALS',str(source),'-o',str(lib)], check=True)
    f = ct.CDLL(str(lib)).cp10_spatial_integral
    f.argtypes = [ct.POINTER(Data),ct.c_int]+[Pointer]*5+[ct.c_int]+[Pointer]*4
    f.restype = ct.c_int
    return f


def integral(f, body, axes, tdot, xi, chi, rule, derivatives):
    d = Data.from_buffer_copy(body.data)
    for c,a in zip('xta',axes):
        setattr(d,'field_'+c,Grid(len(a),pointer(a)))
    out,bad,test = np.zeros(31),np.zeros(17),np.zeros(8)
    cells = np.zeros((len(body.t)-1,len(body.x)-1,len(body.a)-1,2))
    code = f(ct.byref(d),-1 if derivatives else 0,pointer(tdot),pointer(xi),pointer(chi),
             pointer(rule[0]),pointer(rule[1]),len(rule[0]),pointer(out),pointer(bad),pointer(test),pointer(cells))
    if code or bad[13] or not np.isfinite(out).all():
        raise ValueError(f'Ordered-cell evaluation failed: status={code}, witness={bad.tolist()}')
    return dict(energy=out[0],analytic=out[1:7].tolist() if derivatives else None,
                min_Jinv_sampled=bad[15],min_chart_det_sampled=bad[16])


def run(area,lab,output,level,degree):
    output.mkdir(parents=True,exist_ok=False)
    def save():
        tmp=output/'partial.json'
        tmp.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
        tmp.replace(output/'readout.json')
    f=kernel(output)
    old=area/'r4_cp10_labels01/labels.npz'
    ctx=load_context(area,old)
    base=ctx[2]
    with np.load(lab,allow_pickle=False) as data:
        if str(data['source_sha256'])!=CHECKPOINT_SHA:raise ValueError('Wrong checkpoint label field')
        axes=tuple(np.ascontiguousarray(data[c]) for c in 'xta')
        xi=np.ascontiguousarray(data['labels'])
        chi=np.ascontiguousarray(data['order'])
        active=data['active_stencil']
        if level>1:
            minus,plus=(np.ascontiguousarray(data['dxi_dt_'+k]) for k in ('minus','plus'))
        else:
            with np.load(area/'r4_cp10_labels01/label_transport.npz',allow_pickle=False) as rate:
                if str(rate['label_sha256'])!=digest(lab):raise ValueError('Transport and labels disagree')
                minus,plus=(np.ascontiguousarray(rate['dxi_dt'][...,i,:]) for i in (0,1))
    for c,a in zip('xta',axes):
        expected=subdivide(base.buffers['field_'+c],level)
        if not np.allclose(a,expected,rtol=0.,atol=4e-16):raise ValueError('Expected linearly subdivided knots')
    if xi.shape!=chi.shape+(3,) or not np.all(np.isfinite(xi[active])) or not np.all(np.isfinite(chi)):
        raise ValueError('Invalid label transfer')
    loc,weight=np.polynomial.legendre.leggauss(degree)
    rule=(np.ascontiguousarray((1.+loc)/2),np.ascontiguousarray(weight/2))
    t=base.buffers['field_t']
    L=base.data.length
    parent_rate=np.where(t*L<L-.12,(t*L/(L-.12))**25,1.)-t
    tdot=subdivide(parent_rate,level)
    row=dict(checkpoint=CHECKPOINT_SHA,label_sha256=digest(lab),level=level,quadrature_degree=degree,
             native_signed_steps=[4.e-7,8.e-8],modes=MODES,fixed=integral(f,base,axes,tdot,xi,chi,rule,True),probes=[])
    save()
    for h in row['native_signed_steps']:
        pair=[]
        for sign in (-1,1):
            signed=signed_geometry(ctx,0,h*sign)
            coarset=signed.buffers['field_t']
            moved_t=subdivide(coarset,level)
            # Compute the difference before subdivision, to preserve small
            # displacements accurately when t is close to its upper endpoint.
            dt=subdivide(coarset-t,level)
            rate=np.where(dt[None,:,None,None]>=0.,plus,minus)
            product=np.multiply(dt[None,:,None,None],rate,out=np.zeros_like(xi),where=(dt!=0.)[None,:,None,None])
            shifted=np.ascontiguousarray(xi+product)
            if not np.all(np.isfinite(shifted[active])):raise ArithmeticError('Missing signed ordered label tangent')
            value=integral(f,signed,(axes[0],moved_t,axes[2]),tdot,shifted,chi,rule,False)
            pair.append(value)
        row['probes'].append(dict(h=h,signed=pair,force=(pair[0]['energy']-pair[1]['energy'])/(2*h)))
        save()
    # Fluid/order and the exterior crystal are unchanged by crystal label interpolation.
    from mouth_exact_checkpoint_control import read
    reference=read(Path(__file__).resolve().parents[3]/'research/pde_ledger_v3/_measurements/mouth_r4_two_metric_labels.json')
    original=reference['data']['total_shape_forces']
    field=np.array(original['field_analytic'])
    material=np.array(row['fixed']['analytic']); material[0]=row['probes'][-1]['force']
    totals=material+field
    row.update(total=totals.tolist(),matched_crystal=material.tolist(),delta_against_saved=(totals-original['native_saved']).tolist(),
               transport_delta=float(material[0]-row['fixed']['analytic'][0]),
               depth_step_spread=abs(row['probes'][0]['force']-row['probes'][1]['force']),
               derivative='frozen label analytic, except log depth from signed energies with one-sided inverse-label nodal velocities and exact parent grid motion',
               completed=True)
    save()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('area','labels','output'):p.add_argument(arg,type=Path)
    p.add_argument('level',type=int,choices=[1,2,4])
    p.add_argument('degree',type=int)
    args=p.parse_args()
    run(args.area.resolve(),args.labels.resolve(),args.output.resolve(),args.level,args.degree)
