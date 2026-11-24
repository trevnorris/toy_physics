"""Resolve material-plane crossings on unresolved scalar gradients with narrowed, promoted banks."""
import argparse
import ctypes as ct
from pathlib import Path
import re
import subprocess

import numpy as np

from mouth_cp10_label_grid import BITS, get_cp10
from mouth_exact_checkpoint_control import digest, read
from mouth_free_shape_solid import Data, Grid
from mouth_r4_full_label_work import STATE_SHA
from mouth_whole_r4_relax import atomic_json

P=ct.POINTER(ct.c_longdouble)
class HPGrid(ct.Structure):
    _fields_=[('n',ct.c_int),('v',P)]
class HPData(ct.Structure):
    _fields_=[(name,HPGrid if dtype==Grid else dtype if dtype==ct.c_int else ct.c_longdouble if dtype==ct.c_double else P) for name,dtype in Data._fields_]


def promoted(area,quad=False):
    # Only scalar precision is changed. Retain the same radial subcells,
    # registered face roots, three-point quadrature and exact production labels.
    directory=Path(__file__).resolve().parent
    take=lambda name: (directory/name).read_text()
    base=take('mouth_free_shape_solid.cpp').split('struct TriSlice {',1)[0]
    chart=take('mouth_cp10_spatial_crystal.cpp').split('namespace {',1)[1].split('struct Evaluated {',1)[0]
    euler=take('mouth_r4_eulerian_solid.cpp').split('namespace euler {',1)[1].split('extern "C"',1)[0]
    cell=take('mouth_r4_full_label_work.cpp').split('#include "mouth_r4_eulerian_solid.cpp"',1)[1]
    code=base+'\n namespace {\n'+chart+'\n}\nnamespace euler {\n'+euler+'\n'+cell
    code=re.sub(r'\b(?:long\s+)?double\b','Real',code)
    generated=area/'audit.cpp'
    prefix='#include <stdexcept>\nusing Real=long double;\n'
    if quad:
        prefix='#include <quadmath.h>\n#include <stdexcept>\nusing Real=__float128;\n'
        for name in ('cos','sin','exp','log','sqrt','log1p','signbit','abs','isfinite'):
            mapped = {'abs':'fabsq','isfinite':'finiteq'}.get(name,name+'q')
            code=re.sub(r'std::'+name+r'\b',mapped,code)
        code=code.replace('std::numeric_limits<Real>::infinity()', 'HUGE_VALQ')
        code=code.replace('n<52;n++','n<120;n++')
        code+='\n'+take('mouth_r4_label_quad_bridge.hpp')
    generated.write_text(prefix+code)
    library=area/'promoted.so'
    subprocess.run(['g++','-std=gnu++17','-O3','-shared','-fPIC',str(generated),'-o',str(library),'-lquadmath'],check=True)
    lib=ct.CDLL(str(library.resolve()))
    f=lib.high_precision_banks if quad else lib.r4_cell_work
    f.argtypes=[ct.POINTER(HPData),P,P,HPGrid,P]+[ct.c_void_p]*3
    f.restype=ct.c_int
    return f

def run(root,coarse,initial,out,target,quad=False):
    out.mkdir(parents=True,exist_ok=False)
    if np.finfo(np.longdouble).eps>1.1e-19:raise ArithmeticError('Promoted precision unavailable')
    broad=read(initial)
    state=root/'r4_barrier_held01/primary/held_stop.npz'
    if digest(state)!=STATE_SHA or broad['source_sha']!=STATE_SHA:raise ValueError('Changed state')
    b,cp,axes=get_cp10(root)
    data=HPData();buffers={key:np.ascontiguousarray(buf,np.longdouble) for key,buf in b.buffers.items()}
    for key,dtype in data._fields_:
        if key in buffers:
            val=buffers[key];ptr=val.ctypes.data_as(P)
            setattr(data,key,HPGrid(len(val),ptr) if dtype==HPGrid else ptr)
        elif dtype!=P:setattr(data,key,getattr(b.data,key))
    with np.load(state) as original:
        labels=np.ascontiguousarray(original['labels'],np.longdouble)
        order=np.ascontiguousarray(original['field'][2,b.ids].transpose(0,2,1),np.longdouble)
    dims=order.shape
    cells=np.indices(tuple(v-1 for v in dims)).reshape(3,-1).T
    ids=np.stack([np.ravel_multi_index(tuple((cells+v).T),dims) for v in BITS],axis=1)
    active=np.max(order.ravel()[ids],axis=1)>0
    f=promoted(out,quad)
    gp,gw=np.polynomial.legendre.leggauss(3)
    gw=np.ascontiguousarray(gw/2,np.longdouble);gp=np.ascontiguousarray((gp+1)/2,np.longdouble)
    rule=HPGrid(3,gp.ctypes.data_as(P))
    def evaluate(ys,subset):
        enabled=np.ascontiguousarray(subset,np.int32)
        values=np.zeros(len(cells),np.longdouble)
        floor=np.zeros(1,np.longdouble)
        rc=f(ct.byref(data),ys.ctypes.data_as(P),order.ctypes.data_as(P),rule,gw.ctypes.data_as(P),
             enabled.ctypes.data,values.ctypes.data,floor.ctypes.data)
        if rc or floor[0]<=0:raise ArithmeticError(('invalid signed bank',rc,str(floor[0])))
        return values
    # Raw and Richardson values on every width are preserved. The final gate uses
    # the last raw centered width and the independently recorded coarser width.
    steps=2.**(-np.arange(25,39,2) if quad else -np.arange(24,38,2))
    rows=[];sum_covered=0
    for c in range(3):
        with np.load(coarse/f'coordinate_{c}.npz') as raw:
            physical=raw['analytic_crystal'];unit=raw['analytic_gauge'];colors=raw['colors'];r=raw['richardson'];idx=raw['free']
            mask=np.zeros(len(physical),bool)
            err=(r[0]-physical)+1e-8*(r[1]-unit)
            mask[idx]=abs(err[idx])>5e-11
            sum_covered+=len(idx)
            nodes=raw['physical_node']
        identified=np.flatnonzero(mask)
        if not len(identified):
            rows.append(dict(component=c,examined=0));continue
        # Restore exactly the stored registry from archived base labels.
        from mouth_r4_eulerian_relax_labels import LabelRelaxation
        opt=LabelRelaxation(b,np.asarray(labels,float),np.asarray(order,float),1e-8)
        cm=opt._components()[c]
        if not np.array_equal(cm.representative,nodes):raise ValueError('Changed registry')
        group=cm.groups[ids];entire=mask[group].any(axis=1)
        anchor=evaluate(labels,entire)
        deriv=np.zeros((len(steps),2,len(mask)),np.longdouble)
        scale=1e-8
        for color in np.unique(colors[identified]):
            bank=mask & (colors==color)
            touched=bank[group]
            targets=np.max(np.where(touched,group,-1),axis=1)
            subset=targets>=0
            if np.any(touched & (group!=targets[:,None])):raise ArithmeticError('Same cell different coordinate')
            node_mask=bank[cm.groups].reshape(dims)
            for i,h in enumerate(steps):
                for s,sign in enumerate((-1,1)):
                    trial=labels.copy();trial[...,c][node_mask]+=sign*h
                    value=evaluate(trial,subset)
                    local=np.where(active,1.,scale)*(value-anchor)/(sign*h)
                    np.add.at(deriv[i,s],targets[subset],local[subset])
        analytic=physical+scale*unit
        # Compare the unit gauge separately from the existing static audit. If
        # any narrow probe has nonzero gauge, recompute both strengths directly.
        if np.any(unit[identified]!=0):raise ArithmeticError('Gauge-coupled narrow-slice coordinate needs separate unit-gauge table')
        central=deriv.mean(axis=1)
        diff=central[:,identified]-analytic[identified]
        peaks=np.argmax(abs(diff),axis=1)
        worst=identified[peaks[-1]]
        disagreements=abs(central[-1,identified]-central[-2,identified])
        row=dict(component=c,examined=len(identified),widths=steps.tolist(),
                 raw_max_abs_errors=np.max(abs(diff),axis=1).astype(float).tolist(),
                 peak_nodes=[list(map(int,np.unravel_index(nodes[identified[i]],dims))) for i in peaks],
                 peak_final_node=list(map(int,np.unravel_index(nodes[worst],dims))),
                 signed_one_sided_minus_analytic_at_final_peak=(deriv[:,:,worst]-analytic[worst]).astype(float).tolist(),
                 previous_to_final_max_change=float(disagreements.max()),
                 untouched_max_error=float(np.max(abs(err[idx[~mask[idx]]]),initial=0)))
        np.savez(out/f'coordinate_{c}.npz',group=identified,analytic=analytic[identified],signed=deriv[:,:,identified],widths=steps)
        rows.append(row)
    result=dict(state_sha=STATE_SHA,raw_measurement_sha=digest(initial),guarded_input_directory=str(coarse),
                full_free_coordinate_count=sum_covered,components=rows,
                note='Only unresolved or marginal signed components were re-evaluated. Both final widths are raw central energy differences. All sites requiring an additional narrow pair had exactly zero gauge gradient.',
                pass_both_widths=all(v['examined']==0 or max(v['raw_max_abs_errors'][-2:])<=1e-10 for v in rows))
    atomic_json(target,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','coarse','initial','out','target'):p.add_argument(k,type=Path)
    p.add_argument('--quad',action='store_true',help='use binary128 and 120 face-root bisections')
    run(**vars(p.parse_args()))
