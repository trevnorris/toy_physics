"""Frozen, local-stencil tests of the R4 joint active/melt Newton product.

Promote the four-dimensional physical gradient, including the exact moving
material integration cuts, before taking signed differences. Unaffected cell
contributions cancel identically. The smooth determinant/gradient gauge is
differenced separately so that the melt-only block is tested at its own scale.
"""
import argparse
import ctypes as ct
from pathlib import Path
import json
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

import numpy as np
from scipy.ndimage import minimum_filter, maximum_filter

from mouth_exact_checkpoint_control import digest, read
from mouth_r4_eulerian_relax_labels import LabelRelaxation, InfeasibleLabels
from mouth_r4_eulerian_release import ShapeTrial, axes_for_depth, order_grid
from mouth_r4_eulerian_solid import SpatialCrystal, InvalidSpatialCrystal
from mouth_r4_full_label_precision import HPData, HPGrid, P
from mouth_r4_newton_anchor import load_state
from mouth_r4_projected_fields import location
from mouth_r4_slaved_melt import FLOW
from mouth_r4_smooth_fields import read_fields
from mouth_r4_smooth_held import SmoothHeld
from mouth_r4_smooth_static import Gauge, BITS, LD
from mouth_whole_r4_relax import atomic_json

WIDTHS = (2.**-16,2.**-18,2.**-20)
SECTORS = ('ordered','front','melted')
FULL_ERROR_TOL = 1.e-4
HALVING_TOL = 1.e-5
SYMMETRY_TOL = 1.e-5


class PromotedGradient:
    """Use the analytic production formula on the cells incident to a difference."""
    def __init__(self, body, directory):
        sources=Path(__file__).resolve().parent
        def code(name):
            return (sources/name).read_text()
        base=code('mouth_free_shape_solid.cpp').split('struct TriSlice {',1)[0]
        chart=code('mouth_cp10_spatial_crystal.cpp').split('namespace {',1)[1].split('struct Evaluated {',1)[0]
        calc=code('mouth_r4_eulerian_solid.cpp').split('namespace euler {',1)[1]
        signature='double *diagnostics,double *intact){'
        loop='Q1 q;bool has_order=false;int ids[2][2][2];'
        if calc.count(signature)!=1 or calc.count(loop)!=1:
            raise ValueError('Cannot identify the current archived gradient kernel structure')
        calc=calc.replace(signature, 'double *diagnostics,double *intact,const int *selected){')
        calc=calc.replace(loop,'if (!selected[(i*(nt-1)+j)*(na-1)+k]) continue;\n'+loop)
        src=re.sub(r'\b(?:long\s+)?double\b','Real',base+'\nnamespace {\n'+chart+'\n}\nnamespace euler {\n'+calc)
        # Resolve diagonal registry cuts to the promoted mantissa.
        src=src.replace('n<52;n++','n<72;n++')
        path=directory/'stencil.cpp'
        path.write_text('#include <stdexcept>\nusing Real=long double;\n'+src)
        so=directory/'stencil.so'
        subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC',str(path),'-o',str(so)],check=True)
        self.function=ct.CDLL(str(so)).r4_eulerian_material
        self.function.argtypes=[ct.POINTER(HPData),P,P,HPGrid]+[P]*6+[ct.c_void_p]
        self.function.restype=ct.c_int
        self.data=HPData()
        self.buffers={name:np.ascontiguousarray(v,LD) for name,v in body.buffers.items()}
        for name,typ in self.data._fields_:
            if name in self.buffers:
                v=self.buffers[name]
                ptr=v.ctypes.data_as(P)
                setattr(self.data,name,HPGrid(len(v),ptr) if typ==HPGrid else ptr)
            elif typ!=P:
                setattr(self.data,name,getattr(body.data,name))
        nodes,weights=np.polynomial.legendre.leggauss(3)
        self.nodes=np.ascontiguousarray((nodes+1)/2,LD)
        self.weights=np.ascontiguousarray(weights/2,LD)
        self.rule=HPGrid(3,self.nodes.ctypes.data_as(P))

    def gradient(self,labels,order,selected):
        def ptr(a):
            return a.ctypes.data_as(P)
        sel=np.ascontiguousarray(selected,np.int32)
        grad=np.zeros_like(labels)
        if not sel.any():return grad
        loads=np.zeros_like(order)
        diag=np.zeros(23,LD)
        intact=np.full((order.shape[1]-1,4),-1.,LD)
        energy=np.zeros(1,LD)
        rc=self.function(ct.byref(self.data),ptr(labels),ptr(order),self.rule,ptr(self.weights),ptr(energy),
                         ptr(grad),ptr(loads),ptr(diag),ptr(intact),sel.ctypes.data)
        if rc or diag[0]<=0 or not np.all(np.isfinite(grad)):
            raise ArithmeticError(f'Infeasible promoted stencil, return code {rc}, Jinv {diag[0]}')
        return grad


def reduce_gradient(components,gradient):
    result=[]
    for c,cm in enumerate(components):
        sums=np.zeros(len(cm.held),gradient.dtype)
        np.add.at(sums,cm.groups,gradient[...,c].ravel())
        result.append(sums[~cm.held])
    return np.concatenate(result)


def directions(trial,opt,order):
    # Use the same spatially interior three-component supports for both test
    # vectors at each site, so transpose work measures actual cross-coupling.
    components=opt._components()
    mins=minimum_filter(order,size=3,mode='nearest')
    maxs=maximum_filter(order,size=3,mode='nearest')
    idx=np.indices(order.shape)
    inside=np.all((idx>1)&(idx<np.array(order.shape)[:,None,None,None]-2),axis=0)
    offsets=np.cumsum([0]+[int(np.sum(~cm.held)) for cm in components])
    index=[]
    for c,cm in enumerate(components):
        rank=np.cumsum(~cm.held)-1+offsets[c]
        index.append(rank[cm.groups].reshape(order.shape))
        inside&=(~cm.held[cm.groups]&(np.bincount(cm.groups)[cm.groups]==1)).reshape(order.shape)
    # Nearest eligible point to each prescribed chart coordinate, independent of
    # derivative agreement, conditioning, or whether the benchmark is faster.
    coord=np.stack(np.meshgrid(*opt.axes,indexing='ij'),axis=-1)
    masks=(mins>0.05,(mins==0)&(maxs>0)&(order>0),maxs==0)
    targets=((.55,.08,.75),(.70,.27,.75),(.6,.57,.75))
    for sector, mask, target in zip(SECTORS,masks,targets):
        ids=np.flatnonzero(mask&inside)
        if not len(ids):raise ValueError(f'No freely movable {sector} stencil')
        dist=np.sum((coord.reshape(-1,3)[ids]-target)**2,axis=1)
        node=np.unravel_index(ids[np.argmin(dist)],order.shape)
        widths=np.array([min(ax[n]-ax[n-1],ax[n+1]-ax[n]) for ax,n in zip(opt.axes,node)])
        spec=dict(**location(trial,opt.axes,node),chi=float(order[node]),patch_order_min=float(mins[node]),
                  patch_order_max=float(maxs[node]),registered_widths=widths.tolist())
        dirs=[]
        for mix in ((1.,-.4,.75),(.65,.85,-.5)):
            lifted=np.zeros(order.shape+(3,))
            lifted[node]=widths*mix
            packed=np.zeros(offsets[-1])
            for c in range(3):
                packed[index[c][node]]=lifted[node+(c,)]
            dirs.append((packed,lifted))
        yield sector,spec,dirs


def make_oracle(study):
    folder=study/'oracle'
    folder.mkdir()
    with TemporaryDirectory(dir=folder) as work:
        work=Path(work)
        driver=SmoothHeld(study.parent,work,1e-8)
        cp=load_state(driver,study/'paused')
        t=ShapeTrial(driver.base,driver.prep,driver.x,8,24,6,3,2,precondition=False,
                     axes=axes_for_depth(driver.shape(driver.x).length,driver.reference_depth))
        order=order_grid(t,driver.fields)
        gauge=Gauge(t,driver.reference_depth,order)
        opt=LabelRelaxation(t.b,driver.labels,order,driver.epsilon)
        solid=PromotedGradient(t.b,work)
        labels=driver.labels.astype(LD)
        order_ld=np.ascontiguousarray(order,LD)
        active=gauge.order.max(axis=1)>0.
        arrays={}
        cases={}
        for sector,spec,dirs in directions(t,opt,order):
            cases[sector]=spec
            for key,(packed,lifted) in zip(('u','v'),dirs):
                stem=sector+'_'+key
                arrays[stem]=packed
                selected=np.any(np.any(lifted!=0.,axis=-1).ravel()[gauge.ids],axis=1)
                arms=[]
                for width in WIDTHS:
                    g=[]
                    for sign in (-1,1):
                        trial_labels=labels+LD(sign*width)*lifted.astype(LD)
                        material=solid.gradient(trial_labels,order_ld,selected&active)
                        _,numerical=gauge.evaluate(trial_labels,np.zeros(6,LD),enabled=selected,gradient=True)
                        g.append(reduce_gradient(opt.components,material+LD(driver.epsilon)*numerical))
                    arms.append(np.asarray((g[1]-g[0])/(2*LD(width)),float))
                arrays[stem+'_Hv']=np.stack(arms)
        path=folder/'differences.npz'
        np.savez_compressed(path,**arrays)
        atomic_json(folder/'record.json',dict(state_sha256=cp['sha256'],local_gradient_bank=digest(path),
            mantissa_bits=int(np.finfo(LD).nmant+1),relative_steps=WIDTHS,sites=cases,
            note='Centered difference of analytic quadrature gradients; both physical and gauge partials are promoted before differencing. Independent of the Newton product.'))


def norm(v):
    return float(np.linalg.norm(v))


def ratio(v,w):
    n=norm(w)
    if n==0 or not np.isfinite(n):raise ArithmeticError('Zero/invalid curvature in test reference')
    return norm(v)/n


def orientation_controls(joint,labels,opt,probes):
    from mouth_cp10_label_gauge import gradient_basis
    from mouth_r4_barrier_held import HeldStop
    records={}
    for name in ('ordered','melted'):
        cell=np.array(probes[name]['node'])
        coords=tuple((cell+BITS).T)
        trial=labels.copy()
        trial[coords]=labels[coords][np.arange(8)^4]
        widths=np.array([np.diff(ax)[ix] for ax,ix in zip(opt.axes,cell)])
        D=trial[coords].T@gradient_basis((.5,.5,.5))/widths[None,:]
        determinant=float(np.linalg.det(D))
        if determinant>=0:raise ValueError('Failed to construct an inverted case')
        try:
            joint.evaluate(trial)
        except (ValueError,HeldStop) as exc:
            records[name]=dict(cell=cell.tolist(),chart_det_center=determinant,reason=str(exc))
        else:raise ArithmeticError(f'Invalid {name} map admitted by candidate evaluator')
    return records


def audit(study,output):
    from mouth_r4_newton_crystal import JointCrystal
    from mouth_r4_stage2_outcome import verify_guard
    for guard in ('baseline_guard','oracle_guard_serial'):
        verify_guard(study/guard,0)
    baseline=read(study/'baseline/record.json')
    ref=read(study/'oracle/record.json')
    bank=study/'oracle/differences.npz'
    if digest(bank)!=ref['local_gradient_bank']:raise ValueError('Changed promoted gradient bank')
    result=dict(command=sys.orig_argv,source_sha256=digest(Path(__file__)),
          newton_source_sha256=digest(Path(__file__).with_name('mouth_r4_newton_crystal.py')),
          base_solver='393aec91',limits=dict(relative_product=FULL_ERROR_TOL,
          relative_reference_convergence=HALVING_TOL,relative_transpose=SYMMETRY_TOL),sectors={})
    with TemporaryDirectory(dir=study) as work,np.load(bank,allow_pickle=False) as data:
        driver=SmoothHeld(study.parent,Path(work),1e-8)
        cp=load_state(driver,study/'paused')
        result['state_sha256']=cp['sha256']
        if cp['sha256']!=baseline['state_sha256'] or cp['sha256']!=ref['state_sha256']:
            raise ValueError('Readouts used different held maps')
        t=driver.trial()
        solid=SpatialCrystal(t.b)
        order=order_grid(t,driver.fields)
        joint=JointCrystal(driver,t,solid,driver.labels,order)
        energy,gradient,material=joint.evaluate(driver.labels)
        packed_active=np.concatenate([cm.free[~cm.held] for cm in joint.components])
        baseline_g=np.load(study/'baseline/active_gradient.npy',allow_pickle=False)
        physical,load,field=read_fields(driver,t,solid,driver.labels,driver.fields)
        current=t.measure(driver.fields,load)
        material_delta={k:material[k]-v for k,v in baseline['material'].items()}
        field_delta={k:field[k]-baseline['field'][k] for k in
                    ('order_raw','order_kkt','order_projected','amplitude_raw','phase_raw','carrier_metric')}
        current_delta={k:current[k]-baseline['current'][k] for k in FLOW}
        gdelta=float(np.max(np.abs(gradient[packed_active]-baseline_g)))
        energy_delta=physical+material['gauge_energy']-baseline['total_energy']
        sum_delta=energy-(baseline['material']['elastic_energy']+baseline['material']['gauge_energy'])
        all_delta=[v for group in (material_delta,field_delta,current_delta) for v in group.values()]
        all_delta.extend([gdelta,energy_delta,sum_delta])
        result['invariance']=dict(material=material_delta,field=field_delta,current=current_delta,
            active_gradient_max_abs=gdelta,total_energy=energy_delta,crystal_plus_gauge=sum_delta,
            passed=bool(max(map(abs,all_delta))<=8*np.finfo(float).eps))
        result['positivity_controls']=orientation_controls(joint,driver.labels,
                                    joint.opt,ref['sites'])
        for sector,spec,dirs in directions(t,joint.opt,order):
            row=dict(location=spec,directions={})
            native=[]
            wanted=[]
            ok=True
            for k,(packed,lifted) in zip(('u','v'),dirs):
                stem=sector+'_'+k
                if not np.array_equal(packed,data[stem]) or not np.array_equal(joint.lift(packed),lifted):
                    raise ValueError(f'Degree-of-freedom or direction mismatch in {stem}')
                reference=data[stem+'_Hv']
                exact=reference[-1]
                got=joint.hvp(driver.labels,packed)
                native.append(got)
                wanted.append(exact)
                sequence=[ratio(a-b,b) for a,b in zip(reference,reference[1:])]
                discrepancy=ratio(got-exact,exact)
                row['directions'][k]=dict(reference_norm_l2=norm(exact),native_norm_l2=norm(got),
                    reference_widths=WIDTHS,relative_successive_differences=sequence,
                    native_minus_reference_linf=float(np.max(np.abs(got-exact))),
                    native_minus_reference_relative_l2=discrepancy)
                ok &= discrepancy<=FULL_ERROR_TOL and sequence[-1]<=HALVING_TOL
            u,v=(x[0] for x in dirs)
            for engine,products in [('production',native),('promoted',wanted)]:
                hu,hv=products
                uv=float(u@hv)
                vu=float(v@hu)
                scale=max(abs(uv),abs(vu))
                if scale==0:raise ArithmeticError('Bilinear test has zero coupled work')
                rel=abs(uv-vu)/scale
                row[engine+'_symmetry']=dict(uHv=uv,vHu=vu,difference=uv-vu,relative=rel)
                ok &= rel<=SYMMETRY_TOL
            row['passed']=bool(ok)
            result['sectors'][sector]=row
            atomic_json(output,result)
        result['baseline']=baseline
        result['passed']=result['invariance']['passed'] and all(v['passed'] for v in result['sectors'].values())
        atomic_json(output,result)
    if not result['passed']:raise ArithmeticError('Read-only Newton audit failed; see signed results')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('operation',choices=('oracle','audit'))
    p.add_argument('study',type=Path)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.operation=='oracle':make_oracle(a.study.resolve())
    else:
        if a.output is None:p.error('--output is required for the gate record')
        audit(a.study.resolve(),a.output.resolve())
