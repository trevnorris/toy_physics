"""Central energy differences of every independent free scalar label on the held stop.

Q1 elements are local. Color the node registry so a cell touches at most one
varied degree of freedom in each signed evaluation. Its signed energy belongs
unambiguously to that scalar (including axis and melted angular ties).
"""
import argparse
import ctypes as ct
import json
import subprocess
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import BITS, get_cp10
from mouth_exact_checkpoint_control import digest
from mouth_free_shape_solid import Data, Grid, Pointer, pointer
from mouth_r4_eulerian_relax_labels import LabelRelaxation
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_melt_barrier import MeltBarrier
from mouth_whole_r4_relax import atomic_json

STATE_SHA='418376c3eb1922f6335543f43dfee11e531cebbb6958d37afd73e844571b653a'


def coloring(cm,ids):
    group=cm.groups[ids]
    neighbors=[set() for _ in range(len(cm.held))]
    for vertices in group:
        live=set(vertices[~cm.held[vertices]])
        for v in live:
            neighbors[v].update(live)
    color=np.full(len(cm.held),-1,int)
    for v in sorted(np.flatnonzero(~cm.held),key=lambda i:-len(neighbors[i])):
        used={color[k] for k in neighbors[v]}
        c=0
        while c in used:c+=1
        color[v]=c
    return color


def run(root,out,record):
    out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).with_suffix('.cpp')
    library=out/'cell.so'
    subprocess.run(['g++','-std=c++17','-O3','-fPIC','-shared',str(source),'-o',str(library)],check=True)
    func=ct.CDLL(str(library.resolve())).r4_cell_work
    func.argtypes=[ct.POINTER(Data),Pointer,Pointer,Grid,Pointer]+[ct.c_void_p]*3
    func.restype=ct.c_int
    state=root/'r4_barrier_held01/primary/held_stop.npz'
    if digest(state)!=STATE_SHA:raise ValueError('Wrong held stop')
    b,cp,axes=get_cp10(root)
    with np.load(state) as f:
        if not np.array_equal(f['x'],cp['x']) or int(f['radius'])!=4 or int(f['winding'])!=1:raise ValueError('Wrong geometry')
        labels=f['labels'];order=np.ascontiguousarray(f['field'][2,b.ids].transpose(0,2,1))
    opt=LabelRelaxation(b,labels,order,1e-8)
    barrier=MeltBarrier(b,opt,out)
    solid=SpatialCrystal(b)
    energy,physical,_,baseinfo=solid.evaluate(labels,order)
    unit,gauge,_,_=barrier.sample(labels)
    dims=order.shape
    cells=np.indices(opt.active.shape).reshape(3,-1).T
    ids=np.stack([np.ravel_multi_index(tuple((cells+v).T),dims) for v in BITS],axis=-1)
    active=opt.active.ravel()

    def evaluate(label,enabled):
        enabled=np.ascontiguousarray(enabled,dtype='i4')
        value=np.zeros(len(cells),dtype=np.longdouble);minimum=np.zeros(1)
        if any(label[...,c].min()<ax[0]-2e-13 or label[...,c].max()>ax[-1]+2e-13 for c,ax in enumerate(axes)):
            raise ValueError('Infeasible reference domain')
        rc=func(ct.byref(b.data),pointer(label),pointer(order),solid.rule,pointer(solid.weights),
                enabled.ctypes.data,value.ctypes.data,minimum.ctypes.data)
        if rc:raise ValueError(f'Infeasible paired state {rc}')
        return value,minimum[0]

    baseline,min_orient=evaluate(labels,np.ones(len(cells),bool))
    cryst_delta=float(baseline[active].sum()-energy);unit_delta=float(baseline[~active].sum()-unit)
    if abs(cryst_delta)>1e-13 or abs(unit_delta)>1e-11:raise ArithmeticError((cryst_delta,unit_delta))
    result=dict(source_sha=STATE_SHA,baseline_crystal_energy=energy,raw_unit_gauge_energy=float(unit),
                cell_sum_minus_native=dict(physical=cryst_delta,unit_gauge=unit_delta),
                max_abs_required=1e-10,baseline_orientation=baseinfo['min_chart_Jinv_bound'],
                convention='Chart-component unit steps at fixed chi, n, phase, q on the same stopped state. Eight-cell axis ties vary together. Every separately free scalar is tested. Crystal and gauge work are summed separately over their frozen, disjoint cell supports. Refinement widths are dyadic.',components=[])
    for c,cm in enumerate(barrier.components):
        group=cm.groups[ids]
        color=coloring(cm,ids)
        independent=np.flatnonzero(~cm.held)
        pg=cm.sum(physical[...,c]);gg=cm.sum(gauge[...,c])
        # Local centered banks avoid cancellation against label-independent energies.
        h0=2.**-19
        steps=[h0,h0/2,h0/4]
        fd=np.zeros((len(steps),2,len(cm.held)),dtype=np.longdouble)
        for code in range(color.max()+1):
            matching=color[group]==code
            targets=np.max(np.where(matching,group,-1),axis=1)
            enabled=targets>=0
            if np.any(matching & (group!=targets[:,None])):raise AssertionError('Overlapping signed stencils')
            nodes=(color[cm.groups]==code).reshape(dims)
            for idx,h in enumerate(steps):
                trial=labels.copy();trial[...,c][nodes]+=h
                plus,pm=evaluate(trial,enabled)
                trial=labels.copy();trial[...,c][nodes]-=h
                minus,mm=evaluate(trial,enabled)
                min_orient=min(min_orient,pm,mm)
                diff=(plus-minus)/(2*np.longdouble(h))
                for term,mask in enumerate((active,~active)):
                    selected=mask & enabled
                    np.add.at(fd[idx,term],targets[selected],diff[selected])
        rich=(4*fd[-1]-fd[-2])/3
        coarse_rich=(4*fd[-2]-fd[-3])/3
        row=dict(component=c,colors=int(color.max()+1),active_count=int(cm.free.sum()),melted_count=int(cm.dead.sum()),
                 tied_free_count=int(np.count_nonzero((np.bincount(cm.groups)>1)& ~cm.held)),
                 widths=steps)
        def maxima(err):
            return {key:float(np.max(np.abs(err[mask]),initial=0)) for key,mask in [('active',cm.free),('melted',cm.dead)]}
        def peak(err,mask):
            groups=np.flatnonzero(mask);idx=groups[np.argmax(abs(err[groups]))]
            return dict(node=np.unravel_index(cm.representative[idx],dims),group=int(idx),signed_error=float(err[idx]),
                        analytic_physical=float(pg[idx]),analytic_unit_gauge=float(gg[idx]),fd_physical=float(rich[0,idx]),
                        fd_unit_gauge=float(rich[1,idx]))
        for eps in (1e-8,1e-10):
            analytic=pg+eps*gg
            entry={}
            for label,estimate in [(f'h{i}',v) for i,v in enumerate(fd)]+[('coarse_richardson',coarse_rich),('fine_richardson',rich)]:
                err=estimate[0]+eps*estimate[1]-analytic
                entry[label]=maxima(err)
            err=rich[0]+eps*rich[1]-analytic
            entry['peak_active']=peak(err,cm.free);entry['peak_melted']=peak(err,cm.dead)
            for item in ('peak_active','peak_melted'):
                entry[item]['node']=list(map(int,entry[item]['node']))
            row[str(eps)]=entry
        row['unit_gauge_error']=maxima(rich[1]-gg)
        row['crystal_error']=maxima(rich[0]-pg)
        np.savez(out/f'coordinate_{c}.npz',free=independent,fd=fd,analytic_crystal=pg,analytic_gauge=gg,colors=color,
                 richardson=rich,physical_node=cm.representative,steps=steps)
        result['components'].append(row)
        atomic_json(out/'progress.json',result)
    result['signed_banks_min_chart_orientation_bound']=float(min_orient)
    result['pass']=all(row[str(e)]['fine_richardson'][k]<=1e-10 for row in result['components']
                       for e in (1e-8,1e-10) for k in ('active','melted'))
    atomic_json(record,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('scratch',type=Path);p.add_argument('output',type=Path);p.add_argument('measurement',type=Path)
    args=p.parse_args();run(args.scratch,args.output,args.measurement)
