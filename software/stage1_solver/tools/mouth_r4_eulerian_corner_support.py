"""Identify the actual ordered-support constraints on the limiting melted Q1 cell.

Only read the registered labels and order. A nonzero derivative at a free vertex
is a local first-order controllability statement, not a re-extension of the map.
"""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def measure(sc,case,out):
    area=sc/'r4_eulerian_stage2_01'; pin=area/'pin/software/stage1_solver/tools'
    hashes=json.loads((area/'sources.json').read_bytes())
    if digest(pin/'mouth_cp10_label_grid.py')!=hashes['mouth_cp10_label_grid.py']:
        raise ValueError('Mismatch in pinned port coordinates')
    sys.path.insert(0,str(pin.resolve()))
    from mouth_cp10_label_grid import get_cp10
    b,cp,axes=get_cp10(sc)
    rec=json.loads((area/case/'label_stall.json').read_bytes()); path=area/case/rec['file']
    if digest(path)!=rec['sha256']:raise ValueError('Mismatch in accepted endpoint')
    with np.load(path,allow_pickle=False) as f:
        if not np.array_equal(f['x'],cp['x']) or f['radius']!=4 or f['winding']!=1:
            raise ValueError('Expected fixed radius-four CP10 shape')
        c=f['field'][2,b.ids].transpose(0,2,1); labels=f['labels']; epsilon=float(f['epsilon'])
    bad=rec['label']['last_rejection']; cell=np.array(bad['cell']); raw=bad['local'] if 'code' in bad else bad['chart_fractions']
    lim=np.array(raw,int)
    if not (bad.get('code')==2 or bad.get('chart_lower',1.)<=0.):raise ValueError('Unknown orientation rejection')
    if np.any(np.asarray(raw)!=lim) or np.any((lim!=0)&(lim!=1)):raise ValueError('Expected a vertex')
    support=np.lib.stride_tricks.sliding_window_view(c,(2,2,2)).max(axis=(-3,-2,-1))>0
    nodes=np.stack(np.meshgrid(*[np.arange(i,i+2) for i in cell],indexing='ij'),axis=-1)
    corners=labels[tuple(np.moveaxis(nodes,-1,0))];width=np.array([ax[i+1]-ax[i] for ax,i in zip(axes,cell)])
    q=np.empty((3,3)); center=tuple(lim)
    for d in range(3):
        other=lim.copy();other[d]=1-other[d]; q[:,d]=(corners[center]-corners[tuple(other)])*(2*lim[d]-1)/width[d]
    cof=np.array([np.cross(q[:,(d+1)%3],q[:,(d+2)%3]) for d in range(3)])
    rows=[]
    for bits in product(range(2),repeat=3):
        p=cell+bits; adjacent=[]
        for v in product(range(2),repeat=3):
            i=p-v
            if np.all((i>=0)&(i<support.shape)) and support[tuple(i)]:adjacent.append(i.tolist())
        dif=np.flatnonzero(np.asarray(bits)!=lim); grad=np.zeros(3)
        if not len(dif):grad=np.sum((cof*((2*lim-1)/width)[:,None]),axis=0)
        if len(dif)==1:
            d=dif[0];grad=cof[d]*(1-2*lim[d])/width[d]
        # All components are registered on the wall, entrance and terminal
        # disk. On the polar and equatorial planes only the angular component
        # is pinned. Neither endpoint in this audit touches the collapsed axis.
        if p[0]==0:
            raise ValueError('Separate axis-equivalence audit required')
        held = np.zeros(3,bool)
        if p[0]==c.shape[0]-1 or p[1] in (0,c.shape[1]-1):held[:]=True
        if p[2] in (0,c.shape[2]-1):held[2]=True
        free = ~held if not adjacent else np.zeros(3,bool)
        rows.append(dict(node=p.tolist(),chi=float(c[tuple(p)]),incident_active_cells=adjacent,
                         independent_melted_label=bool(free.any()),melted_free_components=free.tolist(),
                         d_chart_det_d_vertex_label=grad.tolist(),derivative_norm=float(np.linalg.norm(grad)),
                         melt_only_derivative_norm=float(np.linalg.norm(grad[free]))))
    res=dict(case=case,epsilon=epsilon,state_sha256=rec['sha256'],order_grid_sha256=hashlib.sha256(c.tobytes()).hexdigest(),
             limiting_cell=cell.tolist(),vertex=(cell+lim).tolist(),vertices=rows)
    with out.open('x') as f:json.dump(res,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(res,allow_nan=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('sc',type=Path);p.add_argument('case',choices=('primary','secondary'));p.add_argument('out',type=Path)
    a=p.parse_args();measure(a.sc.resolve(),a.case,a.out.resolve())
