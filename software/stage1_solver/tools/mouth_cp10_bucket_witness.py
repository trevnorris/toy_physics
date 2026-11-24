"""Read back the ten largest signed-cell contributors with an independent strain calculation."""
import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from mouth_cp10_label_grid import get_cp10, prisms
from mouth_cp10_label_orientation import native_state
from mouth_eulerian_r4_spatial_cell import current, roots
from mouth_exact_checkpoint_control import CHECKPOINT_SHA


def cell_coordinates(axes, point):
    ijk = tuple(int(np.clip(axis.searchsorted(p,side='right')-1,0,len(axis)-2))
                for axis,p in zip(axes,point))
    left = np.array([axis[i] for axis,i in zip(axes,ijk)])
    widths = np.array([axis[i+1]-axis[i] for axis,i in zip(axes,ijk)])
    return ijk,(point-left)/widths,widths


def interpolate(axes, values, point):
    ijk,frac,widths = cell_coordinates(axes,point)
    val = np.zeros(values.shape[3:])
    grad = np.zeros(val.shape+(3,))
    for bits in product((0,1),repeat=3):
        node = values[tuple(i+b for i,b in zip(ijk,bits))]
        basis = np.where(bits,frac,1-frac)
        val += np.prod(basis)*node
        for k in range(3):
            grad[...,k] += (2*bits[k]-1)*basis[(k+1)%3]*basis[(k+2)%3]*node/widths[k]
    return val,grad


def inspect_points(root, observations, output):
    body,state,axes = get_cp10(root)
    mesh = prisms(body,state['solid'])
    key_index = {key:i for i,key in enumerate(mesh[0])}
    time_axis = axes[1]
    length = body.data.length
    rate = np.where(time_axis*length < length-.12,
                    (time_axis*length/(length-.12))**25,1.)-time_axis
    with np.load(root/'r4_nested_labels01/4/labels.npz') as inp:
        fine_axes = [inp[k] for k in 'xta']
        xi = inp['labels']
    native_order = state['field'][2,body.ids].transpose(0,2,1)
    def element(m):
        (i,j,k),frac,_ = cell_coordinates((body.x,body.t,body.a),m)
        return j,i,k,int(frac[0]+frac[1]>1.)
    rows = []
    for item in observations:
        # The physical chart can be recovered from the cylindrical location in
        # the report. The radius map is strictly increasing in its first
        # coordinate. Do not use the reported native material coordinates.
        r,w,a = [item[k] for k in ('r','w','theta')]
        def radial(g):
            phys,_ = current(body,(g,w/length,a))
            return np.hypot(*phys[:2])-r
        point = np.array([brentq(radial,0.,1.,xtol=5e-16),w/length,a])
        hits = roots(mesh,point)
        if len(hits)!=1:
            raise ValueError(f'Expected a unique ordered inverse at {point}, found {len(hits)}')
        material = hits[0]
        j,i,k,tri = element(material)
        e = key_index[(j,i,k,tri)]
        ref,placed = mesh[1][e],mesh[2][e]
        ang = (material[2]-body.a[k])/(body.a[k+1]-body.a[k])
        rc = ref[0]*(1-ang)+ref[1]*ang
        pc = placed[0]*(1-ang)+placed[1]*ang
        edges = (rc[1:,:2]-rc[0,:2]).T
        b = np.linalg.solve(edges,material[:2]-rc[0,:2])
        bary = np.r_[1-b.sum(),b]
        native_der = np.column_stack(((pc[1:]-pc[0]).T@np.linalg.inv(edges),
                                     bary@(placed[1]-placed[0])/(body.a[k+1]-body.a[k])))
        label,spatial_der = interpolate(fine_axes,xi,point)
        mat = native_state(body,material,point,np.linalg.inv(native_der))
        grid = native_state(body,label,point,spatial_der)
        chi,dchi = interpolate(axes,native_order,point)
        minus_dchi = dchi[1]*np.interp(point[1],time_axis,rate)
        physical_readback = dict(chi=float(chi),minus_dchi=float(minus_dchi),J_L=mat['J'],
                                 J_E=grid['J'],W_L=mat['W'],W_E=grid['W'])
        deltas = {k:abs(v-item[k]) for k,v in physical_readback.items()}
        if max(deltas.values())>1e-9:
            raise ArithmeticError(f'Pointwise strain or order disagrees: {deltas}')
        cell=cell_coordinates(fine_axes,point)[0]
        if list(cell)!=item['fine_cell']:
            raise ValueError('Witness is outside its ranked cell')
        verts = xi[tuple(slice(v,v+2) for v in cell)].reshape(-1,3)
        source_prisms = {element(m) for m in verts}
        rows.append(dict(fine_cell=cell,native_prism=[j,i,k,tri],source_prisms_at_nodes=len(source_prisms),
                         weld_straddled_by_spatial_labels=bool(verts[:,0].min()<body.cut<verts[:,0].max()),
                         native_label=material.tolist(),projected_label=label.tolist(),
                         gradient_delta_frobenius=float(np.linalg.norm(np.linalg.inv(native_der)-spatial_der)),
                         independent=physical_readback,max_abs_differences=deltas))
    record = dict(checkpoint_sha256=CHECKPOINT_SHA,rows=rows,
                  maximum_absolute_discrepancies={k:max(x['max_abs_differences'][k] for x in rows)
                                                  for k in physical_readback})
    with output.open('x') as stream:
        json.dump(record,stream,indent=2,allow_nan=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','witnesses','output'):
        p.add_argument(key,type=Path)
    args=p.parse_args()
    observed=json.loads(args.witnesses.read_bytes())
    if len(observed['top10'])!=10 or observed['checkpoint_sha256']!=CHECKPOINT_SHA:
        raise ValueError('Wrong witness record')
    inspect_points(args.root.resolve(),observed['top10'],args.output)
