"""Prolong checkpoint-10 order on nested spatial grids and invert its same welded label map.

At zero gauge only the closure of the ordered cells is needed. Outside it we
retain the prolongation of the saved harmonic chart as an unevaluated extension.
The rate arrays distinguish the two material derivatives at an element boundary.
"""
import argparse
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_cp10_label_grid import BITS, allnode_inverse, get_cp10, prisms
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_conforming_prisms import TRIANGLES


def nested_axis(old, factor):
    fractions = np.arange(factor, dtype=float) / factor
    interior = old[:-1, None] + np.diff(old)[:, None] * fractions
    return np.r_[interior.reshape(-1), old[-1]]


def depth_derivatives(body, solid, label, mask, axes):
    """Invert the material prism tangent and choose by sign of the depth step."""
    indices = np.column_stack(np.nonzero(mask))
    coordinate = label[mask]
    grid = (body.x, body.t, body.a)
    candidates = []
    for v, axis in zip(coordinate.T, grid):
        left = np.searchsorted(axis, v - 5e-13, side='left') - 1
        right = np.searchsorted(axis, v + 5e-13, side='right') - 1
        candidates.append((left.clip(0,len(axis)-2), right.clip(0,len(axis)-2)))
    count = np.zeros((len(indices),2), int)
    derivatives = np.zeros((len(indices),2,3))
    max_spread = 0.
    for sides in BITS:
        cells = np.stack([c[s] for s,c in zip(sides,candidates)],axis=-1)
        used = np.ones(len(indices), dtype=bool)
        for s, (low, high) in zip(sides,candidates):
            if s: used &= low != high
        ids = np.flatnonzero(used)
        if not len(ids): continue
        i,j,k = cells[ids].T
        widths = np.stack([ax[idx+1]-ax[idx] for ax,idx in zip(grid,(i,j,k))],axis=1)
        uvw = (coordinate[ids] - np.column_stack((grid[0][i],grid[1][j],grid[2][k])))/widths
        for tri, offsets in enumerate(TRIANGLES):
            sig = 1 - 2 * tri
            rx,zt,y = uvw.T
            bx,bt = (rx,zt) if tri == 0 else (1-rx,1-zt)
            bary = np.column_stack((1-bx-bt,bx,bt))
            hit = bary.min(axis=1) >= -8e-11
            if not np.any(hit): continue
            h = ids[hit]
            y0 = y[hit]
            bar = bary[hit]
            ds = widths[hit]
            corners = [np.stack([solid[(j+dj)[hit],(i+di)[hit],(k+l)[hit]] for dj,di in offsets],axis=1) for l in (0,1)]
            delta = corners[1]-corners[0]
            displ = corners[0]+y0[:,None,None]*delta
            grad = np.empty((len(h),3,3))
            grad[:,:,0] = sig*(displ[:,1]-displ[:,0])/ds[:,:1]
            grad[:,:,1] = sig*(displ[:,2]-displ[:,0])/ds[:,1:2]
            grad[:,:,2] = np.einsum('pk,pkc->pc',bar,delta)/ds[:,2:]
            grad += np.eye(3)
            valid = abs(np.linalg.det(grad)) >= 1e-12
            if not np.all(valid):
                h,ds,bar,y0,grad=h[valid],ds[valid],bar[valid],y0[valid],grad[valid]
            if not len(h): continue
            rhs = np.zeros((len(h),3,1))
            rhs[:,1,0]=1.
            vel = np.linalg.solve(grad,rhs)[:,:,0]
            vx,vt,va = (vel/ds).T
            planes = np.column_stack((bar,y0,1-y0))
            speed = np.column_stack((-sig*(vx+vt),sig*vx,sig*vt,va,-va))
            for n,sign in enumerate((-1,1)):
                selected = np.all((planes>=8e-11) | (sign*speed>=-1e-7),axis=1)
                hh = h[selected]; vv = vel[selected]
                already = count[hh,n] > 0
                if np.any(already):
                    max_spread = max(max_spread,float(np.max(abs(derivatives[hh[already],n]-vv[already]))))
                derivatives[hh[~already],n] = vv[~already]
                count[hh,n] += 1
    required = np.column_stack((indices[:,1]>0, indices[:,1]<len(axes[1])-1))
    missing = (count == 0) & required
    if np.any(missing) or max_spread>1e-7:
        first = np.flatnonzero(missing.any(axis=1))[:5]
        raise ArithmeticError(f'Unresolved signed tangent at {indices[first].tolist()}; spread {max_spread}')
    info=dict(max_same_sign_candidate_spread=max_spread,
              max_two_sided_jump=float(np.max(abs(np.diff(derivatives[count.min(axis=1)>0],axis=1)),initial=0.)))
    out=[]
    for n in (0,1):
        arr=np.zeros((*mask.shape,3))
        arr[mask] = derivatives[:,n,:]
        out.append(arr)
    return out, info


def main(root, dest, factor):
    case = dest / str(factor)
    case.mkdir(parents=True,exist_ok=False)
    base = root/'r4_cp10_labels01'
    orig = read(base/'reconstruction.json')
    if orig['label_sha256'] != digest(base/'labels.npz') or orig['source'] != CHECKPOINT_SHA:
        raise ValueError('Original label or source hash changed')
    body, state, old = get_cp10(root)
    axes = [nested_axis(v,factor) for v in old]
    coords = np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
    with np.load(base/'labels.npz',allow_pickle=False) as native:
        order = RegularGridInterpolator(old,native['order'])(coords)
        if not np.array_equal(order[::factor,::factor,::factor],native['order']):
            raise ArithmeticError('Nested knots did not preserve original order')
        ext = RegularGridInterpolator(old,native['labels'])(coords)
    active_cells = np.any(np.lib.stride_tricks.sliding_window_view(order,(2,2,2))>0,axis=(-3,-2,-1))
    mask = np.zeros(order.shape,bool)
    for bits in BITS:
        mask[tuple(slice(b,b+d) for b,d in zip(bits,active_cells.shape))] |= active_cells
    mesh = prisms(body,state['solid'])
    inv,counts,duplicates,roundtrip = allnode_inverse(axes,coords,mesh,mask)
    report=dict(factor=factor,source=CHECKPOINT_SHA,source_labels=orig['label_sha256'],
                axes_shape=list(order.shape),active_cells=int(active_cells.sum()),required_nodes=int(mask.sum()),
                max_chart_roundtrip=roundtrip,missing=np.argwhere(mask & (counts==0)).tolist(),
                multiple=[dict(node=list(map(int,np.unravel_index(i,order.shape))),roots=[v.tolist() for v in vs]) for i,vs in duplicates.items()])
    if np.any(counts[mask]!=1):
        (case/'issue.json').write_text(json.dumps(report,indent=2)+'\n')
        raise ArithmeticError('Active-cell closure has missing or multiple distinct reference coordinates')
    ext[mask] = inv[mask]
    # Enforce exactly the same registered components as the CP10 chart.
    fixed=np.zeros_like(ext,dtype=bool)
    fixed[:,[0,-1]]=True
    fixed[-1]=True
    fixed[0,:,:,0]=True
    fixed[:,:,[0,-1],2]=True
    check=fixed&mask[...,None]
    snap = float(np.max(abs(ext[check]-coords[check])))
    if snap>2e-9:raise ArithmeticError(f'Registered label error {snap}')
    ext[fixed]=coords[fixed]
    derivatives, report['transport'] = depth_derivatives(body,state['solid'],ext,mask,axes)
    with np.load(base/'labels.npz',allow_pickle=False) as native:
        overlap=mask[::factor,::factor,::factor] & native['active_stencil']
        common=ext[::factor,::factor,::factor][overlap]-native['labels'][overlap]
        report['max_active_nested_knot_delta']=float(np.max(abs(common),initial=0.))
    if report['max_active_nested_knot_delta']>3e-12:
        raise ArithmeticError('Nested held labels disagree at prior active knots')
    report['max_registered_snap']=snap
    temp=case/'labels.partial.npz'
    np.savez_compressed(temp,x=axes[0],t=axes[1],a=axes[2],labels=ext,order=order,
                        active_stencil=mask,melted_cells=~active_cells,source_sha256=CHECKPOINT_SHA,
                        dxi_dt_minus=derivatives[0],dxi_dt_plus=derivatives[1])
    file=case/'labels.npz'
    temp.rename(file)
    report['sha256']=digest(file)
    (case/'result.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('dest',type=Path);p.add_argument('factor',type=int,choices=(2,4))
    a=p.parse_args();main(a.root.resolve(),a.dest.resolve(),a.factor)
