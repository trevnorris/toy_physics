"""Construct the unfolded comparison displacement on checkpoint 10's original mouth."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_cp10_label_gauge import chart_volume_bounds
from mouth_cp10_label_grid import BITS, allnode_inverse, get_cp10, harmonic, matrix_for, prisms
from mouth_cp10_nested_labels import depth_derivatives, nested_axis
from mouth_cp10_spatial_crystal import load_context
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read


def write(path,data):
    with path.open('x') as file:
        json.dump(data,file,indent=2,allow_nan=False)
        file.write('\n')


def determinant_extrema(body,displacement):
    # The normalized determinant of P1(triangle) × P1(angle) placement is
    # affine in triangle barycentrics and quadratic in the polar fraction.
    keys,reference,corners,_,_=prisms(body,np.asarray(displacement,dtype=np.longdouble))
    ref=reference.astype(np.longdouble)
    volume=np.einsum('ei,ei->e',np.cross(ref[:,0,1]-ref[:,0,0],ref[:,0,2]-ref[:,0,0]),ref[:,1,0]-ref[:,0,0])
    if np.any(volume<=0):raise ArithmeticError('Native material prism orientation')
    def values(y):
        v=corners[:,0]*(1-y)+corners[:,1]*y
        normal=np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0])
        return np.einsum('ei,evi->ev',normal,corners[:,1]-corners[:,0])/volume[:,None]
    low,mid,high=[values(y) for y in (0,.5,1)]
    curvature=2*((high-mid)-(mid-low))
    slope=high-low-curvature
    critical=np.divide(-slope,2*curvature,out=np.zeros_like(slope),where=curvature!=0)
    interior=(critical>0)&(critical<1)
    turning=low+critical*(slope+critical*curvature)
    lower=np.minimum(low,high)
    upper=np.maximum(low,high)
    lower[interior]=np.minimum(lower[interior],turning[interior])
    upper[interior]=np.maximum(upper[interior],turning[interior])
    e,vertex=np.unravel_index(lower.argmin(),lower.shape)
    return dict(min_normalized_chart_J=float(lower[e,vertex]),max_normalized_chart_J=float(upper.max()),
                min_prism=keys[e],min_barycentric_vertex=int(vertex),prisms=len(keys),
                nonpositive_prisms=int((lower.min(axis=1)<=0).sum()))


def prepare(root,area):
    area.mkdir(parents=True,exist_ok=True)
    target=area/'control.npz'
    if target.exists():raise FileExistsError(target)
    body,s,_=get_cp10(root)
    with s:
        fields={k:s[k] for k in s.files}
    scale=.5
    old=fields['solid']
    solid=scale*old
    packed=body.pack(solid)
    if not np.array_equal(body.unpack(packed),solid) or not np.array_equal(packed,scale*fields['q']):
        raise ArithmeticError('Displacement failed a registered wall, angular or shared-axis condition')
    cert=determinant_extrema(body,solid)
    if cert['min_normalized_chart_J']<=0:raise ArithmeticError(str(cert))
    # Check saved predecessors with the same unclipped prism-polynomial test.
    history={}
    for index in range(6,11):
        name=f'state_{index:05d}.npz'
        saved=root/'r4_melted_survival01/coarse'/name
        with np.load(saved) as prior:
            if int(prior['radius'])!=4 or int(prior['winding'])!=1:raise ValueError('Wrong predecessor')
            history[str(index)]=determinant_extrema(body,prior['solid'])
    fields.update(solid=solid,q=packed,parent_sha256=CHECKPOINT_SHA,scale_of_cp10_displacement=scale)
    np.savez_compressed(target,**fields)
    with np.load(root/'r4_melted_survival01/coarse/state_00010.npz') as cp, np.load(target) as control:
        passthrough=[k for k in cp.files if k not in ('solid','q')]
        if not all(np.array_equal(cp[k],control[k]) for k in passthrough):
            raise ValueError('Changed noncrystal state')
    report=dict(parent_sha256=CHECKPOINT_SHA,control_sha256=digest(target),factor=scale,
                rule='Half the displacement from the identity material chart to the CP10 map, at every original node. Native piecewise-smooth prism interpolation and existing reference metrics retained.',
                native_orientation=cert,native_prior_chart_orientation=history,
                unchanged_checkpoint_arrays=passthrough,max_displacement=float(abs(solid).max()),
                entrance_wall_terminal_offset_max=float(max(abs(solid[[0,-1]]).max(),abs(solid[:,-1]).max())),
                axis_radial_offset_max=float(abs(solid[:,0,:,0]).max()),
                polar_reflection_max=float(abs(solid[:,:,[0,-1],2]).max()))
    write(area/'control.json',report)
    print(json.dumps(report,allow_nan=False),flush=True)


def load_control_context(root,area):
    """Signed-geometry context for the verified comparison; preserves the CP10 physical chart."""
    check=read(area/'control.json')
    if check['control_sha256']!=digest(area/'control.npz') or check['parent_sha256']!=CHECKPOINT_SHA:
        raise ValueError('Unfolded preparation changed')
    ctx=load_context(root,root/'r4_cp10_labels01/labels.npz')
    labels=area/'1/labels.npz'
    if digest(labels)!=read(area/'1/result.json')['sha256']: raise ValueError('Baseline label table changed')
    with np.load(labels,allow_pickle=False) as data:
        if str(data['source_sha256'])!=check['control_sha256'] or not np.array_equal(data['order'],ctx[-1]):
            raise ValueError('Unfolded model context or order differs')
        return (*ctx[:-2],data['labels'],data['order'])


def orientation(axes,labels,melt):
    minimum={'active_cells':np.inf,'melted_cells':np.inf}
    for i in range(0,len(axes[0])-1,4):
        j=min(i+4,len(axes[0])-1)
        section=melt[i:j]
        test=chart_volume_bounds((axes[0][i:j+1],*axes[1:]),labels[i:j+1],section)
        for key in minimum:
            minimum[key]=min(minimum[key],test[key]['bernstein_lower_diagnostic'])
    if min(minimum.values())<=0:raise ArithmeticError(f'Nonpositive orientation bound: {minimum}')
    return minimum


def generate(root,area,level):
    directory=area/str(level)
    directory.mkdir(exist_ok=False)
    meta=read(area/'control.json');sha=digest(area/'control.npz')
    if sha!=meta['control_sha256'] or meta['parent_sha256']!=CHECKPOINT_SHA:
        raise ValueError('Unfolded input provenance')
    body,cp,coarse=get_cp10(root)
    cp.close()
    with np.load(area/'control.npz') as control:
        solid=control['solid']
        raw_order=control['field'][2,body.ids].transpose(0,2,1)
    axes=[nested_axis(x,level) for x in coarse]
    coords=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
    order=raw_order if level==1 else RegularGridInterpolator(coarse,raw_order)(coords)
    old_path=(root/'r4_cp10_labels01/labels.npz') if level==1 else (root/f'r4_nested_labels01/{level}/labels.npz')
    with np.load(old_path) as old:
        if not np.array_equal(order,old['order']):raise ValueError('Frozen order mismatch')
        mask=old['active_stencil']
        melt=old['melted_cells']
    active=np.any(np.lib.stride_tricks.sliding_window_view(order,(2,2,2))>0,axis=(-1,-2,-3))
    expected=np.zeros_like(mask)
    for bits in BITS:
        expected[tuple(slice(b,b+n) for b,n in zip(bits,active.shape))]|=active
    if not np.array_equal(expected,mask) or not np.array_equal(~active,melt):
        raise ValueError('Frozen order support differs')
    mesh=prisms(body,solid)
    inv,count,dups,roundtrip=allnode_inverse(axes,coords,mesh,mask)
    if np.any(count[mask]!=1) or dups:raise ArithmeticError('Ambiguous unfolded active labels')
    if level==1:
        labels,_,_,gauge=harmonic(matrix_for(body,axes,melt),coords,inv,mask)
    else:
        with np.load(area/'1/labels.npz') as native:
            if str(native['source_sha256'])!=sha:raise ValueError('Coarse label identity')
            labels=RegularGridInterpolator(coarse,native['labels'])(coords)
        labels[mask]=inv[mask]
        gauge={'extension':'Native harmonic extension prolonged only outside positive-order cell closure'}
    fixed=np.zeros_like(labels,dtype=bool)
    fixed[:,[0,-1]]=True
    fixed[-1]=True
    fixed[0,:,:,0]=True
    fixed[:,:,[0,-1],2]=True
    snap=float(np.max(abs(labels[fixed]-coords[fixed])))
    if snap>2e-9:raise ArithmeticError(f'Boundary trace {snap}')
    labels[fixed]=coords[fixed]
    tangents,transport=depth_derivatives(body,solid,labels,mask,axes)
    det=orientation(axes,labels,melt)
    path=directory/'labels.npz'
    np.savez_compressed(path,x=axes[0],t=axes[1],a=axes[2],order=order,labels=labels,
                        active_stencil=mask,melted_cells=melt,source_sha256=sha,
                        parent_sha256=CHECKPOINT_SHA,dxi_dt_minus=tangents[0],dxi_dt_plus=tangents[1])
    result=dict(level=level,source_sha256=sha,sha256=digest(path),native_cp10_labels_sha256=digest(old_path),
                active_nodes=int(mask.sum()),missing=0,duplicate=0,max_landing_residual=roundtrip,
                max_difference_before_registered_snap=snap,transport=transport,gauge=gauge,
                lower_bernstein_bound_chart_Jinv=det)
    write(directory/'result.json',result)
    print(json.dumps(result,allow_nan=False),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('area',type=Path)
    p.add_argument('action',choices=('prepare','1','2','4'))
    args=p.parse_args()
    root,area=args.root.resolve(),args.area.resolve()
    if args.action=='prepare':prepare(root,area)
    else:generate(root,area,int(args.action))
