"""Measure one-sided inverse-label gradients at internal native material faces."""
import argparse
from collections import defaultdict
from itertools import product
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from mouth_cp10_unfolded_depth import control_context,write_new
from mouth_exact_checkpoint_control import CHECKPOINT_SHA,digest
from mouth_r4_conforming_prisms import TRIANGLES
from mouth_joined_angular_work import bump


def q1_grad(axes,labels,point):
    idx=np.column_stack([(np.searchsorted(axis,p,side='right')-1).clip(0,len(axis)-2)
                         for axis,p in zip(axes,point.T)])
    widths=np.column_stack([ax[i+1]-ax[i] for ax,i in zip(axes,idx.T)])
    uv=(point-np.column_stack([ax[i] for ax,i in zip(axes,idx.T)]))/widths
    out=np.zeros((len(point),3,3))
    for bits in product((0,1),repeat=3):
        weights=np.where(bits,uv,1-uv)
        vertex=labels[tuple((idx+bits).T)]
        for derivative in range(3):
            dirs=[i for i in range(3) if i!=derivative]
            scal=(2*bits[derivative]-1)*weights[:,dirs[0]]*weights[:,dirs[1]]/widths[:,derivative]
            out[:,:,derivative]+=vertex*scal[:,None]
    return out


def measure(root,output):
    old=root/'r4_bucket_unfolded01'
    ctx,solid=control_context(root,old)
    b=ctx[2]
    t,x,a=np.meshgrid(b.t,b.x,b.a,indexing='ij')
    y=np.stack([x,t,a],axis=-1)+solid
    with np.load(old/'control.npz') as f:
        order=f['field'][2,b.ids].transpose(0,2,1)
    chi=RegularGridInterpolator(tuple(b.buffers['field_'+c] for c in 'xta'),order,bounds_error=False,fill_value=None)
    # Square/triangular faces use the center in the common original material coordinates.
    faces=( (0,1,2), (3,4,5), (0,1,3,4), (1,2,4,5), (0,2,3,5) )
    pending={}
    groups=defaultdict(list);checked=0;max_landing_gap=0.
    for j in range(len(b.t)-1):
        for i in range(len(b.x)-1):
            for k in range(len(b.a)-1):
                for half,off in enumerate(TRIANGLES):
                    ids=[(j+dj,i+di,k+l) for l in (0,1) for dj,di in off]
                    corners=np.array([y[q] for q in ids]).reshape(2,3,3)
                    cvals=chi(corners.reshape(6,3))
                    for face in faces:
                        indices=np.asarray(face);key=tuple(sorted(ids[n] for n in face))
                        angle=np.mean(indices//3);bary=np.bincount(indices%3,minlength=3)/len(face)
                        local=corners[0]*(1-angle)+corners[1]*angle
                        pos=bary@local
                        sign=1-2*half
                        der=np.column_stack([sign*(local[1]-local[0])/(b.x[i+1]-b.x[i]),
                                             sign*(local[2]-local[0])/(b.t[j+1]-b.t[j]),
                                             (bary@(corners[1]-corners[0]))/(b.a[k+1]-b.a[k])])
                        inv=np.linalg.inv(der)
                        curr=(pos,inv,float(cvals.min()),float(cvals.max()),(j,i,k,half))
                        if key not in pending:
                            pending[key]=curr;continue
                        prev=pending.pop(key);checked+=1
                        max_landing_gap=max(max_landing_gap,float(np.max(abs(prev[0]-pos))))
                        if min(prev[2],curr[2])==0 and max(prev[3],curr[3])>0:grp='mixed_support_vertices'
                        elif max(prev[3],curr[3])==0:grp='fully_melted_vertices'
                        elif min(prev[2],curr[2])>0:grp='all_ordered_vertices'
                        else:raise ValueError('Unexpected signed order')
                        jump=float(np.linalg.norm(prev[1]-inv))
                        groups[grp].append((jump,pos,(prev[-1],curr[-1]),prev[1],inv,min(prev[2],curr[2]),max(prev[3],curr[3])))
    if max_landing_gap>1e-14:raise ArithmeticError('Welded mesh has face position mismatch')
    result=dict(parent_sha256=CHECKPOINT_SHA,control_sha256=digest(old/'control.npz'),internal_faces=checked,
                max_chart_continuity_error=max_landing_gap,
                statistic='Unweighted face centers, Frobenius norm of one-sided inverse-label gradient difference in normalized chart coordinates. Support classes refer to both incident elements\' material vertices sampled in the exact order field.',groups={})
    for name,values in groups.items():
        jumps=np.array([x[0] for x in values])
        pos=np.array([x[1] for x in values]);sensitivity={}
        for lvl in (1,2,4):
            grad=[]
            for kind in ('unshifted','shifted'):
                path=output.parent/f'labels_{lvl}'/(kind+'.npz')
                with np.load(path) as f:
                    grad.append(q1_grad([f[k] for k in ('x','zero_label_t','a')],f['zero_labels'],pos))
            change=np.linalg.norm(grad[1]-grad[0],axis=(1,2))
            sensitivity[str(lvl)]=dict(rms=float(np.sqrt(np.mean(change**2))),mean=float(change.mean()))
        def witness(idx):
            val,chart,keys,left,right,low,high=values[idx]
            x,t,a=chart; w=t*b.sh.length
            r=x*b.sh.radius(w)[0]+x**3*b.epsilon*bump(w,b.sh.length)[0]*(3*np.cos(a)**2-1)/2
            return dict(native_gradient_jump=val,grid_shift_change_in_gradient=float(change[idx]),
                        positive_order_in_elements=bool(high>0),min_vertex_order=low,max_vertex_order=high,
                        chi_at_face=float(chi(chart[None,:])[0]),physical_r_w_theta=[float(r),float(w),float(a)],
                        forward_chart_determinants=[float(1/np.linalg.det(m)) for m in (left,right)],
                        native_left_inverse_gradient=left.tolist(),native_right_inverse_gradient=right.tolist(),
                        unshifted_4x_gradient=grad[0][idx].tolist(),shifted_4x_gradient=grad[1][idx].tolist(),
                        material_prisms=keys)
        result['groups'][name]=dict(faces=len(values),nonzero_above_1e_10=int(sum(jumps>1e-10)),
             median=float(np.median(jumps)),rms=float(np.sqrt(np.mean(jumps*jumps))),grid_shift_gradient_change=sensitivity,
             maximum_native_jump=witness(jumps.argmax()),strongest_shift_at_4x=witness(change.argmax()))
    write_new(output,result)
    print(result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();measure(a.root.resolve(),a.output.resolve())
