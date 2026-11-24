"""One-sided physical state at the limiting corner of an accepted held R4 label map.

Read saved binary64 nodal coefficients as exact numbers for the chart derivative.
No label, order, carrier, or shape equation is advanced. The recorded failed
trial only supplies the corner address, never its accepted determinant.
"""
import argparse
import hashlib
import json
import sys
from itertools import product
from pathlib import Path

import mpmath as mp
import numpy as np


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def num(x):
    return mp.mpf(float(x))


def profile(body, t, ref=False):
    """Evaluate the exact stored binary spline and shape coefficients in high precision."""
    n=body.data.spline_n
    i=min(max(int(t*n),0),n-1); h=t-mp.mpf(i)/n
    a,b,c,d=map(num,body.buffers['spline'][i]); off=4 if ref else 0
    e,m,q,k=map(num,body.buffers['shape'][off:off+4]); v=t*(1-t)
    r=4*((((a*h+b)*h+c)*h+d)+e*t+4*v*(m+q*(2*t-1))+16*k*v*v)
    rt=4*((3*a*h+2*b)*h+c+e+4*m*(1-2*t)+4*q*(-6*t*t+6*t-1)+32*k*v*(1-2*t))
    return r,rt


def point(body, axes, labels, order, cell, fractions, collar):
    low=[num(ax[i]) for ax,i in zip(axes,cell)]
    widths=[num(ax[i+1])-v for ax,i,v in zip(axes,cell,low)]
    f=list(map(num,fractions)); pos=[v+w*s for v,w,s in zip(low,widths,f)]
    xi=mp.matrix(3,1); Q=mp.matrix(3,3); chi=mp.mpf(0); vertex_order=[]
    for corner in product((0,1),repeat=3):
        node=tuple(j+k for j,k in zip(cell,corner)); co=num(order[node]);vertex_order.append(float(co))
        v=list(map(num,labels[node])); b=[s if k else 1-s for k,s in zip(corner,f)]; w=mp.fprod(b)
        chi+=co*w
        for i in range(3):
            xi[i]+=w*v[i]
            for j in range(3):
                Q[i,j]+=v[i]*(2*corner[j]-1)*b[(j+1)%3]*b[(j+2)%3]/widths[j]
    x,t,theta=pos; L=num(body.data.length); xx,tt,aa=xi
    r,rt=profile(body,t); y=(t-num(.06))/num(.56); b=bt=0
    if 0<y<1:
        b=16*y*y*(1-y)**2*num(body.epsilon); bt=32*y*(1-y)*(1-2*y)*num(body.epsilon)/num(.56)
    p=(3*mp.cos(theta)**2-1)/2; R=x*r+x**3*b*p
    B=mp.matrix([[r+3*x*x*b*p,x*rt+x**3*bt*p,-3*x**3*b*mp.sin(theta)*mp.cos(theta)],
                 [0,L,0],[0,0,R]])
    if xx<=num(body.cut):
        s,st=profile(body,tt,True); ratio=num(body.data.core_ratio)
        C=mp.matrix([[s*ratio,xx*st*ratio,0],[0,num(body.data.ref_length),0],[0,0,xx*s*ratio]])
        side='core'
    else:
        C=collar(body,xx,tt);side='collar'
    grad=Q*(B**-1); inv=C*grad
    ih=C[2,2]*mp.sin(aa)/(R*mp.sin(theta)); Jchart=mp.det(Q)
    Ji=Jchart*abs(mp.det(C))*ih/mp.det(B)
    if Ji<=0: raise ArithmeticError('Saved corner has nonpositive physical orientation')
    J=1/Ji; F=inv**-1; W=(sum(v*v for v in F)+1/ih**2-4)/2-mp.log(J)+(J-1)**2/4
    sv=mp.svd(grad,compute_uv=False); principal=list(mp.svd(inv,compute_uv=False))+[abs(ih)]
    minimum,maximum=min(principal),max(principal)
    return dict(cell=list(map(int,cell)),local=list(map(float,f)),material_side=side,
                physical_chart=list(map(float,pos)),physical=dict(radius=float(R),theta=float(theta),
                cylindrical_radius=float(R*mp.sin(theta)),axial_brane_coordinate=float(R*mp.cos(theta)),w=float(t*L)),
                label=list(map(float,xi)),distance_from_weld_in_label=float(xx-num(body.cut)),
                chi=float(chi),max_vertex_chi=max(vertex_order),positive_vertices=sum(v>0 for v in vertex_order),
                chart_det=float(Jchart),Jinv=float(Ji),J=float(J),W_formal=float(W),chi_W=float(chi*W),
                chi_Jinv_W=float(chi*Ji*W),grad_xi_frobenius=float(mp.norm(grad)),grad_xi_spectral=float(sv[0]),
                condition_F=float(maximum/minimum),min_singular_Finv=float(minimum),max_singular_Finv=float(maximum))


def nearest(body, axes, order, node, physical):
    loc=np.argwhere(order>0);x,t,a=(ax[loc[:,d]] for d,ax in enumerate(axes));L=body.data.length
    rm=body.sh.radius(t*L)[0]; y=(t-.06)/.56
    b=np.where((y>0)&(y<1),16*y*y*(1-y)**2*body.epsilon,0.)
    r=x*rm+x**3*b*(3*np.cos(a)**2-1)/2
    coords=np.column_stack((r*np.sin(a),r*np.cos(a),t*L))
    delta=coords-np.array([physical['cylindrical_radius'],physical['axial_brane_coordinate'],physical['w']]); dist=np.linalg.norm(delta,axis=1)
    i=int(np.argmin(dist)); at=tuple(loc[i]); col=np.flatnonzero(order[node[0],:,node[2]]>0)
    windows=np.lib.stride_tricks.sliding_window_view(order,(2,2,2)); cells=np.argwhere(windows.max(axis=(-3,-2,-1))>0)
    separation=np.maximum.reduce((cells-node,node-cells-1,np.zeros_like(cells))).max(axis=1)
    return dict(nearest_strictly_ordered_vertex=dict(node=loc[i].tolist(),chi=float(order[at]),
                physical_distance=float(dist[i]),physical_cylindrical_axial_depth=coords[i].tolist()),
                nearest_active_cell_index_separation=int(separation.min()),nearest_active_cell_count=int((separation==separation.min()).sum()),
                last_positive_vertex_same_radial_angular_column=dict(depth_index=int(col[-1]),
                w=float(axes[1][col[-1]]*L),chi=float(order[node[0],col[-1],node[2]])) if col.size else None)


def main(root,area,case,output):
    pin=area/'pin/software/stage1_solver/tools'; frozen=json.loads((area/'sources.json').read_bytes())
    modules=('mouth_cp10_label_grid.py','mouth_r4_refinement_review.py','mouth_r4_eulerian_release.py',
             'mouth_free_shape_solid.py','mouth_free_shape_solid.cpp','mouth_r4_eulerian_solid.cpp')
    pinned={m:digest(pin/m) for m in modules}
    if any(pinned[m]!=frozen[m] for m in modules):raise ValueError('Pinned production source differs')
    sys.path.insert(0,str(pin.resolve()))
    from mouth_cp10_label_grid import get_cp10
    from mouth_r4_refinement_review import _collar_reference
    b,cp,axes=get_cp10(root)
    doc=json.loads((area/case/'label_stall.json').read_bytes());state=area/case/doc['file'];g=doc['label']
    if digest(state)!=doc['sha256']:raise ValueError('State digest differs')
    with np.load(state,allow_pickle=False) as st:
        if st['radius']!=4 or st['winding']!=1 or st['iteration']!=10 or not np.array_equal(st['x'],cp['x']):
            raise ValueError('Expected held CP10 geometry')
        labels=st['labels'];order=st['field'][2,b.ids].transpose(0,2,1)
        if st['epsilon']!=g['epsilon']:raise ValueError('Gauge does not match')
    if b.sh.length!=doc['depth']:raise ValueError('Reconstructed physical depth differs from guarded state')
    bad=g['last_rejection']; local=np.array(bad['local'] if 'code' in bad else bad['chart_fractions'])
    cell=np.array(bad['cell']); node=cell+local.astype(int)
    if not np.all((local==0)|(local==1)):raise ValueError('Rejection witness is not a corner')
    if not (bad.get('code')==2 or bad.get('chart_lower',1.)<=0.):raise ValueError('Unknown orientation rejection')
    rows=[]
    with mp.workdps(80):
        for vertex in product(range(2),repeat=3):
            c=node-vertex
            if np.all((c>=0)&(c<np.array(order.shape)-1)):
                rows.append(point(b,axes,labels,order,c,vertex,_collar_reference))
        worst=g['worst_global']; q=point(b,axes,labels,order,worst['cell'],worst['local'],_collar_reference)
    delta={k:q[v]-worst[k] for k,v in (('Jinv','Jinv'),('chi','chi'),('W','W_formal'),('grad_xi_norm','grad_xi_frobenius'))}
    if any(abs(v)>3e-12*max(1,abs(worst[k])) for k,v in delta.items()):raise ArithmeticError(f'Physical Gauss cross-check: {delta}')
    limiting=next(r for r in rows if r['cell']==list(cell))
    if limiting['chart_det']<g['min_chart_Jinv_bound']-8e-14:raise ArithmeticError('Reported bound does not enclose limiting saved corner')
    with mp.workdps(48):
        control=point(b,axes,labels,order,cell,local,_collar_reference)
    agreement={k:control[k]-limiting[k] for k in ('chart_det','Jinv','W_formal','grad_xi_frobenius')}
    result=dict(case=case,epsilon=g['epsilon'],state_sha256=doc['sha256'],source_hashes=pinned,
        norm='Frobenius with undimensionalized labels and orthonormal physical frame as in the production kernel; spectral norm also supplied.',
        interpretation='Physical Jinv contains both reference and physical volumes and the azimuthal hoop. In pure melt W is a formal extension of the crystal law only. The physical crystal energy density there is exactly zero.',
        node=node.tolist(),accepted_incident_cells=rows,accepted_limiting_cell=limiting,
        rejected_trial_only=bad,production_accepted_chart_lower=g['min_chart_Jinv_bound'],
        production_accepted_active_chart_lower=g['min_chart_Jinv_active_bound'],
        accepted_worst_quadrature=worst,accepted_worst_active_quadrature=g['worst_active'],
        quadrature_readback_delta=delta,
        decimal_precision_48_vs_80=agreement,**nearest(b,axes,order,node,limiting['physical']))
    with output.open('x') as f: json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(dict(case=case,node=node.tolist(),limiting=limiting,incident_active=sum(r['positive_vertices']>0 for r in rows),
                          nearest=result['nearest_strictly_ordered_vertex'],quadrature_delta=delta),allow_nan=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('area',type=Path)
    p.add_argument('case',choices=('primary','secondary'));p.add_argument('output',type=Path)
    args=p.parse_args();main(args.root.resolve(),args.area.resolve(),args.case,args.output.resolve())
