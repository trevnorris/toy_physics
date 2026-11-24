"""One-sided source-label depth transport, two-metric welded ring and orientation on CP10."""
from argparse import ArgumentParser
from pathlib import Path
import json
from itertools import product

import numpy as np
from scipy.optimize import root as find_zero
from mouth_cp10_label_grid import BITS,get_cp10,prisms
from mouth_eulerian_r4_spatial_cell import current,roots
from mouth_exact_checkpoint_control import CHECKPOINT_SHA,digest
from mouth_joined_angular_work import bump
from mouth_r4_conforming_prisms import TRIANGLES,_geometry
from mouth_supercore_shape import preparation


def source_t(body,solid,labels,stencil,axes):
    grids=(body.x,body.t,body.a);rate=np.zeros((*stencil.shape,2,3));count=np.zeros((*stencil.shape,2),int)
    turn=0.; reverse=0.; missing=[]
    for node in zip(*np.nonzero(stencil)):
        l=labels[node]
        cells=[np.flatnonzero((v[:-1]<=q+5e-13)&(v[1:]>=q-5e-13)) for v,q in zip(grids,l)]
        vals=[[],[]]
        for i,j,k in product(*cells):
            ds=np.array([np.diff(v)[n] for v,n in zip(grids,(i,j,k))])
            uvw=(l-np.array([v[n] for v,n in zip(grids,(i,j,k))]))/ds
            x,t,y=uvw
            for tri,off in enumerate(TRIANGLES):
                sig=1-2*tri; bx=x if not tri else 1-x;bt=t if not tri else 1-t
                bary=np.array([1-bx-bt,bx,bt])
                if min(bary)<-8e-11:continue
                verts=np.array([[solid[j+dj,i+di,k+z] for dj,di in off] for z in (0,1)])
                p=(1-y)*verts[0]+y*verts[1]
                grad=np.eye(3)+np.column_stack((sig*(p[1]-p[0])/ds[0],sig*(p[2]-p[0])/ds[1],bary@(verts[1]-verts[0])/ds[2]))
                if abs(np.linalg.det(grad))<1e-12:continue
                vel=np.linalg.solve(grad,[0.,1.,0.])
                vb=sig*vel[:2]/ds[:2]
                for n,side in enumerate((-1,1)):
                    test=np.r_[bary,y,1-y]
                    slope=side*np.r_[-sum(vb),vb,vel[2]/ds[2],-vel[2]/ds[2]]
                    if np.min(slope[test<8e-11],initial=0.)< -1e-7:continue
                    vals[n].append(vel)
        for n in range(2):
            if vals[n]:
                rate[node+(n,)]=vals[n][0]
                count[node+(n,)]=len(vals[n])
                turn=max(turn,float(np.max(np.abs(np.asarray(vals[n])-vals[n][0]))))
            elif not (node[1]==0 and n==0 or node[1]==len(axes[1])-1 and n==1):
                missing.append((node,n))
        if vals[0] and vals[1]:reverse=max(reverse,float(max(abs(rate[node+(0,)]-rate[node+(1,)]))))
    return rate,count,dict(unresolved_nodes=[(tuple(map(int,nd)),side) for nd,side in missing],
        max_same_side_candidate_spread=turn,max_distinct_signed_side_difference=reverse)


def bond(body,solid,axes,labels,order):
    key=(0,body.cut_node,2,0);re,xc,_=_geometry(body,solid,key)
    mid=lambda arr: np.asarray((arr[0,0]+arr[0,2]+arr[1,0]+arr[1,2])/4,float)
    ref=mid(re); exact=mid(xc)
    init=preparation(4)[0]; core=body.cut*body.data.core_ratio*float(init.radius(ref[1]*init.length)[0])
    re_t=body.buffers['reg_t'];re_e=body.buffers['reg_e']
    collar=4+np.interp(ref[1],re_t,body.buffers['labels'].reshape(len(re_t),len(re_e),2)[:,-1,0])
    phys,_=current(body,exact); exact_lever=phys[0]/np.sin(ref[2])
    mesh=prisms(body,solid);out=[]
    for level in [1,2,4,8,32]:
        cache={}; observed=[]
        def interpolate(p,values=False):
            indices=[];h=[];local=[]
            for a,z in zip(axes,p):
                ci=np.clip(np.searchsorted(a,z,side='right')-1,0,len(a)-2)
                d=(a[ci+1]-a[ci])/level
                n=np.clip(int((z-a[ci])/d),0,level-1)
                indices.append((ci,n));h.append(d);local.append((z-(a[ci]+n*d))/d)
            weights=np.array([[1-y,y] for y in local]);total=np.zeros(3);chis=0.
            for bit in BITS:
                point=[];parent_weight=[];parent_offset=[]
                for a,(ci,n),d,b in zip(axes,indices,h,bit):
                    point.append(a[ci]+(n+b)*d)
                    parent_weight.append((n+b)/level);parent_offset.append(ci)
                if level==1:
                    vi=tuple(c+d for c,d in zip(parent_offset,bit))
                    xi=labels[vi]; c=order[vi]
                else:
                    idx=tuple(point)
                    if idx not in cache:
                        hits=roots(mesh,np.array(point))
                        if len(hits)!=1:raise ValueError(f'Ambiguous local ring inverse {idx}')
                        c=0.
                        for br in BITS:
                            wt=np.prod(np.where(br,parent_weight,1-np.array(parent_weight)))
                            c+=wt*order[tuple(p+b for p,b in zip(parent_offset,br))]
                        cache[idx]=(hits[0],c)
                    xi,c=cache[idx]
                observed.append(c)
                w=np.prod(weights[np.arange(3),bit]);total+=w*xi;chis+=w*c
            return (total,chis) if values else total
        root=find_zero(lambda p:interpolate(p)-ref,exact,tol=1e-12)
        if max(abs(root.fun))>3e-13:raise ArithmeticError('Failed to locate welded material point')
        P,B=current(body,root.x); lever=P[0]/np.sin(ref[2])
        out.append(dict(refinement=level, minimum_nodal_order=float(min(observed)),landed_chart=root.x.tolist(),
              chi=float(interpolate(root.x,True)[1]),core_hoop=lever/core,collar_hoop=lever/collar,
              hoop_error_core=(lever-exact_lever)/core, hoop_error_collar=(lever-exact_lever)/collar,
              physical_position_error=float(np.linalg.norm(P-phys)), common_label_residual=float(max(abs(root.fun))),
              reconstructed_radial_label_error_at_old_bond=float(interpolate(exact)[0]-body.cut)))
    return dict(reference=ref.tolist(),exact_placement=exact.tolist(),exact_hoops=[exact_lever/core,exact_lever/collar],levels=out)


def main():
    p=ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();b,s,axes=get_cp10(a.root.resolve()); out=a.out.resolve()
    data=out/'labels.npz'
    if (out/'label_transport.npz').exists():raise FileExistsError('Preserve completed readback')
    with np.load(data,allow_pickle=False) as d:
        labels=d['labels'];order=d['order'];stencil=d['active_stencil']
    rate,counts,meta=source_t(b,s['solid'],labels,stencil,axes)
    np.savez_compressed(out/'label_transport.npz',dxi_dt=rate,counts=counts,sides=np.array([-1,1]),
                        source_sha256=CHECKPOINT_SHA,label_sha256=digest(data))
    r=dict(bond=bond(b,s['solid'],axes,labels,order), transport=meta,
           transport_sha256=digest(out/'label_transport.npz'),label_sha256=digest(data))
    (out/'readback.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
    print(json.dumps(r,allow_nan=False),flush=True)


if __name__=='__main__':main()
