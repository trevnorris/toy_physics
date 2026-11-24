"""A continuous inverse material chart on the checkpoint-10 physical Q1 grid.

All active cells use inverse labels on their closure. Only strictly melted
cells enter the scalar gauge. The three nondimensional gauge components are
xi minus the identity physical chart (radius fraction, depth fraction, angle).
"""
from pathlib import Path
import argparse
from itertools import product
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import brentq
from scipy.interpolate import RegularGridInterpolator

from mouth_eulerian_r4_spatial_cell import current, digest, roots
from mouth_exact_checkpoint_control import CHECKPOINT_SHA,read
from mouth_free_shape_solid import Solid
from mouth_joined_angular_work import bump
from mouth_r4_conforming_prisms import TRIANGLES
from mouth_r4_fold_order import axes_model
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_shape import preparation

BITS=list(product(range(2),repeat=3))


def get_cp10(root):
    path=root/'r4_melted_survival01/coarse/state_00010.npz'
    prep=root/'supercore_joined_front01/adopt309/state.npz'
    if digest(path)!=CHECKPOINT_SHA:
        raise ValueError('Checkpoint 10 hash mismatch')
    s=np.load(path,allow_pickle=False)
    if str(s['state_sha'])!=digest(prep) or int(s['radius'])!=4 or int(s['winding'])!=1:
        raise ValueError('Requires the held radius-four wound CP10 state')
    shape=Meridian(4.,s['shape_params'],float(s['shape_extra']))
    depth=read(root/'r4_whole_mouth01/coarse/endpoint.json')['depth']
    b=Solid(prep,shape,float(s['x'][5]),8,24,6,axes_model(shape,depth),degree=3,subdivisions=2)
    axes=[b.buffers['field_'+c] for c in ('x','t','a')]
    return b,s,axes


def prisms(body,solid):
    T,X,A=np.meshgrid(body.t,body.x,body.a,indexing='ij')
    P=np.stack((X,T,A),axis=-1); Y=P+solid
    ref=[]; spatial=[]; keys=[]
    for j,i,k,tri in product(range(len(body.t)-1),range(len(body.x)-1),range(len(body.a)-1),range(2)):
        keys.append((j,i,k,tri))
        ix=[[j+dj for l in range(2) for dj,di in TRIANGLES[tri]],
            [i+di for l in range(2) for dj,di in TRIANGLES[tri]], [k+l for l in range(2) for _ in range(3)]]
        ref.append(P[tuple(ix)].reshape(2,3,3));spatial.append(Y[tuple(ix)].reshape(2,3,3))
    landed=np.asarray(spatial)
    return keys,np.asarray(ref),landed,landed.min(axis=(1,2)),landed.max(axis=(1,2))


def allnode_inverse(axes,P,mesh,need):
    """Solve all candidate prism inverses in batches, validating landed residuals.

    The azimuth fraction is a root of a cubic. Degree reduction is relative
    to its coefficients on the unit interval. Every retained root still
    has to land back on the actual multilinear prism.
    """
    _,ref,land,lo,hi=mesh
    shape=need.shape; pair_elem=[]; pair_node=[]
    for e,(L,H) in enumerate(zip(lo,hi)):
        slices=tuple(slice(np.searchsorted(a,l-2e-13),np.searchsorted(a,h+2e-13,side='right'))
                         for a,l,h in zip(axes,L,H))
        inds=np.nonzero(need[slices])
        if len(inds[0]):
            ids=np.ravel_multi_index(tuple(idx+sl.start for idx,sl in zip(inds,slices)),shape)
            pair_elem.append(np.full(len(ids),e));pair_node.append(ids)
    all_e=np.concatenate(pair_elem);all_n=np.concatenate(pair_node)
    labels=np.full((np.prod(shape),3),np.nan)
    max_err=0.; duplicates={}; found=np.zeros(np.prod(shape),int)
    flat=P.reshape(-1,3)
    for offset in range(0,len(all_e),12000):
        es,ns=all_e[offset:offset+12000],all_n[offset:offset+12000]
        low=land[es,0];diff=land[es,1]-low
        l1,l2=low[:,1]-low[:,0],low[:,2]-low[:,0]
        d1,d2=diff[:,1]-diff[:,0],diff[:,2]-diff[:,0]
        t,dt=flat[ns]-low[:,0],-diff[:,0]
        n0=np.cross(l1,l2);n1=np.cross(l1,d2)+np.cross(d1,l2);n2=np.cross(d1,d2)
        def dot(a,b):return np.einsum('ij,ij->i',a,b)
        c=np.column_stack((dot(t,n0),dot(dt,n0)+dot(t,n1),dot(dt,n1)+dot(t,n2),dot(dt,n2)))
        scale=np.max(np.abs(c),axis=1);rank=np.sum(np.maximum.accumulate((abs(c)>scale[:,None]*2e-14)[:,::-1],axis=1),axis=1)-1
        for n in (1,2,3):
            idx=np.flatnonzero(rank==n)
            if not len(idx):continue
            cc=c[idx,:n]/c[idx,n,None]
            if n==1: roots_y=-cc
            else:
                comp=np.zeros((len(idx),n,n));comp[:,np.arange(1,n),np.arange(n-1)]=1
                comp[:,:,-1]=-cc
                roots_y=np.linalg.eigvals(comp)
            for ys in roots_y.T:
                keep=(abs(ys.imag)<1e-9)&(ys.real>=-2e-10)&(ys.real<=1+2e-10)
                jj=idx[keep]; y=ys[keep].real
                # Newton refinement of the scalar polynomial corrects mild
                # conditioning of the companion matrix before barycentric work.
                for _ in range(2):
                    a,b,cr,d=c[jj].T;dy=(3*d*y+2*cr)*y+b
                    yd=np.divide(((d*y+cr)*y+b)*y+a,dy,out=np.zeros_like(y),where=abs(dy)>1e-30)
                    y=y-yd
                v=l1[jj]+y[:,None]*d1[jj];w=l2[jj]+y[:,None]*d2[jj]
                dst=t[jj]+y[:,None]*dt[jj]; norm=np.cross(v,w)
                den=dot(norm,norm)
                good=den>0
                bc1=np.divide(dot(np.cross(dst,w),norm),den,out=np.full_like(y,np.nan),where=good)
                bc2=np.divide(dot(np.cross(v,dst),norm),den,out=np.full_like(y,np.nan),where=good)
                er=np.max(abs(bc1[:,None]*v+bc2[:,None]*w-dst),axis=1)
                keep=(er<=2e-13)&(np.minimum.reduce([bc1,bc2,1-bc1-bc2,y,1-y])>=-3e-10)
                if not np.any(keep):continue
                jjj=jj[keep]; yy=y[keep]; ee=es[jjj];nn=ns[jjj]
                bcs=np.column_stack((1-bc1[keep]-bc2[keep],bc1[keep],bc2[keep]))
                pts=ref[ee,0]+yy[:,None,None]*(ref[ee,1]-ref[ee,0])
                hits=np.einsum('ij,ijk->ik',bcs,pts)
                max_err=max(max_err,float(np.max(er[keep])))
                for node,h in zip(nn,hits):
                    if found[node]==0:
                        labels[node]=h;found[node]=1
                    elif np.max(abs(labels[node]-h))>1e-9:
                        options=duplicates.setdefault(node,[labels[node]])
                        if all(np.max(abs(val-h))>1e-9 for val in options):
                            options.append(h);found[node]+=1
    # Any degenerate projection requires a separate inverse query. Never
    # insert an extrapolated label into an ordered cell.
    missing=np.flatnonzero(need.ravel()&(found==0))
    for node in missing:
        r=roots(mesh,flat[node])
        if len(r)==1:labels[node]=r[0];found[node]=1
    return labels.reshape((*shape,3)), found.reshape(shape),duplicates, max_err


def matrix_for(body,axes,melt):
    """Physical four-volume Dirichlet integral of xi minus the identity chart."""
    X,T,A=axes; dims=tuple(map(len,axes))
    ijk=np.array(np.nonzero(melt)).T
    if not len(ijk):return sparse.csr_matrix((np.prod(dims),)*2)
    ids=np.array([np.ravel_multi_index(tuple((ijk+bits).T),dims) for bits in BITS]).T
    width=np.column_stack([np.diff(v)[ijk[:,i]] for i,v in enumerate(axes)])
    bottom=np.column_stack([v[ijk[:,i]] for i,v in enumerate(axes)])
    blocks=np.zeros((len(ijk),8,8))
    for g in product((.5-.5/np.sqrt(3.),.5+.5/np.sqrt(3.)), repeat=3):
        x,t,a=(bottom+width*g).T
        z=t*body.sh.length
        h,hz=bump(z,body.sh.length)
        r,rz=body.sh.radius(z)
        p=(3*np.cos(a)**2-1)/2;dp=-3*np.sin(a)*np.cos(a)
        R=x*r+body.epsilon*h*x**3*p
        Rx=r+3*body.epsilon*h*x*x*p
        Rz=x*rz+body.epsilon*hz*x**3*p
        Ra=body.epsilon*h*x**3*dp
        dv=4*np.pi*R**2*np.sin(a)*Rx*body.sh.length*np.prod(width,axis=1)/8
        D=[]
        for bits in BITS:
            val=np.where(bits,g,1-np.asarray(g))
            d=(2*np.asarray(bits)-1)*np.prod(val)/val/width
            u=d[:,0]/Rx
            D.append(np.column_stack((u,d[:,1]/body.sh.length-u*Rz,(d[:,2]-u*Ra)/R)))
        D=np.stack(D,axis=1)
        blocks+=np.einsum('eic,ejc,e->eij',D,D,dv,optimize=True)
    rows=np.broadcast_to(ids[:,:,None],blocks.shape)
    cols=np.broadcast_to(ids[:,None,:],blocks.shape)
    return sparse.coo_matrix((blocks.ravel(),(rows.ravel(),cols.ravel())),shape=(np.prod(dims),)*2).tocsr()


def harmonic(K,p,labels,stencil):
    # Use the registered chart on the clamped wall, entrance and terminal disk.
    # Normal flux elsewhere follows from the energy. On the axis the depth
    # correction is shared across angles; polar-limit coefficients are evenly
    # continued from the neighboring radial row in the strictly melted region.
    dims=stencil.shape
    ids=np.arange(stencil.size).reshape(dims)
    constraints=np.zeros_like(labels,dtype=bool)
    fixed_chart=np.zeros_like(constraints)
    fixed_chart[:,[0,-1]]=True;fixed_chart[-1]=True
    fixed_chart[0,:,:,0]=True
    fixed_chart[:,:,[0,-1],2]=True
    good=stencil[...,None]&fixed_chart
    err=float(np.max(np.abs((labels-p)[good])))
    if err>2e-9:raise ArithmeticError(f'Clamped inverse error {err}')
    initial=np.where(stencil[...,None],labels-p,0)
    initial[fixed_chart]=0.
    residual=[]; mappings=[]
    for c in range(3):
        groups=ids.copy()
        if c==1:groups[0]=ids[0,:,:1]
        if c==2:
            # Exact ordered limits are set by inversion, not reimposed from
            # a finite spatial-cell spacing. Only dead axis nodes use the
            # regular continuation as a constraint.
            for j in range(dims[1]):
                for k in range(dims[2]):
                    if not stencil[0,j,k] and not fixed_chart[0,j,k,c]:groups[0,j,k]=ids[1,j,k]
        _,ix=np.unique(groups,return_inverse=True);ix=ix.ravel();nm=ix.max()+1
        C=sparse.coo_matrix((np.ones(stencil.size),(np.arange(stencil.size),ix)),shape=(stencil.size,nm)).tocsr()
        km=(C.T@K@C).tocsr()
        prescribed=stencil|fixed_chart[...,c]
        fixed=np.unique(ix[prescribed.ravel()]);vals=np.zeros(nm)
        vals[ix[prescribed.ravel()]]=initial[...,c].ravel()[prescribed.ravel()]
        if np.max(abs(vals[ix[prescribed.ravel()]]-initial[...,c].ravel()[prescribed.ravel()]))>2e-9:
            raise ArithmeticError('Conflicting prescribed regularity on the crystalline axis')
        free=np.ones(nm,bool);free[fixed]=False
        vals[free]=spsolve(km[free][:,free],-(km[free][:,fixed]@vals[fixed]))
        initial[...,c]=vals[ix].reshape(dims)
        constraints[...,c]=(~free[ix]).reshape(dims)
        residual.append(float(np.max(abs((km@vals)[free]),initial=0.)))
        mappings.append(ix.reshape(dims))
    labels=p+initial
    d=initial.reshape(-1,3); raw=K@d
    return labels,constraints,np.stack(mappings,axis=-1),dict(max_projected_gradient=residual,
        max_unscaled_active_trace_reaction=np.max(abs(raw[stencil.ravel()]),axis=0).tolist(),
        unscaled_energy=float(.5*np.sum(raw*d)),max_snap_in_clamped_chart=err)


def main():
    pp=argparse.ArgumentParser(description=__doc__);pp.add_argument('root',type=Path);pp.add_argument('out',type=Path)
    args=pp.parse_args();root=args.root.resolve();out=args.out.resolve();out.mkdir(exist_ok=True,parents=True)
    if (out/'labels.npz').exists():raise FileExistsError('Preserve prior labels')
    b,s,axes=get_cp10(root)
    P=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1)
    order=s['field'][2,b.ids].transpose(0,2,1)
    dims=order.shape
    neighbors=np.lib.stride_tricks.sliding_window_view(order,(2,2,2))
    active=np.any(neighbors>0,axis=(-3,-2,-1));melt=~active
    stencil=np.zeros(dims,bool)
    for bits in BITS:
        ss=tuple(slice(bit, d-1+bit) for bit,d in zip(bits,dims));stencil[ss]|=active
    mesh=prisms(b,s['solid'])
    lab,counts,dup,r0=allnode_inverse(axes,P,mesh,stencil)
    summary=dict(source=CHECKPOINT_SHA,axes_shape=dims,n_active_cells=int(sum(active.ravel())),n_melt_cells=int(sum(melt.ravel())),
             required_nodal_inverse=int(stencil.sum()),strictly_ordered_nodes=int((order>0).sum()),
             missing_indices=np.argwhere(stencil&(counts==0)).tolist(),multiple_indices=np.argwhere(counts>1).tolist(),
             multiple_node_labels={str(np.unravel_index(i,dims)):[v.tolist() for v in vs] for i,vs in dup.items()},
             max_inverse_chart_roundtrip=r0)
    def emit(): (out/'reconstruction.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    emit()
    if np.any(stencil&(counts!=1)):
        raise ArithmeticError('Ambiguous inverse on active support closure; inspect roots')
    K=matrix_for(b,axes,melt)
    labels,prescribed,groups,extra=harmonic(K,P,lab,stencil)
    sparse.save_npz(out/'melt_stiffness.npz',K)
    np.savez_compressed(out/'labels.npz',x=axes[0],t=axes[1],a=axes[2],labels=labels,order=order,active_stencil=stencil,
          melted_cells=melt,prescribed=prescribed,constraint_groups=groups,source_sha256=CHECKPOINT_SHA)
    summary.update(harmonic=extra,label_sha256=digest(out/'labels.npz'),stiffness_sha256=digest(out/'melt_stiffness.npz'))
    emit()
    print(json.dumps(summary,allow_nan=False),flush=True)


if __name__=='__main__':main()
