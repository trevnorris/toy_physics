"""Whole-grid signed orientation of the checkpoint-10 spatial Q1 material chart.

The chart determinant is tri-quadratic. Its tensor Bernstein coefficients
bound the polynomial throughout each cell, supplemented by one 4D physical
constitutive readout at the least positive sampled Gauss determinant.
"""
import argparse,json
from itertools import product
from pathlib import Path

import numpy as np

from mouth_cp10_label_grid import BITS,get_cp10
from mouth_eulerian_r4_spatial_cell import current
from mouth_exact_checkpoint_control import digest
from mouth_r4_refinement_review import _collar_reference
from mouth_supercore_shape import preparation
import mpmath as mp


def basis_and_grad(labels,axes,g):
    vals=np.zeros(tuple(len(a)-1 for a in axes)+(3,))
    jac=np.zeros(vals.shape+(3,));q=np.array(g)
    dx=[np.diff(a).reshape([len(a)-1 if d==j else 1 for d in range(3)]) for j,a in enumerate(axes)]
    for bits in BITS:
        block=labels[tuple(slice(n,n+len(ax)-1) for n,ax in zip(bits,axes))]
        weight=np.where(bits,q,1-q)
        vals+=np.prod(weight)*block
        for j in range(3):
            c=(2*bits[j]-1)*weight[(j+1)%3]*weight[(j+2)%3]/dx[j]
            jac[...,j]+=block*c[...,None]
    return vals,jac


def native_state(b,xi,pos,G):
    x,t,a=pos
    _,B=current(b,pos)
    # Orthogonal (radial, depth, polar arc) basis for both saved references.
    xx,tt,aa=xi
    if xx<b.cut:
        ref=preparation(4)[0]
        rr,rz=ref.radius(tt*ref.length)
        R=xx*rr*b.data.core_ratio
        C=np.array([[rr*b.data.core_ratio,xx*rz*ref.length*b.data.core_ratio,0],[0,ref.length,0],[0,0,R]])
    else:
        with mp.workdps(48):
            C=np.array(_collar_reference(b,mp.mpf(float(xx)),mp.mpf(float(tt))).tolist(),float)
        R=C[2,2]
    phys,_=current(b,pos)
    hoop=R*np.sin(aa)/phys[0]
    inv=C@G@np.linalg.inv(B)
    # |det(C)| and |det(B)| account for reference frame handedness.
    jinv=np.sign(np.linalg.det(G))*abs(np.linalg.det(inv)*hoop)
    sv=np.r_[np.linalg.svd(inv,compute_uv=False),abs(hoop)]
    J=1/jinv
    W=(sum(1/sv**2)-4)/2-np.log(J)+(J-1)**2/4 if J>0 else None
    return dict(J_inverse=jinv,J=J,grad_xi_spectral_norm=float(np.linalg.norm(G@np.linalg.inv(B),2)),W=W)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('out',type=Path)
    args=p.parse_args();out=args.out; b,s,axes=get_cp10(args.root);path=out/'labels.npz';f=np.load(path,allow_pickle=False)
    labels,order,melt=f['labels'],f['order'],f['melted_cells'];shape=melt.shape
    all_bern=np.empty(shape+(3,3,3)); corners=[]
    for idx in product(range(3),repeat=3):
        g=np.array(idx)/2; val,D=basis_and_grad(labels,axes,g)
        all_bern[(...,*idx)]=np.linalg.det(D)
    samples=all_bern.copy()
    for j in (-1,-2,-3):
        c=all_bern.swapaxes(-1,j)
        c[...,1]=2*c[...,1]-.5*(c[...,0]+c[...,2])
    lower=np.min(all_bern,axis=(-3,-2,-1)); minimum=samples.min(axis=(-3,-2,-1))
    info=dict(label_sha=digest(path),bernstein_chart_J_cells=dict(
        active_nonpositive_bernstein=int(np.sum((lower<=0)&~melt)),melt_nonpositive_bernstein=int(np.sum((lower<=0)&melt)),
        active_min_coeff=lower[~melt].min(),melt_min_coeff=lower[melt].min(),
        active_min_node_midpoint=minimum[~melt].min(),melt_min_node_midpoint=minimum[melt].min()))
    readings={}
    ref=preparation(4)[0]
    for g in product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)),repeat=3):
        val,D=basis_and_grad(labels,axes,g)
        psi,_=basis_and_grad(order[...,None]*np.ones(3),axes,g)
        pos=np.stack(np.broadcast_arrays(*np.ix_(*[a[:-1]+np.diff(a)*xi for xi,a in zip(g,axes)])),axis=-1)
        x,t,a=pos.transpose(3,0,1,2)
        rm,rz=b.sh.radius(t*b.sh.length)
        band=16*((t-.06)/.56)**2*(1-(t-.06)/.56)**2*b.epsilon
        band=np.where((t>.06)&(t<.62),band,0.)
        P=(3*np.cos(a)**2-1)/2
        R=x*rm+x**3*band*P; Rx=rm+3*x*x*band*P
        xx,tt,aa=val.transpose(3,0,1,2)
        # Natural reference hoop and meridional area from the two embeddings.
        rr=ref.radius(tt*ref.length)[0]; S=xx*rr*b.data.core_ratio; area=rr*b.data.core_ratio*ref.length
        collar=xx>=b.cut
        e=-np.log(xx[collar])/b.data.pressure;tr=b.buffers['reg_t'];er=b.buffers['reg_e']
        jj=np.clip(np.searchsorted(tr,tt[collar],side='right')-1,0,len(tr)-2);ii=np.clip(np.searchsorted(er,e,side='right')-1,0,len(er)-2)
        dt=np.diff(tr)[jj];de=np.diff(er)[ii];u=(e-er[ii])/de;v=(tt[collar]-tr[jj])/dt
        hi=(u+v>1).astype(int);sig=1-2*hi;u=np.where(hi,1-u,u);v=np.where(hi,1-v,v)
        lab=b.buffers['labels'].reshape(len(tr),len(er),2)
        c=lab[jj+hi,ii+hi];c1=lab[jj+hi,ii+1-hi]-c;c2=lab[jj+1-hi,ii+hi]-c
        S[collar]=4+c[:,0]+u*c1[:,0]+v*c2[:,0]
        area[collar]=abs((c1[:,0]*c2[:,1]-c1[:,1]*c2[:,0])/(dt*de*b.data.pressure*xx[collar]))
        ji=area*S**2*np.sin(aa)*np.linalg.det(D)/(R**2*Rx*b.sh.length*np.sin(a))
        for name,mask in [('active',~melt),('melt',melt)]:
            idx=np.unravel_index(np.argmin(np.where(mask,ji,np.inf)),shape)
            if name not in readings or ji[idx]<readings[name]['J_inverse']:
                row=native_state(b,val[idx],pos[idx],D[idx])
                if not np.isclose(row['J_inverse'],ji[idx],rtol=2e-13,atol=0):
                    raise ValueError('Inconsistent independent physical determinant')
                readings[name]=dict(**row,chi=float(psi[idx][0]),cell=list(map(int,idx)),gaussian_location=g,physical_chart=pos[idx].tolist())
    info['min_J_inverse_at_gauss']=readings
    (out/'orientation.json').write_text(json.dumps(info,indent=2,allow_nan=False)+'\n')
    print(json.dumps(info,allow_nan=False))


if __name__=='__main__':main()
