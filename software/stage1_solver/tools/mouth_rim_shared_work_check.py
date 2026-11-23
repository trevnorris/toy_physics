"""Compare signed ring work with the independent change of the deposited source."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_bonded_front_adjoint import angular_weights
from mouth_bonded_material_motion import Bonded
from mouth_last_ring_solid import interpolation
from mouth_rim_shared_crystal import TerminalQuadrature, deposit
from mouth_supercore_joined_front import check_guard,file_sha


def native_load(state,quad):
    b=Bonded(state,11,32,shape=quad.sh,full_depth=True)
    D=quad.d.reshape(-1,2)[b.id]
    pos=np.einsum('qi,qij->qj',b.basis,D)
    B=np.eye(2)+np.stack([np.einsum('qi,qij->qj',g,D) for g in (b.dr,b.dw)],axis=2)
    r,z=b.r+pos[:,0],b.w+pos[:,1]
    j=b.j*np.linalg.det(B)*(r/b.r)**2
    K=B@b.K@B.transpose(0,2,1)
    hoop=b.h*r/b.r
    W=.5*(K[:,0,0]+K[:,1,1]+2*hoop**2-4)-np.log(j)+.25*(j-1)**2
    plane=np.zeros((len(quad.ax),len(quad.z)))
    deposit(plane,quad.ax,quad.z,r/quad.sh.radius(z)[0],z,b.measure*W)
    return plane


def test(state,checkpoint,load,field):
    for p in (load,field):check_guard(p)
    with np.load(load) as l, np.load(field) as f:
        if (not np.array_equal(l['ax'],f['inside_x']) or not np.array_equal(l['z'],f['duct_z'])
                or str(f['checkpoint_sha'])!=str(l['checkpoint_sha']) or str(l['checkpoint_sha'])!=file_sha(checkpoint)):
            raise ValueError('Unmatched grids or material state')
        ax,z=l['ax'],l['z']
        mean=np.einsum('raj,a->jr',f['duct_dc'],angular_weights(f['angles']))
        plane=l['plane']
    quad=TerminalQuadrature(state,checkpoint,ax,z,64)
    native=native_load(state,quad)
    free=z>quad.active_start+1.e-12
    diff=(plane-native)[:-1,free]
    tests=[mean,np.sin(7*ax[None,:])*np.sin(15*z[:,None]/quad.L)]
    tests[-1][z<quad.active_start+1e-12,:]=0
    tests[-1][:,-1]=0
    h=1e-6;p=np.zeros((2,3));mixed=np.zeros((2,len(tests)));deposits=np.zeros((2,2,len(ax),len(z)))
    for v in quad.blocks():
        ids=v['ids']==31*len(quad.xt)+32
        sel=ids.any(axis=1)
        if not sel.any():continue
        take=lambda key:v[key][sel]
        ids=ids[sel].astype(float)
        u=np.sum(take('bary')*ids,axis=1)
        G=np.stack([np.sum(take(k)*ids,axis=1) for k in ('gx','gt')],axis=1)
        G=np.einsum('qi,qij->qj',G,np.linalg.inv(take('B')))
        K,J,V,hoop,W=[take(k) for k in ('K','J','measure','hoop','W')]
        r,zs,x=[take(k) for k in ('r','z','x')]
        b,slope=quad.sh.radius(zs)
        for direction in (0,1):
            H=u/r if direction==0 else 0.
            P=np.sum(G*K[:,direction,:],axis=1)+G[:,direction]*(-1+.5*J*(J-1))+2*H*(hoop**2-1+.5*J*(J-1))
            for k,psi in enumerate([*tests,quad.inner[2]]):
                az,ar=((z,ax) if k<2 else (quad.inner[0]*quad.L,quad.inner[1]))
                chi,cx,cz=interpolation(az,ar,psi,zs,x)
                dchi=(cx/b,cz-x*slope*cx/b)[direction]*u
                derivative=V@(W*dchi+chi*P)
                if k<len(tests):mixed[direction,k]+=derivative
                if k!=1:p[direction,0]+=derivative
            for n,sgn in enumerate((-h,h)):
                rr=r+sgn*u*(direction==0);zz=zs+sgn*u*(direction==1)
                xx=rr/quad.sh.radius(zz)[0]
                B=np.broadcast_to(np.eye(2),K.shape).copy()
                B[:,direction,:]+=sgn*G
                j=J*(1+sgn*G[:,direction])*(rr/r)**2
                K2=B@K@B.transpose(0,2,1)
                W2=.5*(K2[:,0,0]+K2[:,1,1]+2*(hoop*rr/r)**2-4)-np.log(j)+.25*(j-1)**2
                c=interpolation(z,ax,mean,zz,xx)[0]+interpolation(quad.inner[0]*quad.L,quad.inner[1],quad.inner[2],zz,xx)[0]
                c0=interpolation(z,ax,mean,zs,x)[0]+interpolation(quad.inner[0]*quad.L,quad.inner[1],quad.inner[2],zs,x)[0]
                p[direction,1+n]+=V@(c*W2-c0*W)
                deposit(deposits[direction,n],ax,z,xx,zz,V*W2)
    dot=np.einsum('drz,ztr->dt',(deposits[:,1]-deposits[:,0])/(2*h),np.stack(tests,axis=1))
    e=float(np.max(abs(dot-mixed)));fd=(p[:,2]-p[:,1])/(2*h)
    if e>1.e-4 or np.max(abs(fd-p[:,0]))>1.e-4:
        raise ArithmeticError('Mixed variation failed')
    return dict(checkpoint=file_sha(checkpoint),field=file_sha(field),source=file_sha(load),
                deposited_source_native_free_norm=float(abs(native[:-1,free]).sum()),
                corrected_minus_native_free_L1=float(abs(diff).sum()),corrected_minus_native_free_max=float(abs(diff).max()),
                signed_step=h,reassembled_load_mixed_work=dot.tolist(),
                direct_material_mixed_work=mixed.tolist(),max_mixed_work_error=e,
                signed_total_radial_axial_force=(-fd).tolist(),analytic_total_radial_axial_force=(-p[:,0]).tolist())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','checkpoint','load','field'):p.add_argument(name,type=Path)
    a=p.parse_args()
    print(json.dumps(test(a.state,a.checkpoint,a.load,a.field),allow_nan=False))
