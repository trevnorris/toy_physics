"""Independent tangent and signed work of a fixed material hat at the R4 exit.

The original six P1 triangles are integrated locally. Composite quadrature
preserves the original material displacement, collar registry and physical
hat exactly. The analytic Q1-order tangent is element-interior; finite signed
steps independently test occupancy crossing its piecewise-linear chart.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from mouth_bonded_front_adjoint import angular_weights
from mouth_bonded_material_motion import Bonded
from mouth_sealed_shape_relax import Meridian
from mouth_internal_collar_check import read
from mouth_supercore_joined_front import emit, file_sha, check_guard


def q1(grid_x, grid_z, arr, x, z):
    i = np.searchsorted(grid_x, x, side='right')-1
    j = np.searchsorted(grid_z, z, side='right')-1
    i = np.clip(i,0,len(grid_x)-2)
    j = np.clip(j,0,len(grid_z)-2)
    dx, dz = grid_x[i+1]-grid_x[i],grid_z[j+1]-grid_z[j]
    a,b=(x-grid_x[i])/dx,(z-grid_z[j])/dz
    v00,v01,v10,v11=arr[j,i],arr[j,i+1],arr[j+1,i],arr[j+1,i+1]
    mix=v11-v10-v01+v00
    return (v00+a*(v01-v00)+b*(v10-v00)+a*b*mix,
            (v01-v00+b*mix)/dx,(v10-v00+a*mix)/dz,mix/(dx*dz))


def material_points(body,disp,idx,scale,order):
    """Assemble the six original test triangles and exact P1 interpolants."""
    select=np.arange((len(body.t)-1)*(len(body.x)-1)).reshape(len(body.t)-1,-1)
    iz,ix=idx
    nodes=np.arange(len(body.t)*len(body.x)).reshape(len(body.t),-1)
    coords=np.stack(np.meshgrid(body.x,body.t),axis=-1)
    peak=nodes[iz,ix]
    q,g=leggauss(order)
    a,b=np.meshgrid((q+1)/2,(q+1)/2,indexing='ij')
    uv=np.stack((a.ravel(),(b*(1-a)).ravel()))
    weights=np.outer(g,g).ravel()*(1-a).ravel()/4
    # Split each original triangle on the integer barycentric lattice.
    bary=[]; ws=[]
    for i in range(scale):
        for j in range(scale-i):
            bary.append((np.array([[i],[j]])+uv)/scale)
            ws.append(weights/scale**2)
            if i+j<scale-1:
                bary.append((np.array([[i+1],[j+1]])-uv)/scale)
                ws.append(weights/scale**2)
    uv=np.concatenate(bary,axis=1)
    weights=np.concatenate(ws)
    lam=np.column_stack((1-uv[0]-uv[1],*uv))
    for j,i in np.ndindex(select.shape):
        corners=nodes[j:j+2,i:i+2].ravel()
        if peak not in corners:
            continue
        for tri in (corners[[0,1,2]],corners[[3,2,1]]):
            if peak not in tri:continue
            tri_coords=coords.reshape(-1,2)[tri]
            v=disp.reshape(-1,2)[tri]
            h=(tri==peak).astype(float)
            dX=(tri_coords[1:]-tri_coords[0]).T
            data=np.column_stack((v,h))
            grad=np.linalg.solve(dX.T,data[1:]-data[0]).T
            yield lam@tri_coords, lam@data, grad, abs(np.linalg.det(dX))*weights


def point_terms(sh, chart, labels, sites, grid, live, step, direction):
    coord,val,grad,wt=sites
    x,t=coord.T
    br, brp, _=sh.profile(t)
    P=np.exp(chart.co[0])
    r0,z0=br*x,sh.length*t
    raw_eta=-np.log(x)/P
    u=(raw_eta+.5)*2*chart.na
    v=t*chart.nt
    i=np.minimum(u.astype(int),2*chart.na-1)
    j=np.minimum(v.astype(int),chart.nt-1)
    a,b=u-i,v-j
    upper=a+b>1
    ind=chart.triangles[2*(j*(2*chart.na)+i)+upper]
    vertex=labels.reshape(-1,2)[ind]
    a=np.where(upper,1-a,a); b=np.where(upper,1-b,b)
    start=sh.neck+vertex[:,0,0]+a*(vertex[:,1,0]-vertex[:,0,0])+b*(vertex[:,2,0]-vertex[:,0,0])
    inv=np.stack(((vertex[:,1]-vertex[:,0])*2*chart.na,(vertex[:,2]-vertex[:,0])*chart.nt),axis=-1)
    inv*=np.where(upper,-1,1)[:,None,None]
    A=np.zeros_like(inv)
    A[:,0,0]=-P*r0;A[:,0,1]=x*sh.length*brp;A[:,1,1]=sh.length
    F=A@np.linalg.inv(inv)
    hoop0=r0/start
    J0=abs(np.linalg.det(F))*hoop0**2
    volume=4*np.pi*r0**2*br*sh.length*wt/J0
    position=(np.column_stack((r0,z0))+val[:,:2])
    B=np.broadcast_to(np.eye(2),(len(x),2,2)).copy()
    b0=grad[:2,0]/br[:,None]
    B[:,:,0]+=b0
    B[:,:,1]+=grad[:2,1]/sh.length - b0*(x*brp)[:,None]
    dv=val[:,2,None]*direction
    dr,dz=dv.T
    Dh=grad[2,0]/br
    DW=grad[2,1]/sh.length-Dh*x*brp
    U=direction[None,:,None]*np.stack((Dh,DW),axis=-1)[:,None,:]
    T=B@F; V=U@F
    r,z=position.T
    b,bp,bpp=sh.profile(z/sh.length)
    xx=r/b
    dx=(dr-xx*bp*dz)/b
    ddx=-(2*bp*dz*dx+xx*bpp*dz**2)/b
    C,Cx,Cz,Cxz=(u+v for u,v in zip(q1(*grid,xx,z),q1(*live,xx,z)))
    Cp=Cx*dx+Cz*dz
    Cpp=Cx*ddx+2*Cxz*dx*dz
    hoop=hoop0*r/r0
    J=J0*np.linalg.det(B)*(r/r0)**2
    invU=np.linalg.solve(B,U)
    dh=dr/r
    logp=np.trace(invU,axis1=1,axis2=2)+2*dh
    logpp=-np.trace(invU@invU,axis1=1,axis2=2)-2*dh*dh
    def energy(T,hoop,J):
        return .5*((T*T).sum(axis=(1,2))+2*hoop**2-4)-np.log(J)+.25*(J-1)**2
    W=energy(T,hoop,J)
    pcoef=-1+.5*J*(J-1)
    Wp=(T*V).sum(axis=(1,2))+2*hoop**2*dh+pcoef*logp
    Wpp=(V*V).sum(axis=(1,2))+2*(hoop*dh)**2+pcoef*logpp+.5*J*(2*J-1)*logp**2
    deriv=volume@(Cp*W+C*Wp)
    curv=volume@(Cpp*W+2*Cp*Wp+C*Wpp)
    en=[]
    stretches=[]
    for s in [-step,0.,step]:
        dest_r=r+s*dr;dest_z=z+s*dz
        bx=sh.radius(dest_z)[0];ex=dest_r/bx
        ct=q1(*grid,ex,dest_z)[0]+q1(*live,ex,dest_z)[0]
        Bt=B+s*U
        jac=J0*np.linalg.det(Bt)*(dest_r/r0)**2
        H=hoop0*dest_r/r0;Dt=T+s*V
        if min(jac)<=0: raise ValueError('Folded material site')
        value=energy(Dt,H,jac)
        en.append(volume@((ct-C)*value+C*(value-W)))
        sing=np.linalg.svd(Dt,compute_uv=False)
        stretches.append(float(min(H.min(),sing[:,-1].min())))
    return deriv,curv,np.asarray(en),np.asarray(stretches)


def main(a):
    check_guard(a.state)
    outcome=json.loads((a.failed.parent/'outcome.json').read_text())
    limit=json.loads((a.failed.parent/'limit-validation.json').read_text())
    if not outcome['limitsVerified'] or outcome['childOutcome']['exitCode']!=1 or not limit['verified'] or int(limit['actual']['memory.swap.max']):
        raise ValueError('Wrong guarded failed shape')
    with np.load(a.failed) as f:
        sh=Meridian(4.,f['shape_params'],extra=float(f['shape_extra']))
        if file_sha(a.state)!=f['state_sha'].item() or file_sha(a.failed)!='d564e2ab9b6d4c25cf3050bac01ab42ae954f4019b059e1e4eac8e715789e1b8':
            raise ValueError('Incorrect rejected material generation')
        d=f['displacement']; field=f['field']
    b=Bonded(a.state,11,32,shape=sh,full_depth=True)
    chart,labels,*_=read(a.state)
    if a.fields:
        check_guard(a.fields)
        with np.load(a.fields) as f:
            if f['checkpoint_sha'].item()!=file_sha(a.failed):raise ValueError('Mismatched receiving field')
            dims=tuple(int(i) for i in f['grid'])
            rx,zz,angles,dc=(f[key] for key in ('inside_x','duct_z','angles','duct_dc'))
    else:
        from mouth_polar_joined import PolarCoupled
        dims=(13,10,14)
        m=PolarCoupled(sh,b.grid,1,dims,precondition=False)
        rx,zz,angles,dc=m.rr[0],m.z[1],m.a,field[2,m.ids[1,0]]
    correction=np.einsum('raj,a->jr',dc,angular_weights(angles))
    grid=(np.linspace(0,1,b.grid.shape[1]),sh.length*np.linspace(0,1,b.grid.shape[0]),b.grid)
    live=(rx,zz,correction)
    emit(kind='source',shape=file_sha(a.failed),solid=file_sha(a.state),mode_site=[31,32],field_sha=file_sha(a.fields) if a.fields else None,grid=dims)
    for refinement in a.subdivide:
        for step in a.steps:
            work,k,hleft,hright,kfd,stretch=[],[],[],[],[],[]
            for direction in ((1.,0.),(0.,1.),(1.,1.)):
                A,B,en,mins=0.,0.,np.zeros(3),np.full(3,np.inf)
                for sites in material_points(b,d,(31,32),refinement,4):
                    result=point_terms(sh,chart,labels,sites,grid,live,step,np.asarray(direction))
                    A+=result[0];B+=result[1];en+=result[2];mins=np.minimum(mins,result[3])
                work.append(float((en[2]-en[0])/(2*step)));k.append(float(B))
                hleft.append(float(-en[0]/step));hright.append(float(en[2]/step))
                kfd.append(float((en[2]-2*en[1]+en[0])/step**2))
                stretch.append(mins.tolist())
                if not np.isfinite(A+B) or abs(work[-1]-A)>max(.004,step**2*2e6):
                    raise ArithmeticError(f'Signed work and element tangent differ: {A}, {work[-1]}')
            def tensor(h):
                h12=(h[2]-h[0]-h[1])/2
                return [[h[0],h12],[h12,h[1]]]
            emit(kind='independent_patch',subdivide=refinement,step=step,
                 force=(-np.array(work[:2])).tolist(),signed_negative_slope=hleft,signed_positive_slope=hright,
                 analytic_in_cell_H=tensor(k),finite_step_H=tensor(kfd),
                 signed_eigenvalues=np.linalg.eigvalsh(tensor(kfd)).tolist(),min_stretch_minus_zero_plus=stretch)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','failed'):
        p.add_argument(name,type=Path)
    p.add_argument('--fields',type=Path)
    p.add_argument('--subdivide',nargs='+',type=int,default=[1,2,4,8,16])
    p.add_argument('--steps',nargs='+',type=float,default=[1e-5,2.5e-6])
    main(p.parse_args())
