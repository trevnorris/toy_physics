"""Fixed-footprint forces and two-direction stiffness at saved downstream R4 rings.

An original triangular nodal hat and the saved displacement are prolonged
exactly. Subdivision adds integration sites on the same material map. The
crystalline order is sampled in space after each displacement. The coherent
field can be supplied from an independently matched fluid discretization.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from mouth_bonded_front_adjoint import angular_weights
from mouth_internal_collar_check import read
from mouth_joined_front_audit import Front
from mouth_joint_polar_shape import Trial
from mouth_sealed_shape_relax import Meridian
from mouth_supercore_joined_front import check_guard, emit, file_sha


def source(state, point):
    check_guard(state)
    meta = json.loads((point.parent/'outcome.json').read_text())
    limits = json.loads((point.parent/'limit-validation.json').read_text())
    history = list(map(json.loads, (point.parent/'stdout').read_text().splitlines()))
    if (not limits['verified'] or int(limits['actual']['memory.swap.max']) or not meta['limitsVerified']
            or meta['childOutcome']['exitCode'] != 1
            or file_sha(point) != 'd564e2ab9b6d4c25cf3050bac01ab42ae954f4019b059e1e4eac8e715789e1b8'):
        raise ValueError('Expected the saved second complete sweep of the rejected R4 wall')
    with np.load(point) as row:
        if row['state_sha'].item()!=file_sha(state) or json.loads(row['readout'].item())['energy']!=history[-1]['energy']:
            raise ValueError('Unmatched material state')
        return (Meridian(4.,row['shape_params'],extra=float(row['shape_extra'])),
                row['displacement'],row['field'],row['x'])


def order_tables(state, point, source_data, archive):
    sh, d, f, x = source_data
    _,_,inner,_ = read(state)
    if archive is None:
        base=Front(state,(13,10,14),1)
        model=Trial(base,state,x,(d.shape[1]-1)//3,d.shape[0]-1,False).model
        rx,z,mu,dc=model.rr[0],model.z[1],model.a,f[2,model.ids[1,0]]
    else:
        check_guard(archive)
        with np.load(archive) as a:
            if a['checkpoint_sha'].item()!=file_sha(point):
                raise ValueError('Carrier is on a different prepared wall')
            rx,z,mu,dc=[a[key] for key in ('inside_x','duct_z','angles','duct_dc')]
    correction=np.einsum('raj,a->jr',dc,angular_weights(mu))
    return [inner,(z/sh.length,rx,correction)]


def interpolation(t_axis,x_axis,values,t,x):
    kt=np.clip(np.searchsorted(t_axis,t)-1,0,len(t_axis)-2)
    ix=np.clip(np.searchsorted(x_axis,x)-1,0,len(x_axis)-2)
    dt=t_axis[kt+1]-t_axis[kt]; dx=x_axis[ix+1]-x_axis[ix]
    a=(x-x_axis[ix])/dx; b=(t-t_axis[kt])/dt
    lo,hi=values[kt,ix],values[kt+1,ix]
    dl,dh=values[kt,ix+1]-lo,values[kt+1,ix+1]-hi
    return (1-b)*(lo+a*dl)+b*(hi+a*dh),((1-b)*dl+b*dh)/dx,(hi-lo+a*(dh-dl))/dt


def tri(xt,tt,nx,nt,x,t):
    # Both triangles use the same local orientation as the bonded material law.
    i=np.clip(np.searchsorted(xt,x)-1,0,nx-2)
    j=np.clip(np.searchsorted(tt,t)-1,0,nt-2)
    dx,dt=xt[i+1]-xt[i],tt[j+1]-tt[j]
    a,b=(x-xt[i])/dx,(t-tt[j])/dt
    upper=a+b>1.
    a,b=np.where(upper,1-a,a),np.where(upper,1-b,b)
    ids=np.stack([j*nx+i+np.where(upper,nx+1,0),j*nx+i+np.where(upper,nx,1),
                  j*nx+i+np.where(upper,1,nx)],axis=1)
    return ids,np.stack((1-a-b,a,b),axis=1),np.where(upper,-1.,1.)/dx,np.where(upper,-1.,1.)/dt


def seam_samples(xt,tt,j,i,p,chart,sh,d,order,refine):
    """Cut quadrature on both P1 charts and the current spatial order tables."""
    q,w=leggauss(5)
    q,w=(q+1)/2,w/2
    tx=np.linspace(tt[j-1],tt[j+1],2*refine+1)
    tx=np.unique(np.r_[tx,chart.t[(chart.t>tx[0])&(chart.t<tx[-1])]])
    xlo,xhi=xt[i-1],xt[i+1]
    xa=np.unique(np.concatenate([o[1] for o in order]))
    za=np.unique(np.concatenate([o[0] for o in order]))
    pieces=[]
    for tlo,thi in zip(tx[:-1],tx[1:]):
        for t,wt in zip(tlo+q*(thi-tlo),w*(thi-tlo)):
            k=np.searchsorted(tt,t)-1
            v=(t-tt[k])/(tt[k+1]-tt[k])
            grid=xt[i-1:i+2]
            rays=np.sort(np.r_[grid,grid[:-1]+np.diff(grid)*(1-v)])
            m=np.searchsorted(chart.t,t)-1
            b=(t-chart.t[m])/(chart.t[m+1]-chart.t[m])
            cracks=np.exp(-p*np.r_[chart.a,chart.a[:-1]+np.diff(chart.a)*(1-b)])
            cuts=list(rays)+list(cracks[(cracks>xlo)&(cracks<xhi)])
            # Position is affine in x on each native solid triangle.
            ids,bary,_,_=tri(xt,tt,len(xt),len(tt),rays,np.full(len(rays),t))
            off=np.einsum('qi,qik->qk',bary,d.reshape(-1,2)[ids])
            r=rays*sh.radius(t*sh.length)[0]+off[:,0]
            z=t*sh.length+off[:,1]
            xphys=r/sh.radius(z)[0]
            for n,(lo,hi) in enumerate(zip(rays[:-1],rays[1:])):
                r0,r1=r[n:n+2];z0,z1=z[n:n+2]
                a,b=xphys[n:n+2]
                for knot in xa[(xa>min(a,b))&(xa<max(a,b))]:
                    u=brentq(lambda u:(r0+u*(r1-r0))/sh.radius(z0+u*(z1-z0))[0]-knot,0.,1.)
                    cuts.append(lo+u*(hi-lo))
                for depth in za[(za*sh.length>min(z0,z1))&(za*sh.length<max(z0,z1))]:
                    cuts.append(lo+(hi-lo)*(depth*sh.length-z0)/(z1-z0))
            cuts=np.unique(cuts)
            dx=np.diff(cuts);x=cuts[:-1,None]+dx[:,None]*q
            pieces.append((x.ravel(),np.full(x.size,t),(dx[:,None]*w*wt).ravel()))
    return tuple(np.concatenate(a) for a in zip(*pieces))


class Ring:
    def __init__(self,state,sh,d,order,node,refine,native_vertex=False,quadrature=4,seams=False):
        self.sh,self.order,self.node=sh,order,node
        chart,labels,_,_=read(state)
        pressure=np.exp(chart.co[0]); cut=np.exp(-pressure/2)
        k,n=d.shape[0]-1,(d.shape[1]-1)//3
        xt=np.r_[np.linspace(0,cut,2*n+1)[:-1],np.linspace(cut,1,n+1)]
        tt=np.linspace(0,1,k+1)
        j,i=node
        base_dx=xt[i+1]-xt[i]; base_dt=tt[j+1]-tt[j]
        stretch=1/refine if native_vertex else 1.
        # Common-footprint studies keep the original hat scale. A separate
        # native-vertex readout explicitly shrinks the integration footprint.
        self.vcell=(base_dx*stretch,base_dt*stretch)
        nside=1 if native_vertex else refine
        X=np.linspace(-stretch,stretch,2*nside+1)*base_dx+xt[i]
        T=np.linspace(-stretch,stretch,2*nside+1)*base_dt+tt[j]
        q,w=leggauss(quadrature);q=(q+1)/2
        aa,bb=np.meshgrid(q,q,indexing='ij')
        wf=(np.outer(w,w)*(1-aa)/4).ravel()
        aa=aa.ravel();bb=(bb*(1-q[:,None])).ravel()
        xx,zz=np.meshgrid(X[:-1],T[:-1],indexing='ij')
        dx,dt=X[1]-X[0],T[1]-T[0]
        xx,zz=xx.ravel()[:,None],zz.ravel()[:,None]
        x=np.concatenate([xx+dx*aa,xx+dx*(1-aa)],axis=1).ravel()
        t=np.concatenate([zz+dt*bb,zz+dt*(1-bb)],axis=1).ravel()
        quadrature_weight=np.tile(wf,2*len(xx))*dx*dt
        if seams:
            if native_vertex:raise ValueError('The seam comparison uses fixed physical support')
            x,t,quadrature_weight=seam_samples(xt,tt,j,i,pressure,chart,sh,d,order,refine)
        # Unit triangular hat given by interpolation through the original node.
        u=(x-xt[i])/(base_dx*stretch); v=(t-tt[j])/(base_dt*stretch)
        distance=np.stack((u,-u,v,-v,u+v,-u-v))
        h=1-np.max(distance,axis=0);valid=h>0
        face=np.argmax(distance[:,valid],axis=0)
        self.phi=h[valid]
        px=-np.array([1,-1,0,0,1,-1])[face]/(base_dx*stretch)
        pt=-np.array([0,0,1,-1,1,-1])[face]/(base_dt*stretch)
        x,t,quadrature_weight=x[valid],t[valid],quadrature_weight[valid]
        radii,slopes=sh.radius(t*sh.length)
        self.rbase=x*radii;self.zbase=t*sh.length
        self.gr=px/radii;self.gz=pt/sh.length-x*slopes*self.gr
        ids,bary,gx,gt=tri(xt,tt,len(xt),len(tt),x,t)
        D=d.reshape(-1,2)[ids]
        offs=np.einsum('qi,qik->qk',bary,D)
        dr=(D[:,1]-D[:,0])*(gx/radii)[:,None]
        dz=(D[:,2]-D[:,0])*(gt/sh.length)[:,None]-(x*slopes)[:,None]*dr
        self.B=np.eye(2)+np.stack([dr,dz],axis=2)
        self.r,self.z=self.rbase+offs[:,0],self.zbase+offs[:,1]
        if min(self.z)/sh.length<.62:
            raise ValueError('Probe touches active internal angular mode')
        # The original collar registration is piecewise linear on its own
        # triangle grid (e,t), not the computational material grid (x,t).
        e=-np.log(x)/pressure
        ids,B,gx,gt=tri(chart.a,chart.t,len(chart.a),len(chart.t),e,t)
        ell=labels.reshape(-1,2)[ids]
        ref_radius=sh.neck+np.einsum('qi,qi->q',B,ell[:,:,0])
        inv=np.stack([(ell[:,1]-ell[:,0])*gx[:,None],
                      (ell[:,2]-ell[:,0])*gt[:,None]],axis=2)
        spatial=np.zeros_like(inv)
        spatial[:,0,0]=-pressure*self.rbase
        spatial[:,0,1]=x*sh.length*slopes
        spatial[:,1,1]=sh.length
        F=spatial@np.linalg.inv(inv)
        self.K=F@F.transpose((0,2,1))
        self.refhoop=self.rbase/ref_radius
        self.j0=abs(np.linalg.det(F))*self.refhoop**2
        self.w=4*np.pi*self.rbase**2*radii*sh.length*quadrature_weight/self.j0
        self.site=[xt[i]*sh.radius(tt[j]*sh.length)[0]+d[j,i,0],tt[j]*sh.length+d[j,i,1]]
        self.reference_volume=float(self.w.sum())
        self.baseline=None

    def sample(self,radial,axial,stats=False):
        r=self.r+self.phi*radial; z=self.z+self.phi*axial
        b,bw=self.sh.radius(z); x=r/b; t=z/self.sh.length
        B=self.B.copy()
        for component,amount in enumerate((radial,axial)):
            B[:,component,0]+=self.gr*amount
            B[:,component,1]+=self.gz*amount
        scaled=r/self.rbase
        J=self.j0*np.linalg.det(B)*scaled**2
        hoop=self.refhoop*scaled
        if min(J)<=0 or min(scaled)<=0 or max(x)>1 or max(t)>1:
            raise ValueError('Invalid signed material placement')
        K=B@self.K@B.transpose((0,2,1))
        W=(K[:,0,0]+K[:,1,1]+2*hoop**2-4)/2-np.log(J)+.25*(J-1)**2
        c,dx,dt=map(sum,zip(*(interpolation(*grid,t,x) for grid in self.order)))
        density=c*W
        if self.baseline is None:self.baseline=density.copy()
        Delta=self.w@(density-self.baseline)
        A=B@self.K+(-1+.5*J*(J-1))[:,None,None]*np.linalg.inv(B).transpose(0,2,1)
        dW=A[:,:,0]*self.gr[:,None]+A[:,:,1]*self.gz[:,None]
        dW[:,0]+=2*(hoop**2-1+.5*J*(J-1))*self.phi/r
        dchi=np.stack([dx/b,dt/self.sh.length-x*bw*dx/b],axis=1)*self.phi[:,None]
        g=self.w@(c[:,None]*dW+W[:,None]*dchi)
        if not stats:return float(Delta),g
        singular=np.sqrt(np.linalg.eigvalsh(K))
        stretches=np.column_stack([singular,hoop])
        stress=c[:,None]*(stretches**2-1+.5*J[:,None]*(J[:,None]-1))/J[:,None]
        return float(Delta),g,dict(strain_force=(-self.w@(c[:,None]*dW)).tolist(),
                                  order_position_force=(-self.w@(W[:,None]*dchi)).tolist(),
                                  min_stretch=float(min(singular[:,0].min(),hoop.min())),
                                  max_stretch=float(max(singular[:,1].max(),hoop.max())),
                                  min_J=float(J.min()),max_J=float(J.max()),min_order=float(c.min()),
                                  max_abs_principal_elastic_stress=float(np.max(abs(stress))))


def signed(ring,h):
    _,slope,base=ring.sample(0,0,True)
    er,ez=np.eye(2)*h
    vals={}
    steps=(er,-er,ez,-ez,er+ez,er-ez,-er+ez,-er-ez)
    for step in steps:
        de,grad,stats=ring.sample(*step,True)
        vals[tuple(np.rint(step/h).astype(int))]=(de,grad,stats)
    E=lambda i,j:vals[i,j][0]
    Krr=(E(1,0)+E(-1,0))/h**2;Kzz=(E(0,1)+E(0,-1))/h**2
    Krz=(E(1,1)-E(1,-1)-E(-1,1)+E(-1,-1))/(4*h*h)
    K=np.array([[Krr,Krz],[Krz,Kzz]])
    grad_sec=np.column_stack([(vals[1,0][1]-vals[-1,0][1])/(2*h),(vals[0,1][1]-vals[0,-1][1])/(2*h)])
    return dict(step=h,analytic_force=(-slope).tolist(),
                energy_force=[(E(-1,0)-E(1,0))/(2*h),(E(0,-1)-E(0,1))/(2*h)],
                curvature=K.tolist(),principal_curvatures=np.linalg.eigvalsh(K).tolist(),
                gradient_secant=grad_sec.tolist(),base=base,
                signed_steps=[dict(direction=[int(q) for q in v],energy_change=vals[v][0],**vals[v][2]) for v in vals])


def main(a):
    data=source(a.state,a.snapshot)
    sh,d,_,_=data;order=order_tables(a.state,a.snapshot,data,a.fluid)
    emit(kind='source',state=file_sha(a.state),snapshot=file_sha(a.snapshot),
         native_terminal_field=a.fluid is None,fluid_sha=None if a.fluid is None else file_sha(a.fluid),
         profile=[sh.length,sh.exit],fluid_order_grid=[len(order[1][1])-1,len(order[1][0])-1])
    for j,i in a.nodes:
        for factor in a.factor:
            for shrink in (False,True) if a.shrinking and factor>1 else (False,):
                ring=Ring(a.state,sh,d,order,(j,i),factor,shrink,seams=a.seams)
                if a.partition:
                    _,slope,terms=ring.sample(0,0,True)
                    emit(kind='partition',node=[j,i],refinement=factor,seam_fitted=a.seams,
                         site=ring.site,force=(-slope).tolist(),**terms)
                    continue
                for h in a.steps:
                    row=signed(ring,h)
                    emit(kind='ring',node=[j,i],refinement=factor,shrinking=shrink,seam_fitted=a.seams,site=ring.site,
                         footprint_cell=list(ring.vcell),reference_volume=ring.reference_volume,**row)
                    if j==31 and i==32 and factor==1 and not a.seams and a.fluid is None and abs(row['analytic_force'][0]+.0018658717522196423)>2e-8:
                        raise ArithmeticError('Native force is not reproduced')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path);p.add_argument('snapshot',type=Path)
    p.add_argument('--fluid',type=Path)
    p.add_argument('--factor',nargs='+',type=int,default=[1,2,4,8,16])
    p.add_argument('--nodes',nargs=2,action='append',type=int)
    p.add_argument('--shrinking',action='store_true')
    p.add_argument('--seams',action='store_true')
    p.add_argument('--partition',action='store_true',
                   help='Only the first-variation terms from elastic energy and motion through the order gradient.')
    p.add_argument('--steps',nargs='+',type=float,default=[2e-4,1e-4,5e-5])
    a=p.parse_args();a.nodes=a.nodes or [(31,32),(30,32),(31,31)]
    main(a)
