"""Single neutral carrier between a lateral feeding slab and bulk, through a rounded lip.

Axisymmetric in three brane dimensions. Quadratic meridional elements cover
both reservoirs and the smooth semicircular end of the reflecting separator.
All physical lengths are in attachment radii. The aperture has no boundary data.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.integrate import simpson
from scipy.sparse.linalg import spsolve


RATE = .25
K_FRONT = .4
RIM_SLOPE = RATE/(4*np.pi*K_FRONT)


def quadratic(x):
    return np.stack((x*(x-1)/2,1-x*x,x*(x+1)/2),axis=-1), np.stack((x-.5,-2*x,x+.5),axis=-1)


def tensor(a,b):
    return (a[..., :,None]*b[...,None,:]).reshape(*a.shape[:-1],9)


def divided(edges):
    values=np.empty(2*len(edges)-1)
    values[::2]=edges
    values[1::2]=(edges[1:]+edges[:-1])/2
    return values


class Domain:
    def __init__(self, n, far, rounding=.5, slope=RIM_SLOPE):
        self.far,self.rounding,self.slope=far,rounding,slope
        angle=np.arctan(slope)
        self.rc=1+rounding*np.sin(angle)
        self.lip=rounding*(1-np.cos(angle))
        self.beta_lip=angle-np.pi/2
        self.bottom_turn=-np.arctan((1+rounding)/self.rc)
        self.top_turn=np.arctan((far+rounding)/self.rc)
        splits=[-np.pi/2,self.beta_lip,self.bottom_turn,0,self.top_turn,np.pi/2]
        splits=np.unique(np.round(splits,14))
        angular=np.r_[np.concatenate([np.linspace(a,b,max(1,int(np.ceil((b-a)*n))),endpoint=False)
                                     for a,b in zip(splits[:-1],splits[1:],strict=True)]),splits[-1]]
        stretch=np.sinh(2.4*np.linspace(0,1,n+1))/np.sinh(2.4)
        radial=self.rc+(far-self.rc)*stretch
        self.blocks={}; all_nodes=[]; self.dofs=[]
        coeff=[]; interp_rows=[]; interp_cols=[]; interp_values=[]; weights=[]; q_offset=0
        roots,w=leggauss(4)
        p_loc,v_loc=np.meshgrid(roots,roots,indexing='ij'); p_loc=p_loc.ravel(); v_loc=v_loc.ravel()
        p0,dp0=quadratic(p_loc); v0,dv0=quadratic(v_loc)
        val=tensor(p0,v0); derivp=tensor(dp0,v0); derivv=tensor(p0,dv0)
        node_offset=0
        for block,pgrid in enumerate((angular,radial,radial)):
            vgrid=stretch
            pp,vv=np.meshgrid(divided(pgrid),divided(vgrid),indexing='ij')
            r,z,*_=self.chart(block,pp,vv)
            all_nodes.append(np.c_[r.ravel(),z.ravel()])
            ids=np.arange(pp.size).reshape(pp.shape)+node_offset; node_offset+=pp.size
            i,j=np.indices((len(pgrid)-1,len(vgrid)-1)); i=i.ravel(); j=j.ravel()
            dofs=(2*i[:,None,None]+np.arange(3)[None,:,None])*ids.shape[1]+2*j[:,None,None]+np.arange(3)[None,None,:]+ids[0,0]
            dofs=dofs.reshape(-1,9); self.dofs.append(dofs)
            self.blocks[block]=(pgrid,vgrid,dofs)
            hp=np.diff(pgrid)[i,None]/2; hv=np.diff(vgrid)[j,None]/2
            p=pgrid[i,None]+hp*(1+p_loc); v=vgrid[j,None]+hv*(1+v_loc)
            r,z,rp,rv,zp,zv=self.chart(block,p,v)
            det=rp*zv-zp*rv
            gp=derivp[None,:,:]/hp[:,:,None]; gv=derivv[None,:,:]/hv[:,:,None]
            gr=(zv[:,:,None]*gp-zp[:,:,None]*gv)/det[:,:,None]
            gz=(-rv[:,:,None]*gp+rp[:,:,None]*gv)/det[:,:,None]
            measure=4*np.pi*r*r*abs(det)*hp*hv*(w[:,None]*w).ravel()
            coeff.append(np.einsum('eq,eqi,eqj->eij',measure,gr,gr)+np.einsum('eq,eqi,eqj->eij',measure,gz,gz))
            weights.append(measure.ravel()); idx=np.arange(measure.size).reshape(measure.shape)+q_offset
            q_offset+=measure.size
            interp_rows.append(np.broadcast_to(idx[:,:,None],gr.shape).ravel())
            interp_cols.append(np.broadcast_to(dofs[:,None,:],gr.shape).ravel())
            interp_values.append(np.broadcast_to(val[None,:,:],gr.shape).ravel())
        points,mapping=np.unique(np.round(np.concatenate(all_nodes),11),axis=0,return_inverse=True)
        self.points=points
        for block,(pgrid,vgrid,dofs) in self.blocks.items():
            self.blocks[block]=(pgrid,vgrid,mapping[dofs])
        dofs=mapping[np.concatenate(self.dofs)]; count=len(points)
        self.stiff=sparse.coo_matrix((np.concatenate(coeff).ravel(),(np.repeat(dofs,9,axis=1).ravel(),
            np.tile(dofs,(1,9)).ravel())),shape=(count,count)).tocsr()
        self.N=sparse.coo_matrix((np.concatenate(interp_values),(np.concatenate(interp_rows),mapping[np.concatenate(interp_cols)])),shape=(q_offset,count)).tocsr()
        self.weight=np.concatenate(weights)
        self.feed=(points[:,0]==far)&(points[:,1]<=0)
        self.bulk=((points[:,0]==far)&(points[:,1]>=2*rounding)) | (points[:,1]==2*rounding+far)
        self.free=~(self.feed|self.bulk)
        assert np.any(self.feed) and np.any(self.bulk) and np.max(abs(self.stiff@np.ones(count)))<1e-8

    def reach(self,beta):
        c,s=np.cos(beta),np.sin(beta)
        middle=(beta>=self.bottom_turn)&(beta<=self.top_turn)
        m=np.where(middle,self.rc/np.maximum(c,1e-12),
                   np.where(s<0,-(1+self.rounding),self.far+self.rounding)/np.where(abs(s)>1e-10,s,1))
        dm=np.where(middle,m*np.tan(beta),-m*c/np.where(abs(s)>1e-10,s,1))
        return m,dm

    def chart(self,block,p,v):
        one=np.ones_like(p); zero=np.zeros_like(p)
        if block:
            a=-1 if block==1 else self.far
            return p,(0 if block==1 else 2*self.rounding)+a*v,one,zero,zero,a*one
        m,dm=self.reach(p); c,s=np.cos(p),np.sin(p)
        t=self.rounding+(m-self.rounding)*v
        return (self.rc-t*c, self.rounding+t*s, t*s-dm*v*c, -(m-self.rounding)*c,
                t*c+dm*v*s, (m-self.rounding)*s)

    def read(self,field,r,z):
        # The aperture and lip lie wholly within the rounded-end chart.
        p=np.arctan2(z-self.rounding,self.rc-r)
        t=np.hypot(z-self.rounding,self.rc-r)
        m,_=self.reach(p); v=(t-self.rounding)/(m-self.rounding)
        pg,vg,dofs=self.blocks[0]
        i=np.clip(np.searchsorted(pg,p)-1,0,len(pg)-2); j=np.clip(np.searchsorted(vg,v)-1,0,len(vg)-2)
        hp=np.diff(pg)[i]/2; hv=np.diff(vg)[j]/2
        ap,dap=quadratic((p-pg[i])/hp-1); av,dav=quadratic((v-vg[j])/hv-1)
        values=field[dofs[i*(len(vg)-1)+j]]
        value=np.sum(tensor(ap,av)*values,axis=-1)
        fp=np.sum(tensor(dap,av)*values,axis=-1)/hp; fv=np.sum(tensor(ap,dav)*values,axis=-1)/hv
        _,_,rp,rv,zp,zv=self.chart(0,p,v)
        det=rp*zv-zp*rv
        return value,(zv*fp-zp*fv)/det,(-rv*fp+rp*fv)/det

    def current(self,field,r,z):
        f,fr,fz=self.read(field,r,z)
        return (f.conj()*fr).imag,(f.conj()*fz).imag

    def solve(self,rate):
        # Only the relative phase of the distant baths is adjusted. Reuse the
        # stationary complex-carrier equation and bath normalization of Stage 1.
        K=self.stiff*.5; I=self.free; B=self.bulk
        A=K[I][:,I]; N=self.N; Nf=N[:,I]; wt=self.weight
        harmonic=B.astype(float); harmonic[I]=spsolve(A,-(K@harmonic)[I])
        conductance=float((self.stiff@harmonic)[B].sum()); delta=rate/conductance
        psi=np.exp(1j*delta*harmonic); prev=None
        for outer in range(8):
            psi[B]=np.exp(1j*delta)
            for step in range(12):
                f=N@psi; n=abs(f)**2; local=wt*(n**4-1)/4
                residual=K@(psi-1)+N.T@(local*f)
                error=float(np.max(abs(residual[I])))
                if error<3e-11:
                    break
                slope=wt*2*n**3
                def mass(v):
                    return Nf.T@Nf.multiply(v[:,None])
                xy=mass(slope*f.real*f.imag)
                jac=sparse.bmat([[A+mass(local+slope*f.real**2),xy],
                                 [xy,A+mass(local+slope*f.imag**2)]],format='csc')
                correction=spsolve(jac,-np.r_[residual[I].real,residual[I].imag]); dx,dy=np.split(correction,2)
                psi[I]+=dx+1j*dy
            else:
                raise ArithmeticError(f'Carrier equation failed to converge at bias {delta}: residual {error}')
            out=2*np.sum((psi[B].conj()*residual[B]).imag)
            feed=-2*np.sum((psi[self.feed].conj()*residual[self.feed]).imag)
            if abs(out-rate)<3e-11:
                if abs(feed-out)>3e-9 or not np.all(n>.5):
                    raise ArithmeticError('Signed bath flux or finite-density check failed')
                return psi,dict(bias=delta,feed=feed,bulk=out,n_min=n.min(),n_max=n.max(),
                    residual=np.max(abs(residual[I])),conductance=conductance)
            derivative=conductance if prev is None else (out-prev[1])/(delta-prev[0])
            prev=(delta,out); delta+=(rate-out)/derivative
        raise ArithmeticError('Bath calibration failed')


def cap(domain,field,r,factor):
    slope=factor*domain.slope
    return domain.lip+slope*(1-r*r)/2,-slope*r


def readings(domain,field,r,h,hp):
    jr,jw=domain.current(field,r,h)
    q=jw-hp*jr
    return dict(total=simpson(4*np.pi*r*r*q,x=r),center_q=q[0],q_half=q[len(r)//2],
        q_at_09=q[round(.9*(len(r)-1))],q_at_099=q[round(.99*(len(r)-1))],
        min_to_099=q[:round(.99*(len(r)-1))+1].min(),min_to_0999=q[:round(.999*(len(r)-1))+1].min(),
        min_all=q.min(),raw_lip_q=q[-1],lip_jr=jr[-1],lip_jw=jw[-1],lip_wall_error=jw[-1]+domain.slope*jr[-1],
        max_depth=h.max(), slope_at_lip=-hp[-1], min_wall_clearance=(np.hypot(domain.rc-r,h-domain.rounding)-domain.rounding).min())


def integrated_radial_current(domain,field,r,z,panels):
    # Integrate from the reflecting bottom to the requested surface. By
    # continuity this gives the axial flux inside radius r divided by -4πr².
    roots,weights=VERTICAL_RULE
    u=((np.arange(panels)[:,None]+(roots+1)/2)/panels).ravel()
    weights=np.tile(weights/(2*panels),panels)
    r=np.atleast_1d(r); z=np.atleast_1d(z)
    values=[]
    for start in range(0,len(r),50):
        rs=r[start:start+50,None]; zs=z[start:start+50,None]
        heights=-1+(zs+1)*u
        jr,_=domain.current(field,np.broadcast_to(rs,heights.shape),heights)
        values.extend((zs[:,0]+1)*(jr@weights))
    return np.array(values)


VERTICAL_RULE=leggauss(4)


def front(domain,field,K,steps,panels):
    # In spherical symmetry -K Δh=q integrates once. Use independently
    # integrated radial flux, which also avoids differentiating element jumps.
    from scipy.interpolate import CubicHermiteSpline
    r=np.linspace(1,0,steps+1); h=np.zeros_like(r); h[0]=domain.lip
    def slope(radius,height):
        return integrated_radial_current(domain,field,radius,height,panels)[0]/K
    for i,dr in enumerate(np.diff(r)):
        a=slope(r[i],h[i]); b=slope(r[i]+dr/2,h[i]+dr*a/2)
        c=slope(r[i]+dr/2,h[i]+dr*b/2); d=slope(r[i+1],h[i]+dr*c)
        h[i+1]=h[i]+dr*(a+2*b+2*c+d)/6
    hp=integrated_radial_current(domain,field,r,h,panels)/K
    # A regular axisymmetric continuum solution has zero radial current at r=0.
    hp[-1]=0
    return CubicHermiteSpline(r[::-1],h[::-1],hp[::-1])


def run(args):
    dom=Domain(args.n,args.far)
    f,state=dom.solve(args.rate)
    print(json.dumps(dict(test='stationary',n=args.n,far=args.far,rate=args.rate,**state)),flush=True)
    r=np.linspace(0,1,2001)
    for fac in (0.,1.):
        h,hp=cap(dom,f,r,fac); result=readings(dom,f,r,h,hp)
        print(json.dumps(dict(test='cap',factor=fac,**result)),flush=True)
    data={'r':r,'psi':f,'points':dom.points}
    for k in (.4,.8):
        solution=front(dom,f,k,args.steps,args.panels); h,hp=solution(r),solution(r,1)
        result=readings(dom,f,r,h,hp)
        shell_r=np.r_[np.linspace(0,.98,50),.99,1.]
        cylinder=-4*np.pi*shell_r**2*integrated_radial_current(dom,f,shell_r,solution(shell_r),2*args.panels)
        shell_increments=np.diff(cylinder)
        wall_curvature=(1+dom.slope**2)**1.5/dom.rounding
        rim_derivative=result['lip_jr']*(wall_curvature-2*dom.slope) if k==K_FRONT and args.rate==RATE else None
        print(json.dumps(dict(test='kinetic_front',K=k,steps=args.steps,vertical_panels=args.panels,
             independent_last_shell=shell_increments[-1],min_shell=shell_increments.min(),
             cylinder_total=cylinder[-1],grazing_qprime_from_slip=rim_derivative,**result)),flush=True)
        data[f'h_{k}']=h; jr,jw=dom.current(f,r,h); data[f'q_{k}']=jw-hp*jr
        if abs(result['total']-args.rate)>max(.0009,.005*args.rate) or not np.isfinite(h).all():
            raise ArithmeticError('Crossing integral disagrees with reservoirs')
    if args.output:
        np.savez_compressed(args.output,**data)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,required=True); parser.add_argument('--far',type=float,required=True)
    parser.add_argument('--rate',type=float,default=RATE); parser.add_argument('--output',type=Path)
    parser.add_argument('--steps',type=int,default=40); parser.add_argument('--panels',type=int,default=24)
    run(parser.parse_args())
