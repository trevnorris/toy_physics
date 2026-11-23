"""Two-sided coherent drain with a straight, externally prepared phase vortex.

The complex phase is N phi plus an unconstrained meridional phase. Spherical
branches in the three spatial brane directions acquire a physical polar angle.
The straight vortex is factored out of the amplitude to integrate excess energy
over a fixed box without subtracting the much larger homogeneous vortex energy.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.integrate import solve_bvp
from scipy.sparse.linalg import spsolve

from mouth_exchange_carrier import BIAS, REGIONS, assemble, mapped_nodes
from mouth_exit_witness import inputs
from mouth_live_waist_fields import WIDTH, ORDER_INERTIA
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian


def straight_vortex(count, extent):
    """Regular solution for U(n)=(n**5-5*n+4)/20 with unit asymptotic density."""
    n=abs(count)
    if n==0:
        return lambda s:(np.ones_like(s),np.zeros_like(s))
    # At s->0 f~a s^n, at infinity f=1-n²/(4s²)+O(s^-4).
    lo=1.e-5
    hi=max(80.,2*extent)
    x=np.unique(np.r_[np.geomspace(lo,1.,60),np.linspace(1.,12.,121),np.linspace(12.,hi,121)])
    def rhs(s,y):
        f,fp=y
        return np.stack((fp, -fp/s+n*n*f/s**2+.5*(f**8-1)*f))
    def bc(a,b):
        return np.array([lo*a[1]-n*a[0],b[0]-1+n*n/(4*hi**2)])
    f=(x/np.hypot(x,1.4*n))**n
    init=np.stack((f,f*n*(1-x*x/(x*x+(1.4*n)**2))/x))
    sol=solve_bvp(rhs,bc,x,init,tol=2.e-7,max_nodes=8000)
    if not sol.success or np.any(sol.y[0]<=0) or np.min(sol.y[1]) < -1e-6:
        raise ArithmeticError('Regular homogeneous vortex boundary value failed')
    def evaluate(s):
        f,fp=sol.sol(np.maximum(s,lo))
        inside=s<lo
        return np.where(inside,sol.y[0,0]*(s/lo)**n,f),np.where(inside,sol.y[1,0]*(s/lo)**(n-1),fp)
    return evaluate


def prepared(radius):
    trial=Trial(radius,72,bonded_core=True)
    path=ROOT/f'_scratch/s11c/mouth_exit_release01/r{radius}n72/endpoint.npz'
    shifts,chi,_,sha,status,*_=inputs(path,radius,72,trial)
    prior=json.loads((ROOT/'research/pde_ledger_v3/_measurements/mouth_breathing_exit_support.json').read_bytes())
    if sha!=prior['physical_preparations'][f'r{radius}n72'] or status!='balanced':
        raise ValueError('Unmatched or unsupported held preparation')
    return (trial.x+shifts[:4]),chi,sha


def order_slopes(chi,x,t,radius,slope,length):
    nt,nr=np.subtract(chi.shape,1)
    i=np.minimum((x*nr).astype(int),nr-1)
    j=np.minimum((t*nt).astype(int),nt-1)
    a,b=x*nr-i,t*nt-j
    v0,v1,v2,v3=chi[j,i],chi[j,i+1],chi[j+1,i],chi[j+1,i+1]
    gr=nr*((1-b)*(v1-v0)+b*(v3-v2))/radius
    gw=nt*((1-a)*(v2-v0)+a*(v3-v1))/length-x*slope*gr
    return gr,gw


def local_matrix(ids,values,size):
    rows=np.broadcast_to(ids[:,:,None],values.shape)
    cols=np.broadcast_to(ids[:,None,:],values.shape)
    return sparse.coo_matrix((values.ravel(),(rows.ravel(),cols.ravel())),shape=(size,size)).tocsr()


def basis(a,b,c):
    factors=np.array([[1-a,a],[1-b,b],[1-c,c]])
    bits=np.array([[i,j,k] for i in (0,1) for j in (0,1) for k in (0,1)])
    terms=np.array([factors[i][bits[:,i]] for i in range(3)])
    return (terms.prod(axis=0),*((2*bits[:,i]-1)*terms[(i+1)%3]*terms[(i+2)%3] for i in range(3)))


class Wound:
    def __init__(self,shape,chi,count,n,nt,outer,roof):
        self.shape,self.count,self.n,self.nt=shape,count,n,nt
        self.vortex=straight_vortex(count,outer)
        self.radial,self.axial=mapped_nodes(n,shape.length,outer,roof)
        self.angles=np.pi/2*np.linspace(0.,1.,nt+1)**1.35
        _,_,grids,_,_,_,_=assemble(0.,n,WIDTH/12,outer,roof,lambda _:shape,True,True)
        radial_axis=np.unique(np.concatenate([grids[k,0][:,0] for k in range(3)]))
        ids=np.arange((max(g.max() for g in grids.values())+1)*(nt+1)).reshape(-1,nt+1)
        ids[radial_axis]=ids[radial_axis,:1]
        _,ids=np.unique(ids,return_inverse=True)
        ids=ids.reshape(-1,nt+1)
        self.grid={key:ids[grid] for key,grid in grids.items()}
        density=np.arange(ids.max()+1)
        density[self.grid[2,0][0].ravel()]=self.grid[1,0][-1].ravel()
        _,self.did=np.unique(density,return_inverse=True)
        self.nf=max(self.did)+1
        self.np=len(density)
        self.feed=np.unique(self.grid[0,1][:,-1])
        self.bulk=np.unique(np.concatenate([self.grid[1,1][:,-1].ravel(),self.grid[2,1][:,-1].ravel(),
                    self.grid[2,1][-1].ravel(),self.grid[2,0][-1].ravel()]))
        self.free=np.ones(self.np,bool); self.free[np.r_[self.feed,self.bulk]]=False
        self.ffree=np.ones(self.nf,bool);self.ffree[self.did[~self.free]]=False
        self.cells=[]; self.cap=None; self.wall=[]
        quad,wt=leggauss(3);quad=(quad+1)/2;wt=wt/2
        self.geometry=[]
        for slab,zs in enumerate(self.axial):
            for side,xs in enumerate(self.radial):
                grid=self.grid[slab,side]
                # Order of indices within each cube = radial, axial, polar bits.
                cid=np.stack([grid[j:n+j,i:n+i,k:nt+k] for i in (0,1) for j in (0,1) for k in (0,1)],axis=-1).reshape(-1,8)
                cells=len(cid); K=np.zeros((cells,8,8));P=np.zeros_like(K)
                vol=np.zeros(cells); rest=np.zeros(cells); breg=np.zeros(cells)
                x,t,th=np.meshgrid((xs[1:]+xs[:-1])/2,(zs[1:]+zs[:-1])/2,(self.angles[1:]+self.angles[:-1])/2,indexing='xy')
                x=x.ravel();w=t.ravel();th=th.ravel()
                dx,dw,dt=[a.ravel() for a in np.meshgrid(np.diff(xs),np.diff(zs),np.diff(self.angles),indexing='xy')]
                rad,slope=shape.radius(w); width=outer-rad if side else rad
                gr,gw=(order_slopes(chi,x,w/shape.length,rad,slope,shape.length) if (slab,side)==(1,0) else (np.zeros(cells),np.zeros(cells)))
                self.geometry.append((x,w,th,dx,dw,dt,gr,gw))
                for a,wa in zip(quad,wt):
                    xx=x+(a-.5)*dx
                    for b,wb in zip(quad,wt):
                        ww=w+(b-.5)*dw
                        radius,slope=shape.radius(ww);width=outer-radius if side else radius
                        rr=radius+xx*width if side else xx*radius
                        for c,wc in zip(quad,wt):
                            ang=th+(c-.5)*dt; sa=np.sin(ang)
                            _,gx,gz,gt=basis(a,b,c)
                            dr=gx/(width*dx)[:,None]
                            dz=gz/dw[:,None]-(slope*((1-xx) if side else xx))[:,None]*dr
                            da=gt/(rr*dt)[:,None]
                            measure=4*np.pi*rr**2*sa*width*dx*dw*dt*(wa*wb*wc)
                            f0,_=self.vortex(rr*sa)
                            weight=measure*f0*f0
                            vol+=measure; rest+=weight; breg+=measure*f0**10
                            for gradient in (dr,dz,da):
                                K+=weight[:,None,None]*(gradient[:,:,None]*gradient[:,None,:])
                            g=gr[:,None]*dr+gw[:,None]*dz
                            P+=weight[:,None,None]*(g[:,:,None]*g[:,None,:])
                self.cells.append((cid,K,P,rest/vol,breg,ORDER_INERTIA*(gr*gr+gw*gw)))
        self.cid=np.concatenate([a[0] for a in self.cells]);self.fid=self.did[self.cid]
        self.K=np.concatenate([a[1] for a in self.cells]);self.P=np.concatenate([a[2] for a in self.cells])
        self.n0=np.concatenate([a[3] for a in self.cells]); self.volU=np.concatenate([a[4] for a in self.cells])
        self.q=np.concatenate([a[5] for a in self.cells])
        self.A0=local_matrix(self.fid,self.K,self.nf)
        wallmat=sparse.csr_matrix((self.nf,self.nf))
        w0=self.axial[1][:-1,None];dw=np.diff(self.axial[1])[:,None]
        angle=self.angles[:-1][None,:];dt=np.diff(self.angles)[None,:]
        for side in (0,1):
            grid=self.grid[1,side][:,0 if side else -1]
            face=self.did[np.stack([grid[j:n+j,k:nt+k] for j in (0,1) for k in (0,1)],-1).reshape(-1,4)]
            W=np.zeros((len(face),4,4))
            for b,wb in zip(quad,wt):
                r,_=shape.radius(w0+b*dw)
                for c,wc in zip(quad,wt):
                    sa=np.sin(angle+c*dt);f0,df=self.vortex(r*sa)
                    weight=(1-2*side)*4*np.pi*r*r*sa**2*dw*dt*wb*wc*f0*df
                    vec=np.kron([1-b,b],[1-c,c])
                    W+=weight.reshape(-1,1,1)*np.outer(vec,vec)
            self.wall.append((face,W))
            wallmat=wallmat+local_matrix(face,W,self.nf)
        self.wallmat=wallmat
        # Density is continuous across the half-layer. The two terminal residual
        # phase traces are free; no winding jumps across this surface.
        cm=self.grid[1,0][-1]; cp=self.grid[2,0][0]
        def footprint(grid):
            return np.stack([grid[i:n+i,k:nt+k] for i in (0,1) for k in (0,1)],-1).reshape(-1,4)
        self.minus,self.plus=footprint(cm),footprint(cp)
        self.cids=np.column_stack((self.minus,self.plus));self.cf=self.did[self.minus]
        area=np.zeros((n,nt)); self.C=np.zeros((n*nt,8,8))
        xs=self.radial[0][:-1,None];dx=np.diff(self.radial[0])[:,None]
        for a,wa in zip(quad,wt):
            r=shape.exit*(xs+a*dx)
            for c,wc in zip(quad,wt):
                sa=np.sin(angle+c*dt); f0,_=self.vortex(r*sa)
                vec=np.kron([1-a,a],[1-c,c]); jump=np.r_[-vec,vec]
                measure=4*np.pi*r*r*sa*shape.exit*dx*dt*wa*wc
                area+=measure
                self.C+=(measure*f0**4/(WIDTH/12)).reshape(-1,1,1)*np.outer(jump,jump)
        self.area=area
        self.outer=outer
        allids=np.r_[self.cid,self.cids]
        row=np.broadcast_to(allids[:,:,None],(len(allids),8,8)).ravel()
        col=np.broadcast_to(allids[:,None,:],(len(allids),8,8)).ravel()
        self.keys=row*self.np+col
        self.keys,self.inverse=np.unique(self.keys,return_inverse=True)
        self.indices=self.keys%self.np
        self.ptr=np.r_[0,np.cumsum(np.bincount(self.keys//self.np,minlength=self.np))]

    def phase(self,v,bias):
        t=v[self.fid].mean(axis=1)
        den=self.n0*t*t+self.q
        k=t[:,None,None]**2*(self.K-(ORDER_INERTIA/den)[:,None,None]*self.P)
        c=v[self.cf].mean(axis=1)
        cap=c[:,None,None]**4*self.C
        local=np.concatenate([k.ravel(),cap.ravel()])
        mat=sparse.csr_matrix((np.bincount(self.inverse,weights=local),self.indices,self.ptr),shape=(self.np,self.np))
        phase=np.zeros(self.np);phase[self.bulk]=bias
        if bias!=0:
            phase[self.free]=spsolve(mat[self.free][:,self.free],-(mat@phase)[self.free])
        trace=mat@phase
        current=trace[self.bulk].sum()
        capj=np.sum((cap@phase[self.cids][...,None])[:,4:,0],axis=1)
        if max(abs(capj.sum()-current),abs(current+trace[self.feed].sum()))>2e-7*(1+abs(current)):
            raise ArithmeticError('Both reservoirs and cap do not carry the same current')
        return phase, current, capj

    def amplitude(self,v,phase,hessian=False):
        t=v[self.fid].mean(axis=1); x=t*t
        c=v[self.cf].mean(axis=1)
        q=phase[self.cid]-phase[self.cid[:,:1]]
        p=phase[self.cids]-phase[self.cids[:,:1]]
        ek=.5*np.einsum('bi,bij,bj->b',q,self.K,q)
        ep=.5*np.einsum('bi,bij,bj->b',q,self.P,q)
        ec=.5*np.einsum('bi,bij,bj->b',p,self.C,p)
        denom=self.n0*x+self.q
        one=ORDER_INERTIA/denom
        two=ORDER_INERTIA*self.q/denom**2
        phases=x*(ek-one*ep)
        defect=(x-1)**2*(x**3+2*x*x+3*x+4)
        interaction=self.volU*defect/20
        delta=v-1.
        elastic=.5*(delta@(self.A0@delta))
        boundary=.5*(delta@(self.wallmat@(v+1.)))
        # Exact reference area terms on the two side faces cancel for v=1.
        E=elastic+boundary+interaction.sum()+phases.sum()+(c**4*ec).sum()
        fder=2*t*(ek-two*ep)+.5*self.volU*t*(x**4-1)
        g=self.A0@delta+self.wallmat@v
        np.add.at(g,self.fid,np.broadcast_to(fder[:,None]/8,self.fid.shape))
        np.add.at(g,self.cf,np.broadcast_to((4*c**3*ec)[:,None]/4,self.cf.shape))
        if not hessian: return E,g
        fd2=2*(ek-two*ep)+8*x*self.n0*two*ep/denom+.5*self.volU*(9*x**4-1)
        h=(self.A0+self.wallmat+local_matrix(self.fid,np.broadcast_to(fd2[:,None,None]/64,self.K.shape),self.nf)
           +local_matrix(self.cf,np.broadcast_to((12*c*c*ec)[:,None,None]/16,(len(c),4,4)),self.nf))
        fields=dict(phase_and_entrainment=float(phases.sum()),terminal=float((c**4*ec).sum()),potential=float(interaction.sum()),
                    excess_amplitude_gradient=float(elastic),two_sided_boundary=float(boundary))
        return E,g,h,fields

    def solve(self,bias,start=None):
        v=np.ones(self.nf) if start is None else start.copy();free=self.ffree
        for iteration in range(60):
            phase,current,capj=self.phase(v,bias)
            E,g,H,terms=self.amplitude(v,phase,True)
            change=spsolve(H[free][:,free],-g[free])
            err=np.max(abs(change))
            if iteration in (0,1,2,5,10,20,35,50,59):
                emit(kind='convergence',iteration=iteration,energy=E,min_ratio=v.min(),max_ratio=v.max(),
                     max_newton_step=err,max_gradient=float(np.max(abs(g[free]))),current=current)
            if err<1.e-8:
                return v,phase,dict(excess_energy=E,total_current=current,capj=capj,terms=terms,
                                    min_amplitude_ratio=v.min(),max_amplitude_ratio=v.max(),residual=err,iterations=iteration)
            step=1.; slope=g[free]@change
            if slope>=0:
                raise ArithmeticError('Unstable local material Hessian')
            while step>1.e-7:
                nxt=v.copy();nxt[free]+=step*change
                if nxt.min()>0:
                    new,_=self.amplitude(nxt,phase)
                    if new<=E+1.e-4*step*slope+2.e-12*(1+abs(E)):
                        break
                step*=.5
            else:raise ArithmeticError('No positive descending amplitude step')
            v=nxt
        raise ArithmeticError('Unresolved amplitude/phase stationarity')

    def center_fields(self,v,phase,info):
        """Disk quadrature and gradient of the held core energy w.r.t. native order slopes."""
        n=self.n; nt=self.nt
        grid=self.grid[1,0]
        ids=np.stack([grid[j:n+j,i:n+i,k:nt+k] for i in (0,1) for j in (0,1) for k in (0,1)],-1).reshape(-1,8)
        x,w,th,dx,dw,dt,gr,gw=self.geometry[2]
        matv=v[self.did[ids]].mean(axis=-1)
        order_n=self.cells[2][3]*matv**2
        den=order_n+ORDER_INERTIA*(gr*gr+gw*gw)
        dgr=np.zeros(len(gr));dgw=np.zeros(len(gr))
        quad,weights=leggauss(3);quad=(quad+1)/2;weights/=2
        caprows=[]
        p=phase[ids]
        for a,wa in zip(quad,weights):
            xx=x+(a-.5)*dx
            for b,wb in [*zip(quad,weights),(1.,0.)]:
                ww=w+(b-.5)*dw
                br,slope=self.shape.radius(ww);rr=xx*br
                for c,wc in zip(quad,weights):
                    ang=th+(c-.5)*dt; axis_distance=rr*np.sin(ang)
                    f0,_=self.vortex(axis_distance)
                    val,gx,gz,ga=basis(a,b,c)
                    dpr=(p@gx)/(dx*br)
                    dpw=(p@gz)/dw-xx*slope*dpr
                    gpd=gr*dpr+gw*dpw
                    vol=4*np.pi*rr*rr*np.sin(ang)*br*dx*dw*dt*wa*wb*wc
                    weight=-vol*f0*f0*matv*matv*ORDER_INERTIA/den
                    dgr+=weight*(gpd*dpr-ORDER_INERTIA*gr*gpd*gpd/den)
                    dgw+=weight*(gpd*dpw-ORDER_INERTIA*gw*gpd*gpd/den)
                    if b!=1.: continue
                    last=ww>=self.shape.length-1.e-12
                    density=f0[last]**2*(v[self.did[ids[last]]]@val)**2
                    cminus=phase[self.minus];cplus=phase[self.plus]
                    capbasis=np.kron([1-a,a],[1-c,c])
                    # Match the assembled terminal current cell by cell.
                    vr=v[self.cf].mean(axis=1)
                    j=vr**4*f0[last]**4*((cplus-cminus)@capbasis)/(WIDTH/12)
                    ur=(dpr-gr*ORDER_INERTIA*gpd/den)[last]
                    utan=(p@ga/(dt*rr))[last]
                    uspin=self.count/axis_distance[last]
                    area=4*np.pi*rr[last]**2*np.sin(ang[last])*br[last]*dx[last]*dt[last]*wa*wc
                    caprows.append((area,xx[last],ang[last],axis_distance[last],j,density,ur,utan,uspin))
        samples=np.column_stack([np.concatenate(samples) for samples in zip(*caprows)])
        if abs(samples[:,0]@samples[:,4]-info['total_current'])>2.e-7*(1+abs(info['total_current'])):
            raise ArithmeticError('Native and independently sampled terminal current disagree')
        return dict(disk=samples,core_x=x,core_t=w/self.shape.length,core_theta=th,core_order_gr_work=dgr,
                    core_order_gw_work=dgw,core_mean_density=order_n,core_g_r=gr,core_g_w=gw,
                    core_amp=v[self.did[grid]],core_phase=phase[grid],core_x_nodes=self.radial[0],
                    core_height_nodes=self.axial[1]/self.shape.length,core_angle_nodes=self.angles,
                    amplitude_ratio=v,phase=phase,mesh=self.n,polar=self.nt,outer=self.outer,roof=self.axial[2][-1])


def run(args):
    cfg,chi,sha=prepared(args.radius)
    b=Meridian(args.radius,cfg).exit
    emit(kind='input',radius=args.radius,winding=args.wind,circulation=float(2*np.pi*args.wind),
         mesh=args.mesh,polar_cells=args.polar,reference_exit=b,outer=args.outer,roof=args.roof,
         material_sha256=sha,half_offsets=args.half_offsets,mode=args.mode,
         shape_order='saved bonded material, fixed order labels',reservoir='matched straight vortex N phi on both leads, north-south even domain')
    sizes=[float(t) for t in args.half_offsets.split(',')]
    offsets=[0.] if args.mode=='center' else [v for h in sizes for v in (-h,h)]
    results={}
    for offset in offsets:
        changed=cfg.copy();changed[1]+=np.log1p(offset)
        sh=Meridian(args.radius,changed)
        model=Wound(sh,chi,args.wind,args.mesh,args.polar,args.outer,args.roof)
        v,phase,info=model.solve(BIAS)
        capj=info.pop('capj');info['capj']=capj
        if offset==0:
            fields=model.center_fields(v,phase,info)
            disk=fields['disk']
            area,_,_,_,j,dens,ur,ut,up=disk.T
            speed=np.sqrt((j/dens)**2+ur**2+ut**2+up**2)
            ma=speed/dens**2
            details=dict(density_min=float(dens.min()),density_max=float(dens.max()),
                         sampled_total_mach_min=float(ma.min()),sampled_total_mach_max=float(ma.max()),
                         supersonic_area_fraction=float(area[ma>=1].sum()/area.sum()),
                         supersonic_current_fraction=float(area[ma>=1]@j[ma>=1]/(area@j)),
                         meridional_max_mach=float(np.max(np.sqrt((j/dens)**2+ur**2+ut**2)/dens**2)),
                         equator_rim_mach=float(ma[np.argmax(disk[:,3])]),
                         rim_equator_s=float(max(disk[:,3])))
            if args.save is not None:
                np.savez_compressed(args.save,**fields,material_sha256=sha,winding=args.wind,radius=args.radius,
                        disk_columns=['area','r/b','polar_angle','s','j_w','density','v_r','v_polar','v_phi'])
            emit(kind='center_readout',**details,shape_exit=sh.exit,depth=sh.length)
        del info['capj']
        emit(kind='geometry',relative_exit=offset,**info)
        results[offset]=info
    if args.mode=='signed':
        for h in sizes:
            m,p=results[-h],results[h]
            F=-(p['excess_energy']-m['excess_energy'])/(2*h*b)
            components={k:-(p['terms'][k]-m['terms'][k])/(2*h*b) for k in p['terms']}
            emit(kind='signed_work',radius=args.radius,winding=args.wind,half_offset=h,force_on_exit=F,
                 components=components)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('radius',type=int,choices=(1,2))
    parser.add_argument('wind',type=int,choices=(-2,-1,0,1,2))
    parser.add_argument('--mode',choices=('center','signed'),default='signed')
    parser.add_argument('--half-offsets',default='.006')
    parser.add_argument('--mesh',type=int,default=8)
    parser.add_argument('--polar',type=int,default=7)
    parser.add_argument('--outer',type=float,default=28.)
    parser.add_argument('--roof',type=float,default=22.)
    parser.add_argument('--save',type=Path)
    run(parser.parse_args())
