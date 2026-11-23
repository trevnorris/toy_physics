"""Axisymmetric coherent carrier on the two sides of a moving radius-two meridian.

Phase and density-amplitude have independent spatial degrees of freedom. The thin-link
reference takes the limit of a translating Route A kink. Its two phase traces can slip;
sqrt(n) is continuous at that disk. Real-coordinate assembly supplies independent
reaction/flux and signed geometry controls. Harmonic convention exp(-i omega t).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import splu

from mouth_exchange_carrier import REGIONS
from mouth_live_waist_fields import WIDTH

def gram(x,y):
    return x[..., :,None]*y[...,None,:]

def sparse_blocks(ids,data,size):
    return sparse.coo_matrix((data.ravel(),(np.broadcast_to(ids[:,:,None],data.shape).ravel(),
        np.broadcast_to(ids[:,None,:],data.shape).ravel())),shape=(size,size)).tocsr()

def power_matrix(W,omega):
    return -.5*omega*(W-W.conj().T)/(2j)


def mean_front(j, kt=0., extent=.48, refinement=2):
    """Planar intrinsic front conditioned on a subsonic mean current in the disk."""
    from mouth_curved_cap_ports import port
    def sheet(omega,r):
        matrix=port(omega,j,kt,extent,refinement)
        return np.broadcast_to(matrix,(len(r),5,5))
    return sheet


class CurvedSound:
    def __init__(self,shape,exit_scale,n=12,box=7.,pml=4.,sheet=None,base=None,
                 stretch=2.1,advection=True):
        self.shape,self.exit_scale,self.n=shape,exit_scale,n
        self.sheet,self.base,self.advection=sheet,base,advection
        self.box,self.pml,self.stretch=box,pml,stretch
        self._build()
        self._dc=None

    def coordinates(self,x):
        h=np.maximum(x-self.box,0)/self.pml
        return x+1j*self.stretch*self.pml*h**3/3,1+1j*self.stretch*h**2

    def _fields(self,slab,side,r,w):
        if self.base is None:
            return (np.ones_like(r),)+(np.zeros_like(r),)*4
        f=self.base.details(REGIONS[2*slab+side],r,w)
        if self.advection:
            return f[:5]
        return f[0],0*f[1],0*f[2],f[3],f[4]

    def _build(self):
        s,n=self.shape,self.n
        u=np.linspace(0.,1.,n+1)
        # The mesh topology follows the native six welded patches. Near the lip the
        # horizontal grids on both sides cluster toward it.
        xi=1-(1-u)**1.45
        nz=n; nf=max(5,(n+1)//2)
        nr=max(n,int(np.ceil((self.box+self.pml-s.neck)*n/1.7)))
        nw=max(n,int(np.ceil((self.box+self.pml-s.length)*n/1.7)))
        xout=np.linspace(0,1,nr+1)**1.22
        axes=(np.linspace(-1.,0.,nf+1),np.linspace(0.,s.length,nz+1),
              s.length+(self.box+self.pml-s.length)*np.linspace(0,1,nw+1)**1.15)
        coords=(xi,xout)
        self.axes,self.coords=axes,coords
        gridids={}; nodes={}; self.patches={}
        for slab,zs in enumerate(axes):
            for side,xs in enumerate(coords):
                grid=np.empty((len(zs),len(xs)),int)
                for k,i in np.ndindex(grid.shape):
                    key=((0,nf,nf+nz)[slab]+k,side*n+i,int(slab==0 or (slab==1 and side==0)))
                    if key not in nodes: nodes[key]=len(nodes)
                    grid[k,i]=nodes[key]
                gridids[slab,side]=grid
        self.grid=gridids
        total=len(nodes); density=np.arange(total)
        if self.sheet is None:
            density[gridids[2,0][0]]=gridids[1,0][-1]
        _,density=np.unique(density,return_inverse=True)
        self.size=total+max(density)+1
        self.density=total+density
        self.dirichlet=np.unique(np.r_[*(gridids[i,1][:,-1] for i in range(3)),
                                  gridids[2,0][-1],gridids[2,1][-1]])
        self.free=np.setdiff1d(np.arange(self.size),np.r_[self.dirichlet,self.density[self.dirichlet]])
        self.K=[];self.G=[]
        self.B0=np.zeros((self.size,3));self.Bv=np.zeros_like(self.B0)
        gp,gw=leggauss(3);gp=(gp+1)/2;gw/=2
        for slab,zs in enumerate(axes):
            for side,xs in enumerate(coords):
                grid=gridids[slab,side]
                ids=np.stack([grid[:-1,:-1],grid[:-1,1:],grid[1:,:-1],grid[1:,1:]],-1).reshape(-1,4)
                gi=np.concatenate([ids,self.density[ids]],axis=1)
                dims=(len(zs)-1,len(xs)-1,4,4)
                ktt=np.zeros(dims,complex);kaa=np.zeros_like(ktt);kta=np.zeros_like(ktt);gt=np.zeros_like(ktt)
                # Real, unstretched FE cut at each physical port for independent flux extraction.
                dx=np.diff(xs)[None,:]; dw=np.diff(zs)[:,None]
                for b,wb in zip(gp,gw):
                    w=zs[:-1,None]+b*dw; radius,slope=s.radius(w)
                    width=(self.box+self.pml-radius) if side else radius
                    for a,wa in zip(gp,gw):
                        y=xs[None,:-1]+a*dx
                        r=radius+width*y if side else width*y
                        rc,sr=self.coordinates(r); _,sz=self.coordinates(w)
                        N=np.array([(1-a)*(1-b),a*(1-b),(1-a)*b,a*b])
                        dr=np.array([b-1,1-b,-b,b])/(dx*width)[...,None]
                        dz=np.array([a-1,-a,1-a,a])/dw[...,None]-(((1-y) if side else y)*slope)[...,None]*dr
                        dr=dr/sr[...,None]; dz=dz/sz[...,None]
                        vol=4*np.pi*wa*wb*rc**2*width*dx*dw*sr*sz
                        rho,ur,uz,_,_=self._fields(slab,side,r,w)
                        coupling=2*np.sqrt(rho)[...,None]*(ur[...,None]*dr+uz[...,None]*dz)
                        v=vol[...,None,None]
                        grad=gram(dr,dr)+gram(dz,dz)
                        ktt+=(rho*vol)[...,None,None]*grad
                        # Hessian of coherent grand Hamiltonian at bath chemical frequency zero.
                        local=ur**2+uz**2+(rho**4-1)/2+4*rho**4
                        kaa+=v*(grad+local[...,None,None]*gram(N,N))
                        kta+=v*gram(coupling,N)
                        gt+=(2*np.sqrt(rho)*vol)[...,None,None]*gram(N,N)
                block=np.concatenate([np.concatenate([ktt,kta],axis=-1),
                                     np.concatenate([kta.swapaxes(-2,-1),kaa],axis=-1)],axis=-2)
                z=np.zeros_like(gt)
                gblock=np.concatenate([np.concatenate([z,gt],axis=-1),np.concatenate([-gt,z],axis=-1)],axis=-2)
                self.K.append(sparse_blocks(gi,block.reshape(-1,8,8),self.size))
                self.G.append(sparse_blocks(gi,gblock.reshape(-1,8,8),self.size))
                # Flux test at an interior w/r cut: use the nodal reaction of outer
                # elements projected onto boundary nodes, not the PML anti-Hermitian form.
                rmid=s.radius((zs[:-1,None]+zs[1:,None])/2)[0]
                rc=(xs[None,:-1]+xs[None,1:])/2
                rmid=rmid+(self.box+self.pml-rmid)*rc if side else rmid*rc
                out=((rmid>self.box-1.0)|(((zs[:-1]+zs[1:])/2)[:,None]>self.box-1.0)).ravel()
                if np.any(out):
                    self.patches[slab,side]=(sparse_blocks(gi[out],block.reshape(-1,8,8)[out],self.size),
                                            sparse_blocks(gi[out],gblock.reshape(-1,8,8)[out],self.size),gi[out],gi[~out])
        self.wall_data=[]
        zs=axes[1]; dw=np.diff(zs)
        for side in (0,1):
            # The moving side bounds BOTH domains. The same wall motion gives opposite
            # signed volume and base-flow shape work on its two traces.
            grid=gridids[1,side]; ix=0 if side else len(coords[side])-2
            ids=np.stack([grid[:-1,ix],grid[:-1,ix+1],grid[1:,ix],grid[1:,ix+1]],axis=-1)
            x=0 if side else 1
            for t,wt in zip(gp,gw):
                w=zs[:-1]+t*dw; r,slope=s.radius(w)
                tau=w/s.length
                h=np.column_stack([-slope*s.length*tau,self.exit_scale*tau,4*s.neck*tau*(1-tau)])
                dx=coords[side][ix+1]-coords[side][ix]
                width=self.box+self.pml-r if side else r
                N=np.broadcast_to([(1-x)*(1-t),x*(1-t),(1-x)*t,x*t],ids.shape)
                dr=np.array([t-1,1-t,-t,t])/(dx*width[:,None])
                dz=np.array([x-1,-x,1-x,x])/dw[:,None]-slope[:,None]*dr
                rho,ur,uz,fr,fw=self._fields(1,side,r,w)
                area=(1-2*side)*4*np.pi*r**2*wt*dw
                vol=area[:,None]*h
                np.add.at(self.Bv,ids,gram(rho[:,None]*N,vol))
                np.add.at(self.B0,ids,-gram(ur[:,None]*rho[:,None]*dr+uz[:,None]*rho[:,None]*dz,vol))
                pot=np.sqrt(rho)*(ur**2+uz**2+(rho**4-1)/2)
                np.add.at(self.B0,self.density[ids],-gram(fr[:,None]*dr+fw[:,None]*dz+pot[:,None]*N,vol))
                if np.isclose(t,.5):
                    self.wall_data.append((ids,N,rho,ur,uz,dr,dz))
        a=gridids[1,0][-1]; b=gridids[2,0][0]; self.capi=np.column_stack([a[:-1],a[1:],b[:-1],b[1:]])
        self.radial_nodes=xi*s.exit
        self.radii=(self.radial_nodes[:-1]+self.radial_nodes[1:])/2
        self.disk_weights=4*np.pi*np.diff(self.radial_nodes**3)/3
        self.cap_quad=[]
        for t,wt in zip(gp,gw):
            rad=self.radial_nodes[:-1]+t*np.diff(self.radial_nodes)
            self.cap_quad.append((rad,4*np.pi*rad**2*wt*np.diff(self.radial_nodes),np.array([1-t,t])))
        self.k=sum(self.K);self.g=sum(self.G)
        self.outs=[]
        for subset in ((0,),(1,2)):
            Ks=[];Gs=[];outside=[];inside=[]
            for slab in subset:
                for side in (0,1):
                    if (slab,side) in self.patches:
                        k,g,outer,inner=self.patches[slab,side]
                        Ks.append(k);Gs.append(g);outside.extend(outer.ravel());inside.extend(inner.ravel())
                    else:
                        ids=gridids[slab,side].ravel()
                        inside.extend(np.r_[ids,self.density[ids]])
            cut=np.intersect1d(outside,inside)
            self.outs.append((sum(Ks),sum(Gs),cut))

    def _cap(self,w):
        # External forces on the layer: theta rows = (-j,+j) upward actual current j.
        # Thin moving kink: theta_plus-theta_minus=ell*(j-V), V=-i omega X.
        n=self.n
        ids=self.capi
        allids=np.concatenate([ids,self.density[ids]],axis=1)
        kp=np.zeros((n,8,8),complex);ky=np.zeros((self.size,3),complex);yk=np.zeros((3,self.size),complex)
        xx=0j
        for r,area,N in self.cap_quad:
            P=np.zeros((4,8))
            for i in range(4):P[i,2*i:2*i+2]=N
            if self.sheet is None:
                S=np.zeros((n,5,5),complex)
                g=np.array([-1.,1.]); ell=WIDTH/6
                S[:,:2,:2]=np.outer(g,g)/ell
                S[:,:2,4]=-1j*w*g; S[:,4,:2]=1j*w*g
            else:
                S=self.sheet(w,r)
            kp+=area[:,None,None]*(P.T@S[:,:4,:4]@P)
            np.add.at(ky[:,0],allids,(area*self.shape.length)[:,None]*(S[:,:4,4]@P))
            col=np.zeros(self.size,complex)
            np.add.at(col,allids,(area*self.shape.length)[:,None]*(S[:,4,:4]@P))
            yk[0]+=col
            xx+=np.sum(area*self.shape.length**2*S[:,4,4])
        return sparse_blocks(allids,kp,self.size),ky,yk,xx

    def _solve(self,w):
        cap,col,row,xx=self._cap(w)
        A=self.k+1j*w*self.g+cap
        R=self.B0-1j*w*self.Bv-col
        # Port matrices can be retarded; sheet work must use its actual reciprocal
        # reaction row. The analytic shape source is real before the -i omega factor.
        Rt=self.B0.T+1j*w*self.Bv.T-row
        fac=splu(A[self.free][:,self.free].tocsc())
        X=np.zeros((self.size,3),complex)
        X[self.free]=fac.solve(R[self.free])
        W=-Rt@X;W[0,0]+=xx
        return W,X,R

    def static(self):
        if self._dc is None:
            if not self.B0.any():
                # A translationally neutral planar sheet has zero static X reaction.
                self._dc=np.zeros((3,3))
            else:
                # Use the same finite exterior with real coordinates for static compliance.
                real=CurvedSound(self.shape,self.exit_scale,self.n,self.box,self.pml,
                                 self.sheet,self.base,stretch=0.,advection=self.advection)
                static=real._solve(0.)[0].real
                self._dc=(static+static.T)/2
        return self._dc

    def stiffness(self,omega):
        return self._solve(omega)[0]-self.static()

    def evaluate(self,omega):
        w=float(omega)
        W,X,_=self._solve(w)
        W=W-self.static()
        power=[]; cutpower=[]
        for regs,(k,g,cut) in zip(((0,1),(2,3,4,5)),self.outs):
            op=sum(self.K[i]+1j*w*self.G[i] for i in regs)
            power.append(power_matrix(X.conj().T@(op@X),w))
            reaction=(k+1j*w*g)@X
            cutpower.append(power_matrix(X[cut].conj().T@reaction[cut],w))
        ids=self.capi
        vals=np.stack([X[ids[:,:2]].mean(axis=1),X[ids[:,2:]].mean(axis=1),
                       X[self.density[ids[:,:2]]].mean(axis=1),X[self.density[ids[:,2:]]].mean(axis=1)],axis=1)
        wc=np.full(self.radii.shape,self.shape.length)
        nminus=self._fields(1,0,self.radii,wc)[0]
        nplus=self._fields(2,0,self.radii,wc)[0]
        pminus=2*nminus[:,None]**4.5*vals[:,2]
        pplus=2*nplus[:,None]**4.5*vals[:,3]
        Xcap=np.broadcast_to(np.array([self.shape.length,0.,0.]),(self.n,3))
        if self.sheet is None:
            rel=(vals[:,1]-vals[:,0])/(WIDTH/6)
            absolute=rel-1j*w*Xcap
            rel=absolute+1j*w*nplus[:,None]*Xcap
            normal=1j*w*(vals[:,1]-vals[:,0])
        else:
            S=self.sheet(w,self.radii)
            traces=np.concatenate([vals,Xcap[:,None,:]],axis=1)
            conormal=S@traces
            absolute=conormal[:,1];rel=absolute+1j*w*nplus[:,None]*Xcap;normal=conormal[:,4]
        wall=[]
        for ids,N,rho,ur,uz,dr,dz in self.wall_data:
            phase=np.einsum('ki,kij->kj',N,X[ids]); amp=np.einsum('ki,kij->kj',N,X[self.density[ids]])
            wall.append(dict(phase=phase,thermodynamic_pressure=2*rho[:,None]**4.5*amp))
        return dict(W=W,feed_power=power[0],bulk_power=power[1],feed_cut=cutpower[0],bulk_cut=cutpower[1],
                    radii=self.radii,disk_weights=self.disk_weights,relative_current=rel,absolute_current=absolute,
                    front_force=normal,pressure_minus=pminus,pressure_plus=pplus,wall=wall,dc=self._dc)
