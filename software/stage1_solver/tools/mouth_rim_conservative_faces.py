"""Conservative normal-face reconstruction on the fixed R4 terminal strip.

The paired amplitude/phase/order solve does not impose axial material momentum
balance at its held wall. Project the native configurational stress to shared
face moments, allocating the measured weak residual to a cell source. Boundary
source extrapolation and the relative stress-fit norm are explicit controls.
The saved separate wall virtual works are held out during this reconstruction.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from mouth_bonded_front_adjoint import angular_weights
from mouth_disk_pillbox import dot, physical, stresses
from mouth_exchange_carrier import mapped_nodes
from mouth_last_ring_solid import interpolation
from mouth_rim_disk_work_review import Disk
from mouth_rim_shared_crystal import TerminalQuadrature
from mouth_supercore_joined_front import check_guard

GAUSS=(np.array([-1.,1.])/np.sqrt(3)+1)/2
CORNERS=((0,0),(1,0),(0,1),(1,1))
FACE_VERTICES=((0,2),(1,3),(0,1),(2,3))
SIGNS=np.array([-1,1,-1,1])


class Block:
    """Particle momentum through oriented faces in one material chart."""
    def __init__(self,disk,slab,side,knots):
        m=disk.model
        self.slab,self.side=slab,side
        self.x=m.rr[side]
        self.z=m.z[slab]
        z=self.z
        bottom=np.searchsorted(z,knots[0]+1e-11)-1 if slab==1 else 0
        top=np.searchsorted(z,knots[-1]-1e-11) if slab==2 else len(z)-1
        # Include adjacent cells so the enclosure is interior to numerical cuts.
        self.start=max(0,bottom-1)
        end=min(len(z)-1,top+1)
        if z[self.start]<.62*disk.sh.length-1e-9:
            self.start+=1
        self.nz=end-self.start
        self.nx=len(self.x)-1
        i,k=np.meshgrid(np.arange(self.nx),np.arange(self.start,end),indexing='ij')
        self.i,self.k=i.ravel(),k.ravel()
        self.ids=np.column_stack([m.ids[slab,side][self.i+u,0,self.k+v] for u,v in CORNERS])
        n=len(self.i)
        self.target=np.zeros((n,4))
        self.grad=np.zeros_like(self.target)
        self.lump=np.zeros_like(self.target)
        self.area=np.zeros((n,2))

    def integrate(self,disk):
        a=disk.model.a
        na=len(a)-1
        ix,ia,iz=np.meshgrid(np.arange(self.nx),np.arange(na),
                             np.arange(self.start,self.start+self.nz),indexing='ij')
        ix,ia,iz=(v.ravel() for v in (ix,ia,iz))
        dx,da,dz=np.diff(self.x)[ix],np.diff(a)[ia],np.diff(self.z)[iz]
        def angular_sum(t):
            return t.reshape(self.nx,na,self.nz).sum(axis=1).ravel()
        for u in GAUSS:
            for h in GAUSS:
                for v in GAUSS:
                    x,theta,z=self.x[ix]+u*dx,a[ia]+h*da,self.z[iz]+v*dz
                    b,s=disk.sh.radius(z)
                    width=28-b if self.side else b
                    r=b+width*x if self.side else width*x
                    slope=s*(1-x if self.side else x)
                    F,Fp=disk.reader.vortex(r*np.sin(theta))
                    vol=.125*4*np.pi*width*dx*da*dz*r*r*np.sin(theta)
                    tx,tz=np.zeros(len(x)),np.zeros(len(x))
                    for sign,interp in zip((1,-1),disk.reader.fields[self.slab,self.side]):
                        field=np.zeros((len(x),3))
                        grad=np.zeros((len(x),3,3))
                        for l,j,k in np.ndindex(2,2,2):
                            raw=interp.values[ix+l,ia+j,iz+k]
                            p,q,t=(u if l else 1-u,h if j else 1-h,v if k else 1-v)
                            field+=p*q*t*raw
                            grad[:,:,0]+=((2*l-1)*q*t/(dx*width))[:,None]*raw
                            grad[:,:,1]+=(p*(2*j-1)*t/(r*da))[:,None]*raw
                            grad[:,:,2]+=(p*q*(2*k-1)/dz)[:,None]*raw
                        grad[:,:,2]-=slope[:,None]*grad[:,:,0]
                        state=physical(1,r,theta,field.T,tuple(grad[:,j] for j in range(3)),F,Fp)
                        tr,tw=stresses(state)
                        gf,_,_,n,swirl,_=state
                        gv=grad[:,0]
                        y=field[:,0]**2
                        delta=(.5*(F**2*dot(gv,gv)-dot(gf,gf)-n*swirl)
                               +F**10*(y-1)**2*(y**3+2*y*y+3*y+4)/20-(n**5-5*n+4)/20)
                        tx+=sign*(tr.sum(axis=0)+F**2*gv[:,0]*gv[:,2]-gf[:,0]*gf[:,2])
                        tz+=sign*(tw.sum(axis=0)+F**2*gv[:,2]**2-gf[:,2]**2-delta)
                    ex=angular_sum(vol*(tx-slope*tz)/(dx*width))
                    ez=angular_sum(vol*tz/dz)
                    self.target+=np.column_stack([(1-u)*ex,u*ex,(1-v)*ez,v*ez])
                    self.area+=np.column_stack([angular_sum(vol/(dx*width)),angular_sum(vol/dz)])
                    for j,(l,k) in enumerate(CORNERS):
                        p,q=(u if l else 1-u,v if k else 1-v)
                        self.grad[:,j]+=(2*l-1)*q*ex+(2*k-1)*p*ez
                        self.lump[:,j]+=angular_sum(vol)*p*q
        if self.slab==1:
            self.native_bank(disk)

    def native_bank(self,disk):
        a=disk.model.a
        j,k=np.meshgrid(np.arange(len(a)-1),np.arange(self.start,self.start+self.nz),indexing='ij')
        dz=np.diff(self.z)[k]; da=np.diff(a)[j]
        i=0 if self.side else self.nx-1
        row=np.arange(i*self.nz,(i+1)*self.nz)
        x=0. if self.side else 1.
        for h in GAUSS:
            for v in GAUSS:
                t=a[j]+h*da; z=self.z[k]+v*dz
                b,_=disk.sh.radius(z)
                F,Fp=disk.reader.vortex(b*np.sin(t))
                pts=np.column_stack((np.full(z.size,x),t.ravel(),z.ravel()))
                amps=[f(pts)[:,0].reshape(z.shape) for f in disk.reader.fields[1,self.side]]
                # Surface configurational term is minus the bank energy derivative.
                ez=(-.125*(1-2*self.side)*4*np.pi*b*b*da*np.sin(t)**2*
                    F*Fp*(amps[0]**2-amps[1]**2)).sum(axis=0)
                self.target[row,2]+=(1-v)*ez;self.target[row,3]+=v*ez
                j0=0 if self.side else 1
                self.grad[row,j0]-=ez;self.grad[row,j0+2]+=ez

    def solid(self,disk,state,wall,fields):
        with np.load(fields) as f:
            chi=np.einsum('raj,a->jr',(f['wound'][2]-f['rest'][2])[disk.model.ids[1,0]],
                          angular_weights(disk.model.a))
        rule=TerminalQuadrature(state,wall,self.x,self.z,64)
        for bl in rule.blocks():
            X,Z=bl['x'],bl['z']
            wt=bl['measure']*interpolation(self.z,self.x,chi,Z,X)[0]
            if not wt.any():continue
            J,K=bl['J'],bl['K']
            vr=-wt*K[:,0,1]
            vz=-wt*(K[:,1,1]-1+.5*J*(J-1))
            i=np.clip(self.x.searchsorted(X)-1,0,self.nx-1)
            k=np.clip(self.z.searchsorted(Z)-1,0,len(self.z)-2)
            sel=(k>=self.start)&(k<self.start+self.nz)
            b,s=disk.sh.radius(Z)
            dx=np.diff(self.x)[i];dz=np.diff(self.z)[k]
            u=(X-self.x[i])/dx;v=(Z-self.z[k])/dz
            ex=(vr-s*X*vz)/(b*dx);ez=vz/dz
            row=(i*self.nz+k-self.start)[sel]
            def add(matrix,j,values):
                matrix[:,j]+=np.bincount(row,values[sel],minlength=matrix.shape[0])
            for j,t in enumerate(((1-u)*ex,u*ex,(1-v)*ez,v*ez)):
                add(self.target,j,t)
            for j,(l,h) in enumerate(CORNERS):
                p,q=(u if l else 1-u,v if h else 1-v)
                add(self.grad,j,(2*l-1)*q*ex+(2*h-1)*p*ez)


def readrows(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def enclosure_weight(z,w):
    return np.clip((z-w[0])/(w[1]-w[0]),0,1)*np.clip((w[3]-z)/(w[3]-w[2]),0,1)


def results(a):
    disk=Disk(a.state,a.wall,a.fields)
    data=readrows(a.fields.parent/'stdout')
    bankrows=readrows(a.banks/'stdout')
    check_guard(a.banks/'stdout')
    provenance=next(row for row in bankrows if row.get('kind')=='source')
    if any(provenance[key]!=value for key,value in disk.source_hash.items()):
        raise ValueError('Bank measurements refer to other held fields')
    bankrows=[r for r in bankrows if r.get('kind')=='bank_test' and r['quadrature']==2]
    labels={(r['bank'],r['strip'],r['radial_cells']) for r in bankrows}
    expected={(s,h,j) for s in ('duct','outside') for h in ('near','middle','upstream') for j in (4,8)}
    if labels!=expected or len(bankrows)!=len(expected):
        raise ValueError('Incomplete held-out wall probes')
    zlo=min(disk.windows[1][0],min(r['depth_knots'][0] for r in bankrows))
    blocks=[Block(disk,s,b,[zlo,disk.windows[1][-1]]) for s in (1,2) for b in (0,1)]
    # Assembly runs before reading any target force values.
    for bl in blocks:bl.integrate(disk)
    blocks[0].solid(disk,a.state,a.wall,a.fields)
    def gather(name):return np.concatenate([getattr(bl,name) for bl in blocks])
    vertices=np.unique(gather('ids'))
    nodes=np.searchsorted(vertices,gather('ids'))
    nh,nv=len(nodes),len(vertices)
    edge=np.stack([np.sort(nodes[:,j],axis=1) for j in FACE_VERTICES],axis=1)
    unique,finv=np.unique(edge.reshape(-1,2),axis=0,return_inverse=True)
    faces=finv.reshape(-1,4); nf=len(unique)
    counts=np.bincount(faces.ravel(),minlength=nf)
    if (counts>2).any():raise ValueError('Non-manifold topology')
    interior_counts=np.bincount(faces.ravel(),np.tile(SIGNS,nh),minlength=nf)
    if np.max(abs(interior_counts[counts==2])):raise ValueError('Face orientation mismatch')
    for bl,f in zip(blocks,np.split(faces,np.cumsum([len(b.i) for b in blocks])[:-1])):bl.faces=f
    axis=np.concatenate([bl.faces[bl.i==0,0] for bl in blocks if bl.side==0])
    free=np.ones(nf,dtype=bool);free[axis]=False
    boundary=np.zeros(nv,dtype=bool);boundary[unique[counts==1].ravel()]=True
    lump=gather('lump');grad=gather('grad');volume=lump.sum(axis=1)
    lump_v=np.bincount(nodes.ravel(),lump.ravel(),minlength=nv)
    weak=np.bincount(nodes.ravel(),grad.ravel(),minlength=nv)
    interior=~boundary
    source0=-weak/lump_v
    density=source0.copy();density[boundary]=0
    filled=interior.copy(); extended=density.copy()
    while not filled.all():
        i,j=unique.T
        i,j=np.r_[i,j],np.r_[j,i]
        use=filled[i]&~filled[j]
        weights=np.bincount(j[use],minlength=nv)
        new=weights>0
        if not new.any():raise ValueError('No internal source near free boundary')
        extended[new]=np.bincount(j[use],extended[i[use]],minlength=nv)[new]/weights[new]
        filled|=new
    B=sparse.coo_matrix((np.tile(SIGNS,nh),(np.repeat(np.arange(nh),4),faces.ravel())),shape=(nh,nf)).tocsr()[:,free]
    target=gather('target'); area=gather('area')
    cut=mapped_nodes(13,disk.sh.length,28.,22.)[0][1][4:6]
    output=[]
    for name,den,ratio in [('volume_only',density,1.),('extrapolated',extended,1.),('axial_weight_4',extended,4.)]:
        source=(lump*den[nodes]).sum(axis=1)
        wt=volume[:,None]/area**2;wt[:,1]*=ratio
        weights=np.repeat(wt,2,axis=1)
        partner=faces[:,[1,0,3,2]]
        H=sparse.coo_matrix((np.r_[weights.ravel()/3,weights.ravel()/6],
                (np.tile(faces.ravel(),2),np.r_[faces.ravel(),partner.ravel()])),shape=(nf,nf)).tocsc()[free,:][:,free]
        l=np.bincount(faces.ravel(),(weights*target).ravel(),minlength=nf)
        A=sparse.bmat([[H,B.T],[B,None]],format='csc')
        ans=spsolve(A,np.r_[l[free],source]);q=np.zeros(nf);q[free]=ans[:free.sum()]
        balance=np.max(abs(B@q[free]-source))
        if not np.isfinite(balance) or balance>1e-7:raise ArithmeticError(f'Unresolved balance {balance}')
        probes=[]
        for r in bankrows:
            bl=blocks[int(r['bank']=='outside')]
            side=bl.side; j=0 if side else 1
            mask=(bl.i==(0 if side else bl.nx-1))
            za,zb=bl.z[bl.k[mask]],bl.z[bl.k[mask]+1]
            lo,mid,hi=r['depth_knots']
            def hat(z):return np.maximum(0,np.minimum((z-lo)/(mid-lo),(hi-z)/(hi-mid)))
            react=SIGNS[j]*q[bl.faces[mask,j]]@(.5*(hat(za)+hat(zb)))
            probes.append(dict(strip=r['strip'],bank=r['bank'],radial_cells=r['radial_cells'],
                trace=react,target=r['total_native'],difference=react-r['total_native']))
        surfaces=[]
        for band,w in enumerate(disk.windows,2):
            src=0.; wall=0.; exterior=0.
            for bl,ss in zip(blocks,np.split(source,np.cumsum([len(b.i) for b in blocks])[:-1])):
                i,k=bl.i,bl.k;xa=bl.x[i];xb=bl.x[i+1]
                za,zb=bl.z[k],bl.z[k+1];al,ah=enclosure_weight(za,w),enclosure_weight(zb,w)
                def radial(t):return np.clip((cut[1]-t)/np.diff(cut)[0],0,1) if bl.side else np.ones_like(t)
                rl,rh=radial(xa),radial(xb)
                src+=ss@(.25*(rl+rh)*(al+ah))
                fl=q[bl.faces]
                exterior-=np.sum((fl[:,0]+fl[:,1])*.25*(rh-rl)*(al+ah)+(fl[:,2]+fl[:,3])*.25*(ah-al)*(rh+rl))
                if bl.slab==1:
                    j=0 if bl.side else 1
                    sel=i==(0 if bl.side else bl.nx-1)
                    wall+=SIGNS[j]*fl[sel,j]@(.5*(al[sel]+ah[sel]))
            orig=next(r for r in data if r.get('kind')=='closed_work' and r['axial_band']==band and r['lateral_band']==4)
            defect=wall+exterior-src
            if abs(defect)>2e-7:raise ArithmeticError(f'Enclosure identity {defect}')
            surfaces.append(dict(axial_band=band,wall_recovered=wall,exterior_recovered=exterior,
               recovered_total=wall+exterior,source=src,identity_error=defect,
               wall_plus_archived_exterior=wall+orig['outward_momentum'],archived_exterior=orig['outward_momentum']))
        output.append(dict(closure=name,cell_balance_max=balance,probes=probes,surfaces=surfaces))
    return dict(level=disk.level,hashes=disk.source_hash,grid_cells=nh,faces=nf,
                formulation='Mapped angle-integrated RT0; two-point Q1 native stress, vortex-bank work, matched crystal quadrature; weak coordinate residual lumped to control volumes; zero normal momentum at symmetry axis.',
                precision='Bank comparison uses saved second-order gap-matched native virtual work.',results=output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('state','wall','fields','banks','output'):p.add_argument(arg,type=Path)
    a=p.parse_args()
    with a.output.open('x') as dest:json.dump(results(a),dest,indent=2,allow_nan=False)
