"""Outgoing torsional spectrum of the saved, strained bonded openings."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh, splu

from mouth_exchange_shear import OutgoingSlab, block_matrices
from mouth_exit_witness import inputs
from mouth_live_waist_wave import transfer
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian
from mouth_shallow_puncture_fields import PORT
from mouth_strained_wave import cell_tangents, core_factors
from numpy.polynomial.legendre import leggauss, legval


def form(ids, values, size):
    indices = np.broadcast_arrays(ids[:, :, None], ids[:, None, :])
    return sparse.coo_matrix((values.ravel(), tuple(a.ravel() for a in indices)),
                             shape=(size, size)).tocsc()


def channel_port(points, degree):
    port = OutgoingSlab(points, PORT)
    port.cutoff[0] = 0.   # Exact constant through-thickness eigenfunction.
    port.coeff = np.array([1j**n * math.factorial(degree+n) /
                          (math.factorial(n)*math.factorial(degree-n)*2**n)
                          for n in range(degree+1)])
    return port


class Sector:
    def __init__(self, trial, shifts, order, degree):
        self.degree = degree
        shape = self.shape = Meridian(trial.radius, trial.x+shifts[:4])
        mesh = trial.mesh
        self.coeff = trial.co.copy()
        self.coeff[0] += shifts[4]
        self.reference = trial.core_reference
        self.mesh = mesh
        nt = int(np.ceil(2.77*trial.radius*mesh))
        side = round(mesh*(PORT-1))
        if PORT <= trial.radius:
            raise ValueError('The outgoing port must enclose the aperture')
        x = np.linspace(0.,1.,mesh+1)
        z = np.linspace(0.,1.,nt+1)
        feed = np.linspace(-1.,0.,mesh+1)
        def slab(r,w):
            return r,np.zeros_like(w)
        def duct(x,w):
            radius,slope = shape.radius(w)
            return x*radius,x*slope
        chunks = [block_matrices(nodes,feed,slab,0.,shape.length,degree) for nodes in
                  (trial.radius*x,np.linspace(trial.radius,PORT,side+1))]
        chunks.append(block_matrices(x,shape.length*z,duct,0.,shape.length,degree))
        cloud = np.round(np.vstack([piece[0] for piece in chunks]),10)
        _,first,remap = np.unique(cloud,axis=0,return_index=True,return_inverse=True)
        perm = np.argsort(first)
        index = np.argsort(perm)[remap]
        allpoints = cloud[first[perm]]
        size = len(allpoints)
        K = sparse.csc_matrix((size,size))
        M = K.copy()
        offset = 0
        for pts,cells,k,m in chunks[:-1]:
            ids = index[offset+cells]
            K += form(ids,k,size)
            M += form(ids,m,size)
            offset += len(pts)
        _,cells,k,m = chunks[-1]
        ids = index[offset+cells]
        dummy = SimpleNamespace(mesh=mesh,ids=ids,kc=k,mc=m)
        k,m = cell_tangents(dummy,shape,self.coeff,degree=degree,core_reference=self.reference)
        weights = sparse.kron(transfer(trial.n,x),transfer(2*trial.n,z))@order.T.ravel()
        inside = form(ids,weights[:,None,None]*m,size)
        K += form(ids,weights[:,None,None]*k,size)
        M += inside
        # The unknown is angular displacement divided by r. A regular degree-1
        # rotation has finite amplitude at the axis; higher degrees vanish there.
        live = (M.diagonal()>0) & ((allpoints[:,0]>0) | (degree==1))
        self.allpoints = allpoints
        self.points = allpoints[live]
        self.live = live
        self.ids,self.mc,self.weights = ids,m,weights
        self.K,self.M,self.inside = [v[live][:,live].tocsc() for v in (K,M,inside)]
        self.M.eliminate_zeros()
        self.port = channel_port(self.points,degree)
        ncomp,self.labels = connected_components(self.M,directed=False)
        self.connected = np.unique(self.labels[self.port.ids])
        if self.connected.size!=1:
            raise ValueError('Broken feeding slab')
        self.ncomp = ncomp

    def profile(self,v):
        full = np.zeros(len(self.live),complex)
        full[self.live] = v
        local = full[self.ids]
        cell = self.weights*np.einsum('ci,cij,cj->c',local.conj(),self.mc,local).real
        mass = np.vdot(v,self.M@v).real
        axial = cell.reshape(self.mesh,-1).sum(axis=0)
        nt=len(axial)
        t=(np.arange(nt)+.5)/nt
        widths = self.shape.profile(np.linspace(0.,1.,2001))[0]
        waist = float(np.argmin(widths)/2000)
        # Exact collar/core split with the same cut-cell quadrature as cell_tangents.
        order=self.weights.reshape(self.mesh,nt)
        vec=local.reshape(self.mesh,nt,4)
        x0=np.arange(self.mesh)[:,None]/self.mesh
        x1=x0+1/self.mesh
        core=0.
        nodes,wg=leggauss(3)
        for b,wb in zip((1+nodes)/2,wg/2):
            at=(np.arange(nt)[None,:]+b)/nt
            width=self.shape.profile(at)[0]
            p=np.exp(legval(2*at-1,self.coeff))
            hi=np.clip(np.exp(-p/2),x0,x1)
            stretch,_,axial_stretch=core_factors(self.shape,self.coeff,self.reference,at)
            for a,wa in zip((1+nodes)/2,wg/2):
                x=x0+a*(hi-x0)
                u=(x-x0)*self.mesh
                nodal = ((1-u)*(1-b)*vec[:,:,0]+u*(1-b)*vec[:,:,1]
                         +(1-u)*b*vec[:,:,2]+u*b*vec[:,:,3])
                core += np.sum(order*abs(nodal)**2*wa*wb*(hi-x0)*width*(x*width)**4*
                               self.shape.length/nt/(stretch**3*axial_stretch))
        if abs(sum(cell)-np.vdot(v,self.inside@v).real)>1e-9*mass:
            raise ArithmeticError('Kinetic mass partition failed')
        return dict(cavity_fraction=sum(cell)/mass,collar_fraction=(sum(cell)-core)/mass,
                    cavity_mean_depth_fraction=axial@t/axial.sum(),waist_depth_fraction=waist,
                    cavity_waist_band_fraction=axial[np.abs(t-waist)<=.125].sum()/axial.sum(),
                    cavity_exit_quarter_fraction=axial[t>=.75].sum()/axial.sum())


def closed_seeds(op,count):
    """First clamped-aperture levels as starts. These seeds do not establish a bound state."""
    cavity=(op.points[:,1]>0) & (op.labels==op.connected[0])
    K,M = [a[cavity][:,cavity] for a in (op.K,op.M)]
    val,vec=eigsh(K,k=count,M=M,sigma=0.,which='LM',v0=np.ones(K.shape[0]))
    out=[]
    for i in np.argsort(val):
        v=np.zeros(op.K.shape[0],complex)
        v[cavity]=vec[:,i]
        out.append((complex(np.sqrt(val[i]),-1e-3),v))
    return out


def pole(op,seed):
    om,vector=seed
    # Isolate the occupied connected piece. On a detached component the direct
    # Hermitian eigenproblem supplies the exact zero outgoing field.
    weights=np.array([np.vdot(vector[op.labels==i],(op.M@vector)[op.labels==i]).real
                      for i in range(op.ncomp)])
    comp=int(np.argmax(weights))
    if weights[comp]<.98*sum(weights):
        raise ArithmeticError('No unique occupied shear component')
    use=op.labels==comp
    k,m=op.K[use][:,use],op.M[use][:,use]
    vv=vector[use].copy()
    vv/=np.sqrt(np.vdot(vv,m@vv).real)
    scale=sparse.diags(1/np.sqrt(m.diagonal())).tocsc()
    k,m=(scale@k@scale).tocsc(),(scale@m@scale).tocsc()
    u=vv/scale.diagonal()
    is_open=comp==op.connected[0]
    if not is_open:
        val,q=eigsh(k,k=1,M=m,sigma=om.real**2,which='LM',v0=u.real,tol=1e-12)
        om=complex(np.sqrt(val[0]))
        u=q[:,0].astype(complex)
    else:
        ids=np.flatnonzero(use)
        edge=np.searchsorted(ids,op.port.ids)
        ii,jj=np.repeat(edge,len(edge)),np.tile(edge,len(edge))
        factor=scale.diagonal()[edge]
        dual=(m@u).conj()
        for it in range(45):
            block=op.port.block(om)*np.outer(factor,factor)
            A=k-om**2*m-sparse.coo_matrix((block.ravel(),(ii,jj)),shape=k.shape).tocsc()
            res=A@u
            err=np.linalg.norm(res)/max(np.linalg.norm(k@u),1e-15)
            if err<5e-12 and abs(dual@u-1)<2e-10:
                break
            dB=(op.port.block(om+2e-6)-op.port.block(om-2e-6))/4e-6
            d=(-2*om)*(m@u)
            d[edge]-=factor*(dB@(factor*u[edge]))
            bordered=sparse.bmat([[A,sparse.csc_matrix(d[:,None])],
                                  [sparse.csc_matrix(dual[None,:]),None]],format='csc')
            dx=splu(bordered).solve(-np.r_[res,dual@u-1])
            u+=dx[:-1]
            om+=dx[-1]
            if not np.isfinite(om) or om.real<=0:
                raise ArithmeticError('Left positive-frequency outgoing sheet')
        else:
            raise ArithmeticError(f'Unresolved outgoing pole, last residual {err:.3g}')
    result=np.zeros_like(vector,dtype=complex)
    result[use]=scale@u
    result/=np.sqrt(np.vdot(result,op.M@result).real)
    edge=result[op.port.ids]
    rv=(op.K-om**2*op.M)@result
    rv[op.port.ids]-=op.port.block(om)@edge
    residual=np.linalg.norm(rv)/np.linalg.norm(op.K@result)
    if residual>1e-7 or om.imag>1e-10 or not np.isfinite(residual):
        raise ArithmeticError('Physical pole residual or outgoing sign failed')
    channels=op.port.cutoff<om.real
    phi=op.port.project@edge
    flux=.5*om.real*PORT**4*abs(phi[channels])**2*op.port.traction(om.real)[channels].imag
    E=.25*(abs(om)**2+np.vdot(result,op.K@result).real)
    phase_overlap=abs(np.vdot(vector,op.M@result))**2/np.vdot(vector,op.M@vector).real
    loss=-2*om.imag
    stats=dict(frequency=om.real,loss=loss,Q=(om.real/loss if loss>0 else None),
               power_per_energy=flux.sum()/E,channel_power_per_energy=(flux/E).tolist(),
               depth_thresholds=op.port.cutoff[:4].tolist(),residual=residual,
               seed_overlap=phase_overlap,detached=not is_open,**op.profile(result))
    return (om,result),stats


def run(radius,n,degree,detached):
    tr=Trial(radius,n,bonded_core=True)
    base=ROOT/'_scratch/s11c/mouth_exit_release01'/f'r{radius}n{n}'/'endpoint.npz'
    sh,c,seed,digest,status,*_=inputs(base,radius,n,tr)
    if status!='balanced':
        raise ValueError('Unvalidated reference shape')
    if detached:
        path=ROOT/'_scratch/s11c/mouth_exit_material01'/f'r{radius}n{n}'/'endpoint.npz'
        meta=json.loads((path.parent/'outcome.json').read_bytes())
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if (meta['exitCode'] or meta['childOutcome']['exitCode'] or not meta['limitsVerified']
                or digest not in (path.parent/'stdout').read_text()):
            raise ValueError('Unvalidated detached state')
        with np.load(path) as f:
            if f['parent_sha256'].item()!=hashlib.sha256(base.read_bytes()).hexdigest():
                raise ValueError('Unmatched preparation')
            sh,c,seed=f['shifts'],f['chi'],(complex(f['pole'].item()),f['vector'])
    offsets=(-.01,-.005,0.,.005,.01)
    ops={}
    for h in offsets:
        s=sh.copy()
        s[1]+=h
        ops[h]=Sector(tr,s,c,degree)
    op=ops[0.]
    if not all(np.array_equal(other.live,op.live) for other in ops.values()):
        raise ArithmeticError('A geometric probe changed active material labels')
    emit(kind='setup',radius=radius,mesh=n,wave_mesh=tr.mesh,degree=degree,sha256=digest,
         components=op.ncomp,exit=op.shape.exit,depth_thresholds=op.port.cutoff[:4].tolist(),
         disconnected_challenge=detached)
    seeds=[] if detached else closed_seeds(op,3)
    if degree==4:
        from mouth_live_waist_wave import OccupiedWave
        archived=OccupiedWave(op.shape,n,2*n,tr.mesh,PORT)
        ref=np.zeros(len(op.live),complex)
        ref[op.allpoints[:,0]>0]=seed[1]
        if not np.array_equal(archived.uniform[0],op.allpoints[op.allpoints[:,0]>0]):
            raise ValueError('Archived material node numbering differs')
        seeds.insert(0,(seed[0],ref[op.live]))
    for i,start in enumerate(seeds):
        label=('archived' if degree==4 and i==0 else f'clamped_{i+(degree!=4)}')
        try:
            middle,stats=pole(op,start)
        except (ValueError,ArithmeticError,RuntimeError) as e:
            emit(kind='unresolved_start',label=label,start=start[0].real,reason=str(e))
            continue
        if label=='archived' and (abs(middle[0]-seed[0])>1e-6 or stats['seed_overlap']<.9999):
            raise ArithmeticError('Degree-four sector failed archived pole calibration')
        emit(kind='mode',label=label,start=start[0].real,**stats)
        ends={}
        for h in offsets:
            if h==0:
                continue
            try:
                trial,data=pole(ops[h],middle)
                if data['seed_overlap']<.96:
                    raise ArithmeticError('Mode changed on exit variation')
                ends[h]=trial[0]
                emit(kind='signed',label=label,log_exit_step=h,exit=ops[h].shape.exit,**data)
            except (ArithmeticError,ValueError,RuntimeError) as e:
                emit(kind='unresolved_signed',label=label,log_exit_step=h,reason=str(e))
        for h in (.01,.005):
            if h in ends and -h in ends:
                den=ops[h].shape.exit-ops[-h].shape.exit
                emit(kind='exit_work',label=label,step=h,domega_dexit=(ends[h]-ends[-h]).real/den,
                     dloss_dexit=-2*(ends[h]-ends[-h]).imag/den)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius',choices=(1,2),type=int,required=True)
    p.add_argument('--mesh',choices=(48,72),type=int,required=True)
    p.add_argument('--degree',choices=range(1,7),type=int,required=True)
    p.add_argument('--detached',action='store_true')
    args=p.parse_args()
    run(args.radius,args.mesh,args.degree,args.detached)
