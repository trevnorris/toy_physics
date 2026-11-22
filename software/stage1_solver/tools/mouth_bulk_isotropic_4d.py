"""Compression and shear of a neutral four-dimensional bulk with isotropic pair energies."""
import argparse
import ctypes as ct
import itertools
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.linalg import eigh, helmert
from scipy.optimize import minimize

from neutral_patch_sites import TEMPERATURE, WITHIN, emit, morse


def library(work):
    work.mkdir(parents=True,exist_ok=True)
    so = (work/'bulk.so').resolve()
    subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(Path(__file__).with_suffix('.cpp')),'-o',str(so)],check=True)
    lib = ct.CDLL(str(so))
    array = np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    lib.make_cell.argtypes = [ct.c_int,array,ct.c_double,array]
    lib.make_cell.restype = ct.c_void_p
    lib.free_cell.argtypes = [ct.c_void_p]
    lib.free_cell.restype = None
    lib.evaluate.argtypes = [ct.c_void_p,array,array,array]
    lib.evaluate.restype = ct.c_double
    lib.integrate.argtypes = [ct.c_void_p,ct.c_int,ct.c_double,array,array]
    lib.integrate.restype = ct.c_int
    return lib


class Cell:
    def __init__(self,lib,box,n,rc):
        self.lib,self.box,self.n,self.rc = lib,np.ascontiguousarray(box),n,rc
        self.handle = lib.make_cell(n,self.box,rc,np.array(WITHIN))
        self.volume = self.box.prod()

    def __del__(self):
        if hasattr(self,'handle'):
            self.lib.free_cell(self.handle)

    def read(self,pos):
        pos = np.ascontiguousarray(pos)
        gradient,diag = np.empty_like(pos),np.empty(4)
        energy = self.lib.evaluate(self.handle,pos,gradient,diag)
        return energy,gradient,diag

    def settle(self,pos):
        def cost(flat):
            e,g,_ = self.read(flat)
            return e,g
        opt = minimize(cost,pos.ravel(),jac=True,method='L-BFGS-B',
                       options=dict(maxiter=10000,maxls=60,ftol=2.e-16,gtol=1e-9))
        pos = opt.x.reshape(-1,4)
        # Remove optimizer energy-tolerance effects from tiny signed-force probes.
        for k in range(4):
            e,g,_ = self.read(pos)
            if abs(g).max()<2.e-9:
                break
            h = derivatives(self,pos)[0]
            ev,vec = eigh(h)
            if ev[4]<1.e-5:
                break
            change = vec[:,4:]@((vec[:,4:].T@g.ravel())/ev[4:])
            if np.max(abs(change))>.004:
                break
            pos = pos-change.reshape(-1,4)
        err = np.max(abs(self.read(pos)[1]))
        if err>2.e-7:
            raise ArithmeticError(f'Unresolved equilibrium forces {err}, {opt.message}')
        return pos


def initial(box,seed):
    lattice = np.indices((4,4,4,4)).reshape(4,-1).T
    lattice = lattice[lattice.sum(axis=1)%2==0]*box/4
    if seed==0:
        return lattice
    rng = np.random.default_rng(seed)
    pos=[]
    while len(pos)<len(lattice):
        point = rng.uniform(size=4)*box
        d = np.array(pos).reshape(-1,4)-point
        d -= box*np.rint(d/box)
        if np.all(np.sum(d*d,axis=1)>.60**2):
            pos.append(point)
    return np.array(pos)


def images(cell,pos):
    # Independently enumerate unordered pairs and physical images in Python.
    i,j = np.triu_indices(cell.n,1)
    ids = np.arange(cell.n)
    raw = pos[j]-pos[i]
    raw -= cell.box*np.rint(raw/cell.box)
    extent = np.ceil(cell.rc/cell.box+.5).astype(int)
    all_i,all_j,all_d = [],[],[]
    for tile in itertools.product(*(range(-n,n+1) for n in extent)):
        shift = np.array(tile)*cell.box
        if np.linalg.norm(np.maximum(abs(shift)-.5*cell.box,0.))>cell.rc:
            continue
        d = raw+shift
        mask = np.einsum('ij,ij->i',d,d)<cell.rc**2
        all_i.append(i[mask]); all_j.append(j[mask]); all_d.append(d[mask])
        if tile>(0,0,0,0) and shift@shift<cell.rc**2:
            all_i.append(ids); all_j.append(ids); all_d.append(np.tile(shift,(cell.n,1)))
    return np.concatenate(all_i),np.concatenate(all_j),np.vstack(all_d)


def radial(r,rc):
    depth,alpha,length = WITHIN
    u,du = morse(r,*WITHIN)
    w = np.exp(-alpha*(r-length))
    ddu = 2*alpha**2*depth*w*(2*w-1)
    t = np.clip((rc-r)/.25,0.,1.)
    s = t**3*(10+t*(-15+6*t))
    ds = -120*t*t*(1-t)**2
    dds = 960*t*(1-t)*(1-2*t)
    return u*s,du*s+u*ds,ddu*s+2*du*ds+u*dds


def strain_basis():
    # q=log(V/V0); every shear generator has tr(A²)=1/2.
    basis = [np.eye(4)/4]
    basis.extend(np.diag(a)/np.sqrt(2) for a in helmert(4))
    for i,j in itertools.combinations(range(4),2):
        a = np.zeros((4,4)); a[i,j]=a[j,i]=.5
        basis.append(a)
    return np.array(basis)


def derivatives(cell,pos):
    i,j,d = images(cell,pos)
    r = np.linalg.norm(d,axis=1)
    e,du,ddu = radial(r,cell.rc)
    g = du/r
    k = (ddu-g)/r**2
    grad, native = np.zeros_like(pos),cell.read(pos)
    for a in range(4):
        grad[:,a] = np.bincount(j,g*d[:,a],cell.n)-np.bincount(i,g*d[:,a],cell.n)
    np.testing.assert_allclose([e.sum(),*np.sum(g[:,None]*d*d,axis=0)],[native[0],*native[2]],atol=2e-8)
    np.testing.assert_allclose(grad,native[1],atol=1e-10)
    h = np.empty((cell.n,4,cell.n,4))
    for a,b in itertools.product(range(4),repeat=2):
        v = k*d[:,a]*d[:,b]+(a==b)*g
        edge = np.bincount(i*cell.n+j,v,cell.n**2).reshape(cell.n,cell.n)
        incidence = -edge-edge.T
        incidence[np.diag_indices(cell.n)] += np.bincount(i,v,cell.n)+np.bincount(j,v,cell.n)
        h[:,a,:,b] = incidence
    modes = strain_basis()
    change = d@modes
    dot = np.einsum('qpa,pa->qp',change,d)
    born = (dot*k)@dot.T + 2*np.einsum('qpa,spa,p->qs',change,change,g,optimize=True)
    bmat = np.zeros((cell.n,4,len(modes)))
    for index,v in enumerate(change):
        v = g[:,None]*v+(k*dot[index])[:,None]*d
        for a in range(4):
            bmat[:,a,index] = np.bincount(j,v[:,a],cell.n)-np.bincount(i,v[:,a],cell.n)
    stress = d.T@(g[:,None]*d)/cell.volume
    return h.reshape(4*cell.n,-1),bmat.reshape(4*cell.n,-1),born,stress


def structure(pos,box):
    d = pos[:,None,:]-pos[None,:,:]
    d -= box*np.rint(d/box)
    r = np.linalg.norm(d,axis=2)
    bonds = np.triu((r<1.15)&(r>0),1)
    return bonds,r


def endpoint(args,lib):
    n=128
    box=np.full(4,(n/args.density)**.25)
    cell=Cell(lib,box,n,args.cutoff)
    start=initial(box,args.seed)
    pos=cell.settle(start)
    h,b,born,stress=derivatives(cell,pos)
    ev,v=eigh(h)
    emit(test='equilibrium',seed=args.seed,density=n/cell.volume,cutoff=cell.rc,sites=n,
         max_force=float(abs(cell.read(pos)[1]).max()),energy_per_site=cell.read(pos)[0]/n,
         soft_modes=ev[:7].tolist(),highest_mode=ev[-1],stress=stress.tolist(),pressure=-np.trace(stress)/4,
         initial_neighbors=int(structure(start,box)[0].sum()),final_neighbors=int(structure(pos,box)[0].sum()))
    assert max(abs(ev[:4]))<1e-7 and ev[4]>1e-7
    modes=v[:,4:].T@b
    np.testing.assert_allclose(v[:,:4].T@b,0.,atol=1e-7)
    pressure=-np.trace(stress)/4
    born[0,0]+=pressure*cell.volume
    for t in [0.,.01,.1,1.,args.transit,None]:
        weights = (1 if t is None else -np.expm1(-t*ev[4:]))/ev[4:]
        tangent = (born-(modes.T*weights)@modes)/cell.volume
        emit(test='response',time=t,bulk=tangent[0,0],shear_spectrum=np.linalg.eigvalsh(tangent[1:,1:]).tolist(),
             first_shear=tangent[1,1],bulk_shear_coupling=tangent[0,1:].tolist(),
             least_full_eigenvalue=np.linalg.eigvalsh(tangent)[0])
    for index in (0,1):
        a=np.diag(strain_basis()[index])
        delta=2e-5
        signed=[]
        for s in (delta,-delta):
            stretch=np.exp(s*a)
            tilted=Cell(lib,box*stretch,n,cell.rc)
            q=tilted.settle(pos*stretch)
            signed.append(tilted.read(q)[2]@a/tilted.volume)
        fd=(signed[0]-signed[1])/(2*delta)
        emit(test='signed_relaxation',mode=index,step=delta,slope=fd,analytic=tangent[index,index])
        np.testing.assert_allclose(fd,tangent[index,index],rtol=2e-5,atol=1e-5)
    np.savez(args.work/'state.npz',pos=pos,box=box,seed=args.seed,rc=cell.rc,highest=ev[-1],density=n/cell.volume)


def warm(args,lib):
    source=args.warm_from
    valid=json.loads((source.parent.parent/'outcome.json').read_text())
    assert valid['limitsVerified'] and valid['exitCode']==valid['childOutcome']['exitCode']==0
    with np.load(source) as archive:
        pos,box=archive['pos'].copy(),archive['box'].copy()
        cutoff=float(archive['rc']); seed=int(archive['seed']); highest=float(archive['highest'])
    n=len(pos)
    rng=np.random.default_rng(1701+seed)
    cell=Cell(lib,box,n,cutoff)
    nread=48
    steps=int(np.ceil(args.transit/(nread*.001)))
    dt=args.transit/(nread*steps)
    assert highest*dt<1.5
    def blocks():
        for begin in range(0,steps,128):
            count=min(128,steps-begin)
            noise=rng.normal(size=(count,n,4))*np.sqrt(2*TEMPERATURE*dt)
            yield count,noise
    def advance(mat,p,count,noise):
        if lib.integrate(mat.handle,count,dt,noise,p):
            raise ArithmeticError('Warm integration produced an invalid increment')
    # Each drive branches from one independently thermalized copy of its cold basin.
    for _ in range(nread):
        for count,noise in blocks():
            advance(cell,pos,count,noise)
    a=np.diag(strain_basis()[args.mode])
    amps=[0.,1e-4,-1e-4]
    if args.half:
        amps.extend([5e-5,-5e-5])
    states=[]
    for amplitude in amps:
        scale=np.exp(amplitude*a)
        states.append((Cell(lib,box*scale,n,cutoff),pos*scale))
    bonds,_=structure(pos,box)
    initial_bonds=bonds.copy()
    start=pos.copy()
    emit(test='warm_setup',seed=seed,mode=args.mode,n=n,temperature=TEMPERATURE,dt=dt,density=n/cell.volume,
         cutoff=cutoff,transit=args.transit,amplitudes=amps)
    for j in range(nread+1):
        if j:
            for count,noise in blocks():
                for mat,p in states:
                    advance(mat,p,count,noise)
        reading=[]
        for mat,p in states:
            # Free-energy derivative at fixed T. Ideal pressure is isotropic.
            reading.append(mat.read(p)[2]@a/mat.volume-n*TEMPERATURE*np.sum(a)/mat.volume)
        p=states[0][1]
        fresh,r=structure(p,box)
        bonds=np.triu((bonds&(r<1.30))|fresh,1)
        travel=p-start; travel-=travel.mean(axis=0)
        emit(test='warm_hold',fraction=j/nread,signed_stress=reading,slope=(reading[1]-reading[2])/2e-4,
             half_slope=(reading[3]-reading[4])/1e-4 if args.half else None,original_bonds=int(initial_bonds.sum()),
             lost_original=int((initial_bonds&~bonds).sum()),new_neighbors=int((bonds&~initial_bonds).sum()),
             rms_displacement=float(np.sqrt(np.mean(np.sum(travel**2,axis=1)))))
    np.savez(args.work/'warm_endpoint.npz',pos=states[0][1],box=box)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--density',type=float,default=2.)
    parser.add_argument('--cutoff',type=float,default=2.5)
    parser.add_argument('--warm-from',type=Path)
    parser.add_argument('--mode',type=int,choices=(0,1),default=1)
    parser.add_argument('--half',action='store_true')
    parser.add_argument('--transit',type=float,default=48.88171243122866)
    args=parser.parse_args()
    lib=library(args.work)
    if args.warm_from:
        warm(args,lib)
    else:
        endpoint(args,lib)
