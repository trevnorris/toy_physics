"""Release all nine four-dimensional cell shears at fixed volume, then probe the relaxed shape."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm, eigh
from scipy.optimize import minimize

from mouth_bulk_isotropic_4d import radial, strain_basis, structure
from neutral_patch_sites import emit


class PeriodicMaterial:
    def __init__(self,n,rc):
        self.n,self.rc=n,rc
        self.skin=.32
        self.reference=None

    def pairs(self,fractional,cell):
        if self.reference is not None:
            old,old_cell=self.reference
            deformation=np.linalg.solve(old_cell,cell)
            least=np.linalg.svd(deformation,compute_uv=False)[-1]
            travel=np.linalg.norm((fractional-old)@cell,axis=1).max()
            if least*(self.rc+self.skin)-2*travel>self.rc+1e-12:
                return self.indices
        reach=self.rc+self.skin
        widths=reach*np.linalg.norm(np.linalg.inv(cell),axis=0)
        extent=np.floor(.5+widths).astype(int)
        ids=np.arange(self.n)
        i,j=np.triu_indices(self.n,1)
        raw=fractional[j]-fractional[i]
        winding=-np.rint(raw)
        near=raw+winding
        ii,jj,offset=[],[],[]
        for tile in itertools.product(*(range(-k,k+1) for k in extent)):
            if np.any(np.maximum(np.abs(tile)-.5,0.)>widths):
                continue
            delta=(near+tile)@cell
            keep=np.einsum('ij,ij->i',delta,delta)<reach**2
            ii.append(i[keep]); jj.append(j[keep]);offset.append(winding[keep]+tile)
            if tile>(0,0,0,0) and np.linalg.norm(tile@cell)<reach:
                ii.append(ids);jj.append(ids);offset.append(np.tile(tile,(self.n,1)))
        self.indices=(np.concatenate(ii),np.concatenate(jj),np.vstack(offset))
        self.reference=(fractional.copy(),cell.copy())
        return self.indices

    def read(self,fractional,cell):
        i,j,offset=self.pairs(fractional,cell)
        q=fractional[j]-fractional[i]+offset
        d=q@cell
        r=np.linalg.norm(d,axis=1)
        e,du,_=radial(r,self.rc)
        force=(du/r)[:,None]*d
        grad=np.column_stack([np.bincount(j,f,minlength=self.n)-np.bincount(i,f,minlength=self.n) for f in force.T])
        return e.sum(),grad,q.T@force


def shape(params,L):
    diag=params[:3]
    out=np.diag(np.exp(np.r_[diag,-sum(diag)]))*L
    out[np.triu_indices(4,1)]=L*params[3:]
    return out


def joint(solver,pos,box):
    n=solver.n
    L=box[0]
    # Optimize fractional positions in nominal length units and nine upper
    # triangular shape coordinates. Upper triangular cells remove rigid rotations.
    weight=np.sqrt(n)
    def evaluate(flat,details=False):
        params=flat[-9:]/weight
        cell=shape(params,L)
        y=flat[:-9].reshape(-1,4)/L
        e,g,w=solver.read(y,cell)
        diag=np.diag(w)*np.diag(cell)
        gs=np.r_[diag[:3]-diag[-1],L*w[np.triu_indices(4,1)]]/weight
        if details:
            return y,cell,e,g,gs
        return e,np.r_[((g@cell.T)/L).ravel(),gs]
    x=np.r_[pos.ravel(),np.zeros(9)]
    # Bounds only protect far-off line-search trials. Accepted cells must be interior.
    bounds=[(None,None)]*(4*n)+[(-.55*weight,.55*weight)]*3+[(-.9*weight,.9*weight)]*6
    opt=minimize(evaluate,x,jac=True,method='L-BFGS-B',bounds=bounds,
                 options=dict(maxiter=15000,maxls=70,gtol=1e-9,ftol=2e-16,maxcor=25))
    y,cell,e,g,gs=evaluate(opt.x,True)
    p=opt.x[-9:]/weight
    if max(abs(p[:3]))>=.549 or max(abs(p[3:]))>=.899:
        raise ArithmeticError('Relaxed unit cell is still pressing against a numerical shape bound')
    if max(abs(g).max(),abs(gs).max())>1e-4:
        raise ArithmeticError(f'Unresolved shape equilibrium, {opt.message}, {abs(g).max()}, {abs(gs).max()}')
    return y,cell,p


def constitutive(solver,y,cell):
    # Build the 4N position Hessian from radial second derivatives. The strain
    # stiffness below uses independent central differences of measured stress.
    i,j,offset=solver.pairs(y,cell)
    d=(y[j]-y[i]+offset)@cell
    r=np.linalg.norm(d,axis=1)
    _,du,ddu=radial(r,solver.rc)
    g=du/r
    k=(ddu-g)/r**2
    h=np.zeros((solver.n,4,solver.n,4))
    for a,b in itertools.product(range(4),repeat=2):
        w=k*d[:,a]*d[:,b]+(a==b)*g
        edge=np.bincount(i*solver.n+j,w,minlength=solver.n**2).reshape(solver.n,solver.n)
        block=-edge-edge.T
        block[np.diag_indices(solver.n)]+=np.bincount(i,w,solver.n)+np.bincount(j,w,solver.n)
        h[:,a,:,b]=block
    return h.reshape(4*solver.n,-1)


def read_stress(solver,y,cell,basis):
    e,g,w=solver.read(y,cell)
    stress=(cell.T@w)/np.linalg.det(cell)
    return e,g,stress,np.einsum('ij,qij->q',stress,basis)


def run(args):
    args.work.mkdir(parents=True,exist_ok=True)
    source=args.source
    outcome=json.loads((source.parent.parent/'outcome.json').read_text())
    assert outcome['exitCode']==outcome['childOutcome']['exitCode']==0 and outcome['limitsVerified']
    with np.load(source) as data:
        pos,box=data['pos'].copy(),data['box'].copy()
        seed,rc=int(data['seed']),float(data['rc'])
    solver=PeriodicMaterial(len(pos),rc)
    basis=strain_basis()
    start=read_stress(solver,pos/box,np.diag(box),basis)
    y,cell,params=joint(solver,pos,box)
    # One final Newton polish including all internal and shape variables is
    # unnecessary for the tangent: the residual stress is independently checked.
    e,g,s,base=read_stress(solver,y,cell,basis)
    emit(test='free_shape',seed=seed,cutoff=rc,sites=len(y),density=len(y)/np.linalg.det(cell),
         initial_pressure=-np.trace(start[2])/4,pressure=-np.trace(s)/4,
         volume_ratio=np.linalg.det(cell)/np.prod(box),initial_shear_stress_norm=np.linalg.norm(start[3][1:]),
         final_shear_stress_norm=np.linalg.norm(base[1:]),max_site_force=np.max(abs(g)),
         energy_change_per_site=(e-start[0])/len(y),shape=cell.tolist(),shape_parameters=params.tolist())
    np.testing.assert_allclose(np.linalg.det(cell),np.prod(box),rtol=5e-13)
    assert np.linalg.norm(base[1:])<1e-4
    h=constitutive(solver,y,cell)
    ev,v=eigh(h)
    assert max(abs(ev[:4]))<1e-6 and ev[4]>1e-6
    # Response to a signed affine shape step, with the sites equilibrated
    # independently at each sign. Newton in laboratory coordinates. For these
    # tiny steps the reference Hessian is a convergent fixed preconditioner.
    def relax_at(new_cell):
        moved=y.copy()
        inv=np.linalg.inv(new_cell)
        for _ in range(10):
            _,force,_,value=read_stress(solver,moved,new_cell,basis)
            if np.max(abs(force))<2e-10:
                return value
            delta=v[:,4:]@((v[:,4:].T@force.ravel())/ev[4:])
            moved-=delta.reshape(-1,4)@inv
        raise ArithmeticError('Signed probe did not equilibrate')
    for eps in (1.e-4,5.e-5):
        cols=[]
        for mode in basis:
            pos_sig=relax_at(cell@expm(eps*mode))
            neg_sig=relax_at(cell@expm(-eps*mode))
            cols.append((pos_sig-neg_sig)/(2*eps))
        tangent=np.array(cols).T
        asym=np.max(np.abs(tangent-tangent.T))
        # Isotropic pre-stress drops out of shear; remaining antisymmetry is a
        # diagnostic of shear prestress and finite-difference error.
        shear=(tangent[1:,1:]+tangent[1:,1:].T)/2
        total=(tangent+tangent.T)/2
        shear_modes=np.linalg.eigvalsh(shear)
        emit(test='relaxed_probe',step=eps,bulk=tangent[0,0],min_internal=ev[4],
             shear_spectrum=shear_modes.tolist(),max_asymmetry=asym,bulk_shear_coupling=tangent[0,1:].tolist(),
             least_full_eigenvalue=np.linalg.eigvalsh(total)[0],
             stress_free_shape_bulk=(tangent[0,0]-tangent[0,1:]@np.linalg.solve(shear,tangent[1:,0])))
        if asym>.01:
            raise ArithmeticError('Resolved stress derivatives fail tangent symmetry')
    np.savez(args.work/'shape.npz',y=y,cell=cell,rc=rc,seed=seed)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    run(parser.parse_args())
