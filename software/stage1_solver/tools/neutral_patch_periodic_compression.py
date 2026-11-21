"""Compare cold bilayers with and without in-plane edges under the same compression."""
import argparse
from pathlib import Path
import json

import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigvalsh

from neutral_patch_carrier_compression import CALIBRATION, DENSITIES, carrier_pressure
from neutral_patch_sites import BETWEEN, COUNT, GAP, WITHIN, emit, morse

BASIS = np.array([[1.,0.],[.5,np.sqrt(3.)/2]])
PRIOR = Path('research/pde_ledger_v3/_measurements/neutral_patch_carrier_compression.json')


class Bilayer:
    """Count each bond per simulation cell, including distinct periodic images."""
    def __init__(self, nx, ny, periodic, shell):
        indices = np.indices((nx,ny)).reshape(2,-1).T
        self.count = nx*ny
        self.periodic = periodic
        self.reference = np.tile(indices@BASIS,(2,1))
        edge = np.any((indices==0)|(indices==np.array([nx,ny])-1),axis=1)
        self.mobile = np.ones(2*self.count,bool) if periodic else np.tile(~edge,2)
        pairs, vectors, cross = [], [], []
        extent = int(np.ceil(2*shell/np.sqrt(3)))
        for i in range(-extent,extent+1):
            for j in range(-extent,extent+1):
                shift = np.array([i,j])
                vector = shift@BASIS
                if vector@vector>shell**2+1e-10:
                    continue
                image = indices+shift
                keep = np.ones(self.count,bool) if periodic else np.all((image>=0)&(image<[nx,ny]),axis=1)
                if not np.any(keep):
                    continue
                start = np.flatnonzero(keep)
                dest = image[keep]%[nx,ny]
                end = dest[:,0]*ny+dest[:,1]
                layer_pairs = [(0,1)]
                if i>0 or (i==0 and j>0):
                    layer_pairs.extend([(0,0),(1,1)])
                for a,b in layer_pairs:
                    pairs.append(np.column_stack((start+a*self.count,end+b*self.count)))
                    vectors.append(np.tile(vector,(len(start),1)))
                    cross.append(np.full(len(start),a!=b))
        self.i,self.j = np.concatenate(pairs).T
        self.vectors = np.concatenate(vectors)
        self.cross = np.concatenate(cross)
        self.params = [np.where(self.cross,b,a) for a,b in zip(WITHIN,BETWEEN)]

    def read(self,n,strained_gap,offset,stiffness=False):
        stretch = n**(-1/3)
        gap = GAP*(stretch if strained_gap else 1.)
        d = stretch*self.vectors+offset[self.j]-offset[self.i]
        g2 = self.cross*gap**2
        d2 = np.einsum('bi,bi->b',d,d)
        radius = np.sqrt(d2+g2)
        e,first = morse(radius,*self.params)
        kr = first/radius
        forces = kr[:,None]*d
        size = len(offset)
        grad = np.column_stack([np.bincount(self.j,weights=f,minlength=size)-
                                np.bincount(self.i,weights=f,minlength=size) for f in forces.T])
        if stiffness:
            depth,alpha,rest = self.params
            y = np.exp(-alpha*(radius-rest))
            second = 2*depth*alpha**2*y*(2*y-1)
            hess = np.zeros((2*size,2*size))
            for a in range(2):
                for b in range(2):
                    block = (second-kr)*d[:,a]*d[:,b]/radius**2+(a==b)*kr
                    for u,v,sign in ((self.i,self.i,1),(self.j,self.j,1),
                                     (self.i,self.j,-1),(self.j,self.i,-1)):
                        np.add.at(hess,(2*u+a,2*v+b),sign*block)
            return hess
        virial = kr@(d2+(g2 if strained_gap else 0.))
        return e.sum()/self.count, grad, -n*virial/(3*self.count)

    def equilibrium(self,n,strained_gap,seed=None):
        offset = np.zeros_like(self.reference)
        if seed is not None:
            offset[self.mobile] = seed[self.mobile]
        def objective(flat):
            offset[self.mobile] = flat.reshape(-1,2)
            energy,gradient,_ = self.read(n,strained_gap,offset)
            return energy*self.count,gradient[self.mobile].ravel()
        energy, grad = objective(offset[self.mobile].ravel())
        if max(abs(grad))>1e-9:
            result = minimize(objective,offset[self.mobile].ravel(),jac=True,method='L-BFGS-B',
                              options=dict(ftol=1e-15,gtol=1e-9,maxiter=1000,maxls=40))
            offset[self.mobile] = result.x.reshape(-1,2)
            if self.periodic:
                offset -= offset.mean(axis=0)
            # Finish a converged minimization with a translation-fixed Newton step
            # when energy differences at machine precision obscure the last forces.
            for _ in range(3):
                energy,grad = objective(offset[self.mobile].ravel())
                if max(abs(grad))<=1e-9:
                    break
                hess = self.read(n,strained_gap,offset,True)
                active = np.repeat(self.mobile,2)
                if self.periodic:
                    active[:2] = False
                full_grad = self.read(n,strained_gap,offset)[1]
                offset.flat[active] -= np.linalg.solve(hess[active][:,active],full_grad.ravel()[active])
            if self.periodic:
                offset -= offset.mean(axis=0)
        energy,grad,pressure = self.read(n,strained_gap,offset)
        if np.max(abs(grad[self.mobile]))>2e-8 or np.max(abs(offset))>.40:
            raise ArithmeticError('Unresolved force balance or large departure from the registered branch')
        return energy,pressure,offset

    def stability(self,n,strained_gap,offset):
        hess = self.read(n,strained_gap,offset,True)
        active = np.repeat(self.mobile,2)
        low = eigvalsh(hess[active][:,active],subset_by_index=(0,2))
        if self.periodic:
            np.testing.assert_allclose(low[:2],0.,atol=5e-9)
        minimum = low[2] if self.periodic else low[0]
        if minimum<=1e-7:
            raise ArithmeticError(f'Non-translational instability at density {n}: {low}')
        return float(minimum)


def experiment(nx,ny,periodic):
    cell = Bilayer(nx,ny,periodic,6.)
    shell_control = Bilayer(nx,ny,periodic,8.)
    rng = np.random.default_rng(1200+nx*ny)
    disturbed = rng.normal(0,.025,cell.reference.shape)
    disturbed -= disturbed.mean(axis=0)
    prior = json.loads(PRIOR.read_text())['runs']['neutral_patch_carrier_compression01']['literal_stdout']
    for strained_gap in (False,True):
        stack = 'brane' if strained_gap else 'w'
        former = next(r for r in map(json.loads,prior) if r.get('stack')==stack and r.get('interior')=='relaxed')
        # Archive used total patch work. Here pressure is per two-site cell.
        previous_scale = COUNT*former['energy_per_rest_four_volume']
        _,base,x0 = cell.equilibrium(1.,strained_gap)
        _,alternate,shaken = cell.equilibrium(1.,strained_gap,disturbed)
        np.testing.assert_allclose(shaken,x0,atol=2e-7)
        np.testing.assert_allclose(alternate,base,atol=2e-8)
        rows = []
        for n in DENSITIES:
            energy,p,offset = cell.equilibrium(n,strained_gap,x0)
            stiffness = cell.stability(n,strained_gap,offset)
            large_e,large_p,_ = shell_control.equilibrium(n,strained_gap,offset)
            np.testing.assert_allclose([energy,p],[large_e,large_p],rtol=0,atol=3e-9)
            # Check the virial via independently strained and relaxed energies.
            dv = 1e-5
            hi = cell.equilibrium(1/(1/n+dv),strained_gap,offset)[0]
            lo = cell.equilibrium(1/(1/n-dv),strained_gap,offset)[0]
            fd=(lo-hi)/(2*dv)
            np.testing.assert_allclose(fd,p,rtol=1e-7,atol=4e-8)
            rows.append(dict(n=n,raw_change=p-base,min_stiffness=stiffness,work_error=abs(fd-p)))
        fit=next(row for row in rows if row['n']==CALIBRATION)
        scale = carrier_pressure(CALIBRATION)/fit['raw_change']
        curvature=[]
        for h in (.001,.0005):
            plus,minus=[cell.equilibrium(np.exp(s*h),strained_gap,x0)[1]-base for s in (1,-1)]
            slope=(plus-minus)/(2*h)
            curvature.append((plus+minus)/(h*h*slope))
        for r in rows:
            n=r['n']; target=carrier_pressure(n)
            r.update(target=target,predicted=scale*r['raw_change'],
                     previous_finite_scale_prediction=previous_scale*r['raw_change'],
                     fractional_error=scale*r['raw_change']/target-1 if n!=1 else None)
        emit(test='compression',nx=nx,ny=ny,sites=2*cell.count,mobile=int(sum(cell.mobile)),
             boundary='periodic' if periodic else 'open_perimeter',stack=stack,
             cutoff_reference_planar=[6.,8.],rest_virial=base,
             fit_density=CALIBRATION,energy_per_reference_pair_volume=scale,
             scale_from_finite_patch=previous_scale,curvature_steps=[.001,.0005],
             curvature_ratio=curvature,disturbed_branch_max_difference=float(np.max(abs(shaken-x0))),rows=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('nx',type=int)
    parser.add_argument('ny',type=int)
    a=parser.parse_args()
    for periodic in (True,False):
        experiment(a.nx,a.ny,periodic)
