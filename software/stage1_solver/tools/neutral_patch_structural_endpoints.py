"""Cold structural endpoints of the neutral Morse bilayer at matched density."""
import argparse
from itertools import product

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

from neutral_patch_sites import BETWEEN, GAP, WITHIN, emit, morse


PITCH = .9975140800337228  # independently measured ordered zero-load pitch
PULSE_TIME = 1.           # native unit-mobility patch time


class Cell:
    def __init__(self, side, cutoff):
        self.side = side
        self.n = 2*side**2
        self.box = PITCH*side*np.array([1., np.sqrt(3.)/2])
        self.area = np.prod(self.box)
        self.cutoff = cutoff
        # Expand around each pair's actual nearest rectangular image. Truncate
        # on physical bond radius; preserve self-images with one representative
        # of each opposite pair. No neighbour list refers to the initial packing.
        i, j = np.triu_indices(self.n, 1)
        extent = np.ceil(cutoff/self.box+.6).astype(int)
        tiles = list(product(range(-extent[0], extent[0]+1), range(-extent[1], extent[1]+1)))
        pairs, shifts = [], []
        for a, b in tiles:
            min_distance = np.maximum(np.abs([a,b])-.5, 0) * self.box
            if np.linalg.norm(min_distance) > 1.01*cutoff:
                continue
            ii, jj = i, j
            if a > 0 or (a == 0 and b > 0):
                ii, jj = np.r_[i,np.arange(self.n)], np.r_[j,np.arange(self.n)]
            pairs.append((ii,jj))
            shifts.append(np.tile([a,b], (len(ii),1)))
        self.i, self.j = map(np.concatenate, zip(*pairs))
        self.shifts = np.concatenate(shifts)*self.box

    def bonds(self, positions, strain):
        delta = positions[self.j]-positions[self.i]
        delta += self.shifts - self.box*np.rint(delta/self.box)
        dy = delta[:,1]
        dx = delta[:,0] + strain*dy
        cross = (self.i//(self.n//2) != self.j//(self.n//2))
        r2 = dx**2+dy**2+cross*GAP**2
        keep = r2 < self.cutoff**2
        params = [np.where(cross[keep], b, a) for a,b in zip(WITHIN,BETWEEN)]
        return self.i[keep], self.j[keep], dx[keep], dy[keep], np.sqrt(r2[keep]), params

    def read(self, positions, strain=0., tangent=False):
        i,j,x,y,r,params = self.bonds(positions,strain)
        e,first = morse(r,*params)
        g = first/r
        def incidence(values):
            return np.bincount(j, values, minlength=self.n)-np.bincount(i, values, minlength=self.n)
        # Variables supplied to the minimizer are coordinates in the unsheared cell.
        grad = np.column_stack((incidence(g*x),incidence(g*(y+strain*x))))
        tension = (g*x)@y/self.area
        if not tangent:
            return e.sum(), grad, tension
        if strain != 0.:
            raise ValueError('Tangent is defined at the unstrained preparation')
        depth,alpha,length = params
        z = np.exp(-alpha*(r-length))
        s = (2*depth*alpha**2*z*(2*z-1)-g)/r**2
        xx,xy,yy = g+s*x*x,s*x*y,g+s*y*y
        hess = np.zeros((2*self.n,2*self.n))
        for a,b,block in ((0,0,xx),(0,1,xy),(1,0,xy),(1,1,yy)):
            for u,v,sign in ((i,i,1),(j,j,1),(i,j,-1),(j,i,-1)):
                np.add.at(hess,(2*u+a,2*v+b),sign*block)
        # In unsheared coordinates there is also incidence(g*x) in the y
        # component. It vanishes at force balance; retain it for consistency.
        coupling = np.column_stack((incidence(xx*y),incidence(xy*y+g*x))).ravel()
        return hess, coupling, xx@(y*y)

    def settle(self, seed, strain=0.):
        # Remove common translation for this minimization only.
        def state(flat):
            return np.vstack((seed[:1], flat.reshape(-1,2)))
        def objective(flat):
            energy,grad,_ = self.read(state(flat),strain)
            return energy, grad[1:].ravel()
        fit = minimize(objective, seed[1:].ravel(), jac=True, method='L-BFGS-B',
                       options=dict(ftol=2e-16,gtol=1e-10,maxiter=8000,maxls=70,maxcor=15))
        pos = state(fit.x)
        if strain == 0.:
            for _ in range(3):
                _,grad,_ = self.read(pos)
                if np.max(abs(grad)) < 1e-9:
                    break
                hess = self.read(pos,tangent=True)[0][2:,2:]
                pos[1:] -= np.linalg.solve(hess,grad[1:].ravel()).reshape(-1,2)
        energy,grad,stress = self.read(pos,strain)
        residual = np.max(abs(grad))
        if residual > (1e-7 if strain == 0. else 2e-5):
            raise ArithmeticError(f'Unresolved site forces at shear {strain}: {residual}')
        return pos,energy,stress,residual


def preparation(cell, seed):
    if seed == 0:
        i,j = np.indices((cell.side,cell.side))
        plane = PITCH*np.column_stack(((i+.5*(j%2)).ravel(),(j*np.sqrt(3.)/2).ravel()))
        return np.tile(plane,(2,1))
    rng = np.random.default_rng(seed)
    planes = []
    for _ in range(2):
        plane = []
        while len(plane) < cell.n//2:
            candidate = rng.random(2)*cell.box
            delta = np.asarray(plane).reshape(-1,2)-candidate
            delta -= cell.box*np.rint(delta/cell.box)
            if np.all(np.sum(delta**2,axis=1) > .55**2):
                plane.append(candidate)
        planes.extend(plane)
    return np.array(planes)


def structure(cell, positions):
    moments = []
    six = []
    for plane in np.split(positions,2):
        delta = plane[:,None]-plane
        delta -= cell.box*np.rint(delta/cell.box)
        radii = np.linalg.norm(delta,axis=2)
        near = (radii > 0.) & (radii < 1.35*PITCH)
        count = near.sum(axis=1)
        moment = np.sum(near*np.exp(6j*np.arctan2(delta[:,:,1],delta[:,:,0])),axis=1)
        moments.extend(moment/np.maximum(count,1))
        six.extend(count==6)
    return dict(local_sixfold=float(np.mean(np.abs(moments))),
                coherent_sixfold=float(abs(np.mean(moments))),six_neighbors=float(np.mean(six)))


def run(side, seed, refine):
    cell = Cell(side,6.)
    initial = preparation(cell,seed)
    pos,energy,stress,residual = cell.settle(initial)
    hess,b,born = cell.read(pos,tangent=True)
    eig,modes = eigh(hess)
    if max(abs(eig[:2]))>1e-7 or eig[2]<1e-7:
        raise ArithmeticError(f'State is not a resolved fixed-cell local minimum: {eig[:3]}')
    np.testing.assert_allclose(modes[:,:2].T@b,0,atol=2e-7)
    k = eig[2:]
    weights = (modes[:,2:].T@b)**2/(k*cell.area)
    affine = born/cell.area
    relaxed = affine-weights.sum()
    # Quasistatic strain controls: start on the linear equilibrium branch, then
    # independently allow all sites to move again and differentiate actual stress.
    correction = -(modes[:,2:]@((modes[:,2:].T@b)/k)).reshape(-1,2)
    correction -= correction[0].copy()
    checks = []
    steps = [.0002,.0001,.00005,.000025] if refine else [.0002,.0001]
    for h in steps:
        sides = [cell.settle(pos+s*h*correction,s*h) for s in (-1,1)]
        _,em,sm,_ = sides[0]
        _,ep,sp,_ = sides[1]
        checks.append(dict(strain=h,stress_slope=(sp-sm)/(2*h),
                           energy_curvature=(ep+em-2*energy)/(cell.area*h**2)))
    for check in checks[-2:]:
        np.testing.assert_allclose([check['stress_slope'],check['energy_curvature']],
                                   relaxed,atol=.003,rtol=5e-4)
    large = Cell(side,7.)
    big_e,big_g,big_s = large.read(pos)
    np.testing.assert_allclose([big_e,big_s],[energy,stress],rtol=0,atol=2e-8)
    np.testing.assert_allclose(big_g,0.,rtol=0,atol=2e-7)
    def held(t):
        return relaxed + weights@np.exp(-k*t)
    times = [0.,.5,1.,1.5,2.,4.]
    pulse = [held(t) if t<PULSE_TIME else
             weights@(np.exp(-k*t)-np.exp(-k*(t-PULSE_TIME))) for t in times]
    emit(test='structural_endpoint',side=side,sites=cell.n,seed=seed,pitch=PITCH,
         initial=structure(cell,initial),resting=structure(cell,pos),
         energy_per_site=energy/cell.n,baseline_shear=stress,max_site_force=residual,
         lowest_internal_stiffness=k[0],affine_shear_modulus=affine,
         relaxed_shear_modulus=relaxed,times=times,held=[held(t) for t in times],
         pulse_duration=PULSE_TIME,pulse=pulse,checks=checks)


if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument('side',type=int,choices=[4,6])
    args.add_argument('seed',type=int)
    args.add_argument('--refine-strain',action='store_true',
                      help='check convergence to the tangent with two additional halvings')
    pars=args.parse_args()
    run(pars.side,pars.seed,pars.refine_strain)
