"""Cold excess work for a dilute constituent on a fixed-area flat material sample.

The host sites can relax in-plane. Their substrate plane is held; the added
site has two lateral coordinates and a free height. The coherent surface
potential is the periodic sum of the calibrated cross-sheet Morse bond. The
second law assigns the unregistered fraction the lateral mean of the same
potential. Both retain the same energy zero at infinite separation.
"""
import argparse

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

from mouth_calibrated_flare_coating import read_inputs
from mouth_middle_tag_grip import cubic_trace
from neutral_patch_periodic_compression import BASIS, Bilayer
from neutral_patch_sites import BETWEEN, WITHIN, emit, morse


def sheet_mean(z, pitch):
    """Energy and height derivative over a uniform sheet with one site per cell."""
    d, a, r0 = BETWEEN
    y = np.exp(-a*(z-r0))
    area = abs(np.linalg.det(BASIS))*pitch**2
    factor = 2*np.pi*d/area
    integral = y*y*(z/(2*a)+1/(4*a*a)) - 2*y*(z/a+1/(a*a))
    return factor*integral, -factor*z*y*(y-2)


def add_at_sites(ids, forces, count):
    return np.column_stack([np.bincount(ids, weights=g, minlength=count) for g in forces.T])


class Surface:
    def __init__(self, side, shell, pitch, rigid):
        cell = Bilayer(side, side, True, shell)
        host = (cell.i < cell.count) & (cell.j < cell.count)
        self.i, self.j = cell.i[host], cell.j[host]
        self.bonds = pitch*cell.vectors[host]
        self.bare = morse(np.linalg.norm(self.bonds,axis=1),*WITHIN)[0]
        self.prepared_energy = self.bare.sum()
        self.size, self.pitch, self.rigid = cell.count, pitch, rigid
        span = int(np.ceil((shell+3)/(.86*pitch)))
        ij = np.arange(-span,span+1)
        tiles = np.array(np.meshgrid(ij,ij,indexing='ij')).reshape(2,-1).T
        ids = (tiles+side//2) % side
        self.ids = ids[:,0]*side + ids[:,1]
        self.origins = pitch*(tiles@BASIS)
        self.dim = 3 if rigid else 2*(cell.count-1)+3

    def evaluate(self, state, chi, mean_fraction):
        offsets = np.zeros((self.size,2))
        if not self.rigid:
            offsets[1:] = state[3:].reshape(-1,2)  # remove common translation
        bonds = self.bonds+offsets[self.j]-offsets[self.i]
        lengths = np.linalg.norm(bonds,axis=1)
        energy, deriv = morse(lengths,*WITHIN)
        host_delta = (energy-self.bare).sum()
        lateral = state[:2]-self.origins-offsets[self.ids]
        z = state[2]
        r = np.sqrt((lateral*lateral).sum(axis=1)+z*z)
        pair, k = morse(r,*BETWEEN)
        coherent = pair.sum()
        uniform, normal = sheet_mean(z,self.pitch)
        spread = (1-chi)*float(mean_fraction)
        extra = host_delta+chi*coherent+spread*uniform
        projected = chi*(k/r)
        tag_gradient = np.r_[projected@lateral, z*projected.sum()+spread*normal]
        if self.rigid:
            gradient = tag_gradient
        else:
            forces = bonds*(deriv/lengths)[:,None]
            host_gradient = add_at_sites(self.j,forces,self.size)-add_at_sites(self.i,forces,self.size)
            host_gradient -= add_at_sites(self.ids,lateral*projected[:,None],self.size)
            gradient = np.r_[tag_gradient,host_gradient[1:].ravel()]
        return extra,gradient,coherent-mean_fraction*uniform,host_delta

    def equilibrate(self, chi, mixed, start):
        def objective(state):
            return self.evaluate(state,chi,mixed)[:2]
        bounds = [(None,None)]*len(start)
        bounds[2] = (.25,3.)  # select the exposed surface; no endpoint touches either stop
        sol = minimize(objective,start,method='L-BFGS-B',jac=True,bounds=bounds,
                       options={'gtol':1e-10,'ftol':2e-15,'maxiter':3000,'maxcor':24,'maxls':45})
        e,g,slope,host = self.evaluate(sol.x,chi,mixed)
        if np.max(abs(g))>1.5e-5 or not .26<sol.x[2]<2.99:
            raise ArithmeticError(f'Unresolved bound site: {sol.message}; E={e}; max |g|={np.max(abs(g))}')
        return sol.x, e, slope, dict(chi=chi, height=sol.x[2], excess=e, host_strain_cost=host,
                                     order_derivative=slope, max_force=np.max(abs(g)))

    def stability(self,state,chi,mixed):
        h=2e-5
        steps=np.eye(self.dim)*h
        jac=np.column_stack([(self.evaluate(state+d,chi,mixed)[1]-self.evaluate(state-d,chi,mixed)[1])/(2*h)
                             for d in steps])
        return np.linalg.eigvalsh((jac+jac.T)/2)[0]


def run(args):
    trace,rate = cubic_trace(args.mesh)
    patch,_,_ = read_inputs(args.mesh)
    pitch,scale = patch['pitch'],2*patch['shared_scale']
    inlet=trace.sample(trace.chi,np.array([1.175]))[0].item()
    outlet=.5
    sample=Surface(args.side,args.shell,pitch,args.rigid)
    empty=np.zeros(sample.dim)
    empty[2]=1.
    # Removing the probe leaves an exact equilibrium by inversion symmetry in the
    # unchanged periodic cell; subtract the host at each order on this branch.
    _,d=morse(np.linalg.norm(sample.bonds,axis=1),*WITHIN)
    f=sample.bonds*(d/np.linalg.norm(sample.bonds,axis=1))[:,None]
    residual=np.max(abs(add_at_sites(sample.j,f,sample.size)-add_at_sites(sample.i,f,sample.size)))
    emit(test='sample', mesh=args.mesh, intake=rate, side=args.side, shell=args.shell,
         host_sites=sample.size, rigid_host=args.rigid, host_energy=sample.prepared_energy,
         matched_empty_host_max_force=residual,pitch=pitch,scale=scale,entry_order=inlet,
         exit_order=outlet, half_order_depth=trace.crossing(.5), energy_unit='per participating site mass after scaling')
    if residual>1e-10 or abs(rate-.25)>1e-8:
        raise ArithmeticError('Mismatched reference')
    start=empty.copy()
    start[:2]=pitch*BASIS.sum(axis=0)/3
    start[2]=.75
    for mixed in (False,True):
        label='weaken_all_contacts' if not mixed else 'retain_mean_attraction'
        first,e1,_,row1=sample.equilibrate(inlet,mixed,start)
        last,e2,_,row2=sample.equilibrate(outlet,mixed,first)
        # Reverse relaxation plus a separate nonsymmetric local preparation.
        back,eback,_,_=sample.equilibrate(inlet,mixed,last)
        perturbed=start.copy()
        perturbed[:3]+=[.14,-.07,.06]
        _,alternate,_,_=sample.equilibrate(outlet,mixed,perturbed)
        ev=[sample.stability(s,chi,mixed) for s,chi in [(first,inlet),(last,outlet)]]
        nodes,weights=leggauss(args.order_nodes)
        integral=0.
        warm=first
        for x,weight in zip(nodes,weights,strict=True):
            chi=(inlet+outlet)/2+(outlet-inlet)*x/2
            warm,_,slope,_=sample.equilibrate(chi,mixed,warm)
            integral+=weight*slope*(outlet-inlet)/2
        transfer=(e2-e1)*scale
        emit(test='transfer', law=label, entry=row1, destination=row2,
             min_hessian=ev, energy_cycle=transfer, integrated_work=integral*scale,
             half_order_alternative_difference=scale*(alternate-e2),
             reverse_difference=scale*(eback-e1), order_nodes=args.order_nodes)
        if min(ev)<=0 or max(abs(eback-e1),abs(alternate-e2),abs(integral-(e2-e1)))>2e-5:
            raise ArithmeticError('Cold stable-branch endpoint or work check failed')

    # Read a weak probe over the entire lateral unit cell on the blank, supported
    # host. The registry average has a different normal force for the two laws.
    z=1.1
    ref,normal=sheet_mean(z,pitch)
    ij=np.array(np.meshgrid(np.arange(28)+.5,np.arange(28)+.5,indexing='ij')).reshape(2,-1).T/28
    vals=[]
    for xy in pitch*(ij@BASIS):
        pos=empty.copy()
        pos[:3]=[*xy,z]
        e,g,_,_=sample.evaluate(pos,1.,False)
        vals.append([e,g[2]])
    avg=np.mean(vals,axis=0)
    emit(test='mean_normal_probe',z=z,area_average_sampled=avg.tolist(),area_average_exact=[ref,normal],
         upward_holding_per_specific_scale_at_entry_half=[
             [inlet*normal,outlet*normal],[normal,normal]],
         destination_to_entry_ratio=[outlet/inlet,1.],
         force_unit='multiply tabulated dE/dz by specific energy scale / microscopic length')
    np.testing.assert_allclose(avg,[ref,normal],atol=2e-9,rtol=2e-9)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(26,34),required=True)
    p.add_argument('--side',type=int,required=True)
    p.add_argument('--shell',type=float,default=5.)
    p.add_argument('--order-nodes',type=int,default=6)
    p.add_argument('--rigid',action='store_true')
    a=p.parse_args()
    if a.side<3 or a.shell<4 or a.order_nodes<2:
        p.error('Use at least 3 sites per side, a shell of 4 and 2 integration nodes')
    run(a)
