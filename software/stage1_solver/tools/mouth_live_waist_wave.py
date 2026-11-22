"""Order-dependent outgoing transverse pole and its reciprocal force in the occupied duct."""
import numpy as np
from scipy import sparse

from mouth_exchange_shear import DEGREE, OutgoingSlab, assemble, block_matrices, solve


def transfer(n, edges, derivative=False):
    """Cell means of a linear nodal interpolant or of its coordinate derivative."""
    nodes = np.linspace(0., 1., n+1)
    weights = np.zeros((len(edges)-1,n+1))
    for j, (lo, hi) in enumerate(zip(edges[:-1],edges[1:])):
        split = np.r_[lo, nodes[(nodes>lo)&(nodes<hi)], hi]
        for a,b in zip(split[:-1],split[1:]):
            i=min(n-1,int(np.floor(n*(a+b)/2)))
            if derivative:
                left, right = -n,n
            else:
                right=n*(a+b)/2-i
                left=1-right
            weights[j,i:i+2] += (b-a)/(hi-lo)*np.array([left,right])
    return sparse.csr_matrix(weights)


class OccupiedWave:
    def __init__(self, shape, nr, nz, mesh, radius, seed=None, target=None):
        self.shape,self.mesh,self.radius=shape,mesh,radius
        self.seed=seed
        self.target=target
        # Keep numbering fixed under both ± shape steps for continuing the same material field.
        axial_cells=int(np.ceil(2.77*shape.neck*mesh))
        self.uniform=assemble(shape,mesh,radius,0.,axial_cells)
        points, k,m,_=self.uniform
        self.port=OutgoingSlab(points,radius)
        x=np.linspace(0.,1.,mesh+1)
        z=np.linspace(0.,1.,axial_cells+1)
        def mapping(x,w):
            b,s=shape.radius(w)
            return x*b,x*s
        pp,ic,kc,mc=block_matrices(x,shape.length*z,mapping,0.,shape.length,DEGREE)
        ix={tuple(c):i for i,c in enumerate(np.round(points,10))}
        # Axis has zero degree-four amplitude; use index -1 only in padded vectors.
        nodes=np.array([ix.get(tuple(c),-1) for c in np.round(pp,10)])
        self.ids=nodes[ic]
        self.kc,self.mc=kc,mc
        shape4=kc.shape
        self.rows=np.broadcast_to(self.ids[:,:,None],shape4).ravel()
        self.cols=np.broadcast_to(self.ids[:,None,:],shape4).ravel()
        self.live=(self.rows>=0)&(self.cols>=0)
        self.size=k.shape
        self.project=sparse.kron(transfer(nr,x),transfer(nz,z),format='csr')
        self.nz,self.nr=nz,nr
        self.cut_k=k-self.matrix(kc)
        self.cut_m=m-self.matrix(mc)

    def matrix(self, data):
        ids=self.live
        return sparse.coo_matrix((data.ravel()[ids],(self.rows[ids],self.cols[ids])),
                                 shape=self.size).tocsc()

    def pole(self, weights):
        k=self.matrix(weights[:,None,None]*self.kc)
        m=self.matrix(weights[:,None,None]*self.mc)
        system=(self.uniform[0],self.cut_k+k,self.cut_m+m,m)
        prior,vector=(None,None) if self.seed is None else self.seed
        out,v,mass=solve(self.shape,self.mesh,self.radius,0.,prior=prior,incoming_vector=vector,
                         shape_work=False,supplied_system=system,target_frequency=self.target)
        frequency=complex(out['frequency'],-out['amplitude_loss'])
        overlap=1. if vector is None else abs(np.vdot(vector,mass@v))**2/np.vdot(vector,mass@vector).real
        if overlap < .88:
            raise ArithmeticError(f'Occupied transverse mode correspondence lost: {overlap}')
        self.seed=(frequency,v)
        self.target=out['frequency']
        return out,v,mass,frequency,overlap

    def response(self, order):
        weight=np.asarray(self.project@order.T.ravel()).ravel()
        if min(weight) <= 0:
            raise ValueError('Local shear field lost all stiffness and inertia in a wave element')
        overlap_floor=1.
        if self.seed is None:
            self.pole(np.ones_like(weight))
            for s in (.25,.5,.75):
                *_,ov=self.pole(1+s*(weight-1))
                overlap_floor=min(ov,overlap_floor)
        out,v,mass,omega,ov=self.pole(weight)
        vec=np.r_[v,0j][self.ids]
        kcell=np.einsum('ci,cij,cj->c',vec.conj(),self.kc,vec).real
        mcell=np.einsum('ci,cij,cj->c',vec.conj(),self.mc,vec).real
        ktotal=np.vdot(v,self.cut_k@v).real + weight@kcell
        self.strain_weights=np.asarray(self.project.T@kcell).reshape(self.nr+1,self.nz+1).T/ktotal
        self.inertia_weights=np.asarray(self.project.T@mcell).reshape(self.nr+1,self.nz+1).T/np.vdot(v,mass@v).real
        ds=self.kc-omega**2*self.mc
        diff=np.einsum('ci,cij,cj->c',vec,ds,vec)
        boundary_v=v[self.port.ids]
        derivative=(self.port.block(omega+2e-6)-self.port.block(omega-2e-6))/(4e-6)
        den=2*omega*(v@(mass@v))+boundary_v@(derivative@boundary_v)
        grad=np.asarray(self.project.T@(diff/den).real).reshape(self.nr+1,self.nz+1).T.copy()
        return out['frequency'],grad,dict(wave_frequency=out['frequency'],
                   wave_pocket_fraction=out['confined_kinetic_fraction'],
                   wave_loss_rate=out['outward_power_over_interior_energy'],
                   pole_amplitude_loss=out['amplitude_loss'],mode_overlap=min(overlap_floor,ov),
                   pole_residual=out['pole_residual'],wave_minimum_cell_order=float(min(weight)))
