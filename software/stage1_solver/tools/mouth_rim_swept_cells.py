"""Linear-in-wall-step integration of the slivers between two Q1 field meshes.

This is a quadrature correction for the original mixed kinetic energy. Interior
faces of zero relative normal speed contribute nothing. At zero displacement
the meshes coincide. The sign of the step selects the appropriate one-sided
carrier trace. Energy and all three variations use the same face measure.
"""
from itertools import product

import numpy as np
from scipy import sparse
from scipy.interpolate import RegularGridInterpolator

from mouth_material_order_transfer import BITS, dot, grad, shape_basis
from mouth_supercore_joined_front import ORDER_INERTIA
from mouth_wound_carrier import straight_vortex


def face_basis(model, slab, bank, grids, indices, coords, radius, slope):
    widths = [np.diff(g)[i] for g,i in zip(grids,indices)]
    frac = np.column_stack([(point-g[i])/h for point,g,i,h in zip(coords,grids,indices,widths)])
    value, derivatives = shape_basis(frac)
    x = coords[0]
    wall = model.shape.radius(coords[2])[0]
    width = 28.-wall if bank else wall
    radial = derivatives[0] / (width*widths[0])[:,None]
    derivatives = [radial, derivatives[1] / (radius*widths[1])[:,None],
                   derivatives[2]/widths[2][:,None] - (slope*(1-x if bank else x))[:,None]*radial]
    corners = tuple(i[:,None]+BITS[:,j] for j,i in enumerate(indices))
    ids = model.ids[slab,bank][corners]
    return ids, corners, value, derivatives


def matrix(ids, values, size):
    pointers = np.arange(0,ids.size+1,8,dtype=np.int32)
    return sparse.csr_matrix((values.ravel(),ids.ravel(),pointers),shape=(len(ids),size))


class SweptCells:
    def __init__(self, model, reference, heights, material, carrier):
        self.size = model.size
        legs = [[],[]]
        weights, factors = [], []
        grid2 = tuple(np.linspace(0,1,k) for k in reference.shape)
        inside = RegularGridInterpolator(grid2,reference)
        vortex = straight_vortex(1,28.)
        for slab,bank in product((1,2),(0,1)):
            grids = (model.rr[bank],model.a,model.z[slab])
            X,Z = np.meshgrid(grids[0],grids[2],indexing="ij")
            R = model.shape.radius(Z)[0]
            if bank:
                ref = model.terminal_profile.outside(R+X*(28.-R),Z,slab==2)
            elif slab == 1:
                ref = inside((Z/model.shape.length,X))
            else:
                ref = model.terminal_profile(X,Z-model.shape.length)
            for axis in (0,2):
                ranges = [np.arange(len(g)-1) for g in grids]
                ranges[axis] = np.arange(1,len(grids[axis])-1)
                j = ranges[2]
                ranges[2] = j[(grids[2][j+int(axis!=2)]>heights[0])&(grids[2][j]<heights[-1])]
                indices = tuple(u.ravel() for u in np.meshgrid(*ranges,indexing="ij"))
                for uv in product((.5-.5/np.sqrt(3),.5+.5/np.sqrt(3)),repeat=2):
                    abscissas = iter(uv)
                    coords = []
                    measure = np.ones(len(indices[0]))
                    for d,(g,i) in enumerate(zip(grids,indices)):
                        if d == axis:
                            coords.append(g[i])
                        else:
                            widths = np.diff(g)[i]
                            coords.append(g[i]+next(abscissas)*widths)
                            measure *= widths/2.
                    x,a,z = coords
                    wall,slope = model.shape.radius(z)
                    width = 28.-wall if bank else wall
                    r = wall+width*x if bank else width*x
                    relative = carrier(bank,r,z)[0]-material(bank,r,z)[0]
                    oriented = relative*(-slope*(1-x if bank else x) if axis==0 else width)
                    take = oriented != 0
                    if not take.any():
                        continue
                    weights.append((measure*4*np.pi*r*r*np.sin(a)*oriented)[take])
                    factors.append(vortex((r*np.sin(a))[take])[0])
                    pts = [p[take] for p in coords]
                    for side in (0,1):
                        idx = tuple(i[take]-int(d==axis and side==0) for d,i in enumerate(indices))
                        ids,corners,v,dv = face_basis(model,slab,bank,grids,idx,pts,r[take],slope[take])
                        reference_values = ref[corners[0],corners[2]]
                        reference_gradient = np.column_stack([dot(reference_values,d) for d in dv])
                        legs[side].append((ids,v,*dv,reference_gradient))
        self.weight = np.concatenate(weights)
        self.factor = np.concatenate(factors)
        self.order_matrices = []
        self.reference_gradients = []
        for leg in legs:
            ids,v,*other = [np.concatenate(data) for data in zip(*leg)]
            derivatives,reference_gradient = other[:3],other[3]
            self.order_matrices.append([matrix(ids,dv,self.size) for dv in derivatives])
            self.reference_gradients.append(reference_gradient)
        self.values = matrix(ids,v,self.size)
        self.phase_matrices = {}
        for sign in (-1,1):
            lower = (sign*self.weight > 0)[:,None]
            self.phase_matrices[sign] = [
                (lo.multiply(lower)+hi.multiply(~lower)).tocsr()
                for lo,hi in zip(*self.order_matrices)
            ]

    def evaluate(self,u,step,gradient=None):
        phase_matrices = self.phase_matrices[int(np.sign(step))]
        amplitude = self.values@u[0]
        density = (self.factor*amplitude)**2
        gp = grad(phase_matrices,u[1])
        energy = 0.
        for sign,operators,ref in zip((-1,1),self.order_matrices,self.reference_gradients):
            gc = grad(operators,u[2])+ref
            cc = dot(gc,gc)
            pc = dot(gp,gc)
            coefficient = ORDER_INERTIA*pc/(density+ORDER_INERTIA*cc)
            w = sign*step*self.weight
            energy -= .5*w@(density*coefficient*pc)
            if gradient is not None:
                gradient[0] -= self.values.T@(w*self.factor**2*amplitude*coefficient**2*cc)
                load = w*density*coefficient
                for j in range(3):
                    gradient[1] -= phase_matrices[j].T@(load*gc[:,j])
                    gradient[2] -= operators[j].T@(load*(gp[:,j]-coefficient*gc[:,j]))
        return float(energy)

    def work(self,u,step):
        return self.evaluate(u,step)

    def attach(self,evaluate_volume,step):
        def total(u,read=False):
            if read:
                raise ValueError("Recover flux from the face and volume records separately")
            energy,gradient = evaluate_volume(u)
            return energy+self.evaluate(u,step,gradient),gradient
        return total
