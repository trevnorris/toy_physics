"""Material-attached conversion order with independent carrier trial functions.

The order well, gradient and crystalline load are integrated in the material
chart. The phase/order projection samples the material basis at carrier
integration sites. All force terms are transposes of the same value and
gradient operators that enter the energy. Both spaces coincide at zero shift.
"""
from itertools import product

import numpy as np
from scipy import sparse

from mouth_supercore_joined_front import KAPPA, ORDER_INERTIA, WELL

BITS = np.array(list(product((0, 1), repeat=3)))


def shape_basis(fraction):
    factors = np.where(BITS[None, :, :], fraction[:, None, :], 1-fraction[:, None, :])
    value = factors.prod(axis=2)
    raw = [(2*BITS[:, axis]-1) * factors[:, :, (axis+1)%3] * factors[:, :, (axis+2)%3]
           for axis in range(3)]
    return value, raw


def pull_gradient(G, dr, dz, step):
    inverse = 1/(1+step*dz)
    return [(G[0]+G[2].multiply((-step*dr*inverse)[:, None])).tocsr(), G[1],
            G[2].multiply(inverse[:, None]).tocsr()]


def pull_reference(ref, dr, dz, step):
    moved = ref.copy()
    moved[:, 3] /= 1+step*dz
    moved[:, 1] -= step*dr*moved[:, 3]
    return moved


def grad(G, v):
    return np.stack([M@v for M in G], axis=1)


def dot(x, y):
    return np.einsum('ij,ij->i', x, y)


class MixedOrder:
    """Apply to a model with its carrier volume and vortex surface already moved.

    Callback (side,r,z) returns axial motion and its physical derivatives
    in r,z on the undeformed mesh. Supply the same material map for all
    carrier mesh extensions; it must coincide on every physical boundary.
    """
    def __init__(self, model, reference_order=None):
        self.m = model
        self.evaluate_native = model.evaluate
        self.data, self.base_ref, self.base_G = model.D, model.ref, model.G
        q = (1+np.array([-1.,1.])/np.sqrt(3))/2
        self.gauss_inverse = np.linalg.inv(shape_basis(np.array(list(product(q,repeat=3))))[0])
        self.active = None

    def reset(self):
        self.m.evaluate = self.evaluate_native
        self.active = None

    def sampled_basis(self, slab, side, r, angle, z, first, cells, material_motion, step):
        m = self.m
        # Locate z in the undeformed MATERIAL mesh, preserving physical r.
        source_z = z.copy()
        for _ in range(8):
            displacement, _, dz = material_motion(side, r, source_z)
            increment = (source_z + step*displacement-z)/(1+step*dz)
            source_z -= increment
            if np.max(abs(increment)) < 5e-14:
                break
        displacement, dr, dz = material_motion(side, r, source_z)
        if np.max(abs(source_z+step*displacement-z)) > 3e-13:
            raise ArithmeticError('Unresolved inverse material map')
        R, slope = m.shape.radius(source_z)
        width = 28-R if side else R
        x = (r-R)/width if side else r/R
        grid = (m.rr[side], m.a, m.z[slab])
        xyz = (x, angle, source_z)
        if any(y.min()<g[0]-2e-12 or y.max()>g[-1]+2e-12 for y,g in zip(xyz,grid)):
            raise ValueError('Material point outside the matching slit bank')
        ix = [np.clip(np.searchsorted(g,y,side='right')-1,0,len(g)-2) for y,g in zip(xyz,grid)]
        gaps = [np.diff(g)[i] for g,i in zip(grid,ix)]
        fractions = np.stack([(y-g[i])/h for y,g,i,h in zip(xyz,grid,ix,gaps)],axis=1)
        basis, raw = shape_basis(fractions)
        ids = m.ids[slab,side][tuple(i[:,None]+BITS[:,j] for j,i in enumerate(ix))]
        rr = raw[0]/(gaps[0]*width)[:,None]
        gz = raw[2]/gaps[2][:,None] - (slope*(1-x if side else x))[:,None]*rr
        gz /= (1+step*dz)[:,None]
        derivative = [rr-(step*dr)[:,None]*gz,raw[1]/(gaps[1]*r)[:,None],gz]
        refcell = np.ravel_multi_index(tuple(ix), tuple(len(g)-1 for g in grid))
        local_refs = (self.gauss_inverse @ self.base_ref[first:first+8*cells,0].reshape(8,cells)[:,refcell]).T
        ref = np.column_stack([dot(local_refs,v) for v in [basis,*derivative]])
        def matrix(values):
            indptr = np.arange(0,ids.size+1,8,dtype=np.int32)
            return sparse.csr_matrix((values.ravel(), ids.ravel(), indptr),shape=(len(r),m.size))
        return matrix(basis), [matrix(g) for g in derivative], ref

    def install(self, step, material_motion, carrier_motion, material_load):
        m = self.m
        self.reset()
        m.material_load = material_load
        self.blocks = []
        m.ref = self.base_ref.copy()
        offset = 0
        for (slab,side),(nodes,_) in m.blocks.items():
            start = offset
            offset += 8*len(nodes)
            r, angle, z = self.data[start:offset,1:4].T
            S, Sr, Sz = material_motion(side,r,z)
            V, Vr, Vz = carrier_motion(side,r,z)
            index = np.arange(start,offset)
            # Native operations outside the relative mapping need only the
            # physical pullback of the reference gradients.
            moved = (Vr!=0) | (Vz!=0)
            target = index[moved]
            m.ref[target] = pull_reference(self.base_ref[target], Vr[moved], Vz[moved], step)
            use = (S!=V) | (Sr!=Vr) | (Sz!=Vz)
            if not np.any(use) or step == 0:
                continue
            rows = index[use]
            vb, gb, sampled_ref = self.sampled_basis(
                slab,side,r[use],angle[use],z[use]+step*V[use],start,len(nodes),material_motion,step)
            N = m.V[rows]
            G = [g[rows] for g in m.G]
            Sg = pull_gradient([g[rows] for g in self.base_G], Sr[use], Sz[use],step)
            Sw = self.data[rows,0]*(1+step*Sz[use])
            Sref = pull_reference(self.base_ref[rows],Sr[use],Sz[use],step)
            self.blocks.append((rows,N,G,m.ref[rows],m.hq[rows],Sg,Sw,Sref,vb,gb,sampled_ref))
        self.active = len(self.blocks)
        if self.active:
            m.evaluate = self.evaluate

    def evaluate(self,u,read=False):
        answer = self.evaluate_native(u,read)
        energy = answer[0]
        if read:
            _,chi_out,n_out,curr_out,vel_out,stress_out = answer
        else:
            grads = answer[1]
        a,p,c=u
        for rows,N,G,oldref,w,Sg,sw,Sref,Q,Qg,qref in self.blocks:
            t = N@a
            F = self.m.D[rows,4]
            n = (F*t)**2
            gp = grad(G,p)
            oldchi = oldref[:,0]+N@c
            oldgc = oldref[:,1:]+grad(G,c)
            schi = Sref[:,0]+N@c
            sgc = Sref[:,1:]+grad(Sg,c)
            qgc = qref[:,1:]+grad(Qg,c)
            for sign,weight,ch,gc,ref,operators in ((-1.,w,oldchi,oldgc,oldref,G),(1.,sw,schi,sgc,Sref,Sg)):
                c0 = ref[:,0]
                e = (.5*KAPPA*(dot(gc,gc)-dot(ref[:,1:],ref[:,1:])) +
                     WELL*((ch*(1-ch))**2-(c0*(1-c0))**2))
                energy += sign*(weight@e)
                if not read:
                    grads[2] += sign*(N.T@(weight*2*WELL*ch*(1-ch)*(1-2*ch)))
                    for j, A in enumerate(operators):
                        grads[2] += sign*(A.T@(weight*KAPPA*gc[:,j]))
            for sign,gc,operators in ((-1.,oldgc,G),(1.,qgc,Qg)):
                square = dot(gc,gc)
                gp_dot = dot(gp,gc)
                den = n+ORDER_INERTIA*square
                proj = ORDER_INERTIA*gp_dot/den
                vel = gp-proj[:,None]*gc
                energy -= .5*sign*(w@(n*proj*gp_dot))
                if not read:
                    grads[0] -= sign*(N.T@(w*(F**2*t)*(proj**2*square)))
                    factor = sign*w*n*proj
                    for j in range(3):
                        grads[1] -= G[j].T@(factor*gc[:,j])
                        grads[2] -= operators[j].T@(factor*vel[:,j])
            if read:
                # Sample the physical order on the carrier points for pressure
                # and current, independently of the mesh used for its energy.
                from mouth_supercore_joined_front import KAPPA as stiffness
                qchi = qref[:,0]+Q@c
                n_out[rows] = n
                delta_physical = (.5*stiffness*(dot(qgc,qgc)-dot(oldgc,oldgc))+
                                  WELL*((qchi*(1-qchi))**2-(oldchi*(1-oldchi))**2)+
                                  .5*n*dot(gp,vel-vel_out[rows]))
                stress_out[rows] += (n*(vel[:,2]**2-vel_out[rows,2]**2)+
                                     stiffness*(qgc[:,2]**2-oldgc[:,2]**2)-delta_physical)
                vel_out[rows] = vel
                curr_out[rows] = vel*n[:,None]
                chi_out[rows] = qchi
        return (energy,*answer[1:])
