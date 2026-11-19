"""Free compact links across seven prepared shells and a source-free exterior.

One orienting face, four periodic springs per face loop. Two coherent arms
select the dominant spin rays with visibility proportional to the geometric
mean of their polarized fractions. All energies use one spring as unit.
"""
import argparse

import numpy as np
from scipy import sparse
from scipy.linalg import eigh
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

from mouth_collar_defect import states
from mouth_collar_order_preparation import calibration, predicted
from mouth_material_links import section, spherical_areas


RIM = 64.
MATCH_RADIUS = 16.
ROTATION = [.26, -.19, .17]
SHELLS = 7


def principal(a):
    return (a + np.pi) % (2*np.pi) - np.pi


def W(f):
    return 8*np.sin(f/8)**2


def angular_mesh(level):
    normals, faces = section(level)
    normals = Rotation.from_rotvec(ROTATION).apply(normals)
    i = faces.ravel()
    j = np.roll(faces, -1, axis=1).ravel()
    edges, inv = np.unique(np.sort(np.column_stack((i,j)), axis=1),
                           axis=0, return_inverse=True)
    rows = np.repeat(np.arange(len(faces)), 3)
    B = sparse.csr_matrix((np.sign(j-i), (rows, inv)), shape=(len(faces),len(edges)))
    left,right = edges.T
    D = sparse.csr_matrix((np.tile([-1,1],len(edges)),
                           (np.repeat(np.arange(len(edges)),2), edges.ravel())),
                          shape=(len(edges), len(normals)))
    assert (B @ D).nnz == 0
    areas = spherical_areas(normals, faces)
    centers = normals[faces].sum(axis=1)
    centers /= np.linalg.norm(centers,axis=1)[:,None]
    adj = B.tocsc().indices.reshape(-1, 2)
    arc = lambda u,v: np.arccos(np.clip((u*v).sum(axis=1),-1.,1.))
    conductance = arc(normals[left],normals[right])/arc(centers[adj[:,0]],centers[adj[:,1]])
    # Finite-volume sphere: angular gradient energy on neighboring face cells.
    mass = sparse.diags(areas**(-.5))
    laplacian = mass @ B @ sparse.diags(conductance) @ B.T @ mass
    eigenvalues, modes = eigh(laplacian.toarray())
    assert abs(eigenvalues[0]) < 1e-10 and eigenvalues[1] > 0
    assert np.isclose(areas.sum(), 4*np.pi)
    ell = (np.sqrt(1+4*np.maximum(eigenvalues,0))-1)/2
    transform = modes.T / np.sqrt(areas)
    exterior = transform.T @ ((1/(MATCH_RADIUS*(1+ell)))[:,None]*transform)
    assert np.max(abs(exterior@areas-1/MATCH_RADIUS)) < 1e-10
    return normals, edges, B, D, exterior, areas


class Collar:
    def __init__(self, mesh, eta, retention):
        normals, edges, self.B, D, self.exterior, self.areas = mesh
        self.nf, self.ne = self.B.shape
        self.nv = len(normals)
        self.eta = eta
        shift = sparse.diags([1.,-1.],[0,1],shape=(SHELLS-1,SHELLS))
        self.C = sparse.bmat([
            [sparse.kron(sparse.eye(SHELLS),self.B), None],
            [sparse.kron(shift,sparse.eye(self.ne)), sparse.kron(sparse.eye(SHELLS-1),D)]
        ],format="csr")
        self.Ct = self.C.T.tocsr()
        self.surface = slice((SHELLS-1)*self.nf, SHELLS*self.nf)
        left, right = edges.T
        tangent, radial, lengths, rays = [], [], [], []
        for k in reversed(range(SHELLS)):
            p = predicted(normals,eta,retention,k)
            sizes = np.linalg.norm(p,axis=1)
            assert min(sizes)>0  # rotated vertices avoid a zero between sampled cuts
            z = states(p)
            if rays:
                overlap = (rays[-1].conj()*z).sum(axis=1)
                z *= np.exp(-1j*np.angle(overlap))[:,None]
                radial.append(np.sqrt(lengths[-1]*sizes)*(rays[-1].conj()*z).sum(axis=1))
            tangent.append(np.sqrt(sizes[left]*sizes[right])*(z[left].conj()*z[right]).sum(axis=1))
            rays.append(z)
            lengths.append(sizes)
        bonds = np.concatenate(tangent+radial)
        self.phase = np.angle(bonds)
        self.visibility = abs(bonds)
        self.rim = slice(0,self.ne)
        self.vertex_gradient = sparse.bmat([
            [sparse.kron(sparse.eye(SHELLS),D)],
            [sparse.kron(-shift,sparse.eye(self.nv))]
        ],format="csr")
        assert (self.C @ self.vertex_gradient).nnz == 0

    def contact(self,g):
        stiffness = g*self.visibility
        stiffness[self.rim] = RIM*self.visibility[self.rim]
        return stiffness

    def pieces(self,a,g):
        f = principal(self.C@a)
        tail_gradient = self.exterior@f[self.surface]
        tail = .5*(f[self.surface]@tail_gradient)
        slip = a-self.phase
        weights = self.contact(g)
        contact = weights@(2*np.sin(slip/2)**2)
        df = np.sin(f/4)
        df[self.surface] += tail_gradient
        return np.array([W(f).sum(),contact,tail]), self.Ct@df + weights*np.sin(slip)

    def objective(self,a,g):
        terms, grad = self.pieces(a,g)
        return terms.sum(), grad

    def topology(self,a):
        f = principal(self.C@a)
        assert np.max(abs(f)) < np.pi-1e-9
        angular = f[:SHELLS*self.nf].reshape(SHELLS,self.nf)
        sides = f[SHELLS*self.nf:].reshape(SHELLS-1,self.ne)
        q = angular.sum(axis=1)/(2*np.pi)
        m = (np.diff(angular,axis=0)+(self.B@sides.T).T)/(2*np.pi)
        assert np.max(abs(m-np.rint(m))) < 1e-10 and np.max(abs(q-np.rint(q))) < 1e-10
        assert np.allclose(np.diff(q),m.sum(axis=1))
        return np.rint(q).astype(int), np.sum(abs(np.rint(m)),axis=1).astype(int)

    def seeds(self):
        angular = self.phase[:SHELLS*self.ne].reshape(SHELLS,self.ne)
        radial = np.zeros((SHELLS-1)*self.nv)
        # Radial-parallel gauge. A phase-flat input aligned with the outer ray
        # may have vertex-dependent phases; its outer angular links implement these.
        for name, layers in [('through', [0]*SHELLS), ('preferred', list(range(SHELLS))),
                             ('flat', [SHELLS-1]*SHELLS)]:
            yield name, np.r_[angular[layers].ravel(),radial]

    def checks(self,level):
        rng = np.random.default_rng(1102)
        a = self.phase+rng.normal(0,.009,len(self.phase))
        perturb = rng.normal(size=len(a))
        perturb /= np.linalg.norm(perturb)
        h = 1e-5
        g = .37
        energy, derivative = self.objective(a,g)
        fd = (self.objective(a+h*perturb,g)[0]-self.objective(a-h*perturb,g)[0])/(2*h)
        gradient_error = abs(fd-derivative@perturb)
        twist = self.vertex_gradient@rng.uniform(-np.pi,np.pi,SHELLS*self.nv)
        temp = self.phase
        self.phase = temp+twist
        changed = self.objective(a+twist,g)[0]
        self.phase = temp
        print(f"controls level={level} eta={self.eta:+} derivative_error={gradient_error:.3g} "
              f"gauge_error={abs(energy-changed):.3g}",flush=True)
        assert gradient_error < 1e-7 and abs(energy-changed) < 1e-10


def compare(level,eta,mesh,retention,values):
    collar = Collar(mesh,eta,retention)
    collar.checks(level)
    seeds = dict(collar.seeds())
    for g in values:
        for name,a in seeds.items():
            e0 = collar.objective(a,g)[0]
            fit = minimize(collar.objective,a,args=(g,),method="L-BFGS-B",jac=True,
                           options={"ftol":5e-15,"gtol":2e-7,"maxiter":20000,"maxfun":100000})
            terms,grad = collar.pieces(fit.x,g)
            q, defects = collar.topology(fit.x)
            torque = np.max(abs(grad))
            print(f"relax level={level} eta={eta:+} g={g:g} seed={name} "
                  f"initial={e0:.10g} energy={terms.sum():.10g} collar={terms[0]:.9g} "
                  f"contact={terms[1]:.9g} exterior={terms[2]:.9g} cuts={q.tolist()} "
                  f"absolute_cell_counts={defects.tolist()} max_torque={torque:.3g} "
                  f"iterations={fit.nit} success={fit.success}",flush=True)
            assert fit.fun <= e0+1e-10
            # A local search can stop at the principal-angle cut; report its
            # torque and success, without treating that iterate as stationary.


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--level',type=int,default=2)
    parser.add_argument('--g',nargs='+',type=float,default=[.0001,.03,.3,8.])
    args = parser.parse_args()
    mesh = angular_mesh(args.level)
    retention = calibration()
    for eta in [1,-1]:
        compare(args.level,eta,mesh,retention,args.g)
