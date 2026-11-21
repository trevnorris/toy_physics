"""Measure tangent shear stiffness in a flat curl medium and in ordered bilayers.

The slow alignment is a shift of the upper plane by GAP*s. A weak additional
shift probes the tangent along or across it. The two candidates retain the
neutral-patch interactions and order scaling used before the mouth test.
"""
import argparse
from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize

from mouth_shear_ringdown import digest, prepared, quadratic, tensor
from neutral_patch_periodic_compression import BASIS
from neutral_patch_sites import BETWEEN, GAP, PITCH_UNIT, WITHIN, emit, morse
from neutral_patch_structural_endpoints import PITCH

STRAINS = (1/32, 1/64, 1/128)
DIRECTIONS = np.eye(2)


def collar_orders():
    args = SimpleNamespace(radius=1., mesh=34, state=None, outer=11., feed_spacing=.1)
    domain, chi, intake, path = prepared(args)
    depths, radial, nodes = domain.order_blocks[2]
    def sample(w, v):
        i = np.clip(np.searchsorted(depths, w)-1, 0, len(depths)-2)
        j = np.clip(np.searchsorted(radial, v)-1, 0, len(radial)-2)
        x = quadratic(2*(w-depths[i])/np.diff(depths)[i]-1)[0]
        y = quadratic(2*(v-radial[j])/np.diff(radial)[j]-1)[0]
        return np.einsum('qi,qi->q', tensor(x, y), chi[nodes[i*(len(radial)-1)+j]])
    w = np.array([.60, .85, 1.10])
    orders = sample(w, np.zeros_like(w))
    emit(test='collar', radius=1., intake=intake, depth=w.tolist(), wall_order=orders.tolist(),
         source=str(path), sha256=digest(path))
    if not np.all((orders > 0) & (orders <= 1.000001)):
        raise ValueError('Invalid wall order')
    return orders


def curl_patch(orders):
    # Spatial norm of the material two-form, pulled to the lab by F^{-1}.
    # Base and probe can independently produce either displacement gradient.
    for order in orders:
        for s in (0., *STRAINS, *(-v for v in STRAINS)):
            deform = np.eye(4)
            deform[0, 3] = s
            inv = np.linalg.inv(deform)
            def energy(gradient):
                curl = gradient-gradient.T
                physical = inv.T@curl@inv
                return order*np.linalg.det(deform)*np.sum(physical**2)/4
            results = []
            for axis in (0, 1):
                first = np.zeros((4,4))
                first[axis, 3] = 1.
                for label, gradient in [('tangent_by_depth', first),
                                        ('depth_by_tangent', -first.T),
                                        ('zero_curl', first+first.T)]:
                    step = .0002
                    modulus = (energy(step*gradient)+energy(-step*gradient)-2*energy(0*gradient))/step**2
                    results.append([axis, label, modulus])
            emit(test='curl', order=order, prestrain=s, moduli=results)
            for _, label, modulus in results:
                np.testing.assert_allclose(modulus, 0 if label=='zero_curl' else order, atol=2e-14)


class Crystal:
    def __init__(self, law, shell):
        self.law = law
        self.area = PITCH**2*np.linalg.det(BASIS)
        ij = np.indices((2*shell+1,)*2).reshape(2,-1).T-shell
        xy = PITCH*ij@BASIS
        keep = np.linalg.norm(xy,axis=1) <= PITCH*(shell+.01)
        self.xy, self.ij = xy[keep], ij[keep]

    def bonds(self, displacement, cross=True):
        xy = self.xy+displacement if cross else self.xy[np.any(self.ij!=0,axis=1)]
        square = np.einsum('ij,ij->i', xy, xy)
        if cross and self.law=='layered':
            rad = np.sqrt(GAP**2+square**2/PITCH_UNIT**2)
            p, dp = 2*square/PITCH_UNIT**2, 2/PITCH_UNIT**2
        else:
            rad = np.sqrt(square+(GAP**2 if cross else 0.))
            p, dp = 1., 0.
        pars = BETWEEN if cross else WITHIN
        en, prime = morse(rad, *pars)
        depth, alpha, rest = pars
        y = np.exp(alpha*(rest-rad))
        second = 2*depth*alpha**2*y*(2*y-1)
        k = prime/rad
        radial = (second-k)*p**2/rad**2 + 2*k*dp
        grad = (k*p)[:,None]*xy
        hess = (k*p)[:,None,None]*np.eye(2) + radial[:,None,None]*xy[:,:,None]*xy[:,None,:]
        return en, grad, hess

    def reading(self, shift):
        en, grad, hess = self.bonds(shift)
        return en.sum(), grad.sum(axis=0), hess.sum(axis=0)

    def register(self):
        def objective(u):
            en, grad, _ = self.reading(u)
            return en, grad
        starts = np.array([[0,0],[1/3,1/3],[-1/3,-1/3],[.5,0],[.14,.12],[.38,-.08]])
        rows = []
        for frac in starts:
            res = minimize(objective, PITCH*frac@BASIS, method='BFGS', jac=True,
                           options={'gtol':1e-12})
            e, grad, hess = self.reading(res.x)
            rows.append((e, np.linalg.eigvalsh(hess)[0], np.linalg.norm(grad), res.x))
        stable = [row for row in rows if row[1]>0 and row[2]<1e-8]
        if not stable:
            raise ValueError('No stationary registry with positive shift stiffness')
        lowest = min(row[0] for row in stable)
        # Prefer an exact symmetry site where candidates are energetically tied.
        chosen = next(row for row in stable if row[0] <= lowest+1e-10)
        emit(test='registry', law=self.law, chosen=chosen[3].tolist(),
             trials=[dict(energy=e, least_shift_stiffness=k, gradient=g, offset=x.tolist())
                     for e,k,g,x in rows])
        return chosen[3]

    def stability(self, shift, order, size=16):
        cross = self.bonds(shift)[2]*order/.95
        same = self.bonds(None, cross=False)[2]
        k = 2*np.pi*np.indices((size,size)).reshape(2,-1).T/size
        phases = np.exp(1j*k@self.ij.T)
        onsite = (np.einsum('kb,bij->kij',1-phases[:,np.any(self.ij!=0,axis=1)].real,same)
                  + cross.sum(axis=0))
        mixed = -np.einsum('kb,bij->kij',phases,cross)
        matrix = np.concatenate((np.concatenate((onsite, mixed),axis=2),
                                 np.concatenate((mixed.conj().transpose(0,2,1),onsite),axis=2)),axis=1)
        eig = np.linalg.eigvalsh(matrix)
        np.testing.assert_allclose(eig[0,:2],0,atol=1e-10)
        return min(eig[0,2], eig[1:].min())


def bilayer(law, orders):
    crystal = Crystal(law, 8)
    extended = Crystal(law, 10)
    origin = crystal.register()
    h0 = crystal.reading(origin)[2]
    np.testing.assert_allclose(h0, np.eye(2)*h0.trace()/2, atol=2e-6)
    scale = GAP/crystal.area  # effective test volume = area*gap; dipole shift = gap*strain
    for order in orders:
        factor = scale*order/.95
        for s in (0., *STRAINS, *(-v for v in STRAINS)):
            shift = origin + GAP*s*DIRECTIONS[0]
            e, grad, hess = crystal.reading(shift)
            _, _, wider = extended.reading(shift)
            margin = crystal.stability(shift, order)
            emit(test='tangent', law=law, order=order, prestrain=s, stress=(factor/GAP*grad).tolist(),
                 energy_per_cell=e*order/.95, modulus=(factor*hess).tolist(),
                 least_sampled_internal_stiffness=margin, reciprocal_grid=16)
            np.testing.assert_allclose(hess, wider, atol=3e-10, rtol=0)
            if margin <= 0:
                raise ArithmeticError('A sampled mode is unstable at the stated preparation')
            for probe in (.0002,.0001):
                read = []
                for v in DIRECTIONS:
                    em, gm, _ = crystal.reading(shift-probe*GAP*v)
                    ep, gp, _ = crystal.reading(shift+probe*GAP*v)
                    stress_tangent = v@(gp-gm)/(2*probe*GAP)
                    work_tangent = (ep+em-2*e)/(probe*GAP)**2
                    read.append((factor*stress_tangent,factor*work_tangent))
                emit(test='probe', law=law, order=order, prestrain=s, amplitude=probe,
                     traction_and_energy=read)
                np.testing.assert_allclose(read, np.repeat((factor*hess.diagonal())[:,None],2,axis=1),
                                           atol=3e-6, rtol=4e-6)
        for s in STRAINS:
            positive = crystal.reading(origin+GAP*s*DIRECTIONS[0])[2]
            negative = crystal.reading(origin-GAP*s*DIRECTIONS[0])[2]
            even = ((positive+negative)/2 - h0)/(h0.trace()/2*s*s)
            odd = (positive-negative)/(h0.trace()*s)
            emit(test='strain_response',law=law,order=order,prestrain=s,
                 even_quadratic=even.tolist(),odd_linear=odd.tolist())


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('material',choices=('curl','central','layered'))
    material = parser.parse_args().material
    orders = collar_orders()
    if material=='curl':
        curl_patch(orders)
    else:
        bilayer(material, orders)
