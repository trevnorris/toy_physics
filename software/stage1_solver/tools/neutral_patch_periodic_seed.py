"""Wrap the original 24-site neutral patch onto a 4-by-3 triangular torus.

Particle positions and bonds are unchanged on wrapping, apart from creating the
periodic environment. Assignment to lattice labels only indexes periodic images.
"""
from itertools import product
import json

import numpy as np
from scipy.optimize import linear_sum_assignment, minimize

from neutral_patch_periodic_compression import BASIS, Bilayer
from neutral_patch_carrier_compression import CALIBRATION, DENSITIES, carrier_pressure
from neutral_patch_conversion import NEUTRAL
from neutral_patch_sites import COUNT, emit

BOX = np.diag([4,3]) @ BASIS
SHIFTS = np.array(list(product(range(-2,3),repeat=2))) @ BOX


def label_points(cell, points, stretch):
    offsets = []
    for layer in range(2):
        ref = cell.reference[layer*COUNT:(layer+1)*COUNT]
        xy = points[layer*COUNT:(layer+1)*COUNT]/stretch
        xy = (xy @ np.linalg.inv(BOX) % 1.) @ BOX
        delta = xy[None,:,None,:]-ref[:,None,None,:]+SHIFTS
        dist = np.sum(delta**2,axis=-1)
        near = dist.argmin(axis=-1)
        row,col = linear_sum_assignment(dist.min(axis=-1))
        offsets.extend(delta[row,col,near[row,col]])
    result = np.array(offsets)*stretch
    return result-result.mean(axis=0)


def equilibrate(cell, n, strained_gap, seed):
    # Relabel equivalent sites only before minimization; no pins or restraints.
    def objective(flat):
        offset = flat.reshape(-1,2)
        energy,gradient,_ = cell.read(n,strained_gap,offset)
        return energy*cell.count,gradient.ravel()
    sol = minimize(objective,seed.ravel(),jac=True,method='L-BFGS-B',
                   options=dict(gtol=1e-9,ftol=1e-15,maxiter=3000,maxls=50))
    x = sol.x.reshape(-1,2)
    x -= x.mean(axis=0)
    for _ in range(4):
        energy,g,p = cell.read(n,strained_gap,x)
        if abs(g).max()<=2e-9:
            break
        H = cell.read(n,strained_gap,x,True)
        # Suppress just the two common translations in the numerical correction.
        x[1:] -= np.linalg.solve(H[2:,2:],g[1:].ravel()).reshape(-1,2)
        x -= x.mean(axis=0)
    energy,g,p = cell.read(n,strained_gap,x)
    if np.max(abs(g))>1e-8:
        raise ArithmeticError(f'Periodic seed remains out of equilibrium: {sol.message}')
    return energy,p,x


def run():
    frame = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']['contents']['central']
    original = np.array(frame['rested_xy'])
    original -= original.mean(axis=0)
    cell = Bilayer(4,3,True,8.)
    extended = Bilayer(4,3,True,11.)
    for angle in (0.,np.pi/6):
        rotation = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        initial = label_points(cell,original@rotation,1.)
        rest,reference,x0 = equilibrate(cell,1.,False,initial)
        energy,_,p = cell.read(1.,False,np.zeros_like(initial))
        emit(test='preparation',angle_degrees=np.degrees(angle),original_sites=2*COUNT,
             periodic_sites=2*cell.count,affine_crystal_energy=energy,wrapped_relaxed_energy=rest,
             initial_max_displacement=float(abs(initial).max()),rest_max_displacement=float(abs(x0).max()),
             rest_pressure=reference,crystal_pressure=p)
        for strained_gap in (False,True):
            baseline = cell.read(1.,strained_gap,x0)[2]
            def read(n):
                return equilibrate(cell,n,strained_gap,x0*n**(-1/3))
            result=[]
            for n in DENSITIES:
                e,p,x = read(n)
                big_e,big_p,_ = equilibrate(extended,n,strained_gap,x)
                np.testing.assert_allclose([e,p],[big_e,big_p],atol=3e-9,rtol=0)
                assert np.max(np.linalg.norm(x,axis=1))<2.
                mode = cell.stability(n,strained_gap,x)
                delta_v=1e-5
                hi = read(1/(1/n+delta_v))[0]
                lo = read(1/(1/n-delta_v))[0]
                np.testing.assert_allclose((lo-hi)/(2*delta_v),p,atol=5e-8,rtol=1e-7)
                result.append(dict(n=n,raw_change=p-baseline,minimum_mode=mode))
            scale = carrier_pressure(CALIBRATION)/next(r['raw_change'] for r in result if r['n']==CALIBRATION)
            for r in result:
                n=r['n'];r['pressure_prediction']=scale*r['raw_change']
                r['fractional_error']=r['pressure_prediction']/carrier_pressure(n)-1 if n!=1 else None
            h=.0005
            plus,minus=[read(np.exp(s*h))[1]-baseline for s in (1,-1)]
            ratio=2*(plus+minus)/(h*(plus-minus))
            emit(test='wrapped_pressure',angle_degrees=np.degrees(angle),stack='brane' if strained_gap else 'w',
                 scale=scale,curvature_ratio=ratio,curvature_step=h,rows=result)


if __name__=='__main__':
    run()
