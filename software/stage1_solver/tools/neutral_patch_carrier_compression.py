"""A cold neutral patch under uniform strain along the three brane directions.

Only the 1% density increase sets the pressure scale. Volumes count the fourth
spatial direction explicitly. Each of the two gap orientations is a hypothesis.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import root
from scipy.spatial import ConvexHull

from mouth_segment_contact_hold import PITCH
from neutral_patch_conversion import NEUTRAL
from neutral_patch_sites import COUNT, GAP, Patch, emit, morse
from neutral_patch_turning_neighbor import AREA_PER_SITE

CALIBRATION = 1.01
DENSITIES = (.98, .99, .995, 1., 1.005, CALIBRATION, 1.02, 1.05)


class Cell:
    def __init__(self):
        data = json.loads(NEUTRAL.read_text())['saved_measurements']['frames1321.json']
        self.xy = np.array(data['contents']['central']['rested_xy'])
        self.xy -= self.xy.mean(axis=0)
        wall = []
        for start in (0, COUNT):
            wall.extend(start+ConvexHull(self.xy[start:start+COUNT]).vertices)
        self.wall = np.sort(wall)
        self.free = np.setdiff1d(np.arange(2*COUNT), self.wall)
        self.law = Patch('central')
        self.incidence = np.eye(2*COUNT)[self.law.j]-np.eye(2*COUNT)[self.law.i]
        self.free_incidence = self.incidence[:,self.free]

    def measure(self, density, gap_in_brane, relax, seed=None):
        stretch = density**(-1/3)
        gap = GAP*(stretch if gap_in_brane else 1.)
        xy = self.xy*stretch
        law = self.law

        def forces():
            delta = self.incidence@xy
            r = np.sqrt((delta**2).sum(axis=1)+law.cross*gap**2)
            e,t = morse(r,law.depth,law.alpha,law.length)
            grad = self.incidence.T @ ((t/r)[:,None]*delta)
            return delta,r,e,t,grad

        def gradient(flat):
            xy[self.free] = flat.reshape(-1,2)
            return forces()[-1][self.free].ravel()

        def hessian(flat):
            xy[self.free] = flat.reshape(-1,2)
            delta,r,_,t,_ = forces()
            y = np.exp(-law.alpha*(r-law.length))
            curvature = 2*law.depth*law.alpha**2*y*(2*y-1)
            unit_xy = delta/r[:,None]
            block = ((t/r)[:,None,None]*np.eye(2) +
                     (curvature-t/r)[:,None,None]*unit_xy[:,:,None]*unit_xy[:,None,:])
            b = self.free_incidence
            return np.einsum('bi,bkl,bj->ikjl',b,block,b,optimize=True).reshape(2*len(self.free),-1)

        if relax:
            start = xy[self.free] if seed is None else seed
            sol = root(gradient, start.ravel(), jac=hessian, tol=1e-11)
            free_error = float(np.max(abs(gradient(sol.x))))
            least_mode = float(np.linalg.eigvalsh(hessian(sol.x))[0])
            if free_error>1e-8 or least_mode<=0:
                raise ValueError(f'Unstable or unbalanced inner sites: {density}, {free_error}, {least_mode}')
        else:
            free_error = least_mode = None
        delta,r,e,t,grad = forces()
        gap_work = np.sum((t/r)[law.cross])*gap**2 if gap_in_brane else 0.
        drivers = self.wall if relax else np.arange(2*COUNT)
        applied_work = np.sum(grad[drivers]*xy[drivers])+gap_work
        bond_work = np.sum(t/r*np.sum(delta**2,axis=1))+gap_work
        np.testing.assert_allclose(applied_work,bond_work,atol=5e-8,rtol=1e-9)
        # dV/d(log stretch) = 3 V. Raw pressure is in E_micro / V4_rest.
        return dict(energy=float(e.sum()),raw_pressure=float(-density*applied_work/3),
                    free_error=free_error,least_mode=least_mode),xy[self.free].copy()


def carrier_pressure(n):
    # p - p_rest = n W' - W, for W=(n**5-5*n+4)/20 at rest n=1.
    return np.expm1(5*np.log(n))/5


def experiment(cell, gap_in_brane, relax, volume):
    def read(n):
        return cell.measure(n,gap_in_brane,relax)[0]
    rows = {n:read(n) for n in DENSITIES}
    baseline = rows[1.]['raw_pressure']
    calibration = rows[CALIBRATION]['raw_pressure']-baseline
    if calibration<=0:
        raise ValueError('Compression calibration is not positive')
    scale = carrier_pressure(CALIBRATION)/calibration
    errors = []
    for n,r in rows.items():
        # Differentiate the equilibrated energy in the physical volume, then
        # compare with the force on the moving boundaries, including prestress.
        J = 1/n
        for step in (2e-5,1e-5):
            fd = (read(1/(J-step))['energy']-read(1/(J+step))['energy'])/(2*step)
            errors.append(abs(fd-r['raw_pressure']))
            np.testing.assert_allclose(fd,r['raw_pressure'],atol=2e-6,rtol=1e-7)
        prediction = scale*(r['raw_pressure']-baseline)
        target = carrier_pressure(n)
        r.update(n=n,carrier=float(target),prediction=float(prediction),
                 fractional_error=float(prediction/target-1) if n!=1. else None)
    # Dimensionless curvature for log density. A constant multiplier cannot
    # change it. The carrier has p''(0)/p'(0)=5, with eta=log(n).
    curvature = []
    for h in (.002,.001,.0005):
        plus, minus = [read(np.exp(k*h))['raw_pressure']-baseline for k in (1,-1)]
        slope = (plus-minus)/(2*h)
        curvature.append(dict(log_step=h,log_slope=slope,curvature_ratio=(plus+minus)/(h*h*slope)))
    if relax:
        for n in (DENSITIES[0],DENSITIES[-1]):
            _, start = cell.measure(n,gap_in_brane,True)
            check, back = cell.measure(1.,gap_in_brane,True,seed=start/n**(-1/3))
            _, ref = cell.measure(1.,gap_in_brane,True)
            np.testing.assert_allclose(back,ref,atol=2e-8)
    emit(test='pressure_comparison',stack='brane' if gap_in_brane else 'w',interior='relaxed' if relax else 'affine',
         calibration_density=CALIBRATION, energy_per_rest_four_volume=scale,
         patch_energy_unit=scale*volume, patch_force_unit=scale*volume/PITCH,
         virtual_work_max_abs_error=max(errors), curvature=curvature, rows=list(rows.values()))


if __name__=='__main__':
    cell = Cell()
    # Two triangular layers. Assign one pitch of depth to each plane, and
    # one pitch to the unrepresented direction. All lengths are in mouth radii.
    volume = COUNT*AREA_PER_SITE*2*GAP*PITCH**4
    emit(test='cell',sites=2*COUNT,hull_sites=cell.wall.tolist(),interior_sites=cell.free.tolist(),
         rest_carrier_density=1.,pitch=PITCH,rest_four_volume=volume,represented_number=volume,
         represented_number_per_site=volume/(2*COUNT),assumed_area_in_pitch_units=COUNT*AREA_PER_SITE,
         density_steps=DENSITIES,calibration=CALIBRATION)
    for gap_in_brane in (False,True):
        for relax in (False,True):
            experiment(cell,gap_in_brane,relax,volume)
