"""Deform the reference wall along one fixed compact eight-band variation.

Keep the two wall ends fixed. Re-solve the stationary carrier at fixed intake
and the ground wave in the moved well. Forces are conjugate to the largest
radial displacement, in units of the mouth radius.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import yaml
from scipy.special import expit

from mouth_curved_phase import RATE
from mouth_deep_attachment_traction import BANDS
from mouth_deep_support_calibration import band_work, quadrature
from mouth_flared_receiver import EXIT, FlaredDomain, ROUNDING, emit
from mouth_trapped_wave_ports import ARCHIVE as WAVE_LAW
from mouth_wave_flare_work import Wave, radial_profile

MEASUREMENTS = Path(__file__).resolve().parents[3] / 'research/pde_ledger_v3/_measurements'


def archived_rows(record):
    literal = record['literal_stdout']
    if sha256(literal.encode()).hexdigest() != record['stdout_sha256']:
        raise ValueError('Archived output hash changed')
    return [json.loads(line) for line in literal.splitlines()]


def fixed_coefficients():
    affine = json.loads((MEASUREMENTS / 'mouth_free_equilibrium_inputs.json').read_text())
    radial = json.loads((MEASUREMENTS / 'mouth_occupied_radial_balance.json').read_text())
    run = radial['runs']['occupied_radial_n34_01']
    if not (affine['guard_exit'] == affine['worker_exit'] == run['guarded_exit'] ==
            run['worker_exit'] == 0 and affine['limits_verified'] and run['limits_verified']):
        raise ValueError('Missing successful guarded baseline')
    rows = archived_rows(affine)
    train, = [r for r in rows if r['test']=='affine_work' and r['mesh']==34 and r['method']=='volume']
    geometry, = [r for r in rows if r['test']=='preparation' and r['mesh']==34]
    np.testing.assert_array_equal(geometry['bands'], BANDS)
    wave, = [r for r in archived_rows(run) if r['test']=='wave_state']
    dual = np.array(train['dual_weights'])
    volume = np.array(wave['projected_volumes'])
    coef = dual/volume
    coef /= np.max(abs(coef))
    projected = (coef*volume) @ np.array(geometry['design'])
    emit(test='direction', mesh=34, method='volume', bands=BANDS, coefficients=coef.tolist(),
         original_nuisance_work=projected.tolist(),
         carrier_work=float((coef*volume) @ train['carrier_pressure']))
    np.testing.assert_allclose(projected, 0., atol=2e-13)
    return coef


def bump(w, lo, hi):
    t = np.clip((w-lo)/(hi-lo), 0., 1.)
    return 16*t*t*(1-t)**2, 32*t*(1-t)*(1-2*t)/(hi-lo)


class VariedMouth(FlaredDomain):
    def __init__(self, mesh, amplitude, coef):
        self.amplitude, self.coef = amplitude, coef
        super().__init__(mesh, 'cubic', 10, 9.)

    def assemble(self, grids):
        # Resolve the short wall patches within the physical chart before
        # computing the fixed-intake state. No element straddles a patch edge.
        axial, radial = grids[2]
        patches = [np.linspace(lo, hi, 1+max(4, int(np.ceil((hi-lo)*len(radial)*3))))
                   for lo, hi in BANDS]
        knots = np.unique(np.round(np.concatenate([axial, *patches]), 13))
        grids[2] = knots, radial
        super().assemble(grids)

    def radius(self, w):
        r, slope = super().radius(w)
        for c, (lo, hi) in zip(self.coef, BANDS, strict=True):
            b, db = bump(w, lo, hi)
            r = r+self.amplitude*c*b
            slope = slope+self.amplitude*c*db
        return r, slope


class VariedWell(Wave):
    def __init__(self, step, coef):
        params = yaml.safe_load(WAVE_LAW.read_text())['parameters']
        super().__init__(1.,1.,step,5.,params)
        # Integrate the continuous shape against the linear nodal hats.
        # Division by nodal cell width is the archived axial force projection.
        z = self.w.ravel()
        self.projected = []
        for lo, hi in BANDS:
            profile = np.zeros_like(z)
            w, dw = quadrature(z, lo, hi)
            b, _ = bump(w, lo, hi)
            index = np.searchsorted(z,w)-1
            right = (w-z[index])/step
            np.add.at(profile,index,dw*b*(1-right)/step)
            np.add.at(profile,index+1,dw*b*right/step)
            self.projected.append(profile)
        self.projected = np.array(self.projected)
        self.shape = coef @ self.projected
        self.base = radial_profile(self.w)[0]
        self.width = params['delta_r']/2
        depth_width = params['delta_parallel']/2
        self.axial = expit((self.w-ROUNDING)/depth_width)*expit((EXIT-self.w)/depth_width)

    def gate(self, amplitude):
        radial = expit((self.base+amplitude*self.shape-self.r)/self.width)
        gate = self.axial*radial
        derivative = gate*(1-radial)*self.shape/self.width
        return gate, derivative, radial

    def read_force(self, amplitude):
        omega, u, operator, derivative, radial = self.ground(amplitude)
        residual = np.linalg.norm(operator@u-omega**2*u)
        if residual>1e-7 or omega>=self.p['mu_out']:
            raise ArithmeticError('Ground mode failed')
        nodal = self.contrast*(u.reshape(radial.shape)**2*self.axial*radial*(1-radial)/self.width).sum(axis=0)/(2*omega)
        return omega, self.contrast*((u*u)@derivative.ravel())/(2*omega), self.projected@nodal


def force(mesh, amplitude, coef, wave):
    dom = VariedMouth(mesh, amplitude, coef)
    field, state = dom.solve(RATE)
    density = np.abs(dom.N@field)**2
    energy = (.5*np.vdot(field-1,dom.stiff@(field-1)).real
              + dom.weight@((density-1)**2*(density**3+2*density**2+3*density+4))/20)
    potential = energy-RATE*state['bias']
    area, volume, first_moment = np.zeros(3)
    area_work, volume_work, moment_work, surface_work, carrier_work = np.zeros(5)
    band_loads = []
    for c, (lo,hi) in zip(coef, BANDS, strict=True):
        w, dw = quadrature(dom.blocks[2][0], lo, hi)
        r, slope = dom.radius(w)
        b, _ = bump(w,lo,hi)
        projected_volume = 4*np.pi * dw@(r*r*b)
        curvature, surface, bulk = band_work(dom,field,lo,hi)
        area += 4*np.pi * dw@(r*r*np.sqrt(1+slope**2))
        volume += 4*np.pi/3 * dw@(r**3)
        first_moment += 4*np.pi/3 * dw@(w*r**3)
        area_work += c*projected_volume*curvature
        volume_work += c*projected_volume
        moment_work += 4*np.pi*c*dw@(w*r*r*b)
        carrier_work += c*projected_volume*bulk
        surface_work += c*projected_volume*surface
        band_loads.append(dict(curvature=curvature, carrier_surface=surface, carrier_volume=bulk,
                               mean_depth=(4*np.pi * dw@(w*r*r*b))/projected_volume, volume=projected_volume))
    omega, wave_work, wave_bands = wave.read_force(amplitude)
    for load, work in zip(band_loads,wave_bands,strict=True):
        load["wave_pressure_per_action"] = work/load["volume"]
    row = dict(test='force',mesh=mesh,q=amplitude, feed=state['feed'],bulk=state['bulk'],
               bias=state['bias'],residual=state['residual'], functional=potential,
               area=area,volume=volume,first_moment=first_moment,omega=omega,
               carrier_work=carrier_work,surface_work=surface_work,wave_work=wave_work,
               area_work=area_work,volume_work=volume_work,moment_work=moment_work,
               bands=band_loads)
    emit(**row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=int, required=True)
    parser.add_argument('--step', type=float, required=True)
    parser.add_argument('--amplitude', type=float, default=.0002)
    args = parser.parse_args()
    coef = fixed_coefficients()
    well = VariedWell(args.step,coef)
    delta = args.amplitude
    sample = {q: force(args.mesh,q,coef,well) for q in (0.,delta,-delta,delta/2,-delta/2)}
    center = sample[0.]
    for h in (delta, delta/2):
        def diff(key):
            return (sample[h][key]-sample[-h][key])/(2*h)
        readings = {}
        for observable, primitive, sign in [('carrier_work','functional',-1),('wave_work','omega',-1),
                                            ('area_work','area',1),('volume_work','volume',1),
                                            ('moment_work','first_moment',1)]:
            from_difference = sign*diff(primitive)
            readings[observable] = dict(central=center[observable],finite_difference=from_difference,
                                        directed_slope=diff(observable))
        emit(test='centered_checks',mesh=args.mesh,h=h,readings=readings)
        # Absolute intercept checks are necessary on dual-cancelled weak work.
        np.testing.assert_allclose(readings['carrier_work']['finite_difference'],center['carrier_work'],
                                   rtol=.025,atol=2e-8)
        np.testing.assert_allclose(readings['wave_work']['finite_difference'],center['wave_work'],
                                   rtol=.002,atol=2e-7)
        for key in ('volume_work','moment_work','area_work'):
            # A finite-step area difference contains the cubic term; its
            # convergence is measured explicitly by the half-strength pair.
            np.testing.assert_allclose(readings[key]['finite_difference'],center[key],
                                       atol=2.5e-4 if key=='area_work' else 2e-7)


if __name__=='__main__':
    main()
