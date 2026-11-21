"""Bare passive conversion work on a relaxed, fixed-intake bent 4D mouth.

The carrier and scalar order equilibrate in sequence. The physical energy
readout does not contain the integrating factor of the stationary order PDE.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

from mouth_bend_drive_control import DriveControl
from mouth_carrier_lateral_response import prepare
from mouth_curved_phase_stability import weighted_gram
from mouth_flared_receiver import emit
from mouth_free_order_diffuse import NU, positive_factor
from mouth_order_wave_pulse import digest


class BentOrder:
    def __init__(self, carrier, seed, width):
        self.carrier = carrier
        self.width = width
        self.seed = seed
        self.masks = []
        self.factors = []
        dom = carrier.dom
        phase = np.angle(dom.N @ carrier.base)
        wt = dom.weight * np.exp(-phase / NU)
        chi = dom.N @ seed
        well_curvature = (1 - 6*chi*(1-chi)) / width**2
        stiff = weighted_gram(carrier.dr, wt) + weighted_gram(carrier.dw, wt)
        for degree in range(carrier.modes):
            mask = dom.free & ((dom.points[:, 0] > 0) | (degree == 0))
            mass = weighted_gram(dom.N, wt*(well_curvature + degree*(degree+1)/carrier.r**2))
            factor = positive_factor((stiff + mass)[mask][:, mask])
            if factor is None:
                raise ValueError(f'Nonpositive straight order Hessian at degree {degree}')
            self.masks.append(mask)
            self.factors.append(factor)

    def sample(self, field, center=True):
        s = self.carrier
        # Only the monopole admits a spatial constant. Centering its derivatives
        # at 1 avoids cancellation throughout the uniform feeding material.
        centered = field.copy()
        if center:
            centered[:, 0] -= 1
        return s.dom.N @ field, s.dr @ centered, s.dw @ centered

    def gradient(self, values, radial, axial, cos, poly, dpoly, br, bw, jac):
        ang = (values @ dpoly) / self.carrier.r
        dr = radial @ poly
        dx = (cos*dr + (1-cos*cos)*ang) / jac
        return dr-br*dx, ang, axial@poly-bw*dx

    def evaluate(self, field, phase_field, offset, mobility, measure=False, linearize_at=None):
        s = self.carrier
        br, bw = s.sr @ offset, s.sw @ offset
        c, cr, cw = self.sample(field, center=linearize_at is None)
        psi, pr, pw = phase_field
        vr, vw, value = (np.zeros_like(c) for _ in range(3))
        power = 0.
        force = np.zeros(2)
        strain_energy = 0.
        for cos, weight, y, dy in zip(s.mu,s.weights,s.y,s.dy,strict=True):
            jac = 1+cos*br
            angular = 1-cos*cos
            carrier_value = psi@y
            if np.any(carrier_value.real <= 0) or np.any(jac <= 0):
                raise ValueError('Invalid single-phase or geometry branch')
            order = c@y
            radial,ang,axial = self.gradient(c,cr,cw,cos,y,dy,br,bw,jac)
            if measure:
                vg = self.gradient(psi,pr,pw,cos,y,dy,br,bw,jac)
                ur,ua,uw = (np.imag(v/carrier_value) for v in vg)
                material_rate = ur*radial + angular*ua*ang + uw*axial
                chemical_potential = -material_rate/mobility
                volume = s.dom.weight*weight*jac
                power += volume @ (material_rate**2/mobility)
                force += [volume @ (chemical_potential*(cos*radial+angular*ang)),
                          volume @ (chemical_potential*axial)]
                strain_energy += volume @ ((radial**2+angular*ang**2+axial**2+
                                            (order*(1-order)/self.width)**2)/2)
                continue
            phase_weight = weight * np.exp(-np.angle(carrier_value)/NU)
            h = br*radial + bw*axial
            vr += (phase_weight*(jac*radial-cos*h))[:,None]*y
            vw += (phase_weight*jac*axial)[:,None]*y
            if linearize_at is None:
                well = order*(1-order)*(1-2*order)/self.width**2
            else:
                base_order = linearize_at@y
                well = (1-6*base_order*(1-base_order))*order/self.width**2
            value += ((phase_weight*angular*(jac*ang-h)/s.r)[:,None]*dy +
                      (phase_weight*jac*well)[:,None]*y)
        if measure:
            return dict(passive_loss=power,force_on_carrier_xw=force.tolist(),
                        stored_bare_order_energy=strain_energy*NU/mobility)
        wt = s.dom.weight[:,None]
        residual = s.dr.T@(wt*vr) + s.dw.T@(wt*vw) + s.dom.N.T@(wt*value)
        return residual*(2*np.arange(s.modes)+1)

    def run(self, phase_field, offset, mobility):
        field = np.zeros((len(self.seed),self.carrier.modes))
        field[:,0] = self.seed
        free = np.column_stack(self.masks)

        def precondition(residual):
            update = np.zeros_like(field)
            for degree,(mask,factor) in enumerate(zip(self.masks,self.factors,strict=True)):
                update[mask,degree] = factor.solve(residual[mask,degree])
            return update

        for iteration in range(18):
            residual = self.evaluate(field,phase_field,offset,mobility)
            rhs = -precondition(residual)[free]
            norm = np.linalg.norm(rhs, np.inf)
            if norm < 3e-12:
                return field, norm, iteration
            base_order = self.carrier.dom.N @ field

            def action(vector):
                perturbation = np.zeros_like(field)
                perturbation[free] = vector
                tangent = self.evaluate(perturbation,phase_field,offset,mobility,linearize_at=base_order)
                return precondition(tangent)[free]

            # Right-hand side and operator have both been acted on by the same
            # straight-state inverse. The Jacobian retains the bent order profile.
            matrix = LinearOperator((len(rhs),len(rhs)), matvec=action, dtype=float)
            update,info = gmres(matrix, rhs, rtol=1e-7, atol=2e-14, restart=24, maxiter=10)
            if info!=0:
                raise ArithmeticError(f'Order Newton inner solve failed: {info}, residual {norm}')
            fraction = 1.
            while fraction >= 1/128:
                candidate = field.copy()
                candidate[free] += fraction*update
                test = self.evaluate(candidate,phase_field,offset,mobility)
                if np.max(abs(precondition(test))) < norm:
                    field = candidate
                    break
                fraction /= 2
            else:
                raise ArithmeticError(f'No descent of preconditioned order residual {norm}')
        raise ArithmeticError(f'Order Newton iteration did not converge: residual {norm}')


def run(mesh, amplitudes, modes, output):
    folder = Path('_scratch/s11c')
    source = folder / f'mouth_flared_delayed_n{mesh}'
    seed_path = folder / f'mouth_free_order_delayed_n{mesh}_d04/order_seed1.55.npz'
    archive = Path('research/pde_ledger_v3/_measurements/mouth_free_order_depth.json')
    if digest(seed_path) != json.loads(archive.read_text())['diffuse_orders_sha256'][str(seed_path)]:
        raise ValueError('Deep order source changed')
    meta = json.loads((seed_path.parent/'stdout').read_text().splitlines()[0])
    if meta['source_sha256'] != digest(source/'field.npz'):
        raise ValueError('Unmatched order and carrier')
    dom, base = prepare(SimpleNamespace(source=source))
    with np.load(seed_path) as stored:
        seed = stored['chi']
    s = DriveControl(dom, base, modes)
    conv = BentOrder(s, seed, meta['delta'])
    emit(test='settings',mesh=mesh,modes=modes,order=str(seed_path),order_sha256=digest(seed_path),
         delta=meta['delta'],mobility=meta['mobility'],nu=NU,sigma=meta['sigma'])
    baseline = None
    for amplitude in amplitudes:
        off = np.array([0., amplitude])
        solution = s.solve(off,.25)
        psi = conv.sample(s.last_field)
        chi, error, it = conv.run(psi,off,meta['mobility'])
        result = conv.evaluate(chi,psi,off,meta['mobility'],measure=True)
        max_change = np.max(abs(chi[:,0]-seed))
        order_modes = np.max(abs(chi),axis=0)
        if baseline is None:
            assert amplitude==0
            baseline = result
        emit(test='matched_intake',bend=amplitude,feed=solution['feed'],bulk=solution['bulk'],
             bias=solution['phase_bias'],field_residual=solution['residual'],order_correction=error,
             order_iterations=it,max_monopole_change=max_change,order_mode_norms=order_modes.tolist(),
             **result,saved_power=baseline['passive_loss']-result['passive_loss'],
             stored_energy_change=result['stored_bare_order_energy']-baseline['stored_bare_order_energy'])
        if output:
            output.mkdir(parents=True,exist_ok=True)
            dest=output/f"bend_{amplitude:+.5f}.npz"
            np.savez_compressed(dest,psi=s.last_field,chi=chi)
            emit(test="field",bend=amplitude,path=str(dest),sha256=digest(dest))
        if abs(solution['feed']-.25)>2e-9 or abs(solution['bulk']-.25)>2e-9 or max_change>.05:
            raise ArithmeticError('Unmatched intake or wrong order branch')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesh',type=int,choices=(34,48),required=True)
    p.add_argument('--bends',type=float,nargs='+',default=[0.,.01,-.01,.005,-.005])
    p.add_argument('--modes',type=int,default=3)
    p.add_argument('--fields',type=Path)
    a=p.parse_args()
    run(a.mesh,a.bends,a.modes,a.fields)
