"""Held order surface with two curved and one straight tangential directions in 4D.

Normal coordinate s is in delta_B units and force per three-area in a_B units.
Use the Stage 1 stationary coherent solve on the same normal finite elements.
The neutral action is additive, with unexcited shear, wave, mixing and throat terms.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.special import expit

from mouth_curved_phase import Domain, quadratic
from mouth_curved_phase_stability import weighted_gram


class Collar:
    """Radial Q2 patch, uniform along the straight tangent, per unit pin area."""
    def __init__(self, curvature, elements_per_width):
        extent = 21
        elements = 2*extent*elements_per_width
        self.nodes = np.linspace(-extent, extent, 2*elements+1)
        q, w = leggauss(6)
        h = extent/elements
        self.s = (self.nodes[1::2,None]+h*q).ravel()
        self.measure = np.tile(h*w, elements)
        self.ratio = 1+curvature*self.s/2
        self.k = curvature
        self.weight = self.measure * self.ratio**2
        if np.any(self.ratio<=0):
            raise ValueError('Normal coordinates crossed the cylinder axis')
        rows = np.broadcast_to(np.arange(len(self.s)).reshape(elements,6,1),(elements,6,3))
        cols = np.broadcast_to(2*np.arange(elements)[:,None,None]+np.arange(3),rows.shape)
        def matrix(shape):
            values = np.broadcast_to(shape, rows.shape)
            return sparse.csr_matrix((values.ravel(),(rows.ravel(),cols.ravel())),
                                     shape=(len(self.s),len(self.nodes)))
        shapes, derivs = quadratic(q)
        self.N = matrix(shapes)
        self.grad = matrix(derivs/h)
        self.stiff = weighted_gram(self.grad, self.weight)
        self.feed = self.nodes==self.nodes[0]
        self.bulk = self.nodes==self.nodes[-1]
        self.free = ~(self.feed|self.bulk)

    def order_profile(self):
        chi = expit(-self.nodes)
        chi[0], chi[-1] = 1.,0.
        mobile = self.free & (self.nodes!=0)
        for _ in range(15):
            x = self.N@chi
            potential_derivative = 2*x*(1-x)*(1-2*x)
            residual = 2*(self.stiff@chi)+self.N.T@(self.weight*potential_derivative)
            error = max(abs(residual[mobile]))
            if error<3e-12:
                return chi, error
            hessian = 2*self.stiff + weighted_gram(self.N,self.weight*(2-12*x+12*x*x))
            chi[mobile] -= spsolve(hessian[mobile][:,mobile].tocsc(), residual[mobile])
        raise ArithmeticError(f'Order-sector equilibrium residual {error}')

    def reading(self):
        psi, carrier = Domain.solve(self, 0.)
        chi, error = self.order_profile()
        # Moving the dividing surface in +s by a compact virtual displacement.
        # Both baths stay fixed. Envelope derivative of the minimized energy.
        g = np.maximum(1-(self.s/7.)**2,0.)**2
        dg = -4*self.s/49*np.maximum(1-(self.s/7.)**2,0.)
        x = self.N@chi
        grad = self.grad@chi
        well = (x*(1-x))**2
        tension_normal = well-grad**2
        energy = well+grad**2
        hold = self.measure@(self.k*self.ratio*g*energy + self.ratio**2*dg*tension_normal)
        def shifted_energy(step):
            jac = (1+self.k*(self.s+step*g)/2)**2
            stretch = 1+step*dg
            return self.measure@(jac*(stretch*well+grad**2/stretch))
        step = 1e-4
        by_energy = (shifted_energy(step)-shifted_energy(-step))/(2*step)
        return dict(curvature=self.k, order_residual=error, carrier=carrier,
                    carrier_max_deviation=float(max(abs(psi-1))), half_order=float(chi[len(chi)//2]),
                    hold=hold, energy_difference=by_energy)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--elements-per-width',type=int,required=True)
    parser.add_argument('--spherical-calibration',type=Path,required=True)
    args=parser.parse_args()
    data=args.spherical_calibration.read_bytes()
    meta=json.loads((args.spherical_calibration.parent/'outcome.json').read_text())
    if meta['exitCode'] or meta['childOutcome']['exitCode'] or not meta['limitsVerified']:
        raise ValueError('Spherical calibration failed')
    ref, = [r for line in data.splitlines() if (r:=json.loads(line)).get('test')=='prediction'
            and r['curvature']>0]
    for k in [0., .06, -.06, .03, -.03]:
        row=Collar(k,args.elements_per_width).reading()
        prediction=ref['tension_signed_training']*k
        print(json.dumps(dict(test='cylindrical_patch', source_sha256=sha256(data).hexdigest(),
                              mesh=args.elements_per_width, signed_sphere_prediction=prediction,
                              **row), allow_nan=False),flush=True)
        if row['carrier_max_deviation']>1e-11 or abs(row['half_order']-.5)>1e-12:
            raise ArithmeticError('Coherent rest state or dividing surface changed')
        if abs(row['hold']-row['energy_difference'])>2e-9:
            raise ArithmeticError('Virtual work and energy difference disagree')


if __name__=='__main__':
    main()
