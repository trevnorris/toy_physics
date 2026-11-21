"""Linear carrier response when translation and both flared-wall bends are free.

The far bath surfaces stay fixed. The three coordinate displacements along x
are a common translation and two successive offsets along w. Their common
extension vanishes before either bath boundary. No intrinsic wall mass.
"""
import argparse
import numpy as np
from scipy import sparse
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from mouth_carrier_lateral_inertia import taper
from mouth_carrier_lateral_response import Dipole, prepare, weighted_gram
from mouth_carrier_attachment_tension import area_stiffness
from mouth_empty_bend_stiffness import INTERVALS, quadrature, shapes
from mouth_flared_receiver import emit
from pathlib import Path


def real_pair(values):
    return np.concatenate((values.real, values.imag), axis=0)


class Response:
    def __init__(self, dom, field, extension):
        r, w, gr, gw = quadrature(dom)
        hs, sr, sw = shapes(r, w, extension, dom.far, dom.roof)
        radial, dr = taper(r, 2.0 if extension == 'compact' else 3.4, dom.far-.35)
        depth, dw = taper(w, 3.3 if extension == 'compact' else 4.3, dom.roof-.35)
        shape = np.column_stack((radial*depth, hs))
        sr = np.column_stack((dr*depth, sr))
        sw = np.column_stack((radial*dw, sw))
        f = (dom.N@field)[:, None]
        fr = (gr@(field-1))[:, None]
        fw = (gw@(field-1))[:, None]
        wt = dom.weight[:, None]
        drive = sr*fr + sw*fw
        # One half of the physical Hessian with the common 1/3 angular weight removed.
        bend = .5*(gr.T@(wt*drive) + gw.T@(wt*(sw*fr-sr*fw))
                   + dom.N.T@(wt*(2*drive/r[:, None]-sr*(abs(f)**8-1)*f/2)))
        # The pulled-back Berry term has both Jacobian and transport variations.
        # Integration by parts in time gives the sr*f term in the velocity source.
        motion = .5j*(gr.T@(wt*shape*f)
                      + dom.N.T@(wt*(shape*(2*f/r[:, None]-fr)+sr*f)))
        weight = wt*abs(fr)**2/3
        self.frozen = sr.T@(weight*sr) + sw.T@(weight*sw)
        self.dipole = Dipole(dom, field)
        free = self.dipole.free
        self.bend = real_pair(bend)[free]
        self.motion = real_pair(motion)[free]
        self.operator = self.dipole.operator[free][:, free].tocsc()
        mass = weighted_gram(dom.N, dom.weight)
        skew = sparse.bmat([[None, -mass], [mass, None]], format='csr')
        self.skew = skew[free][:, free].tocsc()
        self.area = np.diag([0., *(area_stiffness(dom, span)[0] for span in INTERVALS)])
        self.cache = {}

    def force_matrix(self, rate):
        if rate not in self.cache:
            rhs = self.bend + rate*self.motion
            system = self.operator - rate*self.skew
            solve = self.dipole.factor.solve if rate == 0 else splu(system.tocsc()).solve
            field = solve(rhs)
            residual = np.max(abs(system@field-rhs))
            output = self.frozen - (2/3)*(self.bend-rate*self.motion).T@field
            if not np.isfinite(output).all() or residual > 3e-8:
                raise ArithmeticError(f'Failed forced field solve at {rate}: {residual}')
            self.cache[rate] = output
        return self.cache[rate]


def main(args):
    dom, field = prepare(args)
    response = Response(dom, field, args.extension)
    static = response.force_matrix(0.)
    emit(test='static', holding=static.tolist(), area=response.area.diagonal().tolist(),
         reciprocity_error=float(np.max(abs(static-static.T))))
    if np.max(abs(static-static.T)) > 1e-8:
        raise ArithmeticError('Static reciprocity failed')
    for omega in (.02, .01, .005):
        z = response.force_matrix(1j*omega)
        emit(test='impedance', omega=omega, real=z.real.tolist(), imag=z.imag.tolist(),
             rigid_inertia=float(-(z[0,0].real-static[0,0])/omega**2),
             hermitian_error=float(np.max(abs(z-z.conj().T))))
        if np.max(abs(z-z.conj().T)) > 1e-8:
            raise ArithmeticError('Lossless spectral force matrix violates Hermitian symmetry')
        for tension in args.tensions:
            load = np.array([1., 0., 0.])  # applied only at the entrance
            matrix = z + tension*response.area
            x = np.linalg.solve(matrix, load)
            acc = -omega**2*np.array([x[0], x.sum()])
            ratio = x.sum()/x[0]
            residual = np.max(abs(matrix@x-load))
            emit(test='push', omega=omega, area_tension=tension,
                 entrance_and_flare_acceleration_real=acc.real.tolist(),
                 entrance_and_flare_acceleration_imag=acc.imag.tolist(),
                 flare_over_entrance=[float(ratio.real), float(ratio.imag)],
                 force_residual=float(residual))
            if residual > 1e-8:
                raise ArithmeticError('Unbalanced free coordinates')
    # A sign change brackets an exponentially growing homogeneous solution;
    # a scan cannot exclude other poles. Real s avoids real-frequency cavity resonances.
    samples = [.01, .025, .05, .10, .15, .22, .32, .45, .65, .9]
    for tension in args.tensions:
        def det(rate):
            return np.linalg.det(response.force_matrix(rate)+tension*response.area)
        values = [det(s) for s in samples]
        emit(test='real_growth_scan', area_tension=tension, rates=samples, det=values)
        for a, b, first, second in zip(samples, samples[1:], values, values[1:]):
            if first*second >= 0:
                continue
            root = brentq(det, a, b, xtol=1e-8)
            matrix = response.force_matrix(root)+tension*response.area
            _, singular, right = np.linalg.svd(matrix)
            mode = right[-1]/np.max(abs(right[-1]))
            emit(test='unforced_mode', area_tension=tension, bracket=[a,b], growth=root,
                 singular_values=singular.tolist(), coordinates=mode.tolist(),
                 relative_balance=float(np.max(abs(matrix@mode))/singular[0]))
            if np.max(abs(matrix@mode))/singular[0] > 1e-7:
                raise ArithmeticError('Unforced mode failed its force balance')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--extension', choices=['compact','wide'], default='compact')
    p.add_argument('--tensions', type=float, nargs='+', default=[0., .02])
    main(p.parse_args())
