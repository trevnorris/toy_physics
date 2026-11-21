"""Entrance force per acceleration with both internal bends free of applied force.

Differentiate the carrier equations at zero Laplace frequency. Compare the result
with finite-frequency forces after subtraction of the small finite-bath spring.
"""
import argparse
from pathlib import Path

import numpy as np

from mouth_carrier_lateral_response import prepare
from mouth_flared_receiver import emit
from mouth_free_bend_drive import Response


def entrance_force(matrix):
    return matrix[0,0] - matrix[0,1:] @ np.linalg.solve(matrix[1:,1:],matrix[1:,0])


def main(args):
    dom, field = prepare(args)
    response = Response(dom, field, args.extension)
    solve = response.dipole.factor.solve
    tangent = solve(response.bend)
    velocity = response.motion+response.skew@tangent
    inertia = (2/3)*velocity.T@solve(velocity)
    cross = response.motion.T@tangent
    gyro = (2/3)*(cross-cross.T-tangent.T@(response.skew@tangent))
    zero = response.force_matrix(0.)
    floor = np.min((abs(dom.N@field)**8-1)/4+1/response.dipole.rq**2)
    for tension in (.02, .04):
        static = zero+tension*response.area
        k = static[1:,1:]
        translation = np.r_[1., -np.linalg.solve(k, static[1:,0])]
        coupling = (gyro@translation)[1:]
        base = translation@inertia@translation
        extra = coupling@np.linalg.solve(k, coupling)
        ev = np.linalg.eigvalsh(static)
        emit(test='entrance_inertia', tension=tension, field_floor=floor,
             static_eigenvalues=ev.tolist(), bath_spring=entrance_force(static),
             static_translation=translation.tolist(), rigid=inertia[0,0],
             projected_inertia=base, bend_flow_increment=extra, total=base+extra)
        if min(ev) <= 0 or floor <= 0 or np.max(abs(gyro+gyro.T))>1e-8:
            raise ArithmeticError('Positive-energy or antisymmetry check failed')
        for omega in (.01, .005, .0025):
            matrix = response.force_matrix(1j*omega)+tension*response.area
            z = entrance_force(matrix)
            measured = (entrance_force(static)-z.real)/omega**2
            emit(test='finite_frequency', tension=tension, omega=omega, apparent_inertia=measured,
                 derivative_limit=base+extra, imaginary_force=float(z.imag),
                 acceleration_per_entrance_force=float((-omega**2/z).real))
        if not np.isclose(measured,base+extra,rtol=.01):
            raise ArithmeticError('Small-frequency impedance failed derivative check')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--extension',choices=['compact','wide'],required=True)
    main(p.parse_args())
