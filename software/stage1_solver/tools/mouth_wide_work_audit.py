"""Independent signed work on a held exit, with both receiving-order samples minimized."""
import argparse
import json

import numpy as np

from mouth_circulation_material import held
from mouth_live_waist_fields import BIAS, Carrier
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import CUT, Order, emit, exterior, source
from mouth_sealed_shape_relax import Meridian


def main(radius, exit_radius, offset, legacy_carrier):
    trial, shape, order, inside, outside, signature = held(radius)
    # Meridian does not retain its input parameter vector. Recover the exact
    # four shape coordinates directly from the archived preparation.
    from mouth_exit_witness import inputs
    from mouth_poly_shallow_release import ROOT
    endpoint = ROOT / f'_scratch/s11c/mouth_exit_release01/r{radius}n72/endpoint.npz'
    shift, _, _, digest, status, *_ = inputs(endpoint, radius, 72, trial)
    if digest != signature or status != 'balanced':
        raise ValueError('Invalid common preparation')
    t = np.linspace(0., 1., order.shape[0])
    x = np.linspace(0., 1., order.shape[1])
    target = shape.exit if exit_radius is None else exit_radius
    emit(kind='input', radius=radius, exit=target, prepared_exit=shape.exit,
         relative_signed_step=offset, material_sha256=digest)
    previous = {}
    for sign in (-1, 1):
        b = target*(1+sign*offset)
        coords = trial.x + shift[:4]
        coords[1] += np.log(b/shape.exit)
        moved = Meridian(radius, coords)
        if abs(moved.exit-b) > 1e-13:
            raise ArithmeticError('Wrong exit variation')
        order_energy = Order(moved, trial.n).energy_gradient(order)[0]
        collar = source(moved, inside, t, x, False)
        core = core_source(moved, inside, *trial.core_reference, t, x)
        outer = exterior(moved, trial.n, CUT, side_order_work=lambda w,u:
                         source(moved, outside, w/moved.length, u, True))
        energies = dict(interior_order=order_energy, inside_collar=order.ravel()@collar.ravel(),
                        core=order.ravel()@core.ravel(), receiver=outer['energy']+outer['static_work'])
        if legacy_carrier:
            energies['legacy_carrier'] = Carrier(moved, trial.n, 2*trial.n, trial.mesh).response(order, BIAS)[0]
        emit(kind='signed', sign=sign, exit=moved.exit, receiver_residual=outer['residual'],
             energies=energies)
        previous[sign] = energies
    forces={k:(v-previous[1][k])/(2*target*offset) for k,v in previous[-1].items()}
    emit(kind='audit_force', radius=radius, exit=target, half_offset=offset,
         material=sum(value for key,value in forces.items() if key!='legacy_carrier'),
         total_with_archived_carrier=sum(forces.values()) if legacy_carrier else None,
         component_forces=forces)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('radius',type=int,choices=[1,2])
    parser.add_argument('--exit',type=float)
    parser.add_argument('--step',type=float,default=.003)
    parser.add_argument('--legacy-carrier',action='store_true')
    args=parser.parse_args()
    main(args.radius,args.exit,args.step,args.legacy_carrier)
