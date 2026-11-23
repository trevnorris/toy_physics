"""Signed material work along the held one-quantum exit family, independent of carrier force."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_circulation_material import collar, held, supported_order_probes
from mouth_poly_bonded_core import core_source
from mouth_poly_shallow_release import CUT, Order, emit, exterior, source
from mouth_sealed_shape_relax import Meridian
from mouth_winding_order_readout import reconstructed


def deformed(reference, exit_radius):
    coordinates = reference.params.copy()
    coordinates[1] += np.log(exit_radius / reference.exit)
    return Meridian(reference.neck, coordinates)


def inside(trial, shape, chi, registration):
    t, x = [np.linspace(0., 1., n) for n in chi.shape]
    load = source(shape, registration, t, x, False)
    load += core_source(shape, registration, *trial.core_reference, t, x)
    owner = Order(shape, trial.n)
    energy, gradient = owner.energy_gradient(chi)
    return owner, energy, load, gradient


def material_energy(trial, reference, chi, co, cr, exit_radius):
    shape = deformed(reference, exit_radius)
    _, order_energy, elastic, _ = inside(trial, shape, chi, co)
    receiving = exterior(shape, trial.n, CUT,
                         side_order_work=lambda w, x: source(shape, cr, w/shape.length, x, True))
    return dict(inner_order=order_energy, inner_elastic=np.sum(chi*elastic),
                outer_order=receiving['energy'], outer_elastic=receiving['static_work']), receiving['residual']


def run(radius, b, steps, field=None):
    trial, reference, chi, co, cr, digest = held(radius)
    if b is None:
        b = reference.exit
    shape = deformed(reference, b)
    emit(kind='material_input', radius=radius, exit=b, reference_exit=reference.exit, order_sha256=digest,
         depth=shape.length, signed_steps=steps, mesh=trial.n)
    owner, _, load, grad = inside(trial, shape, chi, co)
    probes = supported_order_probes(radius, chi, owner)
    targets = owner.free & (chi == 0)
    if np.any(targets):
        ix = np.unravel_index(np.argmin(np.where(targets, grad+load, np.inf)), chi.shape)
        single = np.zeros_like(chi)
        single[ix] = 1.
        probes['least_material_melted'] = single
    emit(kind='held_order_and_collar', inner=collar(shape, co, True), outer=collar(shape, cr, False),
         min_core_order=float(np.min(chi[owner.free])), axial_order_at_waist=owner.axis(chi),
         released_order_zero_count=int(targets.sum()), free_core_nodes=int(owner.free.sum()),
         order_work_before_carrier={k:float(np.sum(v*(grad+load))) for k,v in probes.items()},
         least_melted_node=[int(i) for i in ix] if np.any(targets) else None)
    for h in steps:
        left, res_l = material_energy(trial, reference, chi, co, cr, b*(1-h))
        right, res_r = material_energy(trial, reference, chi, co, cr, b*(1+h))
        components = {k:(left[k]-right[k])/(2*h*b) for k in left}
        emit(kind='material_force', radius=radius, exit=b, step=h, components=components,
             absolute_material_force=sum(components.values()), receiver_residuals=[res_l,res_r],
             left_energy=left, right_energy=right)
    if field is not None:
        field = Path(field)
        guard = json.loads((field.parent/'outcome.json').read_bytes())
        if guard['exitCode'] or guard['childOutcome']['exitCode'] or not guard['limitsVerified']:
            raise ValueError('Unverified central wound solution')
        emit(kind='central_carrier', file_sha256=hashlib.sha256(field.read_bytes()).hexdigest())
        least = np.unravel_index(np.argmin(np.where(owner.free, chi, np.inf)), chi.shape)
        probe = np.zeros_like(chi)
        probe[least] = 1.
        probes['smallest_order_site'] = probe
        emit(kind='minimum_order_site', node=[int(i) for i in least], starting_order=float(chi[least]))
        with np.load(field, allow_pickle=False) as z:
            if int(z['winding']) != 1 or int(z['radius']) != radius or str(z['material_sha256']) != digest:
                raise ValueError('Unmatched carrier and crystalline reference')
            if abs(float(z['exit'])-b)>1.e-9:
                raise ValueError('Wrong central exit')
            for quadrature in (2,3):
                for result in reconstructed(z, shape, chi, probes, quadrature):
                    crystal_work = np.sum(probes[result['probe']]*(load+grad))
                    if abs(result['carrier_work']-result['finite_difference_check'])>2.e-8*(1+abs(result['carrier_work'])):
                        raise ValueError('Carrier order gradient failed its energy check')
                    emit(kind='held_order_total_work', material_work=crystal_work,
                         total_work=crystal_work+result['carrier_work'], **result)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('radius', type=int, choices=(1,2))
    p.add_argument('--exit', type=float)
    p.add_argument('--steps', type=float, nargs='*', default=[.001])
    p.add_argument('--field')
    args = p.parse_args()
    run(args.radius,args.exit,args.steps,args.field)
