"""Evaluate the same outgoing-bath power standard on wave-empty, held-shape breathing modes."""
import argparse
import hashlib
import json

import numpy as np

from mouth_breathing_radiation import calibration
from mouth_exit_witness import inputs
from mouth_poly_shallow_release import ROOT, Trial, emit
from mouth_sealed_shape_relax import Meridian


def run(radius, directional):
    n = 72
    original = ROOT / f'_scratch/s11c/mouth_breathing_matched01/r{radius}n{n}'
    folder = original.with_name(original.name+'_direction') if directional else original
    guard = json.loads((folder / 'outcome.json').read_bytes())
    if not guard['limitsVerified'] or guard['exitCode'] or guard['childOutcome']['exitCode']:
        raise ValueError('Invalid reference spectrum')
    raw = (folder / 'stdout').read_bytes()
    rows = list(map(json.loads, raw.splitlines()))
    source, = [r for r in rows if r['kind'] == 'input']
    expected = .001
    if directional:
        parent = (original / 'stdout').read_bytes()
        if hashlib.sha256(parent).hexdigest() != source['partial_stdout_sha256']:
            raise ValueError('Invalid parent for adopted spectrum')
        expected = source['central_mesh_step']
    center, = [r for r in rows if r['kind'] == 'spectrum' and r['relative_exit'] == 0 and r['shape_step'] == expected]
    modes = [r for r in rows if r['kind'] == 'signed' and r['exit_step'] == .006]
    trial = Trial(radius, n, bonded_core=True)
    path = ROOT / f'_scratch/s11c/mouth_exit_release01/r{radius}n{n}/endpoint.npz'
    shifts, _, _, digest, status, *_ = inputs(path, radius, n, trial)
    if digest != source['material_sha256'] or status != 'balanced':
        raise ValueError('Different material preparation')
    shape = Meridian(radius, trial.x+shifts[:4])
    emit(kind='input',radius=radius,material_sha256=digest,mode_source_sha256=hashlib.sha256(raw).hexdigest(),
         fixed_opening=shape.exit,depth=shape.length,positive_frequency_count=len(modes))
    for mode in modes:
        i = mode['index']
        v = np.array(center['mass_normalized_vectors'][i])
        freq = center['frequency'][i]
        mass = v @ np.array(center['mass']) @ v
        coarse, fine = [calibration(shape,v,source['original_exit_scale'],freq,nq) for nq in (64,128)]
        error = max(abs(a['power_per_squared_modal_speed']-b['power_per_squared_modal_speed'])
                    for a,b in zip(coarse['cases'],fine['cases']))
        emit(kind='empty_mode_power',radius=radius,index=i,frequency=freq,unit_modal_mass=mass,
             source_force=mode['outward_force_per_action'],power_quadrature_error=error,**fine)
        if (abs(mass-1)>1e-9 or abs(fine['unsigned_volume_norm']-abs(mode['volume_change_at_unit_modal_coordinate']))>1e-7
                or error>1e-8 or fine['transformed_volume_errors'][-1]>5e-7):
            raise ArithmeticError('Unresolved absolute source or shape normalization')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('radius',choices=(1,2),type=int)
    p.add_argument('--directional',action='store_true',help='Use a finished directional adoption of the original central spectrum.')
    args = p.parse_args()
    run(args.radius,args.directional)
