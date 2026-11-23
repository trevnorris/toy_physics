"""Two-sided carrier and spatially resolved conversion layer on a recorded radius-four mouth."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mouth_exchange_carrier import BIAS
from mouth_poly_shallow_release import emit
from mouth_sealed_shape_relax import Meridian
from mouth_wider_wound_carrier import summarize_disk
from mouth_winding_resolved_patch import Patch
from mouth_wound_carrier import Wound


def file_sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def material_sha(params, order, extra):
    digest = hashlib.sha256()
    arrays = (params, order, [extra]) if extra != 0. else (params, order)
    for values in arrays:
        values = np.asarray(values, dtype='<f8')
        digest.update(str(values.shape).encode())
        digest.update(np.ascontiguousarray(values).tobytes())
    return digest.hexdigest()


def frozen_extra(fields):
    values = [float(fields[k]) for k in ('shape_extra','extra') if k in fields]
    if values and any(v != values[0] for v in values):
        raise ValueError('Inconsistent extra shape coordinate')
    return values[0] if values else 0.


def check_guard(path):
    outcome = json.loads((path.parent/'outcome.json').read_bytes())
    limits = json.loads((path.parent/'limit-validation.json').read_bytes())
    if (outcome['exitCode'] or outcome['childOutcome']['exitCode']
            or not outcome['limitsVerified'] or not limits['verified']
            or int(limits['actual']['memory.swap.max'])):
        raise ValueError(f'Unvalidated import {path}')


def carrier(source, save, mesh, polar, winding):
    check_guard(source)
    with np.load(source, allow_pickle=False) as archive:
        radius = int(archive['radius'])
        if radius != 4:
            raise ValueError('Unexpected entrance')
        params, order = archive['shape_params'], archive['order']
        extra = frozen_extra(archive)
    shape = Meridian(radius, params, extra=extra)
    digest = material_sha(params, order, extra)
    if not np.isfinite(order).all() or order.min() < 0 or order.max() > 1:
        raise ValueError('Crystalline order outside constitutive domain')
    model = Wound(shape, order, winding, mesh, polar, 28., 22.)
    amplitude, phase, info = model.solve(BIAS)
    fields = model.center_fields(amplitude, phase, info)
    emit(kind='sharp_disk', radius=radius, winding=winding, depth=shape.length, exit=shape.exit, mesh=mesh, polar=polar,
         source_sha256=file_sha(source), material_sha256=digest, shape_extra=extra, **summarize_disk(fields, model, info))
    np.savez_compressed(save, **fields, order=order, shape_params=params, shape_extra=extra,
                        radius=radius, exit=shape.exit, material_sha256=digest, winding=winding,
                        outlet_phase_bias=BIAS)
    emit(kind='carrier_field', file=str(save), sha256=file_sha(save))


def resolved(source, mesh, halfwidth):
    check_guard(source)
    with np.load(source, allow_pickle=False) as fields:
        radius = int(fields['radius'])
        if radius != 4 or int(fields['winding']) not in (0,1) or fields['outlet_phase_bias'] != BIAS:
            raise ValueError('Unmatched entrance, winding or driving phase')
        params = fields['shape_params']
        order = fields['order']
        extra = frozen_extra(fields)
        digest = material_sha(params, order, extra)
        if digest != str(fields['material_sha256']):
            raise ValueError('Geometry and order verification failed')
        b = float(fields['exit'])
    if mesh[2] % 2:
        raise ValueError('Resolved exit must be a mesh plane')
    patch = Patch(source, radius, b, mesh, halfwidth, material_state=(params, order, digest), shape_extra=extra)
    patch.run()


if __name__=='__main__':
    a = argparse.ArgumentParser(description=__doc__)
    modes = a.add_subparsers(dest='action', required=True)
    sharp = modes.add_parser('carrier')
    sharp.add_argument('source', type=Path)
    sharp.add_argument('save', type=Path)
    sharp.add_argument('--mesh', type=int, default=13)
    sharp.add_argument('--polar', type=int, default=10)
    sharp.add_argument('--winding', type=int, default=1, choices=(0,1))
    front = modes.add_parser('layer')
    front.add_argument('source', type=Path)
    front.add_argument('--mesh', type=int, nargs=3, default=[14, 10, 48])
    front.add_argument('--halfwidth', type=float, default=.12)
    args = a.parse_args()
    if args.action == 'carrier':
        carrier(args.source, args.save, args.mesh, args.polar, args.winding)
    else:
        resolved(args.source, args.mesh, args.halfwidth)
