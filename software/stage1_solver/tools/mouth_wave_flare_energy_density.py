"""Predict axial force at two unused radii with a single mode energy per volume.

Occupation is assigned once per radius; each virtual stretch then holds its
action fixed. The depths in this test come from the carrier-only 0PN ratio.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import yaml

from mouth_wave_flare_work import ARCHIVE, Wave, emit, requirements


def verified_lines(folder):
    path=Path(folder)
    outcome=json.loads((path/'outcome.json').read_text())
    limits=json.loads((path/'limit-validation.json').read_text())
    if outcome['exitCode'] or outcome['childOutcome']['exitCode'] or not limits['verified']:
        raise ValueError('Invalid carrier holdout')
    return [json.loads(line) for line in (path/'stdout').read_text().splitlines()]


def holdout(folder):
    row=verified_lines(folder)[-1]
    return row['radius'], (row['depth'],row['required_restoration'],row['pressure_only'],row['tension_only'])


def main(args):
    p=yaml.safe_load(ARCHIVE.read_text())['parameters']
    inputs=requirements()
    for folder in args.holdout:
        radius, data=holdout(folder)
        if radius in inputs:
            raise ValueError('Holdout duplicates the calibration data')
        inputs[radius]=data
    rows=verified_lines(args.wave_reference)
    settings=rows[0]
    if settings['step']!=args.step or settings['padding']!=args.pad:
        raise ValueError('Wave calibration grid differs from holdouts')
    baseline={row['radius']:row for row in rows if row['test']=='occupied_minus_empty'}
    modes={}
    profile={}
    reference=inputs[1.][1]
    for r in sorted(inputs.keys()-{.5}):
        h, required, pressure, tension=inputs[r]
        if r in baseline:
            row=baseline[r]
            if row['depth']!=h:
                raise ValueError('Carrier and wave calibration depths differ')
            omega,force=row['omega'],row['force_per_action']
        else:
            omega, force=Wave(r,h,args.step,args.pad,p).read(sideways=False)
        if omega>=p['mu_out']:
            raise ValueError('Energy-density preparation needs a bound mode')
        modes[r]=(omega,force)
        # A geometrically similar collar-plus-flare has volume proportional
        # to R^3 H. The energy per geometric volume is the input law.
        wave_load=r**3*h*force/omega
        profile[r]=np.array([pressure/reference,tension/reference,-wave_load])
    ref_wave=-profile[1.][2]
    for r in profile:
        profile[r][2]/=ref_wave
    anchors=(1.,2.,4.)
    matrix=np.vstack([profile[r] for r in anchors])
    coeff=np.linalg.solve(matrix,[inputs[r][1] for r in anchors])
    emit(test='energy_density_calibration',radii=anchors,mesh_step=args.step,
         pressure_force_at_one=coeff[0],tension_force_at_one=coeff[1],
         wave_force_at_one=coeff[2],matrix_condition=np.linalg.cond(matrix),
         action_at_one_per_bias_squared=coeff[2]/modes[1.][1])
    for r in sorted(profile):
        predicted=profile[r]@coeff
        h,required,_,_=inputs[r]
        emit(test='energy_density_prediction',radius=r,depth=h,held_out=r not in anchors,
             required=required,predicted=predicted,residual=predicted-required,
             fractional_residual=(predicted-required)/required,
             wave_force= -profile[r][2]*coeff[2],wave_action_per_bias_squared=
             -profile[r][2]*coeff[2]/modes[r][1])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step',type=float,required=True)
    parser.add_argument('--pad',type=float,required=True)
    parser.add_argument('--wave-reference',required=True)
    parser.add_argument('--holdout',nargs='+',required=True)
    main(parser.parse_args())
