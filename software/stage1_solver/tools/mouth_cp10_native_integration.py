"""Refine the integration of the fixed checkpoint-10 Lagrangian map.

This control keeps the native crystal elements, fields and virtual motions.
It measures how much numerical integration in the comparison standard can
contribute to a spatial label refinement residual.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_free_shape_solid import Solid
from mouth_r4_fold_order import axes_model
from mouth_sealed_shape_relax import Meridian


def run(root, level, output):
    checkpoint = root/'r4_melted_survival01/coarse/state_00010.npz'
    prep = root/'supercore_joined_front01/adopt309/state.npz'
    if digest(checkpoint) != CHECKPOINT_SHA:
        raise ValueError('Requires the untouched radius-four checkpoint 10')
    with np.load(checkpoint, allow_pickle=False) as saved:
        if str(saved['state_sha']) != digest(prep) or int(saved['winding']) != 1 or int(saved['radius']) != 4:
            raise ValueError('Mismatched wound R4 preparation')
        pars, waist, polar = saved['shape_params'], float(saved['shape_extra']), float(saved['x'][5])
        q, order = saved['q'], saved['field'][2]
    depth = read(root/'r4_whole_mouth01/coarse/endpoint.json')['depth']

    def energy(mode=None, step=0.):
        p = pars.copy()
        w, a = waist, polar
        if mode == 0:
            p[0] += step
        elif mode in (1,2,3):
            p[mode] += .25*step
        elif mode == 4:
            w += .25*step
        elif mode == 5:
            a += step
        geometry = Meridian(4., p, extra=w)
        b = Solid(prep, geometry, a, 8, 24, 6, axes_model(geometry, depth),
                  degree=level, subdivisions=2)
        b.melted_labels = True
        b.data.exact_melt_support = 1
        b.buffers['base'].fill(0.)
        return b.evaluate(q, order, load=False, gradient=False)[0]

    reference = read(root/'r4_inverse_material01/cp10/result.json')
    h = 1.e-4
    held = energy()
    values = []
    for i in range(6):
        minus,plus = energy(i,-h),energy(i,h)
        values.append(dict(minus=minus,plus=plus,force=(minus-plus)/(2*h)))
    out = dict(state=CHECKPOINT_SHA,degree=level,subdivisions=2,step=h,
               energy=held,energy_minus_saved=held-reference['lagrangian_crystal'],six=values)
    with output.open('x') as f: json.dump(out,f,indent=2,allow_nan=False)
    print(json.dumps(out,allow_nan=False))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('degree',type=int,choices=(4,5))
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    run(args.root,args.degree,args.output)
