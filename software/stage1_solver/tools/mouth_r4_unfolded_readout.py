"""Assemble the independent unfolded-control preparations and six matched-knot derivatives."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_cp10_spatial_crystal import MODES
from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read


def reduce(area, output):
    manifest = {}
    def guarded(filename, guard):
        path = area / filename
        log = area / guard
        out = read(log / 'outcome.json')
        if (out['exitCode'] or not out['limitsVerified'] or out['childOutcome']['exitCode']
                or int(read(log/'effective-limits.json')['memory.swap.max'])):
            raise ValueError(f'Invalid guarded result at {log}')
        manifest[filename] = dict(sha256=digest(path), command=read(log/'invocation.json')['command'])
        return read(path)

    preparation = guarded('control.json', 'control_guard')
    key = digest(area / 'control.npz')
    if preparation['parent_sha256'] != CHECKPOINT_SHA or preparation['control_sha256'] != key:
        raise ValueError('Control provenance mismatch')
    levels = (1, 2, 4)
    labels = [guarded(f'{l}/result.json', f'guard_{l}') for l in levels]
    label_hash = [digest(area / str(l) / 'labels.npz') for l in levels]
    for i, lab in enumerate(labels):
        if (lab['sha256'] != label_hash[i] or lab['source_sha256'] != key or
                lab['duplicate'] or lab['missing'] or min(lab['lower_bernstein_bound_chart_Jinv'].values())<=0):
            raise ValueError('Invalid labels')
    inverse = guarded('depth/prepared.json', 'depth_guard_prepare')
    depth = guarded('depth/result.json', 'depth_guard_combine')
    if depth['source_sha256'] != key or depth['prepared_sha256'] != digest(area/'depth/prepared.json'):
        raise ValueError('Signed depth provenance mismatch')
    for l, hsh in zip(levels, label_hash):
        row = inverse['grids'][str(l)]
        if row['static_labels_sha256'] != hsh or row['signed_sha256'] != digest(area/f'depth/signed_{l}.npz'):
            raise ValueError('Mismatched signed labels')
    sectors = [guarded(f'depth/polar_{k}/result.json', f'depth_guard_{k}') for k in range(6)]
    if any(s['sector'] != i or s['source_sha256'] != key for i,s in enumerate(sectors)):
        raise ValueError('Missing signed sector')
    sources = {}
    results = {}
    for degree in (6, 8):
        shape = guarded(f'shape_{degree}/result.json', f'shape_guard_{degree}')
        if (shape['degree']!=degree or shape['control_sha256']!=key or
                shape['label_sha256']!=label_hash or shape['crystal_columns']!=['energy',*MODES[1:]]):
            raise ValueError('Five-force control has different material state or virtual motions')
        sources[degree] = shape
        values=[]
        raw=[]
        split_energy=[]
        for i, l in enumerate(levels):
            dc=depth['rows'][f'{l}_{degree}']
            delta=shape['spatial_minus_material'][i]
            row=[delta[0],dc['spatial_minus_native'],*delta[1:]]
            fv=np.asarray(shape['total_force'])
            pair=np.array([[shape['crystal_values'][0][0],dc['native_total'],*fv[0]],
                           [shape['crystal_values'][i+1][0],dc['spatial_total'],*fv[i+1]]])
            if np.max(abs(pair[1]-pair[0]-row))>2e-14:
                raise ArithmeticError('Common fields do not cancel between the two formulations')
            values.append(row)
            raw.append(pair.tolist())
            es=np.asarray(dc['signed_energies'])
            if es.shape!=(2,2):
                raise ValueError('Missing signed energy pair')
            # Same two maps on the signed depth quadrature, with all moving metric
            # boundaries split. The mean of the two displaced energies is a
            # near-static cross-check; the fixed-state integral above is primary.
            signed_delta=es[:,1]-es[:,0]
            split_energy.append(dict(signed_delta=signed_delta.tolist(),mean_delta=float(signed_delta.mean()),
                                     max_deviation_from_mean=float(np.ptp(signed_delta)/2)))
        diff=np.diff(np.abs(values),axis=0)
        results[degree]=dict(raw_pairs=raw,spatial_minus_material=values,
                             near_static_energy_on_split_signed_rule=split_energy,
                             change_in_absolute_error=diff.tolist(),
                             monotone_in_absolute_value=np.all(diff<0,axis=0).tolist())
    record=dict(parent_sha256=CHECKPOINT_SHA,control_sha256=key,columns=['stored_crystal_energy',*MODES],
                levels=levels,raw_pair_order=['native material','spatial labels'],epsilon_reg=0,
                phase_drop=0.4510084282857874,fluid_potential='(n^5 - 5n + 4)/20',
                preparation=preparation,labels=labels,exact_signed_inverse=inverse,depth=depth,five_shapes=sources,
                results=results,manifest=manifest)
    with output.open('x') as stream:
        json.dump(record,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(results,allow_nan=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    reduce(args.area.resolve(),args.output.resolve())
