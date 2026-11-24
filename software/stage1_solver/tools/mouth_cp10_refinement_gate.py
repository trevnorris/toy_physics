"""Compare frozen CP10 transfers at the same material points and at matched shape motion."""
import argparse
import json
from pathlib import Path

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read


def guarded(area, name, guard):
    status = read(area/guard/'outcome.json')
    if not status['limitsVerified'] or status['exitCode'] != 0 or status['childOutcome']['exitCode'] != 0:
        raise ArithmeticError(f'Incomplete measurement in {guard}')
    path = area/name
    return read(path), dict(file=str(name), sha256=digest(path),
                            command=read(area/guard/'invocation.json')['command'])


def compare(root, original):
    energy = root/'r4_nested_energy01'
    forces = root/'r4_nested_forces01'
    labels = root/'r4_nested_labels01'
    control = root/'r4_nested_gate01'
    old = read(original)
    d = old['data']
    ref = d['total_shape_forces']
    rows, sources, tests = {}, {}, {}
    first = d['same_material_quadrature']
    full = d['split_registry_total']
    material_depth = full['depth_motion_comparison'][-1]
    native = first['native']
    common = np.array(ref['native_saved'])
    bath = np.array(ref['field_analytic'])
    if old['checkpoint_sha256'] != CHECKPOINT_SHA or len(common) != 6:
        raise ValueError('Incompatible native baseline')
    baseline_forces = np.array(full['analytic'])
    baseline_forces[0] = material_depth['total_material_motion_plus_analytic_bath']
    same = first['by_side_correspondence'][0]
    rows['1'] = dict(energy=first['projected'],delta=first['energy_difference'],
                     deficit_percent=-100*first['energy_difference']/native,
                     same_side_delta=same[2]-same[1],matched_six=baseline_forces.tolist(),
                     residual_six=(baseline_forces-common).tolist())
    for level in (2,4):
        ev,es = guarded(energy,f'{level}/result.json',f'guard{level}')
        lab,ls = guarded(labels,f'{level}/result.json',f'guard{level}')
        if (ev['checkpoint'] != CHECKPOINT_SHA or lab['source'] != CHECKPOINT_SHA or
            ev['label_sha256'] != lab['sha256'] or digest(labels/f'{level}/labels.npz') != lab['sha256'] or
            ev['gauss_points_unchanged'] != sum(int(g[0]) for g in first['by_side_correspondence']) or
            abs(ev['total']['original']-native)>1e-14 or lab['missing'] or lab['multiple']):
            raise ValueError('Mismatch in transfer/source material points')
        v = ev['by_correspondence']
        rows[str(level)] = dict(energy=ev['total']['projected'],delta=ev['total']['difference'],
             deficit_percent=100*ev['signed_relative_deficit'],
             same_side_delta=v['same_material_prism']['difference']+v['other_prism_same_bond_side']['difference'])
        sources[f'energy_{level}']=es
        sources[f'labels_{level}']=ls
        tests[f'energy_{level}']=ev
        tests[f'labels_{level}']=lab

    def force_case(level,deg):
        case=f'{level}_{deg}'
        found=[path.name for path in forces.glob(f'*guard_{case}') if (path/'outcome.json').exists()]
        if len(found)!=1: raise ValueError(f'Expected one admitted measurement for {case}, got {found}')
        result, source=guarded(forces,f'{case}/readout.json',found[0])
        if result['checkpoint']!=CHECKPOINT_SHA or not result['completed']:
            raise ValueError(f'Unverified shapes for {case}')
        expected_label = old['data']['chart']['label_sha256'] if level==1 else tests[f'labels_{level}']['sha256']
        if result['label_sha256']!=expected_label or result['modes']!=ref['modes']:
            raise ValueError('Mismatched label or modal context')
        if not np.array_equal(bath+result['matched_crystal'],result['total']):
            raise ValueError('Total force does not decompose')
        tests[f'force_{case}']=result
        sources[f'force_{case}']=source
        return result

    for level,deg in ((1,8),(2,4),(2,6),(2,8),(4,4),(4,6),(4,8)):
        f=force_case(level,deg)
        if deg==6:
            rows[str(level)].update(matched_six=f['total'],residual_six=f['delta_against_saved'])
    for level in (1,2,4):
        v,src=guarded(energy,f'depth{level}/result.json',f'depth_guard{level}')
        if v['checkpoint']!=CHECKPOINT_SHA or v['active_points']!=5291652:
            raise ValueError('Depth control uses a different checkpoint or measure')
        tests[f'fixed_chart_depth_{level}']=v
        sources[f'fixed_chart_depth_{level}']=src
    for deg in (4,5):
        v,src=guarded(control,f'native_{deg}.json',f'native_guard_{deg}')
        if v['state']!=CHECKPOINT_SHA: raise ValueError('Native integration control source changed')
        v['total_minus_saved']=(bath+np.array([c['force'] for c in v['six']])-common).tolist()
        tests[f'native_{deg}']=v
        sources[f'native_{deg}']=src
    localization,src=guarded(energy,'summary.json','guard_summary')
    sources['regional_energy']=src
    errors=[abs(row['residual_six'][0]) for row in rows.values()]
    losses=[row['deficit_percent'] for row in rows.values()]
    gate=dict(energy_strictly_decreases=bool(np.all(np.diff(losses)<0)),
              transport_corrected_depth_strictly_decreases=bool(np.all(np.diff(errors)<0)),
              reference_degree=6, energy_empirical_orders=localization['successive_empirical_orders'])
    gate['degree8_depth_residuals']=[tests[f'force_{i}_8']['delta_against_saved'][0] for i in (1,2,4)]
    gate['both_primary_trends_pass']=gate['energy_strictly_decreases'] and gate['transport_corrected_depth_strictly_decreases']
    return dict(checkpoint=CHECKPOINT_SHA,input_commit='f40af401',source_measurement_sha256=digest(original),
                phase_drop=0.4510084282857874,reference_energy=native,saved_total_six=common.tolist(),mode_names=ref['modes'],
                rows=rows,gate=gate,localization=localization,checks=tests,sources=sources)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('scratch','original','output'): p.add_argument(key,type=Path)
    args=p.parse_args()
    data=compare(args.scratch,args.original)
    with args.output.open('x') as f: json.dump(data,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(dict(rows=data['rows'],gate=data['gate'])),flush=True)
