"""Collect the guarded matched-knot force measurements and their common-point controls."""
import argparse
import hashlib
import json
from pathlib import Path


def load(path):
    return json.loads(path.read_bytes())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def collect(root, output):
    tools = Path(__file__).resolve().parent
    repo = tools.parents[2]
    previous = repo/'research/pde_ledger_v3/_measurements/mouth_r4_spatial_refinement.json'
    reference = load(previous)
    state_sha = reference['checkpoint']
    if sha(root/'r4_melted_survival01/coarse/state_00010.npz')!=state_sha:
        raise ValueError('Checkpoint 10 changed')
    record = dict(input_commit='7b45bfec',checkpoint_sha256=state_sha,phase_drop=reference['phase_drop'],
                  prior_measurement_sha256=sha(previous),sources={},matched={},order_work={},summary={})
    def measured(area,file,guard):
        result = area/file
        g = area/guard
        outcome = load(g/'outcome.json')
        cmd = load(g/'invocation.json')['command']
        if (outcome['exitCode']!=0 or outcome['childOutcome']['exitCode']!=0 or not outcome['limitsVerified']
                or load(g/'effective-limits.json')['memory.swap.max']!='0'):
            raise ArithmeticError(f'Guard failed: {g}')
        item = load(result)
        if item['checkpoint_sha256']!=state_sha:
            raise ValueError(f'Wrong input at {result}')
        record['sources'][str(result.relative_to(root))]=dict(sha256=sha(result),command=cmd)
        return item
    comparisons = root/'r4_matched_knot01'
    order = root/'r4_order_transport01'
    boundary = root/'r4_matched_boundary01'
    record['boundaries'] = measured(boundary,'boundaries.json','verified_boundaries')
    record['production'] = measured(boundary,'work.json','verified_work')
    record['signed_label_inversion'] = measured(order,'signed_inverses.json','inverse_guard')
    bath = reference['checks']['fixed_chart_depth_4']['analytical_bath']
    saved = reference['saved_total_six'][0]
    record['common_bath_depth_force'] = bath
    record['archived_cp10_depth_force'] = saved
    for degree in (6,8):
        component = measured(order,f'degree_{degree}/result.json',f'guard_{degree}')
        if component['degree']!=degree:
            raise ValueError('Wrong independent order-work quadrature')
        record['order_work'][str(degree)]=component
        for pos,level in enumerate((1,2,4)):
            key = f'{level}_{degree}'
            guard = [p.name for p in comparisons.glob(f'*guard_{key}') if (p/'outcome.json').is_file()]
            if len(guard)!=1:
                raise ValueError(f'Expected one completed job for {key}')
            row = measured(comparisons,f'{key}/readout.json',guard[0])
            if (not row['completed'] or row['degree']!=degree or row['factor']!=level
                    or row['cxx_sha256']!=sha(tools/'mouth_cp10_matched_depth.cpp')
                    or row['labels_sha256']!=component['label_sha256'][pos]
                    or row['labels_sha256']!=record['signed_label_inversion']['grids'][str(level)]['labels_sha256']):
                raise ValueError(f'Incompatible row {key}')
            record['matched'][key]=row
            full,half=row['steps']
            if full['h']!=2*half['h'] or half['h']!=2.e-7:
                raise ValueError(f'Incompatible steps in {key}')
            for step in (full,half):
                if step['min_Jinv_over_positive_samples']<=0 or step['signed_positive_order_evaluations']<=0:
                    raise ValueError(f'Invalid active map {key}')
            native,spatial=half['forces_native_projected']
            delta=spatial-native
            if delta!=half['projected_minus_native']:
                raise ValueError(f'Inconsistent difference for {key}')
            record['summary'][key]=dict(depth_residual=delta,
                full_minus_half=full['projected_minus_native']-delta,
                native_total=native+bath,spatial_total=spatial+bath,
                native_minus_archived=native+bath-saved,
                independent_native_image_stress_delta=component['projected_minus_native'][pos][0],
                independent_native_image_order_delta=component['projected_minus_native'][pos][1])
    with output.open('x') as stream:
        json.dump(record,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(record['summary'],allow_nan=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('output',type=Path)
    arguments=parser.parse_args()
    collect(arguments.root.resolve(),arguments.output.resolve())
