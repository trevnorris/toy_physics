"""Measure force on a held ordered layer and transport across it in the 4D coupon."""
import argparse
from collections import defaultdict
from hashlib import sha256
import json
from math import prod, sqrt
from pathlib import Path
from statistics import mean, stdev


def interval(values):
    n = len(values)
    t975 = {3: 4.30265273, 6: 2.57058184}[n]
    mid = mean(values)
    error = t975 * stdev(values) / sqrt(n)
    return dict(n=n, mean=mid, interval95=[mid-error, mid+error], values=values)


def read_case(record):
    blocks = record['blocks']
    duration = record['measurement_time']
    area = prod(record['box'][:3])
    amounts = {}
    for name in ('sheet', 'middle'):
        for sign in ('pos', 'neg'):
            amounts[f'{name}_{sign}'] = sum(b['ledger'][name]['count_'+sign] for b in blocks)
            amounts[f'{name}_px_{sign}'] = sum(b['ledger'][name]['p_'+sign][0] for b in blocks)
        amounts[name+'_net_current'] = (amounts[name+'_pos']-amounts[name+'_neg'])/(duration*area)
        amounts[name+'_net_px'] = (amounts[name+'_px_pos']-amounts[name+'_px_neg'])/(duration*area)
    for k, letter in ((0,'x'),(3,'w')):
        amounts['whole_force_'+letter] = sum(b['sheet_mobile_force_impulse'][k] for b in blocks)/duration
        amounts['traction_'+letter] = amounts['whole_force_'+letter] / area
        amounts['last_block_traction_'+letter] = blocks[-1]['sheet_mobile_force_impulse'][k]/(record['dt']*blocks[-1]['steps']*area)
        for zone in (1,2):
            samples = [b['drift_zones'][zone] for b in blocks]
            amounts[f'zone_{zone}_v{letter}'] = sum(z['integrated_displacement'][k] for z in samples)/sum(z['residence_site_time'] for z in samples)
    # For each branch the force on the sheet is body + bath - stored momentum rate.
    for key in ('body', 'bath'):
        amounts[key+'_force_x'] = sum(b['ledger']['global'][key][0] for b in blocks)/duration
    amounts['stored_rate_x'] = (blocks[-1]['ledger']['global']['final_p'][0]-blocks[0]['ledger']['global']['initial_p'][0])/duration
    amounts['pww'] = sum(b['steps']*record['dt']*b['pressure_diagonal'][3] for b in blocks)/duration
    for z,side in enumerate(('minus','plus')):
        amounts[f'bath_{side}_heat'] = sum(b['bath_kinetic_energy_delivered_minus_plus'][z] for b in blocks)
    for k,name in enumerate(('held','mobile')):
        end = record['final_structure']
        amounts[name+'_lost_fraction'] = end['starting_contacts_lost'][k]/end['starting_contacts'][k]
    amounts['connected_mobile'] = end['largest_connected_mobile']
    for z,name in enumerate(record['zone_names']):
        for key in ('n','kinetic_T','tangent_bragg_intensity'):
            amounts[f'{name}_{key}'] = end['zones'][z][key]
    return amounts


def analyze(root):
    design = json.loads((root/'design.json').read_text())
    done = json.loads((root/'DONE.json').read_text())
    if set(design['cases']) != set(done['cases']):
        raise ValueError('Incomplete case matrix')
    cases = defaultdict(lambda: defaultdict(dict))
    for task, inp in design['cases'].items():
        raw = (root/task/'result.json').read_bytes()
        if sha256(raw).hexdigest()!=done['cases'][task]:
            raise ValueError(f'Hash mismatch for {task}')
        r = json.loads(raw)
        assert (r['seed'], r['dt'],r['drive'],r['T_minus'])==(inp['seed'],inp['dt'],[inp['ax'],0,0,inp['aw']],inp['t_minus'])
        cases[inp['group']][inp['seed']][inp['arm']] = read_case(r)
    result = {}
    for group, histories in cases.items():
        zeros = [arms['zero'] for arms in histories.values()]
        all_arms = [arm for arms in histories.values() for arm in arms.values()]
        result[group] = data = dict(histories=len(histories),runs=len(all_arms),
            counts={cut: {sg: sum(a[f'{cut}_{sg}'] for a in all_arms) for sg in ('pos','neg')} for cut in ('sheet','middle')},
            zero={k: interval([s[k] for s in zeros]) for k in zeros[0]})
        response = data['signed_response_per_body_acceleration'] = {}
        for axis in ('x','w'):
            for amp,prefix in ((.1,axis),(.05,axis+'h')):
                if prefix+'p' not in next(iter(histories.values())):
                    continue
                # Adjacent normal and lateral drifts are resolved independently from the external drive.
                response[prefix] = {q:interval([(arms[prefix+'p'][q] - arms[prefix+'m'][q])/(2*amp)
                    for arms in histories.values()]) for q in ('whole_force_'+axis,'traction_'+axis,
                        'last_block_traction_'+axis,'zone_1_v'+axis,'zone_2_v'+axis,'sheet_net_current','sheet_net_px',
                        'body_force_x','bath_force_x','stored_rate_x')}
            if axis+'h' in response:
                data[axis+'_full_minus_half_slope'] = {q:interval([a-b for a,b in zip(
                    response[axis][q]['values'],response[axis+'h'][q]['values'])]) for q in response[axis]}
        if 'x' in response:
            partition = {}
            for name, amplitude in (('x',.1),('xh',.05)):
                if name in response:
                    partition[name] = {key:interval([amplitude*v for v in response[name][key]['values']])
                        for key in ('whole_force_x','body_force_x','bath_force_x','stored_rate_x')}
            data['odd_force_budget'] = partition
        if 'pp' in next(iter(histories.values())):
            data['four_sign_mixed_derivative'] = {q:interval([(arms['pp'][q]-arms['pm'][q]-arms['mp'][q]+arms['mm'][q])/.04
                for arms in histories.values()]) for q in ('sheet_net_px','traction_x','middle_net_px')}
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    result = analyze(args.root)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    def compact(d):
        return {k:d[k] for k in ('mean','interval95')}
    for name,group in result.items():
        print(json.dumps(dict(sample=name, runs=group['runs'],counts=group['counts'],
            undriven={k:compact(group['zero'][k]) for k in ('pww','mobile_lost_fraction','connected_mobile','minus_near_tangent_bragg_intensity','plus_near_tangent_bragg_intensity','minus_far_kinetic_T','plus_far_kinetic_T')},
            signed={key:{q:compact(a) for q,a in v.items()} for key,v in group['signed_response_per_body_acceleration'].items()},
            odd_force_budget={k:{term:compact(estimate) for term,estimate in balance.items()} for k,balance in group.get('odd_force_budget',{}).items()},
            mixed={k:compact(a) for k,a in group.get('four_sign_mixed_derivative',{}).items()}),allow_nan=False))
