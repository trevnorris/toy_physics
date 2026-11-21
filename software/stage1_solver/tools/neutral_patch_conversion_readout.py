"""Reduce matched opposite loads and audit the fixed-mouth impulse ledger."""
import argparse
import json
from pathlib import Path

import numpy as np


def emit(**row):
    print(json.dumps(row,allow_nan=False),flush=True)


def load(folder):
    read = lambda name: json.loads((folder/name).read_text())
    outcome, limits = read('outcome.json'), read('limit-validation.json')
    if (outcome['exitCode'] or outcome['childOutcome']['exitCode'] or
            not outcome['limitsVerified'] or not limits['verified'] or
            int(limits['actual']['memory.swap.max'])):
        raise ValueError(f"Guarded measurement incomplete: {folder}")
    rows = [json.loads(line) for line in (folder/'stdout').read_text().splitlines()]
    setup,*gates,loss = rows
    if (setup['test'] != 'preparation' or loss['test'] != 'first_loss_order'
            or len(gates) != 6 or not all(r['test']=='crossing' for r in gates)):
        raise ValueError("Incomplete trajectory")
    return setup,gates,loss


def run(folder):
    setup,gates,loss = load(folder)
    drive = np.array(setup['drive_order'])
    if list(drive) != [0.,.02,-.02,.01,-.01]:
        raise ValueError("Unmatched drives")
    pairs = setup['pairs']
    tag = pairs.index(setup['tag'])
    for row in gates[1:]:
        t = row['t']
        def odd(key,half=False):
            z = np.asarray(row[key])
            ip,im = (3,4) if half else (1,2)
            return (z[:,ip]-z[:,im])/2
        normalized = {k:odd(k)/(t*drive[1]) for k in
                      ('rim_impulse','end_impulse','line_momentum','bath_impulse','port_impulse')}
        ledger = sum(normalized[k] for k in ('rim_impulse','end_impulse','line_momentum','bath_impulse'))
        np.testing.assert_allclose(ledger,1.,atol=1e-8,rtol=0)
        port_ledger = normalized['rim_impulse']+normalized['end_impulse']+normalized['line_momentum']
        np.testing.assert_allclose(port_ledger,normalized['port_impulse'],atol=1e-8,rtol=0)
        groups = row['strict_zero_ids']
        surviving = [len(g) for g in groups]
        alive_tags = np.array(row['loaded_tag_continuous'])
        if [tag in g for g in groups] != list(alive_tags[:,0]):
            raise ValueError("Inconsistent unforced tag record")
        rim = odd('rim')
        fraction = normalized['rim_impulse']
        half_fraction = odd('rim_impulse',True)/(t*drive[3])
        emit(run=folder.name,chi=round(row['chi'],6),t=t,
             unforced_survivors=surviving,distinct_ids=pairs,
             survivor_paths_per_pair=[sum(i in g for g in groups) for i in range(len(pairs))],
             tag_survival_by_drive=alive_tags.sum(axis=0).tolist(),
             neutral_threshold_survivors=row['absolute_zero_count'],
             rim_force_history=rim.tolist(),mean_rim_force=float(rim.mean()),
             mean_tag_force=float(odd('tag_force').mean()),
             rim_fraction_history=fraction.tolist(),mean_rim_fraction=float(fraction.mean()),
             half_fraction_history=half_fraction.tolist(),mean_half_fraction=float(half_fraction.mean()),
             max_full_half_abs_fraction_difference=float(max(abs(fraction-half_fraction))),
             mean_impulse=float(odd('rim_impulse').mean()),
             mean_time_average_force=float(odd('rim_impulse').mean()/t),
             mean_ledger={k:float(v.mean()) for k,v in normalized.items()})
    emit(run=folder.name,tag=setup['tag'],tag_first_loss_chi=[
        history[tag] for history in loss['unforced_by_seed_and_original_pair']])


def combined(folders):
    """Combine disjoint refinements of the same protocol using the saved drives."""
    sets = [load(folder) for folder in folders]
    first = sets[0][0]
    seeds = [s for setup,_,_ in sets for s in setup['seeds']]
    if seeds != [8401,8402,8403,8404]:
        raise ValueError("Missing or duplicate thermal histories")
    for setup,_,_ in sets:
        for key in first.keys()-{'seeds'}:
            if setup[key] != first[key]:
                raise ValueError(f"Different physical or numerical parameter: {key}")
    tag = first['pairs'].index(first['tag'])
    for index in range(1,6):
        rows = [gates[index] for _,gates,_ in sets]
        def gather(key):
            return np.concatenate([row[key] for row in rows])
        t = rows[0]['t']
        impulse = gather('rim_impulse')
        force = gather('rim')
        full = (impulse[:,1]-impulse[:,2])/2
        half = (impulse[:,3]-impulse[:,4])/2
        input_force = first['drive_order'][1]
        groups = [g for row in rows for g in row['strict_zero_ids']]
        emit(test='four_history_refinement',law=first['law'],chi=round(rows[0]['chi'],6),
             t=t,mean_rim_force=float(np.mean(force[:,1]-force[:,2])/2),
             mean_time_average_rim_force=float(full.mean()/t),mean_rim_impulse=float(full.mean()),
             percent_of_applied_impulse=float(100*full.mean()/(t*input_force)),
             half_drive_percent=float(100*half.mean()/(t*first['drive_order'][3])),
             surviving_counts=[len(g) for g in groups],tag_present=[tag in g for g in groups],
             survival_paths_per_id=[sum(i in g for g in groups) for i in range(len(first['pairs']))])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folders',type=Path,nargs='+')
    parser.add_argument('--combine',action='store_true')
    args=parser.parse_args()
    if args.combine:
        combined(args.folders)
    else:
        for folder in args.folders:
            run(folder)
