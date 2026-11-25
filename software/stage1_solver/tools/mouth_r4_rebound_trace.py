"""Read the historical step 90–108 forces and distinguish physical samples from chart certificates."""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import numpy as np

from mouth_r4_edge_descent import edge_data


def read(path):
    return json.loads(path.read_bytes())


def trace(study):
    measures = study.parents[2]/'research/pde_ledger_v3/_measurements'
    prior = read(measures/'mouth_r4_newton_stop_work.json')
    journal = study/'primary_guard/stdout'
    source_sha = hashlib.sha256(journal.read_bytes()).hexdigest()
    if prior['source_sha256']['primary_guard/stdout'] != source_sha:
        raise ValueError('Unverified primary log')
    events = [json.loads(line) for line in journal.read_text().splitlines()]
    steps = {r['outer']:r for r in events if r['kind'] == 'newton_step'}
    updates = {r['outer']:r for r in events if r['kind'] == 'newton_rematched'}
    final = read(study/'primary/last_active.json')
    edge = read(measures/'mouth_r4_newton_edge.json')
    # Only the product of the three chart widths enters the determinant.
    # This saved normalization is rounded to binary64. Cross-check against
    # the independently calculated edge at 108 before using it at 100.
    widths = (Fraction(edge['normalizing_cell_volume']), Fraction(1), Fraction(1))
    saved_edges = {}
    saved_hashes = {}
    for name in ('last_active','rematched'):
        location = study/'primary'/name
        info = read(location.with_suffix('.json'))
        state = location.with_suffix('.npz')
        saved_hashes[name] = hashlib.sha256(state.read_bytes()).hexdigest()
        if saved_hashes[name] != info['sha256']:
            raise ValueError(f'Changed archived labels at {name}')
        with np.load(state, allow_pickle=False) as data:
            _, minimum, _, _ = edge_data(data['labels'], widths)
        saved_edges[info['outer']] = float(minimum)
    if abs(saved_edges[108]-edge['min_on_blocking_edge'])>2e-17:
        raise ArithmeticError('Held edge normalization disagrees with step 108')
    rows = []
    for step in range(90,109):
        row = steps[step]
        # A rebuilt Jacobi records the certificate passed INTO the next step.
        # A reused Jacobi proves only that this certificate was in its reuse window.
        # If the fields were rematched first, the recorded chart is the rematched one.
        if step == final['outer']:
            bound = final['material']['chart']['chart_det_global_lower']
            interval = [bound, bound]
        else:
            jacobi = steps[step+1]['jacobi']
            base = jacobi['chart_at_assembly']
            interval = [base,base] if jacobi['rebuilt'] else [.5*base,2*base]
        update = updates.get(step)
        rows.append(dict(step=step, active_Linf=row['max_force'], joint_L2=row['gradient_l2'],
                         alpha=row['alpha'], cg_iterations=row['cg_iterations'], cg_reason=row['cg_reason'],
                         scaled_diagonal_on_first_nonpositive_curvature=bool(
                             row['cg_iterations']==1 and row['cg_reason']=='nonpositive curvature'),
                         physical_Jinv_min_sampled=row['min_Jinv'], chart_lower_bound_interval=interval,
                         field_update=update, limiting_edge_minimum=saved_edges.get(step)))
    return dict(source_sha256=source_sha, label_hashes=saved_hashes, previous_study='r4_newton_held01',
                meaning='The interval encloses the value of the certified lower bound, not the true physical '
                        'minimum. At a field-update row it applies after that update. The sampled physical '
                        'inverse determinant applies on the accepted crystal step. CG truncation on its first '
                        'iteration selects scaled steepest descent. Additional line-search fallback is not '
                        'logged for each step. The limiting edge is reconstructed at the two stored label maps only.',
                rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', type=Path)
    print(json.dumps(trace(parser.parse_args().study.resolve()), indent=2, allow_nan=False))
