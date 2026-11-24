"""Reduce accepted R4 melted-collar states and their independently checked starting point."""
import argparse
import hashlib
import json
from pathlib import Path


def read(path):
    return json.loads(path.read_text())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def trajectory(folder):
    saved, rows = {}, []
    short_material, short_shape = [], []
    for line in (folder/'stdout').read_text().splitlines():
        row = json.loads(line)
        kind = row.get('kind')
        if kind in ('rejected_material_rematch', 'rejected_shape_rematch'):
            (short_material if kind.startswith('rejected_material') else short_shape).append(row)
        if kind == 'accepted':
            checkpoint = folder/row['file']
            if digest(checkpoint) != row['sha256']:
                raise ValueError(f'Changed accepted state {checkpoint}')
            saved[row['iteration']] = row
        if kind != 'residuals':
            continue
        i = row['iteration']
        state = saved[i]
        domain = state['orientation']
        if not domain['all_active_cells_oriented']:
            raise ValueError(f'Accepted state {i} has unverified active orientation')
        cur, mat, field = row['current'], row['solid'], row['field']
        phase = cur['shaft_order']
        # The native step count precedes field rematching. When the proposal
        # is shortened, distinguish the accepted displacement from this count.
        material_scale = (short_material[-1]['scale']/2 if short_material else 1.) if row['material_steps'] else 0.
        shape_scale = short_shape[-1]['scale']/2 if short_shape else 1.
        rows.append(dict(update=i, sha256=state['sha256'], depth=state['depth'],
            exit=state['exit_radius'], minimum_radius=state['minimum_radius'], pole=state['pole_radius'],
            equator=state['equator_radius'], energy=state['energy'], shape_force=row['shape_force'],
            max_abs_shape_force=max(map(abs,row['shape_force'])), material=mat,
            material_proposal_steps=row['material_steps'], accepted_material_fraction=material_scale,
            prior_shape_step_fraction=shape_scale, material_message=row['material_stop'],
            projected_and_metric_residuals=field['field'], raw_residuals=field['raw_residual'],
            zero_order_check=field['zero_order_check'], domain=domain,
            front_fraction=phase['terminal_from_fraction'], entrance=phase['sections'][0],
            exit_rim=phase['sections'][-2], distinct_nodes=phase['shaft_distinct_nodes'],
            node_count_below=phase['shaft_nodes_at_or_below'],
            flow={k: cur[k] for k in ('feed_current','bulk_current','weak_disk_current','port_error',
                                        'disk_minus_feed_current','incident_slice_current','receiving_slice_current',
                                        'dense_total_subsonic_share')},
            backtracked_material=short_material, backtracked_shape=short_shape))
        short_material, short_shape = [], []
    return rows, saved, short_material, short_shape


def main(area, output):
    folder = area/'coarse'
    outcome = read(folder/'outcome.json')
    limits = read(folder/'limit-validation.json')
    duration = read(folder/'duration-validation.json')
    if (not outcome['limitsVerified'] or not limits['verified'] or limits['actual']['memory.swap.max']!='0'
            or not duration['verified'] or duration['wallDeadlineSeconds'] is not None
            or duration['inactivityDeadlineSeconds'] is not None):
        raise ValueError('Guard settings not verified')
    end = 'DONE' if (area/'DONE.json').is_file() else 'ISSUE'
    receipt = read(area/(end+'.json'))
    if (end=='DONE') != (outcome['exitCode']==0):
        raise ValueError('Disagreement between guard and scientific endpoint')
    inv = read(folder/'invocation.json')
    sources = read(area/'sources.json')
    for pathname, expected in sources.items():
        if digest(area/'pin'/pathname) != expected:
            raise ValueError(f'Changed source {pathname}')
    rows, saved, mat_tail, shape_tail = trajectory(folder)
    extras = {}
    for name in ('saved_state_comparison','endpoint','issue'):
        path = folder/(name+'.json')
        if path.exists():
            data=read(path)
            extras[name] = data if name=='saved_state_comparison' else dict(
                sha256=digest(path), error=data.get('error'), half_force=data.get('half_shape_force'))
    if end=='DONE' and ('endpoint' not in extras or not rows or not saved):
        raise ValueError('Missing coarse equilibrium output')
    if end=='DONE' and receipt['sha256']!=saved[max(saved)]['sha256']:
        raise ValueError('Mismatched coarse endpoint')
    result=dict(scope='24x24 radius-four continuation; sealed carrier sidewall, one winding, continuous order',
        source_sha256=sources, command=inv['command'], start=inv['startedUtc'],
        wall_seconds=outcome['wallSeconds'], exit_code=outcome['exitCode'], note=receipt.get('error'),
        stdout=dict(path=str(folder/'stdout'),bytes=(folder/'stdout').stat().st_size,sha256=digest(folder/'stdout')),
        stderr=(folder/'stderr').read_text()[-3500:], states={str(k):v['sha256'] for k,v in saved.items()},
        trajectory=rows, final_rejected_material=mat_tail, final_rejected_shape=shape_tail, **extras)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('area',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    main(args.area,args.output)
