"""Fixed-state invariance and one production step for the held cubic R4 projector."""
import argparse
from pathlib import Path
import sys

import numpy as np


STEP287 = '5683615a04966ef35ba992026c8990965005c13bfec8e6d440534ade1e4b46c0'
MATCH286 = '530dc77673a9c71ada3fd369bed01922098beb2727adbe9941308aae78886282'


def load_state(study, output, *, with_fields):
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
    from mouth_r4_collar166_state import load_held
    from mouth_r4_cubic_held import approved_reference
    from mouth_r4_eulerian_release import order_grid
    from mouth_r4_stage2_outcome import verify_guard

    verify_guard(study/'primary_guard', 1)
    driver, trial, solid, _ = load_held(study, output, 1e-8, solve_fields=with_fields)
    approved_reference(trial.b)
    folder = study/'primary'
    record = read(folder/'step_00287.json')
    source = folder/record['file']
    if digest(source) != record['sha256'] or record['sha256'] != STEP287 or record['outer'] != 287:
        raise ValueError('Wrong step-287 source')
    with np.load(source, allow_pickle=False) as data:
        if (int(data['radius']) != 4 or int(data['winding']) != 1
                or str(data['parent_sha256']) != CHECKPOINT_SHA or not np.array_equal(data['x'], driver.x)):
            raise ValueError('Held shape or preparation changed')
        driver.labels, driver.fields = data['labels'], data['field']
    return driver, trial, solid, order_grid(trial, driver.fields), record


def force_bank(joint, labels):
    _, gradient, material = joint.evaluate(labels)
    elastic, physical, load, _ = joint.solid.evaluate(labels, joint.order)
    if elastic != material['elastic_energy']:
        raise ArithmeticError('Raw and assembled crystal energies differ')
    return material, dict(physical_nodal=-physical, physical_free=-joint.reduce(physical),
                          total_free=-gradient, physical_order_load=load)


def reference(study, output):
    from mouth_exact_checkpoint_control import digest, read
    from mouth_r4_aligned_crystal import AlignedCrystal
    from mouth_whole_r4_relax import atomic_json

    pinned = study/'pinned'
    tools = pinned/'software/stage1_solver/tools'
    manifest = read(pinned/'manifest.json')
    if Path(sys.modules[AlignedCrystal.__module__].__file__).resolve().parent != tools:
        raise ValueError('Frozen reference must use the original pinned evaluator')
    for path, wanted in manifest.items():
        if digest(tools/path) != wanted:
            raise ValueError(f'Changed pinned evaluator {path}')
    driver, trial, solid, order, _ = load_state(study, output, with_fields=False)
    joint = AlignedCrystal(driver, trial, solid, driver.labels, order)
    material, forces = force_bank(joint, driver.labels)
    target = output/'forces.npz'
    with target.open('xb') as stream:
        np.savez_compressed(stream, **forces)
    atomic_json(output/'result.json', dict(passed=True, source_sha256=STEP287,
                pinned_manifest=digest(pinned/'manifest.json'), force_bank_sha256=digest(target),
                material=material, sizes={k: int(v.size) for k,v in forces.items()}))
    print('frozen step 287', STEP287, 'reference DOFs', forces['total_free'].size, flush=True)


def compare_frozen(study, output, joint, driver, trial, solid, step):
    from mouth_exact_checkpoint_control import digest, read
    from mouth_r4_collar166_state import readout
    from mouth_r4_stage2_outcome import verify_guard
    from mouth_whole_r4_relax import atomic_json

    verify_guard(output.parent/'reference_guard', 0)
    verify_guard(study/'stop_readout_guard', 0)
    path = output.parent/'reference'
    old = read(path/'result.json')
    if (old['source_sha256'] != STEP287 or not old['passed']
            or digest(path/'forces.npz') != old['force_bank_sha256']):
        raise ValueError('Frozen force bank lacks original step-287 provenance')
    material, forces = force_bank(joint, driver.labels)
    archived = read(study/'stop_readout/result.json')
    observed = readout(driver, trial, solid)
    bank_deltas = {}
    with np.load(path/'forces.npz', allow_pickle=False) as before:
        if set(before.files) != set(forces):
            raise ValueError('Incomplete force registry')
        for name, new in forces.items():
            if before[name].shape != new.shape:
                raise ValueError(f'Force registry changed for {name}')
            bank_deltas[name] = float(np.max(np.abs(before[name]-new)))
    anchor = archived['readout']
    compare_keys = ('total_energy', 'physical_energy', 'crystal_energy', 'gauge_energy', 'peak_active')
    scalar_deltas = {k:observed[k]-anchor[k] for k in compare_keys}
    current_deltas = {k:observed['current'][k]-v for k,v in anchor['current'].items()}
    field_deltas = {k:observed['field'][k]-v for k,v in anchor['field'].items() if isinstance(v,float)}
    checks = dict(original_state=archived['source_sha256'] == STEP287,
                  all_force_arrays=max(bank_deltas.values()) == 0, physical_and_gauged_material=material == old['material'],
                  readout=observed == anchor, saved_step_energy=bool(abs(observed['total_energy']-step['total_energy'])
                  <= 2*abs(np.spacing(observed['total_energy']))))
    result = dict(sha256=STEP287,
                  reference_result_sha256=digest(path/'result.json'),
                  archived_stop_result_sha256=digest(study/'stop_readout/result.json'),
                  energy_force_deltas=scalar_deltas, max_force_array_deltas=bank_deltas,
                  max_force_count=old['sizes'], field_deltas=field_deltas, current_deltas=current_deltas,
                  delta_vs_archived_accumulated_total=observed['total_energy']-step['total_energy'],
                  frozen=observed, checks=checks, passed=all(checks.values()))
    atomic_json(output/'invariance.json', result)
    if not result['passed']:
        raise ValueError(f'Frozen invariance failed: {checks}')
    return result


def restore_cache(study, joint, driver, step):
    from mouth_exact_checkpoint_control import digest, read

    folder=study/'primary'
    anchor=read(folder/'matched_00286.json')
    path=folder/anchor['file']
    if digest(path) != anchor['sha256'] or anchor['sha256'] != MATCH286:
        raise ValueError('Wrong step-286 diagonal anchor')
    with np.load(path, allow_pickle=False) as saved:
        if not np.array_equal(driver.fields,saved['field']) or not np.array_equal(driver.x,saved['x']):
            raise ValueError('Fields or geometry changed inside the fixed-order burst')
        _, stats=joint.inverse_diagonal(saved['labels'],anchor['material']['chart']['chart_det_global_lower'])
    search=step['search']
    if stats != search['jacobi'] or search['gradient_l2_before'] != anchor['material']['gradient_l2']:
        raise ValueError('Reconstructed diagonal or CG history differs from step 287')
    joint.initial_force=anchor['material']['max_free_label_force']
    joint.previous_norm=search['gradient_l2_before']
    for rejected in search['rejections']:
        site=rejected['site']
        local=site.get('chart_fractions', site.get('local'))
        if 'cell' in site and local is not None:
            pair=(tuple(site['cell']),tuple(local))
            if pair not in joint.tight:
                joint.tight.append(pair)
    return dict(cache_source=MATCH286, jacobi=stats, initial_force=joint.initial_force,
                previous_norm=joint.previous_norm, retained_sites=len(joint.tight))


def replay(study, output):
    from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
    from mouth_r4_aligned_crystal import AlignedCrystal
    from mouth_r4_cubic_physical_orientation import physical_certificate
    from mouth_whole_r4_relax import atomic_json, save_fields

    driver, trial, solid, order, step = load_state(study, output, with_fields=True)
    joint = AlignedCrystal(driver, trial, solid, driver.labels, order)
    inv=compare_frozen(study,output,joint,driver,trial,solid,step)
    cache=restore_cache(study,joint,driver,step)
    val,g,material=joint.evaluate(driver.labels)
    after,_,_,m,search=joint.step(driver.labels,val,g,material)
    orientation=physical_certificate(trial.b,after,m['chart']['chart_det_global_lower'])
    diag=read(study.parent/'r4_step288_projection01/result.json')['joint']
    expected={(tuple(s['cell']),tuple(s['local'])) for s in diag['constraints']}
    registered={(tuple(c),tuple(l)) for c,l in joint.tight}
    events=search['projection_events']
    protecting={(tuple(s['cell']),tuple(s['local'])) for s in events[-1]['sites']}
    checks=dict(invariance=inv['passed'], diagonal_reconstructed=True,
                direction_fraction=search['alpha'] == diag['alpha'] == .25,
                signed_work=abs(search['energy_change']-diag['total_work']) <= 1e-10,
                peak_force=abs(m['max_free_label_force']-diag['peak_force']) <= 1e-7,
                five_sites=len(expected)==5 and expected.issubset(protecting) and protecting.issubset(registered),
                projected_work=events[-1]['slope']<0 and events[-1]['kkt_norm']<=1e-10,
                oriented=orientation['cell_count']==28080 and orientation['physical_global_lower']>=8.54e-8)
    post=dict(alpha=search['alpha'], work=search['energy_change'],peak_force=m['max_free_label_force'],
              work_error=search['energy_change']-diag['total_work'],
              peak_force_error=m['max_free_label_force']-diag['peak_force'],
              registered_sites=[dict(cell=list(c),local=list(l)) for c,l in sorted(registered)],
              accepted_direction_sites=events[-1]['sites'], projection_fallbacks=search['projection_fallbacks'],
              diagnostic_sha256=digest(study.parent/'r4_step288_projection01/result.json'))
    record=dict(outer=288,initial_sha256=STEP287,epsilon=driver.epsilon,source_sha256=digest(Path(__file__)),
                cache=cache,invariance=inv,orientation=orientation,post_trial=post,checks=checks,
                passed=all(checks.values()))
    if not record['passed']:
        record['candidate_search']=search
        atomic_json(output/'result.json',record)
        raise ArithmeticError(f'Production replay missed a gate: {checks}')
    state=output/'step_00288.npz'
    save_fields(state,labels=after,field=driver.fields,x=driver.x,radius=4,winding=1,epsilon=driver.epsilon,
                outer=288,iteration=10,parent_sha256=CHECKPOINT_SHA,initial_sha256=STEP287)
    state_sha=digest(state)
    info=dict(outer=288,depth=step['depth'],exit_radius=step['exit_radius'],source_sha256=STEP287,
              sha256=state_sha,file=state.name,labels_fixed_order=True,material=m,search=search,
              total_energy=step['total_energy']+search['energy_change'],orientation=orientation)
    atomic_json(output/'step_00288.json',info)
    record.update(state_sha256=state_sha,record_sha256=digest(output/'step_00288.json'))
    atomic_json(output/'result.json',record)
    print('held R4 production step 288',post,'physical lower',orientation['physical_global_lower'],flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('reference','replay'))
    parser.add_argument('study',type=Path)
    parser.add_argument('out',type=Path)
    args=parser.parse_args()
    study,output=args.study.resolve(),args.out.resolve()
    output.mkdir(parents=True,exist_ok=False)
    if args.operation=='reference':
        sys.path.insert(0,str(study/'pinned/software/stage1_solver/tools'))
        reference(study,output)
    else:
        replay(study,output)


if __name__=='__main__':
    main()
