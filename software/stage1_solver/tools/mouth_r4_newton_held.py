"""Joint Newton–CG crystal bursts at the fixed radius-four mouth, starting at step 29."""
import argparse
from pathlib import Path
from time import perf_counter

import numpy as np

from mouth_exact_checkpoint_control import CHECKPOINT_SHA, digest, read
from mouth_r4_eulerian_relax_labels import FORCE_TOL
from mouth_r4_eulerian_release import order_grid
from mouth_r4_eulerian_solid import SpatialCrystal
from mouth_r4_newton_crystal import JointCrystal
from mouth_r4_smooth_endpoint import compare, write_endpoint
from mouth_r4_smooth_fields import coherent, read_fields
from mouth_r4_smooth_held import SmoothHeld
from mouth_whole_r4_relax import atomic_json, log

PAUSED_SHA = 'd2e764212c08bdcacfbc18121812ceb484f2eaffcd0b44f46991a9cf03120202'
FIELD_TARGET = 1e-8
BENCHMARK_STEPS = 5
MAX_BURST = 20
FIELDS = ('order_kkt','order_projected','amplitude_raw','phase_raw','carrier_metric')


def accepted(study, name, state):
    record = read(study / name)
    path = study / state
    if not record['passed'] or digest(path) != record['saved']:
        raise ValueError(f'Unverified initial state at {path}')
    return path


def initial(study, case):
    if case == 'benchmark':
        path = study/'paused/last_active.npz'
        old = read(path.with_suffix('.json'))
        if digest(path) != PAUSED_SHA or old['sha256'] != PAUSED_SHA or old['outer'] != 29:
            raise ValueError('Unexpected source for the step-29 benchmark')
        return path, 29, None
    if case == 'primary':
        gate = read(study/'benchmark_review.json')
        if not gate['passed'] or gate['benchmark_sha256'] != digest(study/'benchmark/result.json'):
            raise ValueError('Fixed-field benchmark requires independent clearance')
        path = accepted(study,'benchmark/result.json','benchmark/last_active.npz')
        state = read(study/'benchmark/result.json')
        return path, state['outer'], state['burst']
    path = accepted(study,'primary/HELD.json','primary/stationary_state.npz')
    return path, 0, None


def build_joint(driver, t, solid, *, registry_update):
    order = order_grid(t,driver.fields)
    joint = JointCrystal(driver,t,solid,driver.labels,order)
    if registry_update:
        driver.labels = joint.compatible(driver.labels)
    return joint, joint.evaluate(driver.labels)


def run(study, case):
    start_time = perf_counter()
    root = study / case
    root.mkdir()
    source, previous_count, resume_burst = initial(study,case)
    epsilon = 1e-10 if case == 'secondary' else 1e-8
    driver = SmoothHeld(study.parent,root,epsilon)
    driver.case, driver.start_sha = case, digest(source)
    with np.load(source,allow_pickle=False) as state:
        if (int(state['radius']) != 4 or int(state['winding']) != 1 or
                str(state['parent_sha256']) != CHECKPOINT_SHA or not np.array_equal(driver.x,state['x'])):
            raise ValueError('Held shape or checkpoint ancestry changed')
        driver.labels, driver.fields = state['labels'].copy(),state['field'].copy()
    t = driver.trial()
    driver.current_trial = t
    solid = SpatialCrystal(t.b)
    joint, (value, gradient, material) = build_joint(driver,t,solid,registry_update=False)
    physical, _, field = read_fields(driver,t,solid,driver.labels,driver.fields)
    total_energy = physical + material['gauge_energy']
    initial_force = material['max_free_label_force']
    force_history = []
    phases = []
    best_alpha = None
    if resume_burst is not None:
        burst = resume_burst.copy()
        best_alpha = read(study/'benchmark/result.json')['minimum_alpha']
    else:
        burst = dict(count=0,start_force=initial_force)
    if case == 'benchmark' and abs(initial_force-read(source.with_suffix('.json'))['material']['max_free_label_force'])>1e-11:
        raise ArithmeticError('Starting force differs from the accepted step-29 force')
    count = previous_count
    setup_seconds = perf_counter() - start_time
    clock = perf_counter()
    initial_energy = total_energy
    log('newton_initial',case=case,source_sha256=driver.start_sha,setup_seconds=setup_seconds,
        energy=initial_energy,force=initial_force,field=field,material=material)
    while True:
        field_norm = max(field[k] for k in FIELDS)
        if material['stationary'] and burst['count'] == 0 and field_norm <= FIELD_TARGET:
            if case != 'benchmark':
                write_endpoint(driver,t,solid,dict(material=material,case=case,starting_sha256=driver.start_sha,
                        smallest_accepted_active_fraction=best_alpha,jointly_matched_energy=total_energy))
                return
            # A native held equilibrium also terminates a short benchmark.
            break
        rematch_due = (burst['count']>=MAX_BURST or
            (burst['count']>0 and material['max_free_label_force']<=max(FORCE_TOL,burst['start_force']/5)) or
            (material['stationary'] and field_norm>FIELD_TARGET))
        if rematch_due:
            before = total_energy
            stamp = perf_counter()
            driver.fields,_,_,field = coherent(driver,t,solid,driver.labels,driver.fields,
                                              field_target=FIELD_TARGET,raw_target=FIELD_TARGET)
            joint,(value,gradient,material)=build_joint(driver,t,solid,registry_update=True)
            physical,_,field=read_fields(driver,t,solid,driver.labels,driver.fields)
            total_energy=physical+material['gauge_energy']
            if total_energy>before+64*abs(np.spacing(before)):
                raise ArithmeticError(f'Joint rematch raised energy: {total_energy-before}')
            if max(field[k] for k in FIELDS)>FIELD_TARGET and burst['count']==0:
                raise ArithmeticError(f'Fixed-label field rematch missed the requested tolerance: {field}')
            burst=dict(count=0,start_force=material['max_free_label_force'])
            record=dict(outer=count,elapsed_seconds=perf_counter()-stamp,total_energy=total_energy,
                        field=field,material=material,burst=burst)
            driver.save_internal('rematched',t,driver.labels,driver.fields,record)
            phases.append(dict(outer=count,wall_seconds=record['elapsed_seconds'],energy_change=total_energy-before))
            log('newton_rematched',outer=count,energy=total_energy,active_force=material['max_free_label_force'],
                field_residual=max(field[k] for k in FIELDS),seconds=record['elapsed_seconds'])
            if material['stationary']:
                continue
        if case=='benchmark' and count-previous_count>=BENCHMARK_STEPS:
            break
        stamp = perf_counter()
        new_labels,new_value,gradient,new_material,search = joint.step(driver.labels,value,gradient,material)
        work = search['energy_change']
        if not np.isfinite(work) or work>=0 or not search['alpha']>1e-4:
            raise ArithmeticError(f'Unresolved accepted step work or fraction: {search}')
        if new_material['min_Jinv_global']<=0 or new_material['chart']['chart_det_global_lower']<=0:
            raise ArithmeticError('Accepted Newton step failed orientation')
        driver.labels = new_labels
        value,material = new_value,new_material
        total_energy += work
        best_alpha = min(search['alpha'],best_alpha or 1.)
        count += 1
        burst['count'] += 1
        row=dict(outer=count,energy=total_energy,work=work,max_force=material['max_free_label_force'],
                 min_Jinv=material['min_Jinv_global'],alpha=search['alpha'],cg_iterations=search['cg_iterations'],
                 cg_reason=search['cg_reason'],eta=search['forcing_eta'],cg_residual=search['cg_relative_residual'],
                 cg_recursive_residual=search['cg_recursive_relative_residual'],forcing_attained=search['forcing_attained'],
                 jacobi=search['jacobi'],
                 gradient_l2=search['gradient_l2_after'],initial_gradient_l2=search['gradient_l2_before'],
                 projections=search['projections'],backtracks=search['backtracks'],timing=search['timing_seconds'],
                 crystal_step_seconds=perf_counter()-stamp)
        force_history.append(row)
        record=dict(**row,burst=burst.copy(),material=material,search=search,
                    cumulative_min_alpha=best_alpha,source_sha256=driver.start_sha)
        driver.save_internal('last_active',t,driver.labels,driver.fields,record)
        log('newton_step',**row)
    if not force_history:
        raise ArithmeticError('No benchmark step from the nonstationary archived state')
    physical,_,field=read_fields(driver,t,solid,driver.labels,driver.fields)
    explicit_energy=physical+material['gauge_energy']
    if abs(total_energy-explicit_energy)>64*abs(np.spacing(total_energy)):
        raise ArithmeticError('Accumulated Newton work disagrees with the endpoint energy')
    # Store the precise exit state even if a field rematch followed the fifth step.
    driver.save_internal('last_active',t,driver.labels,driver.fields,dict(outer=count,burst=burst,
              search=search,material=material,field=field,source_sha256=driver.start_sha,total_energy=explicit_energy))
    prep_free_time = perf_counter()-clock
    total_time = perf_counter()-start_time
    history = read(study.parents[2]/'research/pde_ledger_v3/_measurements/mouth_r4_newton_history.json')
    window = history['windows'][1]
    if window['steps'] != [22,29] or window['step_count'] != 8:
        raise ValueError('Invalid benchmark reference window')
    baseline = window['elapsed_s']/window['step_count']
    final_force = material['max_free_label_force']
    recent = [r['max_force'] for r in force_history[-3:]]
    # A 20% margin distinguishes a useful speedup and residual drop from a
    # marginal change over the fluctuating historical trajectory.
    checks=dict(faster=prep_free_time/len(force_history)<=0.8*baseline,
                lower_force=final_force<=0.8*initial_force and float(np.median(recent))<initial_force,
                lower_l2=material['gradient_l2']<force_history[0]['initial_gradient_l2'],
                lower_energy=explicit_energy<initial_energy,
                positive_chart=material['chart']['chart_det_global_lower']>0)
    result=dict(passed=all(checks.values()),checks=checks,steps=force_history,rematches=phases,
                benchmark_basis='Accepted-step time includes due rematches and checkpoint I/O. Geometry setup is reported separately.',
                setup_seconds=setup_seconds,steady_seconds_per_step=prep_free_time/len(force_history),
                inclusive_seconds_per_step=total_time/len(force_history),baseline_seconds_per_step=baseline,required_time_and_peak_fraction=0.8,
                initial_force=initial_force,final_force=final_force,energy_change=explicit_energy-initial_energy,
                min_Jinv_global=min(r['min_Jinv'] for r in force_history),minimum_alpha=best_alpha,
                outer=count,burst=burst,field=field,material=material,source_sha256=driver.start_sha,
                saved=digest(root/'last_active.npz'))
    atomic_json(root/'result.json',result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',choices=('benchmark','primary','secondary','pair'))
    parser.add_argument('study',type=Path)
    args=parser.parse_args()
    study=args.study.resolve()
    if args.case=='pair':
        primary_start=accepted(study,'benchmark/result.json','benchmark/last_active.npz')
        compare(study,primary_initial=digest(primary_start))
        return
    try:
        run(study,args.case)
    except Exception as exc:
        issue=dict(error=type(exc).__name__,message=str(exc))
        for key in ('diagnostics','details','witness'):
            if hasattr(exc,key):
                issue[key]=getattr(exc,key)
        atomic_json(study/args.case/'ISSUE.json',issue)
        raise


if __name__=='__main__':
    main()
