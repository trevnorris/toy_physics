"""Zero-step tests of the step-382 frozen functional and production CG operator."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import statistics
import sys
from time import perf_counter

import numpy as np


TOOLS = Path(__file__).resolve().parent
SCRATCH = TOOLS.parents[2] / '_scratch/s11c'
SOURCE = SCRATCH / 'r4_projector_repair01'
VECTORS = SCRATCH / 'r4_step383_audit01/cg_vectors.npz'
STEP_SHA = '346fc15d7dc678db51b420bbfde7e07ba8cbffac68aabc9d5cdbd4dd4b149f4c'
VECTOR_SHA = 'b3a633b19af74acfac934330415fe165427f5a319e6cf4c6f3aaa698233d57d0'
HISTORICAL_MEDIAN_SECONDS = 1.0687327338171992


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def store(path, report):
    with path.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')


def open_state(output, *, pinned, fields=False):
    bank = SOURCE / 'pinned/software/stage1_solver/tools' if pinned else TOOLS
    sys.path.insert(0, str(bank))
    from mouth_r4_aligned_crystal import AlignedCrystal
    from mouth_r4_collar166_state import load_held
    from mouth_r4_cubic_held import approved_reference
    from mouth_r4_eulerian_release import order_grid

    if Path(inspect.getfile(AlignedCrystal)).resolve().parent != bank:
        raise RuntimeError('Incorrect source provenance for the material evaluator')
    if pinned:
        for filename, expected in read(SOURCE / 'pinned/manifest.json').items():
            if digest(bank / filename) != expected:
                raise RuntimeError(f'Changed step-382 reference implementation in {filename}')
    driver, trial, solid, _ = load_held(SOURCE, output, 1e-8, solve_fields=fields)
    approved_reference(trial.b)
    record = read(SOURCE / 'primary/step_00382.json')
    path = SOURCE / 'primary' / record['file']
    if not record['outer'] == 382 or not digest(path) == STEP_SHA == record['sha256']:
        raise ValueError('Step-382 provenance failed')
    with np.load(path, allow_pickle=False) as data:
        if int(data['radius']) != 4 or int(data['winding']) != 1 or not np.array_equal(driver.x, data['x']):
            raise ValueError('Frozen radius, winding or six shape coordinates differ')
        driver.labels, driver.fields = data['labels'], data['field']
    driver.labels.flags.writeable = driver.fields.flags.writeable = False
    joint = AlignedCrystal(driver, trial, solid, driver.labels, order_grid(trial, driver.fields))
    return driver, trial, joint, record


def frozen(output, *, pinned):
    driver, trial, joint, record = open_state(output, pinned=pinned, fields=not pinned)
    from mouth_r4_collar166_state import readout
    from mouth_r4_projector_replay import force_bank
    material, bank = force_bank(joint, driver.labels)
    if material != record['material']:
        raise ArithmeticError('Frozen material scalars changed')
    path = output / 'frozen_reference.npz'
    if pinned:
        with path.open('xb') as stream:
            np.savez(stream, **bank)
        store(output / 'reference.json', dict(state_sha256=STEP_SHA, source_manifest=digest(SOURCE/'pinned/manifest.json'),
                                              force_bank_sha256=digest(path)))
        return
    if digest(path) != read(output/'reference.json')['force_bank_sha256'] or digest(VECTORS) != VECTOR_SHA:
        raise ValueError('Changed gradient reference')
    with np.load(path, allow_pickle=False) as archive:
        if set(bank) != set(archive.files) or any(bank[k].shape != archive[k].shape for k in bank):
            raise ValueError('Changed crystal force registry')
        deltas = {k:float(np.max(abs(bank[k] - archive[k]))) for k in bank}
    with np.load(VECTORS, allow_pickle=False) as vectors:
        deltas['total_free_vs_archived_cg_rhs'] = float(np.max(abs(bank['total_free'] + vectors['g'])))
    measured = readout(driver, trial, joint.solid)
    prior = read(SOURCE/'stop_readout/result.json')
    if prior['source_sha256'] != STEP_SHA:
        raise ValueError('Wrong frozen-state measurement')
    old = prior['readout']
    keys = ('total_energy', 'physical_energy', 'crystal_energy', 'gauge_energy', 'peak_active')
    passed = measured == old and all(v == 0 for v in deltas.values())
    store(output/'invariance.json', dict(passed=passed, state_sha256=STEP_SHA,
          reference_sha256=digest(path), force_counts={k:v.size for k,v in bank.items()}, force_max_abs_delta=deltas,
          scalar_deltas={k:measured[k]-old[k] for k in keys}, currents=measured['current'],
          current_deltas={k:v-old['current'][k] for k,v in measured['current'].items()},
          field_deltas={k:v-old['field'][k] for k,v in measured['field'].items() if isinstance(v,float)},
          complete_readout_equal=measured==old, frozen=measured, accepted_steps=0))
    if not passed:
        raise ArithmeticError('Frozen-state invariance failed')


def norm(v):
    return np.sqrt(np.sum(np.square(v, dtype=np.longdouble), dtype=np.longdouble))


def dot(u,v):
    return np.sum(np.multiply(u,v,dtype=np.longdouble), dtype=np.longdouble)


class EndOfLinearSolve(Exception):
    """Raised before the production solver reaches the outer line search."""


def cg(output):
    driver, _, joint, record = open_state(output, pinned=False)
    from mouth_r4_newton_crystal import JointCrystal
    from mouth_r4_eulerian_relax_labels import InfeasibleLabels
    value, g, material = joint.evaluate(driver.labels)
    if material != record['material'] or digest(VECTORS) != VECTOR_SHA:
        raise ValueError('CG replay source differs from step 382')
    with np.load(VECTORS, allow_pickle=False) as bank:
        if not np.array_equal(g,bank['g']) or not np.array_equal(joint.active,bank['active']):
            raise ValueError('Wrong packed residual or row ordering')
        u, v = bank['p30'], bank['p35']
    labels = driver.labels
    actions = [joint.hvp(labels, w) for w in (u,v,u+v,v-u)]
    if joint.last_hvp['relative_step'] != 0.:
        raise ValueError('The physical action still has a finite-difference width')
    hu,hv,hs,hd = (a.astype(np.longdouble) for a in actions)
    errors = dict(add=float(norm(hs-(hu+hv))/norm(hs)), subtract=float(norm(hd-(hv-hu))/norm(hd)),
                  symmetry=float(abs(dot(u,hv)-dot(v,hu))/(norm(u)*norm(hv))))
    algebra = dict(errors=errors, limits=dict.fromkeys(errors,1e-12), passed=max(errors.values())<=1e-12,
                   vector_sha256=VECTOR_SHA, state_sha256=STEP_SHA, implementation=joint.last_hvp, accepted_steps=0)
    store(output/'algebra.json', algebra)
    # The CG replay is itself a zero-step control. It can isolate an algebraic failure without any descent.
    joint.initial_force = read(SOURCE/'primary/matched_00377.json')['material']['max_free_label_force']
    joint.previous_norm = record['search']['gradient_l2_before']
    lines, start = inspect.getsourcelines(JointCrystal.step)
    stop_line = start + next(i for i,line in enumerate(lines) if line.strip()=='tests=[]')
    direct_line = start + next(i for i,line in enumerate(lines) if line.strip()=='if true_rel<=eta:')
    metrics = {}
    checkpoints = []
    durations = []
    jacobi_durations = []
    physical_action = joint.hvp
    def timed(labels, vec):
        before = perf_counter()
        action = physical_action(labels,vec)
        elapsed = perf_counter()-before
        if sys._getframe(1).f_code.co_name == 'inverse_diagonal':
            jacobi_durations.append(elapsed)
        else:
            durations.append(elapsed)
        return action
    joint.hvp = timed
    def observe(frame,event,arg):
        if frame.f_code is not JointCrystal.step.__code__:
            return None
        if event=='line' and frame.f_lineno in (direct_line,stop_line):
            data = frame.f_locals
            if frame.f_lineno==direct_line:
                direct = data['true_residual']
                recurrent = data['residual']
                checkpoints.append(dict(iteration=data['cg'], recursive=data['rel'], direct=data['true_rel'],
                        normalized_vector_gap=float(norm(direct-recurrent)/norm(g))))
            else:
                metrics.update(cg_iterations=data['cg'], eta=data['eta'], reason=data['reason'],
                               replacements=data['replacements'], explicit=data['actual'][0], recursive=data['rel'],
                               jacobi=data['jacobi'])
                raise EndOfLinearSolve()
        return observe
    sys.settrace(observe)
    try:
        joint.step(labels,value,g,material)
    except EndOfLinearSolve:
        pass
    except InfeasibleLabels as exc:
        metrics.update(failure=str(exc),witness=exc.witness)
    else:
        raise AssertionError('CG-only trace stop failed')
    finally:
        sys.settrace(None)
    median = statistics.median(durations)
    passed = (metrics.get('reason') == 'CG forcing attained' and checkpoints
              and all(abs(c['direct']-c['recursive'])<=1e-10 for c in checkpoints)
              and metrics['explicit']<=metrics['eta'])
    store(output/'cg_gate.json',dict(state_sha256=STEP_SHA,accepted_steps=0,passed=bool(passed), metrics=metrics,
          tolerance_checks=checkpoints,operator_calls=len(durations), per_call_seconds=durations,
          preconditioner_calls=len(jacobi_durations), preconditioner_product_seconds=jacobi_durations,
          median_seconds=median, archived_median_seconds=HISTORICAL_MEDIAN_SECONDS,
          timing_ratio=median/HISTORICAL_MEDIAN_SECONDS))
    if not passed or not algebra['passed']:
        raise ArithmeticError('Frozen additivity/symmetry or CG gate failed')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('reference','frozen','cg'))
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if args.operation=='cg':
        cg(args.output)
    else:
        frozen(args.output,pinned=args.operation=='reference')
