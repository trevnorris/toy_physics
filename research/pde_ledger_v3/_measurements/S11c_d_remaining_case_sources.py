#!/usr/bin/env python3
"""Materialize the three remaining native cases from the saved parent operands."""
import argparse
import ast
import contextlib
import json
from pathlib import Path
import resource
import shutil
import signal
import sys
import time

import sympy as sp
from sympy.core.function import AppliedUndef
import S11c_d_continuum_grades as grades
import S11c_d_continuum_boundary as boundary
import S11c_d_reduced_action_source_check as source
import S11c_d_reduced_action_assembly_check as assembly_output

f = grades.f
engine = f.engine
PLAN = f.M/'S11c_d_remaining_case_sources_plan.md'
BASELINE = ('LAB_HELD', 'RHO4_CONSTANT')
CASES = tuple((a, d) for a in ('LAB_HELD', 'MATERIAL_ADVECTED')
              for d in ('RHO4_CONSTANT', 'RHOBR_CONSTANT'))
PREFIX = 'REMAINING_CASE_SOURCES'
SCOPE = ('Complete native reduced sources and local/intact-integral assemblies for all four cases. '
         'The accepted first case is reused. Exact integral identities are reuse eligibility only; '
         'new case scattering, continuum responses, controls, pole diagnostics and exports remain.')


def definition_joins(path, seeds):
    current = {n.name: n for n in ast.parse(engine.HERE.read_text()).body
               if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    previous = {n.name: n for n in ast.parse(path.read_text()).body
                if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    names = set(seeds)
    while True:
        expanded = names | {n.id for name in names for n in ast.walk(current[name])
                            if isinstance(n, ast.Name) and n.id in current}
        if expanded == names:
            break
        names = expanded
    joins = {name: name in previous and ast.dump(current[name]) == ast.dump(previous[name])
             for name in sorted(names)}
    f.require(all(joins.values()), ('unchanged consumed native definition closure', str(path), joins))
    return joins


def load(base):
    saved, scp, spath, scheckpoint = grades.packet(
        'S11c_d_reduced_action_source_checkpoint.json', 'reduced-action.pickle')
    assembled, acp, apath, acheckpoint = grades.packet(
        'S11c_d_reduced_action_assembly_checkpoint.json', 'assembly.pickle')
    actions_path = Path(scp['runDirectory'])/'actions.pickle'
    f.require(f.digest(actions_path) == scp['artifacts']['actions.pickle']['sha256'], 'accepted actions hash')
    actions = f.unpickle(actions_path)
    f.require(actions['reducedActionSha256'] == f.digest(spath)
              and assembled['provenance']['ACTIONS_SHA256'] == f.digest(actions_path)
              and assembled['provenance']['REDUCED_ACTION_SHA256'] == f.digest(spath),
              'original reduction/action/assembly provenance')
    engine_name = str(engine.HERE.relative_to(f.ROOT))
    old_source = Path(scp['runDirectory'])/'source'/engine_name
    old_assembly = Path(acp['runDirectory'])/'source'/engine_name
    f.require(f.digest(old_source) == saved['sourceDigests'][engine_name], 'saved reduction producer')
    joins = {
        'reduction': definition_joins(old_source, {'EdgeReduction', 'ReducedPencil', 'DimensionAnalysis',
                                                  'PhysicalMetadata', 'ChannelInput'}),
        'assembly': definition_joins(old_assembly, {'ReducedActionAssembly', 'ReducedPencil',
                                                  'DimensionAnalysis', 'PhysicalMetadata'})}
    for name, sha in saved['sourceDigests'].items():
        if name != engine_name:
            f.require(f.digest(f.ROOT/name) == sha, ('unchanged consumed physics', name))
    for module, checkpoint in ((source, scp), (assembly_output, acp)):
        name = str(Path(module.__file__).resolve().relative_to(f.ROOT))
        f.require(f.digest(f.ROOT/name) == checkpoint['sourceFiles'][name], ('unchanged persistence/emitter', name))
    f.require(saved['schema'] == 1 and
              {(key, tuple(map(str, case))) for key, case in saved['payloads']} ==
              {(key, BASELINE) for key in engine.CLOSED_KEYS}, 'exact baseline payload census')
    case_keys = {}
    for key in engine.CLOSED_KEYS:
        case_keys[key] = {tuple(map(str, case)): case for case, _ in saved['rows'][key]['value']}
        f.require(set(case_keys[key]) == set(CASES), ('all four original cases', key))
    # Freeze imported local code as well as the actual physical roots. Importing
    # these modules does not call their historical numerical constructors.
    paths = {Path(module.__file__).resolve() for module in tuple(sys.modules.values())
             if getattr(module, '__file__', None) and
             Path(module.__file__).resolve().is_relative_to(f.ROOT) and
             Path(module.__file__).suffix == '.py'}
    paths.update((Path(__file__).resolve(), PLAN, scheckpoint, acheckpoint, engine.HERE,
                  f.M/'S11c_d_variable_profile_development_input.json', f.ACCEPTANCE,
                  f.ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md',
                  f.ROOT/'directives/S11c_d_sympy_build_PROGRAM_BRIEF.md',
                  f.M/'S11c_d_focused_completion_plan.md'))
    paths.update(f.ROOT/name for name in saved['sourceDigests'])
    pins = {str(path.relative_to(f.ROOT)): f.digest(path) for path in sorted(paths)}
    operands = {str(path): f.digest(path) for path in (spath, actions_path, apath, old_source, old_assembly)}
    for name in pins:
        target = base/'source'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/name, target)
    for path, name in ((spath, 'accepted-reduced-action.pickle'),
                       (actions_path, 'accepted-actions.pickle'), (apath, 'accepted-assembly.pickle')):
        shutil.copyfile(path, base/name)
        f.require(f.digest(base/name) == operands[str(path)], ('byte-identical baseline copy', name))
    specification = json.loads((f.M/'S11c_d_variable_profile_development_input.json').read_text())
    f.save(base/'inputs.json', {'sourceFiles': pins, 'inputPackets': operands,
        'nativeHelperJoins': joins, 'cases': CASES, 'baseline': BASELINE,
        'input': specification, 'scope': SCOPE, 'sympyVersion': sp.__version__})
    return saved, actions, assembled, case_keys, specification, pins, operands


def validate_assembly(result, actions):
    f.require({(row['COLUMN'], row['ROW']) for row in result['ROWS']} ==
              {(j, i) for j in range(5) for i in range(5)}, 'complete 25-cell action')
    residuals = [value for row in result['ROWS'] for value in
                 (row['RECONSTRUCTION_RESIDUAL'], row['NONLINEAR_OR_AFFINE_REMAINDER'],
                  *row['DERIVATIVE_EXTRACTION_RESIDUALS'])]
    f.require(all(value == 0 for value in residuals), 'native reconstruction and independent coefficient extraction')
    f.require(all(item['FOREIGN_PROBES'] == item['UNSUBSTITUTED_FIELDS'] == 0
                  for item in actions['census']), 'complete source field/probe join')
    return len(residuals)


def baseline_context(saved, actions, assembled):
    r, dimensions = source.restore_context(saved)
    dimensions.__dict__.update(assembled['dimensionState'])
    values = {key: engine.named(payload, 'VALUE') for (key, _), payload in saved['payloads'].items()}
    pencil = engine.ReducedPencil(*(values[key] for key in engine.CLOSED_KEYS), r)
    f.require(pencil.strong == actions['strong'] and pencil.kernel == actions['kernel']
              and pencil.fields == actions['fields'] and pencil.probes == actions['probes'],
              'accepted full baseline source/action identities')
    f.require(source.action_census(actions['columns'], pencil) == actions['census'], 'accepted action census')
    return r, dimensions, pencil


def exact_integral_matches(result, baseline):
    previous = {integral: i for i, integral in enumerate(baseline['NONLOCAL_INTEGRALS'])}
    matches = []
    for i, integral in enumerate(result['NONLOCAL_INTEGRALS']):
        old_index = previous.get(integral)
        if old_index is not None:
            old_integral = baseline['NONLOCAL_INTEGRALS'][old_index]
            f.require(integral == old_integral and integral.limits == old_integral.limits
                      and integral.function == old_integral.function, 'full integral/field/ordered-limit reuse join')
        matches.append({'index': i, 'baselineIndex': old_index,
                        'status': 'EXACT_SOURCE_IDENTITY' if old_index is not None else 'NEW_SOURCE_OPERAND'})
    # An actual changed ordered limit must not inherit the original address.
    if result['NONLOCAL_INTEGRALS']:
        original = result['NONLOCAL_INTEGRALS'][0]
        limits = list(original.limits)
        v, lower, upper = limits[0]
        changed = sp.Integral(original.function, sp.Tuple(v, lower, sp.S.Zero), *limits[1:])
        f.require(upper != 0 and changed != original and changed not in previous, 'wrong-limit address rejection')
    return matches


def construct_case(base, case, saved, case_keys, baseline, specification, progress):
    destination = base/'cases'/'__'.join(case)
    destination.mkdir(parents=True)
    r, dimensions = source.restore_context(saved)
    payloads, branch_bindings, reduction_counts = {}, {}, {}
    constraint_residual_count = 0
    for ordinal, key in enumerate(engine.CLOSED_KEYS):
        case_key = case_keys[key][case]
        original = dict(saved['rows'][key]['value'])[case_key]
        value = engine.named(original, 'VALUE')
        progress('reducing', case=case, row=key)
        reduced, records = r.value(value)
        f.atomic_pickle(destination/f'reduction-{ordinal}.pickle',
                        {'case': case_key, 'key': key, 'sourcePayload': original,
                         'reduced': reduced, 'records': records})
        f.require(tuple(record[0] for record in records) ==
                  tuple(sorted(value.atoms(sp.Integral), key=sp.default_sort_key)), 'every original action integral')
        replay = engine.memo_xreplace(r.branches(r.strip_local(engine.memo_xreplace(
            value, {record[0]: record[1] for record in records}))), r.normal_map)
        f.require(replay == reduced, 'saved literal native reduction replay')
        for original_integral, image, constraints in records:
            for equations, matrix, rhs, solutions, determinant in constraints:
                solutions = list(solutions)
                f.require(len(solutions) == 1 and determinant != 0, 'unique full tangential constraint solution')
                residual = matrix*sp.Matrix(solutions[0])-rhs
                f.require(all(sp.expand(v) == 0 for v in residual) and matrix.det() == determinant,
                          'actual tangential constraint and Jacobian identities')
                constraint_residual_count += len(residual)+1
        payload = r.payload(original, reduced)
        f.atomic_pickle(destination/f'payload-{ordinal}.pickle', payload)
        f.require([str(k) for k, _ in payload] == ['VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M',
                  'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS'], 'full five-slot payload')
        payloads[(key, case_key)] = payload
        branch_bindings[(key, case_key)] = tuple((eq.lhs, eq.rhs) for eq in engine.named(payload, 'COMPUTED_BRANCH_BINDINGS'))
        reduction_counts[key] = len(records)
        progress('reduction_saved', case=case, row=key, integralRecords=len(records))
    reduced_path = destination/'reduced-action.pickle'
    current_roots = {name: f.digest(f.ROOT/name) for name in saved['sourceDigests']}
    engine.save_reduced_action_cache(reduced_path, saved['rows'], payloads, branch_bindings,
                                    r, dimensions, specification, current_roots)
    # This packet records its actual current producer. The original cache and
    # both frozen native definition closures remain explicit in inputs.json.
    values = {key: engine.named(payload, 'VALUE') for (key, _), payload in payloads.items()}
    pencil = engine.ReducedPencil(*(values[key] for key in engine.CLOSED_KEYS), r)
    columns = pencil.columns()
    slab_payload = next(v for (key, _), v in payloads.items() if key == engine.CLOSED_KEYS[0])
    slab_units = engine.payload_units(pencil.slab, engine.named(slab_payload, 'DIMENSION_L_T_M'))
    row_units = [slab_units[('U', i)] for i in range(3)]+[slab_units[(name,)] for name in ('THETA', 'E_W')]
    actions = {'columns': columns, 'strong': pencil.strong, 'kernel': pencil.kernel,
               'fields': pencil.fields, 'probes': pencil.probes,
               'columnUnits': {(j, i): row_units[i] for j in range(5) for i in range(5)},
               'census': source.action_census(columns, pencil),
               'dimensionState': dict(vars(dimensions)), 'reducedActionSha256': f.digest(reduced_path)}
    f.atomic_pickle(destination/'actions.pickle', actions)
    progress('actions_saved', case=case)
    result = engine.ReducedActionAssembly(pencil, columns).construct()
    provenance = {'CASE': case, 'REDUCED_ACTION_SHA256': f.digest(reduced_path),
                  'ACTIONS_SHA256': f.digest(destination/'actions.pickle'),
                  'INPUT_SHA256': f.digest(base/'inputs.json')}
    f.atomic_pickle(destination/'assembly.pickle', {'result': result,
        'dimensionState': dict(vars(dimensions)), 'provenance': provenance})
    residual_count = validate_assembly(result, actions)
    matches = exact_integral_matches(result, baseline['result'])
    common_orders = set(result['LOCAL_MATRICES']) & set(baseline['result']['LOCAL_MATRICES'])
    differences = {order: result['LOCAL_MATRICES'][order]-baseline['result']['LOCAL_MATRICES'][order]
                   for order in common_orders}
    # Full coefficient operands remain in the assemblies; source reuse does
    # not imply unchanged cell coefficients or a reused scattering response.
    comparison = {'integrals': matches, 'localDifferencesAtCommonOrders': differences,
        'newOrders': tuple(sorted(set(result['LOCAL_MATRICES'])-common_orders)),
        'baselineOnlyOrders': tuple(sorted(set(baseline['result']['LOCAL_MATRICES'])-common_orders)),
        'sourceReductionRecords': reduction_counts, 'constraintResidualScalars': constraint_residual_count,
        'assemblyResidualScalars': residual_count}
    f.atomic_pickle(destination/'comparisons.pickle', comparison)
    inputs = engine.ChannelInput(r, specification)
    live = {r.z, r.regulator, r.symbols['sigma_W']}
    missing = tuple(sorted((s for s in engine.dag_free_symbols(columns)-live
                            if s.name not in inputs.parameters), key=sp.default_sort_key))
    f.require(not missing, ('required physical parameters are missing', case, tuple(map(str, missing))))
    summary = {'case': case, 'sourceReductionRecords': reduction_counts,
        'constraintResidualScalars': constraint_residual_count, 'assemblyResidualScalars': residual_count,
        'localOrders': sorted(map(int, result['LOCAL_MATRICES'])),
        'nonlocalIntegrals': len(matches), 'nonlocalCellTerms': sum(len(row['NONLOCAL']) for row in result['ROWS']),
        'exactBaselineIntegrals': sum(item['baselineIndex'] is not None for item in matches),
        'newIntegralOperands': sum(item['baselineIndex'] is None for item in matches),
        'changedLocalEntries': sum(v != 0 for matrix in differences.values() for v in matrix),
        'artifacts': {p.name: {'bytes': p.stat().st_size, 'sha256': f.digest(p)} for p in sorted(destination.glob('*.pickle'))}}
    f.save(destination/'checks.json', summary)
    progress('case_saved', **summary)
    return summary


def emit_case(base, case, baseline=False):
    location = base if baseline else base/'cases'/'__'.join(case)
    prefix = 'accepted-' if baseline else ''
    saved = f.unpickle(location/(prefix+'reduced-action.pickle'))
    actions = f.unpickle(location/(prefix+'actions.pickle'))
    assembled = f.unpickle(location/(prefix+'assembly.pickle'))
    r, dimensions, pencil = baseline_context(saved, actions, assembled)
    case_prefix = PREFIX+'_'+'_'.join(case)
    for (key, case_key), payload in saved['payloads'].items():
        f.require(tuple(map(str, case_key)) == case, 'emitted case source address')
        original = dict(saved['rows'][key]['value'])[case_key]
        value = engine.named(payload, 'VALUE')
        units = engine.payload_units(engine.named(original, 'VALUE'), engine.named(original, 'DIMENSION_L_T_M'))
        engine.physical(case_prefix+'_PAYLOAD_'+key, payload, operands=value, zero_dimensions=units)
    engine.fingerprinted(case_prefix+'_ACTION_COLUMNS', actions['columns'], actions['columnUnits'])
    assembly_output.PREFIX = case_prefix+'_ASSEMBLY'
    assembly_output.emit_result(assembled['result'], actions, (), assembled['provenance'])
    f.require(not dimensions.constraints, ('case dimension constraints', case))
    unresolved = [str(v) for v, unit in dimensions.unknown.items()
                  if v not in dimensions.known and any(s.free_symbols for s in unit)]
    f.require(not unresolved, ('case unresolved units', case, unresolved))


def output_and_replay(base, catalogue):
    def emit_result():
        for case in CASES:
            emit_case(base, case, case == BASELINE)
        boundary.structural_flags(PREFIX+'_CATALOGUE', catalogue)
    engine.EMISSION_LINES.clear()
    engine.PAYLOAD_ENCODER = grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream, contextlib.redirect_stdout(stream):
        emit_result()
        keys = {tag: 's11cdRemainingCaseSource'+str(i) for i, tag in enumerate(engine.EMISSION_LINES)
                if not tag.startswith('PY_S11CD_METADATA_')}
        boundary.structural_flags(PREFIX+'_WRITE_KEYS', keys)
        index = engine.emission_index(engine.EMISSION_LINES)
        boundary.structural_flags(PREFIX+'_EMISSION_LINES', index)
    entries = {}
    for line in grades.decoded_lines(base/'full.out'):
        tag, _, body = line.rstrip('\n').partition(': ')
        f.require(tag not in entries, 'unique all-case source tags')
        entries[tag] = grades._restore(body)
    previous_emit, seen = engine.emit, set()
    def compare(name, value):
        tag = 'PY_S11CD_'+name
        f.require(tag not in seen and entries.get(tag) == engine.cas(value), ('full emitted payload replay', tag))
        seen.add(tag)
    engine.emit = compare
    try:
        emit_result()
        boundary.structural_flags(PREFIX+'_WRITE_KEYS', keys)
        boundary.structural_flags(PREFIX+'_EMISSION_LINES', index)
    finally:
        engine.emit = previous_emit
    f.require(seen == set(entries) and len(keys) == len(set(keys.values()))
              and not set(keys.values()) & set(engine.IMPORT_KEYS), 'full case/key census')
    paths = 0
    for tag, body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):
            continue
        for path, fields in body:
            info = {str(k): value for k, value in fields}
            unit = info['DIMENSION_L_T_M']
            f.require(isinstance(unit, sp.Tuple) and len(unit) == 3 and
                      all(not v.free_symbols for v in unit) and
                      'MULTIGRADE' in info and 'EPSILON_LAMBDA_SUPPORT' in info, 'complete source metadata')
            paths += 1
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    grades.restore_emission_index({str(k): v for k, v in entries[final]}, list(entries)[:list(entries).index(final)])
    return {'tagCount': len(entries), 'writeKeys': len(keys), 'metadataPaths': paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=('preflight', 'construct'), required=True)
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve()
    base.relative_to(f.STORE)
    base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    started = time.monotonic()
    def timeout(*_):
        raise TimeoutError('remaining-case source budget; retain completed packets')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(900)
    def progress(stage, **values):
        with (base/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'stage': stage, 'wallSeconds': time.monotonic()-started, **values})+'\n')
    saved, actions, assembled, case_keys, specification, pins, operands = load(base)
    baseline_context(saved, actions, assembled)
    baseline_residuals = validate_assembly(assembled['result'], actions)
    baseline_matches = exact_integral_matches(assembled['result'], assembled['result'])
    f.require(all(item['baselineIndex'] == i for i, item in enumerate(baseline_matches)), 'every baseline integral address')
    preflight = {'mode': args.mode, 'baselineResidualScalars': baseline_residuals,
                 'baselineIntegralAddresses': len(baseline_matches), 'cases': CASES,
                 'baselinePacketsCopiedExactly': True, 'newNumericalNodes': 0}
    f.save(base/'preflight.json', preflight)
    progress('preflight_complete')
    summaries, replay = [], {}
    if args.mode == 'construct':
        for case in CASES[1:]:
            summaries.append(construct_case(base, case, saved, case_keys, assembled, specification, progress))
        catalogue = {'cases': summaries, 'baseline': preflight, 'sourceFiles': pins,
                     'inputPackets': operands, 'scope': SCOPE}
        f.atomic_pickle(base/'remaining-case-sources.pickle', catalogue)
        before = {str(p.relative_to(base)): f.digest(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
        f.save(base/'packet-inventory-before-emission.json', before)
        replay = output_and_replay(base, catalogue)
        f.require(all(f.digest(base/name) == sha for name, sha in before.items()), 'all pre/post source packet identities')
    f.require(all(f.digest(f.ROOT/name) == sha and f.digest(base/'source'/name) == sha for name, sha in pins.items()), 'current/frozen sources unchanged')
    f.require(all(f.digest(Path(name)) == sha for name, sha in operands.items()), 'accepted operands unchanged')
    checks = {'runDirectory': str(base), 'mode': args.mode, 'sourceFiles': pins,
        'inputPackets': operands, 'preflight': preflight, 'cases': summaries, **replay,
        'artifacts': {str(p.relative_to(base)): {'bytes': p.stat().st_size, 'sha256': f.digest(p)}
                      for p in sorted(base.rglob('*')) if p.suffix in ('.pickle', '.out')
                      and 'source' not in p.relative_to(base).parts},
        'wallSeconds': time.monotonic()-started,
        'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'scope': SCOPE}
    f.save(base/'checks.json', checks)
    progress('complete')
    signal.alarm(0)
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
