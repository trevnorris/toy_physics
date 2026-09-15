#!/usr/bin/env python3
"""Evaluate source and assembled actions on a bounded numerical test domain."""
import argparse
import ast
import contextlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import time

import numpy as np
import sympy as sp
from sympy.core.function import AppliedUndef

import S11c_d_reduced_action_source_check as source
import S11c_d_variable_profile_input_check as input_check
from S11c_d_reduced_action_source_check import ROOT, STORE, engine, digest, save, atomic_pickle
from S11c_d_output_codec import decoded_lines, restore_emission_index
from ledger_fold import _restore

MEASUREMENTS = ROOT/'_measurements'
ASSEMBLY = MEASUREMENTS/'S11c_d_reduced_action_assembly_checkpoint.json'
INPUT = MEASUREMENTS/'S11c_d_variable_profile_development_input.json'
INPUT_JOIN = MEASUREMENTS/'S11c_d_variable_profile_input_checkpoint.json'
PREFIX = 'NUMERICAL_ACTION_LAB_HELD_RHO4_CONSTANT'
ADDITIONS = ('NumericalReducedAction', 'BoundedActionQuadrature')
SOURCES = (*engine.BUILD_INPUT_PATHS, Path(__file__).resolve(), Path(source.__file__).resolve(),
    Path(input_check.__file__).resolve(), ASSEMBLY, INPUT, INPUT_JOIN,
    MEASUREMENTS/'S11c_d_numerical_action_plan.md',
    MEASUREMENTS/'S11c_d_variable_profile_matching_plan.md',
    *(MEASUREMENTS/('S11c_d_'+name+'.json') for name in (
        'channel_preflight_input', 'variable_profile_development_input_proposal',
        'variable_profile_parameter_checkpoint', 'variable_profile_parameter_inventory',
        'matching_channels_checkpoint')))


def load():
    checkpoint = json.loads(ASSEMBLY.read_text())
    base = Path(checkpoint['runDirectory'])
    for name, item in checkpoint['artifacts'].items():
        if digest(base/name) != item['sha256']:
            raise ValueError(('accepted assembly artifact changed', name))
    publication = checkpoint['publication']
    if digest(ROOT/publication['path']) != publication['sha256']:
        raise ValueError('accepted assembly publication changed')
    for name, sha in checkpoint['sourceFiles'].items():
        if digest(base/'source'/name) != sha:
            raise ValueError(('assembly source snapshot changed', name))
        if name != str(engine.HERE.relative_to(ROOT)) and digest(ROOT/name) != sha:
            raise ValueError(('assembly consumed source changed', name))
    frozen = ast.parse((base/'source'/engine.HERE.relative_to(ROOT)).read_text())
    current = ast.parse(engine.HERE.read_text())
    current.body = [n for n in current.body if getattr(n, 'name', None) not in ADDITIONS]
    if ast.dump(current) != ast.dump(frozen):
        raise ValueError('native engine changed beyond numerical action classes')
    source_path = MEASUREMENTS/'S11c_d_reduced_action_source_checkpoint.json'
    source_record = json.loads(source_path.read_text())
    source_base = Path(source_record['runDirectory'])
    for name, item in source_record['artifacts'].items():
        if digest(source_base/name) != item['sha256']:
            raise ValueError(('source action cache changed', name))
    with (source_base/'reduced-action.pickle').open('rb') as stream:
        packet = pickle.load(stream)
    with (source_base/'actions.pickle').open('rb') as stream:
        actions = pickle.load(stream)
    with (base/'assembly.pickle').open('rb') as stream:
        assembly = pickle.load(stream)
    if (assembly['provenance'] != checkpoint['provenance'] or
            actions['reducedActionSha256'] != digest(source_base/'reduced-action.pickle') or
            checkpoint['provenance']['ACTIONS_SHA256'] != digest(source_base/'actions.pickle')):
        raise ValueError('assembly/action/reduction packet provenance join')
    reduction, dimensions = source.restore_context(packet)
    dimensions.__dict__.update(assembly['dimensionState'])
    values = {key: engine.named(body, 'VALUE') for (key, _), body in packet['payloads'].items()}
    pencil = engine.ReducedPencil(*(values[key] for key in engine.CLOSED_KEYS), reduction)
    if pencil.strong != actions['strong'] or pencil.fields != actions['fields'] or pencil.kernel != actions['kernel']:
        raise ValueError('full reduced source/pencil join')
    selected = input_check.check()
    accepted = json.loads(INPUT_JOIN.read_text())
    if selected != {k: v for k, v in accepted.items() if k != 'instrument'}:
        raise ValueError('approved input projection changed')
    for name, data in selected['addedCoefficientUnits'].items():
        unit = dimensions.measure(reduction.symbols[name])
        if tuple(map(str, unit)) != tuple(data):
            raise ValueError(('approved coefficient unit changed', name))
    return checkpoint, actions, pencil, assembly['result'], selected


def scalar(value):
    return complex(np.asarray(value).item())


def number(value):
    if isinstance(value, np.ndarray):
        return number(value.tolist())
    if isinstance(value, (tuple, list)):
        return sp.Tuple(*(number(v) for v in value))
    if isinstance(value, complex):
        return sp.Float(value.real, 17)+sp.I*sp.Float(value.imag, 17)
    return engine.cas(value)


def progress(base, stage, **fields):
    with (base/'progress.jsonl').open('a') as stream:
        stream.write(json.dumps({'stage': stage, 'wallClock': time.time(), **fields})+'\n')


def construct(base, actions, pencil, assembly):
    adapter = engine.NumericalReducedAction(pencil, assembly, json.loads(INPUT.read_text()))
    r = pencil.r
    local = adapter.local_matrices()
    bound_integrals = tuple(adapter.bind(g) for g in assembly['NONLOCAL_INTEGRALS'])
    live = {r.z, r.regulator}
    all_bound = sp.Tuple(*local.values(), *bound_integrals)
    missing = engine.dag_free_symbols(all_bound)-live
    functions = all_bound.atoms(AppliedUndef)
    foreign = {f.func for f in functions}-set(pencil.probes)
    if missing or foreign or all_bound.has(sp.Limit):
        raise ValueError(('incomplete numerical binding', tuple(map(str, missing)), tuple(map(str, foreign))))
    atomic_pickle(base/'bound-operators.pickle', {'localMatrices': local,
        'nonlocalIntegrals': bound_integrals, 'input': adapter.input.specification,
        'origin': adapter.input.origin, 'profileLimits': adapter.input.limits})
    progress(base, 'bound_operators', localMatrices=len(local), nonlocalIntegrals=len(bound_integrals))
    positions = (-7.0, 0.0, 6.0)
    tests = ((sp.Integer(4), sp.Rational(3, 10)), (sp.Integer(7), -sp.Rational(9, 20)))
    settings = (
        {'momentumBound': 2.0, 'momentumNodes': 12, 'sourceBound': 32.0,
         'sourceNodes': 48, 'profileBound': 10.0, 'profileNodes': 48, 'regulator': 0.2},
        {'momentumBound': 2.0, 'momentumNodes': 16, 'sourceBound': 32.0,
         'sourceNodes': 64, 'profileBound': 10.0, 'profileNodes': 64, 'regulator': 0.2})
    momenta = tuple(r.normal_map[group[2]] for group in r.momentum_groups)
    results = []
    for test_index, (width, momentum) in enumerate(tests):
        # Coefficients in the supplied reference frame; all five field units
        # and all equation units are inherited separately from the source.
        field = lambda z: sp.exp(-(z/width)**2+sp.I*momentum*z)
        direct = adapter.direct_actions(field)
        assembled = adapter.assembled_actions(field)
        atomic_pickle(base/f'test-{test_index}-operands.pickle', {'width': width,
            'momentum': momentum, 'field': field(r.z), 'direct': direct, 'assembled': assembled})
        progress(base, 'test_operands', test=test_index)
        for grid_index, setting in enumerate(settings):
            domains = {k: (-setting['momentumBound'], setting['momentumBound'], setting['momentumNodes'])
                       for k in momenta}
            domains[r.zp] = (-setting['sourceBound'], setting['sourceBound'], setting['sourceNodes'])
            domains[r.xi] = (-setting['profileBound'], setting['profileBound'], setting['profileNodes'])
            assembled_quadrature = engine.BoundedActionQuadrature(domains, contraction='sum')
            direct_quadrature = engine.BoundedActionQuadrature(domains, contraction='dot')
            arrays = {name: np.empty((len(positions), 5, 5), dtype=complex)
                      for name in ('local', 'nonlocal', 'direct', 'residual')}
            contributions = []
            for position_index, position in enumerate(positions):
                environment = {r.z: position, r.regulator: setting['regulator']}
                for record in assembled:
                    j, i = record['COLUMN'], record['ROW']
                    index = (position_index, j, i)
                    local_value = scalar(assembled_quadrature.evaluate(record['LOCAL'], environment))
                    terms = [scalar(assembled_quadrature.evaluate(c, environment))*
                             scalar(assembled_quadrature.evaluate(g, environment))
                             for c, g in record['NONLOCAL']]
                    nonlocal_value = sum(terms)
                    direct_value = scalar(direct_quadrature.evaluate(direct[j][i], environment))
                    arrays['local'][index] = local_value
                    arrays['nonlocal'][index] = nonlocal_value
                    arrays['direct'][index] = direct_value
                    arrays['residual'][index] = local_value+nonlocal_value-direct_value
                    contributions.append({'index': index, 'terms': terms})
                progress(base, 'position_evaluated', test=test_index, grid=grid_index, position=position_index)
            # Mutations act on an actual selected nonlocal integration measure
            # and a local coefficient; both original and modified operands stay.
            candidates = [(abs(v), record['index'], n) for record in contributions
                          for n, v in enumerate(record['terms'])]
            _, index, term_index = max(candidates)
            pi, j, i = index
            environment = {r.z: positions[pi], r.regulator: setting['regulator']}
            record = next(a for a in assembled if (a['COLUMN'], a['ROW']) == (j, i))
            c, g = record['NONLOCAL'][term_index]
            mutated_quadrature = engine.BoundedActionQuadrature(domains, weight_scale=1.001)
            before = scalar(assembled_quadrature.evaluate(c, environment))*scalar(assembled_quadrature.evaluate(g, environment))
            after = scalar(mutated_quadrature.evaluate(c, environment))*scalar(mutated_quadrature.evaluate(g, environment))
            measure_mutation = {'index': index, 'term': term_index, 'weightScale': 1.001,
                                'original': before, 'modified': after, 'difference': after-before}
            local_candidates = []
            for record in assembly['ROWS']:
                for order, coefficient in record['LOCAL'].items():
                    term = adapter.bind(coefficient)*sp.diff(field(r.z), r.z, order)
                    value = scalar(assembled_quadrature.evaluate(term, environment))
                    local_candidates.append((abs(value), record['COLUMN'], record['ROW'],
                                             int(order), coefficient, value))
            _, j, i, order, coefficient, before = max(local_candidates, key=lambda x: x[0])
            modified_coefficient = adapter.bind(coefficient)*sp.Rational(1001, 1000)
            after = scalar(assembled_quadrature.evaluate(modified_coefficient*sp.diff(field(r.z), r.z, order), environment))
            coefficient_mutation = {'index': (pi, j, i), 'order': order,
                'originalCoefficient': adapter.bind(coefficient), 'modifiedCoefficient': modified_coefficient,
                'original': before, 'modified': after, 'difference': after-before}
            item = {'test': test_index, 'grid': grid_index, 'width': width, 'momentum': momentum,
                'positions': positions, 'settings': setting, 'domains': domains, 'arrays': arrays,
                'contributions': contributions, 'measureMutation': measure_mutation,
                'coefficientMutation': coefficient_mutation,
                'assembledEvaluatedIntegrals': tuple(sorted(assembled_quadrature.evaluated_integrals, key=sp.default_sort_key)),
                'directEvaluatedIntegrals': tuple(sorted(direct_quadrature.evaluated_integrals, key=sp.default_sort_key))}
            atomic_pickle(base/f'test-{test_index}-grid-{grid_index}.pickle', item)
            results.append(item)
            progress(base, 'grid_saved', test=test_index, grid=grid_index)
    refinements = [{'test': t, 'difference': results[2*t+1]['arrays']['direct']-results[2*t]['arrays']['direct']}
                   for t in range(len(tests))]
    return {'results': results, 'refinements': refinements, 'localMatrices': local,
            'boundIntegrals': bound_integrals, 'origin': adapter.input.origin,
            'profileLimits': adapter.input.limits, 'profiles': adapter.input.profiles,
            'positions': positions, 'tests': tests}


def emit_result(result, actions, pencil, assembly, provenance):
    dimensions = engine.PHYSICAL_METADATA.dimensions
    units = [actions['columnUnits'][(0, i)] for i in range(5)]
    metadata = engine.FullPencilModes.__new__(engine.FullPencilModes); metadata.r = pencil.r
    def numeric(suffix, value, unit):
        body = number(value)
        engine.emit(PREFIX+'_'+suffix, metadata.compact_fingerprint(body))
        engine.emit('METADATA_'+PREFIX+'_'+suffix, metadata.numeric_metadata(body, unit))
    engine.fingerprinted(PREFIX+'_SYMBOLIC_SOURCE_ACTIONS', actions['columns'],
                         {(j, i): units[i] for j in range(5) for i in range(5)})
    engine.fingerprinted(PREFIX+'_SYMBOLIC_ORDERED_INTEGRALS', engine.cas(assembly['NONLOCAL_INTEGRALS']))
    engine.physical(PREFIX+'_PROVENANCE', provenance)
    engine.physical(PREFIX+'_PROFILE_FORMULAS', result['profiles'])
    parameters = json.loads(INPUT.read_text())['parameters']
    symbols = dict(pencil.r.symbols); symbols.update({s.name: s for s in pencil.r.tangents})
    names = sorted(parameters)
    engine.physical(PREFIX+'_PARAMETER_NAMES', names)
    numeric('PARAMETER_VALUES', [sp.Rational(parameters[n]) for n in names],
            lambda p: dimensions.measure(symbols[names[p[0]]]))
    numeric('HOMOTOPY_ORIGIN', tuple(result['origin'].values()), lambda p: dimensions.zero)
    numeric('TEST_POSITIONS', result['positions'], lambda p: dimensions.measure(pencil.r.z))
    numeric('TEST_WIDTHS_AND_MOMENTA', result['tests'],
            lambda p: tuple((1 if p[1] == 0 else -1)*v for v in dimensions.measure(pencil.r.z)))
    for item in result['results']:
        suffix = str(item['test'])+'_'+str(item['grid'])
        for name, value in item['arrays'].items():
            numeric(name.upper()+'_'+suffix, value, lambda p: units[p[2]])
        domain_names = sorted(item['domains'], key=sp.default_sort_key)
        engine.physical(PREFIX+'_DOMAIN_VARIABLES_'+suffix, domain_names)
        numeric('DOMAIN_BOUNDS_AND_COUNTS_'+suffix, [item['domains'][n] for n in domain_names],
                lambda p: dimensions.zero if p[1] == 2 else dimensions.measure(domain_names[p[0]]))
        numeric('FINITE_ABEL_REGULATOR_'+suffix, item['settings']['regulator'],
                lambda p: dimensions.measure(pencil.r.regulator))
        for name in ('measureMutation', 'coefficientMutation'):
            mutation = item[name]
            numeric(name.upper()+'_'+suffix,
                    tuple(mutation[k] for k in ('original', 'modified', 'difference')),
                    lambda p: units[mutation['index'][2]])
        numeric('NONLOCAL_CONTRIBUTIONS_'+suffix, [r['terms'] for r in item['contributions']],
                lambda p: units[item['contributions'][p[0]]['index'][2]])
    for refinement in result['refinements']:
        numeric('GRID_CHANGE_'+str(refinement['test']), refinement['difference'], lambda p: units[p[2]])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve(); base.relative_to(STORE)
    base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    pins = {str(p.relative_to(ROOT)): digest(p) for p in SOURCES}
    for p in SOURCES:
        target = base/'source'/p.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, target)
    checkpoint, actions, pencil, assembly, input_join = load()
    provenance = {'ASSEMBLY_CHECKPOINT_SHA256': digest(ASSEMBLY),
        'ASSEMBLY_PACKET_SHA256': checkpoint['artifacts']['assembly.pickle']['sha256'],
        'APPROVED_INPUT_SHA256': digest(INPUT), 'INPUT_PROJECTION_SHA256': digest(INPUT_JOIN),
        'BASE_INPUT_SHA256': input_join['files']['base']['sha256']}
    save(base/'preflight.json', {'sourceFiles': pins, 'provenance': provenance})
    result = construct(base, actions, pencil, assembly)
    atomic_pickle(base/'actions.pickle', {'result': result, 'provenance': provenance,
        'dimensionState': dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    before = digest(base/'actions.pickle')
    with (base/'full.out').open('x') as stream, contextlib.redirect_stdout(stream):
        emit_result(result, actions, pencil, assembly, provenance)
        keys = {tag: 's11cd'+''.join(w.title() for w in tag.removeprefix('PY_S11CD_').split('_'))
                for tag in engine.EMISSION_LINES if not tag.startswith('PY_S11CD_METADATA_')}
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        index = engine.emission_index(engine.EMISSION_LINES)
        zero_units = {p: (0, 0, 0) for p, _ in engine.leaves(engine.cas(index))}
        engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=zero_units)
    entries = {}
    for line in decoded_lines(base/'full.out'):
        tag, _, body = line.rstrip('\n').partition(': ')
        if tag in entries:
            raise ValueError('duplicate numerical action emission')
        entries[tag] = _restore(body)
    seen = set(); original = engine.emit
    def compare(name, value):
        tag = 'PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value):
            raise ValueError(('numerical saved/emitted mismatch', tag))
        seen.add(tag)
    engine.emit = compare
    try:
        emit_result(result, actions, pencil, assembly, provenance)
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=zero_units)
    finally:
        engine.emit = original
    if seen != set(entries) or len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('numerical action emission/write-key census')
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    restore_emission_index({str(k): v for k, v in entries[final]}, list(entries)[:list(entries).index(final)])
    metadata_paths = 0
    for tag, body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):
            continue
        for record in body:
            # Numerical metadata groups paths with the same unit/grade;
            # symbolic metadata stores a descriptor for each path.
            if len(record) == 2 and isinstance(record[0], sp.Tuple):
                fields = {str(k): v for k, v in record[1]}
                count = 1
            else:
                fields = {str(k): v for k, v in record}
                count = len(fields['PATHS'])
            unit = fields['DIMENSION_L_T_M']
            if not isinstance(unit, sp.Tuple) or len(unit) != 3 or any(v.free_symbols for v in unit):
                raise ValueError(('unresolved numerical action units', tag))
            if any(k not in fields for k in ('MULTIGRADE', 'EPSILON_LAMBDA_SUPPORT')):
                raise ValueError(('missing numerical action grades', tag))
            metadata_paths += count
    norms = [{'test': item['test'], 'grid': item['grid'],
        'maxResidual': float(np.max(np.abs(item['arrays']['residual']))),
        'maxScaledResidual': float(np.max(np.abs(item['arrays']['residual'])/(1+np.abs(item['arrays']['direct'])))),
        'measureMutationDifference': abs(item['measureMutation']['difference']),
        'coefficientMutationDifference': abs(item['coefficientMutation']['difference']),
        'assembledIntegralCount': len(item['assembledEvaluatedIntegrals']),
        'directIntegralCount': len(item['directEvaluatedIntegrals'])} for item in result['results']]
    summary = {'runDirectory': str(base), 'sourceFiles': pins, 'provenance': provenance,
        'inputProjection': input_join, 'norms': norms,
        'gridChanges': [float(np.max(np.abs(r['difference']))) for r in result['refinements']],
        'actionComponentEvaluations': sum(item['arrays']['direct'].size for item in result['results']),
        'nonlocalTermEvaluations': sum(len(r['terms']) for item in result['results'] for r in item['contributions']),
        'localDerivativeOrders': [int(n) for n in result['localMatrices']],
        'boundIntegralCount': len(result['boundIntegrals']), 'tagCount': len(entries), 'writeKeyCount': len(keys),
        'metadataPaths': metadata_paths,
        'packetSha256BeforeEmission': before, 'packetSha256AfterEmission': digest(base/'actions.pickle'),
        'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts': {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)}
                      for p in sorted(base.iterdir()) if p.suffix in ('.pickle', '.out')},
        'scope': 'Complete five-field source/assembly numerical action comparisons on the recorded finite quadrature domains and positive Abel regulator. '
                 'Shared quadrature nodes with independent sum/dot contractions test assembly; physical tail, infinite-domain, and Abel weak limits remain unresolved. '
                 'No scattering or profile-frequency pole solve.'}
    save(base/'checks.json', summary)
    if pins != {str(p.relative_to(ROOT)): digest(p) for p in SOURCES} or before != digest(base/'actions.pickle'):
        raise ValueError('numerical action source/packet changed during execution')
    if any(n['maxScaledResidual'] > 1e-9 or n['measureMutationDifference'] <= 1e-14 or
           n['coefficientMutationDifference'] <= 1e-14 for n in norms):
        raise ValueError('numerical action residual or mutation sensitivity requires inspection')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
