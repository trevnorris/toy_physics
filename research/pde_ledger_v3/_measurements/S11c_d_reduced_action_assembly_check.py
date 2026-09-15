#!/usr/bin/env python3
"""Construct the complete local-jet/intact-integral assembly from saved actions."""
import argparse
import ast
import contextlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import time

import sympy as sp

import S11c_d_reduced_action_source_check as source
from S11c_d_reduced_action_source_check import ROOT, STORE, engine, digest, save, atomic_pickle
from S11c_d_output_codec import decoded_lines, restore_emission_index
from ledger_fold import _restore

CHECKPOINT = ROOT/'_measurements/S11c_d_reduced_action_source_checkpoint.json'
PREFIX = 'REDUCED_ACTION_ASSEMBLY_LAB_HELD_RHO4_CONSTANT'
SOURCES = (*engine.BUILD_INPUT_PATHS, Path(__file__).resolve(), Path(source.__file__).resolve(),
    CHECKPOINT, source.INPUT, ROOT/'_measurements/S11c_d_variable_profile_matching_plan.md')


def load():
    checkpoint = json.loads(CHECKPOINT.read_text())
    base = Path(checkpoint['runDirectory'])
    for name, item in checkpoint['artifacts'].items():
        if digest(base/name) != item['sha256']:
            raise ValueError(('saved action artifact changed', name))
    for name, item in checkpoint['publications'].items():
        if digest(ROOT/item['path']) != item['sha256']:
            raise ValueError(('source annex payload changed', name))
    for name, sha in checkpoint['sourceFiles'].items():
        if digest(base/'source'/name) != sha:
            raise ValueError(('source snapshot changed', name))
        if name != str(engine.HERE.relative_to(ROOT)) and digest(ROOT/name) != sha:
            raise ValueError(('consumed source changed', name))
    frozen = ast.parse((base/'source'/engine.HERE.relative_to(ROOT)).read_text())
    current = ast.parse(engine.HERE.read_text())
    current.body = [n for n in current.body if getattr(n, 'name', None) != 'ReducedActionAssembly']
    if ast.dump(current) != ast.dump(frozen):
        raise ValueError('native engine changed beyond the new assembly class')
    with (base/'reduced-action.pickle').open('rb') as stream:
        packet = pickle.load(stream)
    with (base/'actions.pickle').open('rb') as stream:
        actions = pickle.load(stream)
    if actions['reducedActionSha256'] != digest(base/'reduced-action.pickle'):
        raise ValueError('action/reduction packet join')
    reduction, dimensions = source.restore_context(packet)
    dimensions.__dict__.update(actions['dimensionState'])
    values = {key: engine.named(body, 'VALUE') for (key, _), body in packet['payloads'].items()}
    pencil = engine.ReducedPencil(*(values[key] for key in engine.CLOSED_KEYS), reduction)
    if pencil.strong != actions['strong'] or pencil.kernel != actions['kernel'] or pencil.fields != actions['fields']:
        raise ValueError('complete reduced source/pencil join')
    return checkpoint, packet, actions, pencil


def emit_result(result, actions, missing, provenance):
    dimensions = engine.PHYSICAL_METADATA.dimensions
    row_units = [actions['columnUnits'][(0, i)] for i in range(5)]
    field_units = [dimensions.known[f] for f in actions['fields']]
    z_unit = dimensions.measure(engine.PHYSICAL_METADATA.reduction.z)
    for order, matrix in result['LOCAL_MATRICES'].items():
        units = {(5*i+j,): tuple(a-b+order*c for a, b, c in zip(row_units[i], field_units[j], z_unit))
                 for i in range(5) for j in range(5)}
        engine.fingerprinted(PREFIX+'_LOCAL_MATRIX_'+str(order), matrix, units)
    for record in result['ROWS']:
        suffix = '_'+str(record['COLUMN'])+'_'+str(record['ROW'])
        row_unit = row_units[record['ROW']]
        operands = sp.Tuple(*(sp.Tuple(g, c) for g, c in zip(record['GENERATORS'], record['COEFFICIENTS'])))
        generator_units = [dimensions.measure(g) for g in record['GENERATORS']]
        coefficient_units = [tuple(a-b for a, b in zip(row_unit, u)) for u in generator_units]
        operand_units = {(i, j): unit for i, pair in enumerate(zip(generator_units, coefficient_units))
                         for j, unit in enumerate(pair)}
        engine.fingerprinted(PREFIX+'_COEFFICIENT_OPERANDS'+suffix, operands, operand_units)
        residuals = sp.Tuple(record['RECONSTRUCTION_RESIDUAL'], record['NONLINEAR_OR_AFFINE_REMAINDER'])
        engine.physical(PREFIX+'_ACTION_RECONSTRUCTION_RESIDUALS'+suffix, residuals,
                        zero_dimensions={(i,): row_unit for i in range(2)})
        engine.physical(PREFIX+'_DERIVATIVE_EXTRACTION_RESIDUALS'+suffix,
                        record['DERIVATIVE_EXTRACTION_RESIDUALS'],
                        zero_dimensions={(i,): u for i, u in enumerate(coefficient_units)})
    engine.fingerprinted(PREFIX+'_ORDERED_NONLOCAL_INTEGRALS', engine.cas(result['NONLOCAL_INTEGRALS']))
    engine.physical(PREFIX+'_UNBOUND_CONSTITUTIVE_PARAMETERS', engine.cas(missing))
    engine.physical(PREFIX+'_SOURCE_PROVENANCE', provenance)


def resume_packet(origin, destination, provenance):
    """Reuse saved construction and transcript after an instrument-only fix."""
    origin = origin.resolve(); origin.relative_to(STORE)
    before = json.loads((origin/'preflight.json').read_text())
    instrument = str(Path(__file__).resolve().relative_to(ROOT))
    for name, sha in before['sourceFiles'].items():
        if digest(origin/'source'/name) != sha:
            raise ValueError(('original assembly snapshot changed', name))
        if name != instrument and digest(ROOT/name) != sha:
            raise ValueError(('assembly input/constructor changed', name))
    frozen = ast.parse((origin/'source'/instrument).read_text())
    current = ast.parse(Path(__file__).read_text())
    joins = {}
    for name in ('load', 'emit_result'):
        a = next(n for n in frozen.body if getattr(n, 'name', None) == name)
        b = next(n for n in current.body if getattr(n, 'name', None) == name)
        joins[name] = ast.dump(a) == ast.dump(b)
    if not all(joins.values()) or before['provenance'] != provenance:
        raise ValueError('assembly source/emitter resume join')
    artifacts = {}
    for name in ('assembly.pickle', 'full.out'):
        path = origin/name
        artifacts[name] = {'bytes': path.stat().st_size, 'sha256': digest(path)}
        shutil.copyfile(path, destination/name)
        if digest(destination/name) != artifacts[name]['sha256']:
            raise ValueError('assembly resume copy mismatch')
    with (destination/'assembly.pickle').open('rb') as stream:
        saved = pickle.load(stream)
    if saved['provenance'] != provenance:
        raise ValueError('assembly packet provenance changed')
    engine.PHYSICAL_METADATA.dimensions.__dict__.update(saved['dimensionState'])
    record = {'runDirectory': str(origin), 'sourceFiles': before['sourceFiles'],
              'artifacts': artifacts, 'consumedHelperAstJoins': joins}
    save(destination/'resumed-inputs.json', record)
    return saved['result'], record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    parser.add_argument('--resume-from', type=Path,
                        help='validate an existing saved packet/transcript without construction')
    args = parser.parse_args()
    base = args.run_directory.resolve(); base.relative_to(STORE)
    base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    pins = {str(p.relative_to(ROOT)): digest(p) for p in SOURCES}
    for p in SOURCES:
        target = base/'source'/p.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, target)
    checkpoint, packet, actions, pencil = load()
    provenance = {'SOURCE_CHECKPOINT_SHA256': digest(CHECKPOINT),
        'REDUCED_ACTION_SHA256': actions['reducedActionSha256'],
        'ACTIONS_SHA256': checkpoint['artifacts']['actions.pickle']['sha256'],
        'INPUT_SHA256': digest(source.INPUT)}
    save(base/'preflight.json', {'sourceFiles': pins, 'provenance': provenance})
    resumed = None
    if args.resume_from:
        result, resumed = resume_packet(args.resume_from, base, provenance)
    else:
        assembly = engine.ReducedActionAssembly(pencil, actions['columns'])
        result = assembly.construct()
        atomic_pickle(base/'assembly.pickle', {'result': result,
            'dimensionState': dict(vars(engine.PHYSICAL_METADATA.dimensions)), 'provenance': provenance})
    free = set(engine.dag_free_symbols(actions['columns']))
    parameters = packet['inputSpecification']['parameters']
    live = {pencil.r.z, pencil.r.regulator, pencil.r.symbols['sigma_W']}
    missing = tuple(sorted((s for s in free-live if s.name not in parameters), key=sp.default_sort_key))
    census = [{'COLUMN': r['COLUMN'], 'ROW': r['ROW'], 'LOCAL_TERMS': len(r['LOCAL']),
               'NONLOCAL_TERMS': len(r['NONLOCAL'])} for r in result['ROWS']]
    zero = engine.PHYSICAL_METADATA.dimensions.zero
    census_units = {p: zero for p, _ in engine.leaves(engine.cas(census))}
    if not args.resume_from:
        with (base/'full.out').open('x') as transcript, contextlib.redirect_stdout(transcript):
            emit_result(result, actions, missing, provenance)
            engine.physical(PREFIX+'_TERM_CENSUS', census, zero_dimensions=census_units)
            keys = {tag: 's11cd'+''.join(w.title() for w in tag.removeprefix('PY_S11CD_').split('_'))
                    for tag in engine.EMISSION_LINES if not tag.startswith('PY_S11CD_METADATA_')}
            engine.physical(PREFIX+'_WRITE_KEYS', keys)
            index = engine.emission_index(engine.EMISSION_LINES)
            index_units = {p: zero for p, _ in engine.leaves(engine.cas(index))}
            engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=index_units)
    entries = {}
    for line in decoded_lines(base/'full.out'):
        tag, _, body = line.rstrip('\n').partition(': ')
        if tag in entries:
            raise ValueError('duplicate assembly emission')
        entries[tag] = _restore(body)
    if args.resume_from:
        keys = {str(k): str(v) for k, v in entries['PY_S11CD_'+PREFIX+'_WRITE_KEYS']}
        index = {str(k): v for k, v in entries['PY_S11CD_'+PREFIX+'_EMISSION_LINES']}
        index_units = {p: zero for p, _ in engine.leaves(engine.cas(index))}
    seen = set(); original = engine.emit
    def compare(name, value):
        tag = 'PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value):
            raise ValueError(('computed/emitted assembly mismatch', tag))
        seen.add(tag)
    engine.emit = compare
    try:
        emit_result(result, actions, missing, provenance)
        engine.physical(PREFIX+'_TERM_CENSUS', census, zero_dimensions=census_units)
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=index_units)
    finally:
        engine.emit = original
    if seen != set(entries) or len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('assembly emission/write-key census')
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    restore_emission_index({str(k): v for k, v in entries[final]}, list(entries)[:list(entries).index(final)])
    metadata_paths = 0
    for tag, body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):
            continue
        for path, fields in body:
            fields = {str(k): v for k, v in fields}
            unit = fields['DIMENSION_L_T_M']
            if not isinstance(unit, sp.Tuple) or len(unit) != 3 or any(s.free_symbols for s in unit):
                raise ValueError(('unresolved assembly dimension', tag, path))
            if any(k not in fields for k in ('MULTIGRADE', 'EPSILON_LAMBDA_SUPPORT')):
                raise ValueError('missing assembly grades')
            metadata_paths += 1
    residuals = [v for r in result['ROWS'] for v in (r['RECONSTRUCTION_RESIDUAL'],
        r['NONLINEAR_OR_AFFINE_REMAINDER'], *r['DERIVATIVE_EXTRACTION_RESIDUALS'])]
    summary = {'runDirectory': str(base), 'sourceFiles': pins, 'provenance': provenance,
        'localDerivativeOrders': [int(n) for n in result['LOCAL_MATRICES']], 'termCensus': census,
        'distinctNonlocalIntegrals': len(result['NONLOCAL_INTEGRALS']),
        'residualScalars': len(residuals), 'nonzeroResidualScalars': sum(v != 0 for v in residuals),
        'unboundConstitutiveParameters': [str(s) for s in missing],
        'tagCount': len(entries), 'metadataPaths': metadata_paths, 'writeKeyCount': len(keys),
        'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts': {n: {'bytes': (base/n).stat().st_size, 'sha256': digest(base/n)} for n in ('full.out', 'assembly.pickle')},
        'scope': 'Exact local-jet and intact-integral assembly with formal reconstruction and derivative-extraction residuals. No numerical quadrature or scattering solution.'}
    if resumed is not None:
        summary['resumedFrom'] = resumed
    save(base/'checks.json', summary)
    if pins != {str(p.relative_to(ROOT)): digest(p) for p in SOURCES}:
        raise ValueError('assembly sources changed during execution')
    if any(v != 0 for v in residuals) or engine.PHYSICAL_METADATA.dimensions.constraints:
        raise ValueError('assembly residual requires inspection; all operands preserved')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
