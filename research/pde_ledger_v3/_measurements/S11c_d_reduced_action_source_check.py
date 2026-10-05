#!/usr/bin/env python3
"""Materialize native reduced rows and probe actions for numerical assembly."""
import argparse
import ast
import contextlib
import hashlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import subprocess
import sys
import time

import sympy as sp
from sympy.core.function import AppliedUndef

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT.parents[1]/'_scratch/s11c'
sys.path.insert(0, str(ROOT/'scripts'))
import S11c_d_mixing_scattering_sympy_audit as engine
from S11c_d_output_codec import decoded_lines, restore_emission_index
from ledger_fold import _restore

CASE = ('LAB_HELD', 'RHO4_CONSTANT')
PREFIX = 'REDUCED_ACTION_SOURCE_'+'_'.join(CASE)
INPUT = ROOT/'_measurements/S11c_d_channel_preflight_input.json'
MATCHING = ROOT/'_measurements/S11c_d_matching_channels_checkpoint.json'
PRODUCER = STORE/'s11c-thickness-coordinate-20260914/d_full/manifest.json'
SOURCES = (*engine.BUILD_INPUT_PATHS, Path(__file__).resolve(), INPUT, MATCHING,
    ROOT/'_measurements/S11c_d_variable_profile_matching_plan.md',
    ROOT/'directives/S11c_d_sympy_build_PROGRAM_BRIEF.md',
    ROOT/'directives/S11c_d_sympy_build_directive.md')


def digest(path):
    with path.open('rb') as stream:
        h = hashlib.sha256()
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def atomic_pickle(path, value):
    temporary = path.with_name(path.name+'.new')
    if path.exists():
        raise ValueError('saved operand already exists')
    with temporary.open('xb') as stream:
        pickle.dump(value, stream, protocol=5)
        stream.flush()
        import os
        os.fsync(stream.fileno())
    temporary.replace(path)


def source_joins():
    matching = json.loads(MATCHING.read_text())
    base = Path(matching['runDirectory'])
    for name, record in matching['artifacts'].items():
        if digest(base/name) != record['sha256']:
            raise ValueError(('matching artifact changed', name))
    publication = matching['publication']
    if digest(ROOT/publication['path']) != publication['sha256']:
        raise ValueError('matching annex payload changed')
    for name, expected in matching['sourceFiles'].items():
        if digest(base/'source'/name) != expected:
            raise ValueError(('matching frozen source changed', name))
    frozen = ast.parse((base/'source'/engine.HERE.relative_to(ROOT)).read_text())
    current = ast.parse(engine.HERE.read_text())
    # Only persistence and its run-time option were added after the accepted
    # matching engine. Every native constructor and helper must still agree.
    frozen.body = [n for n in frozen.body if getattr(n, 'name', None) != 'run']
    current.body = [n for n in current.body if getattr(n, 'name', None) not in
                    ('run', 'save_reduced_action_cache')]
    if ast.dump(frozen) != ast.dump(current):
        raise ValueError('native engine changed beyond reduction-cache instrumentation')
    producer = json.loads(PRODUCER.read_text())
    joins = {}
    for name in ('scripts/S11c_b_exports.py', 'scripts/S11c_c1_exports.py',
                 'scripts/S11c_c2_exports.py', 'scripts/ledger_fold.py',
                 'directives/S11c_d_SHARED_PHYSICS.md'):
        expected = producer['source_hashes_after'][name]
        joins[name] = digest(ROOT/name) == expected
    if not all(joins.values()):
        raise ValueError(('end producer physical inputs changed', joins))
    for end in ('LEFT', 'RIGHT'):
        source = matching['sources'][end]
        if source['inputSha256'] != digest(INPUT) or source['producerManifestSha256'] != digest(PRODUCER):
            raise ValueError(('matching profile/producer identity', end))
    return {'matchingCheckpointSha256': digest(MATCHING), 'inputSha256': digest(INPUT),
        'producerManifestSha256': digest(PRODUCER), 'nativeConstructorAstJoin': True,
        'physicalInputJoins': joins}


def restore_context(packet):
    # These dictionaries are the actual executed constructors' saved state;
    # restoring them avoids recalculating Fourier integrals or guessing symbols.
    reduction = engine.EdgeReduction.__new__(engine.EdgeReduction)
    reduction.__dict__.update(packet['reductionState'])
    dimensions = engine.DimensionAnalysis.__new__(engine.DimensionAnalysis)
    dimensions.__dict__.update(packet['dimensionState'])
    engine.PHYSICAL_METADATA = engine.PhysicalMetadata(dimensions, reduction)
    return reduction, dimensions


def replay_rows(path, packet):
    original = engine.emit
    emitted = {}
    def collect(name, value):
        emitted['PY_S11CD_'+name] = sp.srepr(engine.cas(value))
    engine.emit = collect
    try:
        for (key, case), payload in packet['payloads'].items():
            suffix = '_'.join((key, *map(str, case)))
            source = dict(packet['rows'][key]['value'])[case]
            reduced = engine.named(payload, 'VALUE')
            units = engine.payload_units(engine.named(source, 'VALUE'),
                                         engine.named(source, 'DIMENSION_L_T_M'))
            engine.physical('REDUCED_ACTION_ROWS_'+suffix, reduced, zero_dimensions=units)
            engine.physical('REDUCED_FIVE_SLOT_PAYLOAD_'+suffix, payload,
                            operands=reduced, zero_dimensions=units)
    finally:
        engine.emit = original
    tags, joined, indexed = [], set(), None
    for line in decoded_lines(path):
        tag, _, body = line.rstrip('\n').partition(': ')
        if tag in emitted:
            if emitted[tag] != body:
                raise ValueError(('saved reduction/emission mismatch', tag))
            joined.add(tag)
        if tag == 'PY_S11CD_EMISSION_LINES':
            indexed = restore_emission_index({str(k): v for k, v in _restore(body)}, tags)
        tags.append(tag)
    if joined != set(emitted) or len(tags) != len(set(tags)) or indexed is None:
        raise ValueError('reduction replay/tag-index census')
    if tags[-1] != 'PY_S11CD_PROCESS_COMPLETION':
        raise ValueError('incomplete reduction transcript')
    return {'tagCount': len(tags), 'replayedObjectAndMetadataTags': len(joined),
            'indexedTags': len(indexed)}


def action_census(columns, pencil):
    result = []
    for j, column in enumerate(columns):
        for i, expression in enumerate(column):
            functions = expression.atoms(AppliedUndef)
            result.append({'COLUMN': j, 'ROW': i, 'DAG_NODES': engine.dag_size(expression),
                'TOP_LEVEL_ADDENDS': len(sp.Add.make_args(expression)),
                'DISTINCT_INTEGRALS': len(expression.atoms(sp.Integral)),
                'DISTINCT_DERIVATIVES': len(expression.atoms(sp.Derivative)),
                'PROBE_APPLICATIONS': sum(f.func == pencil.probes[j] for f in functions),
                'FOREIGN_PROBES': sum(f.func in set(pencil.probes)-{pencil.probes[j]} for f in functions),
                'UNSUBSTITUTED_FIELDS': sum(f.func in pencil.fields for f in functions)})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve(); base.relative_to(STORE)
    base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    pins = {str(p.relative_to(ROOT)): digest(p) for p in SOURCES}
    for source in SOURCES:
        target = base/'source'/source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    joins = source_joins()
    save(base/'preflight.json', {'sourceFiles': pins, 'joins': joins})
    cache_path = base/'reduced-action.pickle'
    command = [sys.executable, '-u', str(engine.HERE), '--case', '__'.join(CASE),
        '--dev-stop-after-reduction', '--dev-reduced-action-cache', str(cache_path),
        '--channel-input-file', str(INPUT)]
    with (base/'reduction.out').open('xb') as stdout, (base/'reduction.stderr').open('xb') as stderr:
        child = subprocess.Popen(command, cwd=ROOT, stdout=stdout, stderr=stderr)
        save(base/'reduction-invocation.json', {'command': command, 'pid': child.pid})
        code = child.wait()
    save(base/'reduction-invocation.json', {'command': command, 'pid': child.pid,
        'exitCode': code, 'stderrBytes': (base/'reduction.stderr').stat().st_size})
    if code or (base/'reduction.stderr').stat().st_size:
        raise ValueError('reduction stopped; inspect preserved reduction.stderr')
    with cache_path.open('rb') as stream:
        packet = pickle.load(stream)
    expected = {(key, CASE) for key in engine.CLOSED_KEYS}
    actual = {(key, tuple(map(str, case))) for key, case in packet['payloads']}
    if packet['schema'] != 1 or actual != expected:
        raise ValueError('complete selected-case reduced row pair required')
    if packet['inputSpecification'] != json.loads(INPUT.read_text()):
        raise ValueError('reduced-action profile instance changed')
    if packet['sourceDigests'] != {str(p.relative_to(ROOT)): digest(p) for p in engine.BUILD_INPUT_PATHS}:
        raise ValueError('reduced-action source pins changed')
    for key, payload in packet['payloads'].items():
        if [str(slot) for slot, _ in payload] != ['VALUE', 'MULTIGRADE', 'DIMENSION_L_T_M',
                                                'COMPUTED_BRANCH_BINDINGS', 'FOURIER_PROFILE_BINDINGS']:
            raise ValueError('reduced five-slot census')
        if packet['branchBindings'][key] != tuple((eq.lhs, eq.rhs) for eq in engine.named(payload, 'COMPUTED_BRANCH_BINDINGS')):
            raise ValueError('saved branch binding join')
    reduction, dimensions = restore_context(packet)
    replay = replay_rows(base/'reduction.out', packet)
    save(base/'row-replay.json', replay)
    values = {key: engine.named(payload, 'VALUE') for (key, _), payload in packet['payloads'].items()}
    payloads = {key: payload for (key, _), payload in packet['payloads'].items()}
    pencil = engine.ReducedPencil(*(values[key] for key in engine.CLOSED_KEYS), reduction)
    columns = pencil.columns()
    slab_units = engine.payload_units(pencil.slab,
        engine.named(payloads[engine.CLOSED_KEYS[0]], 'DIMENSION_L_T_M'))
    row_units = [slab_units[('U', i)] for i in range(3)]+[slab_units[(label,)] for label in ('THETA', 'E_W')]
    kernel_units = engine.payload_units(pencil.kernel,
        engine.named(payloads[engine.CLOSED_KEYS[1]], 'DIMENSION_L_T_M'))
    units = {(j, i): unit for j in range(len(columns)) for i, unit in enumerate(row_units)}
    census = action_census(columns, pencil)
    actions = {'columns': columns, 'strong': pencil.strong, 'kernel': pencil.kernel,
        'fields': pencil.fields, 'probes': pencil.probes, 'columnUnits': units,
        'census': census, 'dimensionState': dict(vars(dimensions)),
        'reducedActionSha256': digest(cache_path)}
    atomic_pickle(base/'actions.pickle', actions)
    join_operands = dict(joins, nativeConstructorAstJoin=int(joins['nativeConstructorAstJoin']),
                         physicalInputJoins={k: int(v) for k, v in joins['physicalInputJoins'].items()})
    objects = [('COLUMNS', columns, units, True),
               ('KERNEL', pencil.kernel, kernel_units, True),
               ('ACTION_CENSUS', engine.cas(census), {p: dimensions.zero for p, _ in engine.leaves(engine.cas(census))}, False),
               ('SOURCE_JOINS', engine.cas(join_operands),
                {p: dimensions.zero for p, _ in engine.leaves(engine.cas(join_operands))}, False)]
    def emit_objects():
        for name, value, unit, heavy in objects:
            if heavy:
                engine.fingerprinted(PREFIX+'_'+name, value, unit)
            else:
                engine.physical(PREFIX+'_'+name, value, zero_dimensions=unit)
    with (base/'full.out').open('x') as transcript, contextlib.redirect_stdout(transcript):
        emit_objects()
        names = tuple(engine.EMISSION_LINES)
        keys = {name: 's11cd'+''.join(word.title() for word in name.removeprefix('PY_S11CD_').split('_'))
                for name in names if not name.startswith('PY_S11CD_METADATA_')}
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        index = engine.emission_index(engine.EMISSION_LINES)
        engine.physical(PREFIX+'_EMISSION_LINES', index,
                        zero_dimensions={p: dimensions.zero for p, _ in engine.leaves(engine.cas(index))})
    entries = {}
    for line in decoded_lines(base/'full.out'):
        tag, _, payload = line.rstrip('\n').partition(': ')
        if tag in entries:
            raise ValueError('duplicate action source tag')
        entries[tag] = _restore(payload)
    original = engine.emit; seen = set()
    def compare(name, value):
        tag = 'PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value):
            raise ValueError(('probe action computed/emitted join', tag))
        seen.add(tag)
    engine.emit = compare
    try:
        emit_objects()
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        engine.physical(PREFIX+'_EMISSION_LINES', index,
                        zero_dimensions={p: dimensions.zero for p, _ in engine.leaves(engine.cas(index))})
    finally:
        engine.emit = original
    if seen != set(entries) or len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('action emission/write-key census')
    metadata_paths = 0
    for tag, value in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):
            continue
        for path, fields in value:
            fields = {str(k): v for k, v in fields}
            unit = fields['DIMENSION_L_T_M']
            if not isinstance(unit, sp.Tuple) or len(unit) != 3 or any(v.free_symbols for v in unit):
                raise ValueError(('unresolved action metadata', tag, path))
            if any(key not in fields for key in ('MULTIGRADE', 'EPSILON_LAMBDA_SUPPORT')):
                raise ValueError('missing action grade metadata')
            metadata_paths += 1
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    restore_emission_index({str(k): v for k, v in entries[final]}, list(entries)[:list(entries).index(final)])
    if any(row['FOREIGN_PROBES'] or row['UNSUBSTITUTED_FIELDS'] for row in census):
        raise ValueError('probe action field census; operands preserved')
    if dimensions.constraints:
        raise ValueError('action units introduced inconsistent constraints')
    if pins != {str(p.relative_to(ROOT)): digest(p) for p in SOURCES}:
        raise ValueError('sources changed during action-source construction')
    artifacts = {name: {'bytes': (base/name).stat().st_size, 'sha256': digest(base/name)}
        for name in ('reduction.out', 'reduction.stderr', 'reduced-action.pickle', 'actions.pickle', 'full.out')}
    record = {'runDirectory': str(base), 'sourceFiles': pins, 'sourceJoins': joins,
        'reductionReplay': replay, 'actionCensus': census, 'artifacts': artifacts,
        'emissionReplayTags': len(seen), 'writeKeyCount': len(keys), 'metadataPaths': metadata_paths,
        'wallSeconds': time.monotonic()-started,
        'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'reductionChildPeakRssKiB': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
        'scope': 'Computed reduced source rows and native probe actions. No local/nonlocal numerical quadrature, boundary solve or scattering result.'}
    save(base/'checks.json', record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
