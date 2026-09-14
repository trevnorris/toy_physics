#!/usr/bin/env python3
"""Compute fresh reduced reference and end symbols before spectral regeneration.

This scoped producer executes the engine's reduction and constant-end
constructors. It computes no spectrum and publishes no complete d transcript.
Its source-pinned cache is consumed by the existing current preflight loader.
"""
import argparse
import contextlib
import hashlib
import json
from pathlib import Path
import pickle
import resource
import shutil
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
RUN = None
sys.path.insert(0, str(ROOT / 'scripts'))
import sympy as sp
from ledger_fold import load_model, check_consumer, assert_lookups_equal_manifest
import S11c_d_mixing_scattering_sympy_audit as engine


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def construct():
    fold, audit = load_model(*(str(ROOT / 'scripts' / ('S11c_'+stage+'_exports.py')) for stage in ('b', 'c1', 'c2')))
    closure = check_consumer(fold, engine.IMPORT_KEYS)
    lookups = assert_lookups_equal_manifest(engine.bind, fold, engine.IMPORT_KEYS)
    rows = lookups['result']
    engine.emit('TRACE_REPAIR_END_SYMBOL_IMPORT_FOLD', audit)
    engine.emit('TRACE_REPAIR_END_SYMBOL_IMPORT_LOOKUPS', sorted(lookups['lookups']))
    engine.emit('TRACE_REPAIR_END_SYMBOL_IMPORT_CLOSURE', {key: value for key, value in closure.items() if key != 'resolved_imports'})
    reduction = engine.EdgeReduction(rows)
    dimensions = engine.DimensionAnalysis(rows, reduction)
    engine.PHYSICAL_METADATA = engine.PhysicalMetadata(dimensions, reduction)
    carriers = engine.FourierCarrierReconstruction(reduction)
    reconstruction = engine.EdgeReconstruction(reduction, carriers)
    case = ('LAB_HELD', 'RHO4_CONSTANT')
    reduced = {}
    for root in engine.CLOSED_KEYS:
        payload = next(payload for axes, payload in rows[root]['value'] if tuple(map(str, axes)) == case)
        suffix = '_'.join((root, *case))
        value = engine.named(payload, 'VALUE')
        units = engine.payload_units(value, engine.named(payload, 'DIMENSION_L_T_M'))
        definitions = tuple(engine.named(payload, 'FOURIER_PROFILE_BINDINGS'))
        result, records = reduction.value(value)
        engine.fingerprinted('TRACE_REPAIR_END_SYMBOL_SOURCE_'+suffix, value, units)
        engine.fingerprinted('TRACE_REPAIR_END_SYMBOL_REDUCED_'+suffix, result, units)
        reconstruction.row(value, result, records, suffix, units, definitions)
        reduced_payload = reduction.payload(payload, result)
        for index, equation in enumerate(engine.named(reduced_payload, 'COMPUTED_BRANCH_BINDINGS')):
            engine.physical('REDUCED_BINDING_OPERANDS_'+suffix+'_COMPUTED_BRANCH_BINDINGS_'+str(index),
                            (equation.lhs, equation.rhs), zero_dimensions={(1,): dimensions.known[equation.lhs.func]})
        reduced[root] = result
    pencil = engine.ReducedPencil(reduced[engine.CLOSED_KEYS[0]], reduced[engine.CLOSED_KEYS[1]], reduction)
    ends = engine.ConstantEndPencil(reduction)
    blocks = pencil.sector_blocks()
    curl, gauge, gauge_residual = ends.curl_gauge_operands(pencil)
    weak_units = ends.weak_matrix_units(blocks)
    strong_units = ends.strong_matrix_units(pencil)
    for label,endpoint in (('REFERENCE',None),('LEFT',-sp.oo),('RIGHT',sp.oo)):
        constant_strong = ends.background(pencil.strong, endpoint)
        constant_blocks = ends.background(blocks, endpoint)
        strong = ends.strong_matrix(pencil, constant_strong)
        full = ends.weak_matrix(constant_blocks)
        current = engine.UniformSlabCurrent(reduction, rows['energy_basis_variable'], ends, strong[3, :])
        current.emit(case[0], endpoint, label+'_'+'_'.join(case))
        engine.fingerprinted('TRACE_REPAIR_END_SYMBOL_STRONG_'+label, strong, strong_units)
        engine.fingerprinted('TRACE_REPAIR_END_SYMBOL_WEAK_'+label, full, weak_units)
        cache = RUN / 'symbols' / (label+'_LAB_HELD_RHO4_CONSTANT.pickle')
        cache.parent.mkdir(exist_ok=True)
        cache.write_bytes(pickle.dumps((full,curl,weak_units,dict(dimensions.known),strong,
            rows['energy_basis_variable']['value']),protocol=5))
    engine.physical('TRACE_REPAIR_END_SYMBOL_DIMENSION_CONSTRAINTS', tuple(dimensions.constraints))
    unresolved = tuple((str(atom), units) for atom, units in dimensions.unknown.items()
                       if atom not in dimensions.known and any(unit.free_symbols for unit in units))
    engine.physical('TRACE_REPAIR_END_SYMBOL_DIMENSION_UNRESOLVED', unresolved)
    if dimensions.constraints or unresolved:
        raise ValueError('end-symbol dimensional constraints; see emitted operands')


def run():
    global RUN
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    args=parser.parse_args();RUN=args.run_directory
    RUN.mkdir(parents=True, exist_ok=False)
    paths = [Path(__file__).relative_to(ROOT), Path('scripts/S11c_d_mixing_scattering_sympy_audit.py'),
             Path('scripts/S11c_d_output_codec.py'), Path('scripts/ledger_fold.py'),
             Path('directives/S11c_d_SHARED_PHYSICS.md'),
             Path('directives/S11b_SHARED_PHYSICS.md'),
             Path('directives/S11c_d_sympy_build_directive.md'),
             Path('directives/S11c_d_sympy_build_PROGRAM_BRIEF.md'),
             *(Path('scripts') / ('S11c_'+stage+'_exports.py') for stage in ('b', 'c1', 'c2'))]
    pins = {str(path): digest(ROOT / path) for path in paths}
    for path in paths:
        target = RUN / 'source' / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / path, target)
    started = time.monotonic()
    code = 0
    with (RUN / 'full.out').open('w') as out, (RUN / 'stderr.txt').open('w') as err:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                construct()
            except Exception:
                traceback.print_exc()
                code = 1
    manifest = {'run_directory': str(RUN), 'command': [sys.executable, str(Path(__file__))],
                'scope': 'LAB_HELD/RHO4_CONSTANT reduction and reference/LEFT/RIGHT symbols; no spectrum',
                'exit_code': code, 'source_hashes_before': pins,
                'source_hashes_after': {str(path): digest(ROOT / path) for path in paths},
                'wall_seconds': time.monotonic()-started,
                'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    manifest['artifacts'] = {str(path.relative_to(RUN)): {'bytes': path.stat().st_size, 'sha256': digest(path)}
                             for path in RUN.rglob('*') if path.is_file() and 'source' not in path.relative_to(RUN).parts}
    (RUN / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({key: value for key, value in manifest.items() if key not in ('source_hashes_before', 'source_hashes_after', 'artifacts')}))
    raise SystemExit(code or (0 if pins == manifest['source_hashes_after'] else 1))


if __name__ == '__main__':
    run()
