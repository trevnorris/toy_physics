#!/usr/bin/env python3
"""Construct finite-domain source Fourier factors from accepted native actions."""
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

import S11c_d_numerical_action_check as native
from S11c_d_numerical_action_check import ROOT, STORE, engine, digest, save, atomic_pickle
from S11c_d_output_codec import decoded_lines, restore_emission_index
from ledger_fold import _restore

CHECKPOINT = ROOT/'_measurements/S11c_d_quadrature_domain_checkpoint.json'
PLAN = ROOT/'_measurements/S11c_d_source_fourier_factorization_plan.md'
PREFIX = 'BOUNDED_SOURCE_FOURIER_LAB_HELD_RHO4_CONSTANT'
SOURCES = (*native.SOURCES, Path(native.__file__).resolve(), Path(__file__).resolve(), CHECKPOINT, PLAN)


def load():
    accepted = json.loads(CHECKPOINT.read_text())
    base = Path(accepted['runDirectory'])
    for name, item in accepted['artifacts'].items():
        if digest(base/name) != item['sha256']:
            raise ValueError(('quadrature-domain operand changed', name))
    if digest(ROOT/accepted['publication']['path']) != accepted['publication']['sha256']:
        raise ValueError('quadrature-domain publication changed')
    engine_name = str(engine.HERE.relative_to(ROOT))
    for name, sha in accepted['sourceFiles'].items():
        if digest(base/'source'/name) != sha:
            raise ValueError(('quadrature-domain source snapshot changed', name))
        if name != engine_name and digest(ROOT/name) != sha:
            raise ValueError(('consumed quadrature source changed', name))
    frozen = ast.parse((base/'source'/engine_name).read_text())
    current = ast.parse(engine.HERE.read_text())
    current.body = [n for n in current.body if getattr(n, 'name', None) != 'BoundedSourceFourierAssembly']
    if ast.dump(current) != ast.dump(frozen):
        raise ValueError('native engine changed beyond the new bounded source Fourier class')
    assembly_checkpoint = json.loads(native.ASSEMBLY.read_text())
    assembly_base = Path(assembly_checkpoint['runDirectory'])
    for name, item in assembly_checkpoint['artifacts'].items():
        if digest(assembly_base/name) != item['sha256']:
            raise ValueError(('assembly operand changed', name))
    source_path = ROOT/'_measurements/S11c_d_reduced_action_source_checkpoint.json'
    source_checkpoint = json.loads(source_path.read_text())
    source_base = Path(source_checkpoint['runDirectory'])
    for name, item in source_checkpoint['artifacts'].items():
        if digest(source_base/name) != item['sha256']:
            raise ValueError(('source action operand changed', name))
    with (source_base/'reduced-action.pickle').open('rb') as stream:
        packet = pickle.load(stream)
    with (source_base/'actions.pickle').open('rb') as stream:
        actions = pickle.load(stream)
    with (assembly_base/'assembly.pickle').open('rb') as stream:
        assembly = pickle.load(stream)
    if (actions['reducedActionSha256'] != digest(source_base/'reduced-action.pickle') or
            assembly['provenance'] != assembly_checkpoint['provenance'] or
            assembly['provenance']['SOURCE_CHECKPOINT_SHA256'] != digest(source_path) or
            assembly['provenance']['ACTIONS_SHA256'] != digest(source_base/'actions.pickle')):
        raise ValueError('source action/assembly packet join')
    reduction, dimensions = native.source.restore_context(packet)
    dimensions.__dict__.update(assembly['dimensionState'])
    values = {key: engine.named(body, 'VALUE') for (key, _), body in packet['payloads'].items()}
    pencil = engine.ReducedPencil(*(values[key] for key in engine.CLOSED_KEYS), reduction)
    if pencil.strong != actions['strong'] or pencil.kernel != actions['kernel'] or pencil.fields != actions['fields']:
        raise ValueError('full reduced source row/field/kernel join')
    projection = native.input_check.check()
    return pencil, assembly['result'], {'QUADRATURE_CHECKPOINT_SHA256': digest(CHECKPOINT),
        'ASSEMBLY_CHECKPOINT_SHA256': digest(native.ASSEMBLY),
        'SOURCE_ACTION_PACKET_SHA256': digest(source_base/'actions.pickle'),
        'APPROVED_INPUT_SHA256': projection['files']['input']['sha256']}


def emit_result(result, pencil, provenance):
    dimensions = engine.PHYSICAL_METADATA.dimensions
    zero = dimensions.zero
    zp_unit = dimensions.measure(pencil.r.zp)
    engine.physical(PREFIX+'_PROVENANCE', provenance)
    engine.physical(PREFIX+'_FINITE_CUTOFF_VARIABLES', result['CUTOFFS'])
    for row in result['ROWS']:
        suffix = str(row['INDEX'])
        unit = dimensions.measure(row['ORIGINAL'])
        integrand_unit = dimensions.measure(row['BOUNDED'].function)
        engine.fingerprinted(PREFIX+'_ORIGINAL_AND_BOUNDED_'+suffix,
                             sp.Tuple(row['ORIGINAL'], row['BOUNDED']), {(0,): unit, (1,): unit})
        engine.fingerprinted(PREFIX+'_SOURCE_FIRST_BOUNDED_'+suffix, row['SOURCE_FIRST_BOUNDED'], {(): unit})
        engine.physical(PREFIX+'_INTEGRAND_RECONSTRUCTION_RESIDUAL_'+suffix,
                        row['RECONSTRUCTION_RESIDUAL'], zero_dimensions={(): integrand_unit})
        operands, residuals, operand_units, residual_units = [], [], {}, {}
        for index, factor in enumerate(row['FACTORS']):
            source_unit = dimensions.measure(factor['SOURCE'])
            coefficient_unit = tuple(a-b for a, b in zip(integrand_unit, source_unit))
            units = (source_unit, coefficient_unit, zero, tuple(-v for v in zp_unit), source_unit,
                     tuple(a+b for a, b in zip(source_unit, zp_unit)))
            operands.append(sp.Tuple(*(factor[k] for k in
                ('SOURCE', 'COEFFICIENT', 'CHARACTER', 'FREQUENCY', 'AMPLITUDE', 'SOURCE_INTEGRAL'))))
            operand_units.update({(index, j): u for j, u in enumerate(units)})
            residuals.append(sp.Tuple(*(factor[k] for k in ('AMPLITUDE_RECONSTRUCTION_RESIDUAL',
                'CHARACTER_NORMALIZATION_RESIDUAL', 'CHARACTER_EQUATION_RESIDUAL'))))
            residual_units.update({(index, 0): source_unit, (index, 1): zero,
                                   (index, 2): tuple(-v for v in zp_unit)})
        engine.fingerprinted(PREFIX+'_SOURCE_FACTOR_OPERANDS_'+suffix, sp.Tuple(*operands), operand_units)
        engine.physical(PREFIX+'_SOURCE_CHARACTER_RESIDUALS_'+suffix, sp.Tuple(*residuals),
                        zero_dimensions=residual_units)
    phase_operands, phase_units = [], {}
    for i, (original, record) in enumerate(result['PHASES']):
        phase_operands.append(sp.Tuple(original, *(record[k] for k in ('SOURCE_EXPONENT', 'OTHER_EXPONENT',
                                      'EXPONENT_RESIDUAL', 'SECOND_SOURCE_DERIVATIVE'))))
        phase_units.update({(i, j): zero if j != 4 else tuple(-2*v for v in zp_unit) for j in range(5)})
    engine.fingerprinted(PREFIX+'_PHASE_SPLIT_OPERANDS', sp.Tuple(*phase_operands), phase_units)
    engine.fingerprinted(PREFIX+'_BOUNDED_SOURCE_INTEGRALS', sp.Tuple(*result['SOURCE_INTEGRALS']))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    args = parser.parse_args()
    base = args.run_directory.resolve(); base.relative_to(STORE); base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    paths = tuple(dict.fromkeys(SOURCES))
    pins = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    for path in paths:
        target = base/'source'/path.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    pencil, assembly, provenance = load()
    save(base/'preflight.json', {'sourceFiles': pins, 'provenance': provenance})
    builder = engine.BoundedSourceFourierAssembly(pencil.r, assembly['NONLOCAL_INTEGRALS'])
    result = builder.construct()
    atomic_pickle(base/'factorization.pickle', {'result': result, 'provenance': provenance,
        'dimensionState': dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    before = digest(base/'factorization.pickle')
    with (base/'full.out').open('x') as stream, contextlib.redirect_stdout(stream):
        emit_result(result, pencil, provenance)
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
            raise ValueError('duplicate bounded source Fourier emission')
        entries[tag] = _restore(body)
    seen = set(); original = engine.emit
    def compare(name, value):
        tag = 'PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value):
            raise ValueError(('bounded source Fourier emission mismatch', tag))
        seen.add(tag)
    engine.emit = compare
    try:
        emit_result(result, pencil, provenance)
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=zero_units)
    finally:
        engine.emit = original
    if seen != set(entries) or len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('bounded source Fourier emission/write-key census')
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    restore_emission_index({str(k): v for k, v in entries[final]}, list(entries)[:list(entries).index(final)])
    metadata_paths = 0
    for tag, body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'):
            continue
        for path, descriptor in body:
            fields = {str(k): v for k, v in descriptor}
            unit = fields['DIMENSION_L_T_M']
            if not isinstance(unit, sp.Tuple) or len(unit) != 3 or any(s.free_symbols for s in unit):
                raise ValueError(('unresolved bounded source Fourier unit', tag, path))
            if any(k not in fields for k in ('MULTIGRADE', 'EPSILON_LAMBDA_SUPPORT')):
                raise ValueError(('missing bounded source Fourier grades', tag))
            metadata_paths += 1
    residuals = [r['RECONSTRUCTION_RESIDUAL'] for r in result['ROWS']]
    residuals.extend(f[k] for row in result['ROWS'] for f in row['FACTORS'] for k in
                     ('AMPLITUDE_RECONSTRUCTION_RESIDUAL', 'CHARACTER_NORMALIZATION_RESIDUAL', 'CHARACTER_EQUATION_RESIDUAL'))
    residuals.extend(p[k] for _, p in result['PHASES'] for k in ('EXPONENT_RESIDUAL', 'SECOND_SOURCE_DERIVATIVE'))
    summary = {'runDirectory': str(base), 'sourceFiles': pins, 'provenance': provenance,
        'originalIntegralCount': len(assembly['NONLOCAL_INTEGRALS']), 'factorizedIntegralCount': len(result['ROWS']),
        'factorCounts': [len(r['FACTORS']) for r in result['ROWS']],
        'distinctBoundedSourceIntegrals': len(result['SOURCE_INTEGRALS']), 'phaseCount': len(result['PHASES']),
        'cutoffs': {str(v): str(c) for v, c in result['CUTOFFS'].items()},
        'residualScalars': len(residuals), 'nonzeroResidualScalars': sum(r != 0 for r in residuals),
        'tagCount': len(entries), 'writeKeyCount': len(keys), 'metadataPaths': metadata_paths,
        'packetSha256BeforeEmission': before, 'packetSha256AfterEmission': digest(base/'factorization.pickle'),
        'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts': {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)} for p in sorted(base.iterdir())
                      if p.suffix in ('.pickle', '.out')},
        'scope': 'Source-position and phase factorization of the full native nonlocal integrands on explicit finite domains. '
                 'Source-first iterated integrals require the recorded smooth-test/positive-regulator/nonsingular-domain conditions. '
                 'No equality of unbounded distributional iterated integrals, numerical action limit or scattering solution is asserted.'}
    save(base/'checks.json', summary)
    if pins != {str(p.relative_to(ROOT)): digest(p) for p in paths} or before != digest(base/'factorization.pickle'):
        raise ValueError('bounded source Fourier source/packet changed')
    if any(r != 0 for r in residuals) or engine.PHYSICAL_METADATA.dimensions.constraints:
        raise ValueError('bounded source Fourier residual or dimension requires inspection')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
