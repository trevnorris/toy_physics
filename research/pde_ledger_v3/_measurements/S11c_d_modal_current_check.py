#!/usr/bin/env python3
"""Pinned one-case continuation of the native two-frequency current builder."""
import argparse
import ast
import faulthandler
import hashlib
import json
from pathlib import Path
import pickle
import resource
import time

import sympy as sp

from S11c_d_joint_sheet_check import load, engine


def digest(path):
    with path.open('rb') as stream:
        result = hashlib.sha256()
        for block in iter(lambda:stream.read(2**20), b''):
            result.update(block)
    return result.hexdigest()


def source_node(source, name):
    node = next(node for node in ast.parse(source).body if getattr(node, 'name', None) == name)
    return ast.dump(node, include_attributes=False)


def build(args):
    modes, strong, units, bindings, inputs, provenance = load(args)
    manifest = json.loads(args.manifest.read_text())
    cache = Path(manifest['run_directory'])/'symbols'/'REFERENCE_LAB_HELD_RHO4_CONSTANT.pickle'
    energy = pickle.loads(cache.read_bytes())[5]
    r = modes.r
    r.x = tuple(r.symbols['s11cc2X'+str(i)] for i in (1, 2, 3))
    r.t = r.symbols['s11cc2Time']
    r.z = r.symbols['s11cdNormalPosition']
    ends = engine.ConstantEndPencil.__new__(engine.ConstantEndPencil)
    ends.r, ends.kn = r, modes.k
    current = engine.UniformSlabCurrent(r, {'value':energy}, ends, strong[3, :])
    balance = engine.SlabEnergyBalance(current)
    acoustic = engine.ClosedAcousticEnergy(balance, modes, strong)
    if args.current_manifest:
        saved_manifest = json.loads(args.current_manifest.read_text())
        base = Path(saved_manifest['run_directory'])
        if saved_manifest.get('exit_code') != 0 or saved_manifest['source_hashes_before'] != saved_manifest['source_hashes_after']:
            raise ValueError('stable completed current producer required')
        if saved_manifest['reference_manifest_sha256_after'] != digest(args.manifest):
            raise ValueError('current cache uses a different reference producer')
        objects = base/'objects.pickle'
        if digest(objects) != saved_manifest['artifacts']['objects.pickle']['sha256']:
            raise ValueError('current cache payload mismatch')
        frozen = base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
        if digest(frozen) != saved_manifest['source_hashes_after']['scripts/S11c_d_mixing_scattering_sympy_audit.py']:
            raise ValueError('current cache source snapshot mismatch')
        before, now = frozen.read_text(), Path(engine.__file__).read_text()
        constructors = ('UniformSlabCurrent', 'SlabEnergyBalance', 'ClosedAcousticEnergy')
        for name in constructors:
            if source_node(before, name) != source_node(now, name):
                raise ValueError(('cached constructor changed', name))
        old_ends = next(n for n in ast.parse(before).body if getattr(n, 'name', None) == 'ConstantEndPencil')
        old_phase = next(n for n in old_ends.body if getattr(n, 'name', None) == 'phase_terms')
        new_ends = next(n for n in ast.parse(now).body if getattr(n, 'name', None) == 'ConstantEndPencil')
        new_phase = next(n for n in new_ends.body if getattr(n, 'name', None) == 'phase_terms')
        if ast.dump(old_phase, include_attributes=False) != ast.dump(new_phase, include_attributes=False):
            raise ValueError('cached phase reduction changed')
        slab, conservative, bulk = pickle.loads(objects.read_bytes())
        current.construct = lambda anchoring, endpoint:conservative
        balance.construct = lambda anchoring, endpoint:slab
        acoustic.construct = lambda anchoring, endpoint:bulk
        provenance['currentCacheSha256'] = digest(objects)
        provenance['currentManifestSha256'] = digest(args.current_manifest)
        provenance['cachedConstructors'] = constructors+('ConstantEndPencil.phase_terms',)
    provenance['instrumentSha256'] = digest(Path(__file__))
    return engine.ClosedCurrentPairing(acoustic), inputs, provenance


def verified_checks_cache(base, provenance, derivatives):
    """Reuse completed coefficient proofs when only emission code changes."""
    records = [json.loads(line) for line in (base/'progress.jsonl').read_text().splitlines()]
    start = next(record for record in records if record.get('stage') == 'construction')
    finish = next(record for record in records if record.get('stage') == 'checks_constructed')
    report = json.loads((base/'checks.json').read_text())
    if any(value['nonzeroScalars'] for key, value in report['records'].items() if key.endswith('_RESIDUAL')):
        raise ValueError('cached verification contains nonzero residuals')
    for key in ('producerManifestSha256', 'cacheSha256', 'currentCacheSha256', 'currentManifestSha256', 'inputSha256'):
        if start['provenance'][key] != provenance[key]:
            raise ValueError(('verification cache input mismatch', key))
    source = base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
    helper = base/'source/_measurements/S11c_d_modal_current_check.py'
    if digest(source) != start['provenance']['engineSha256'] or digest(helper) != start['provenance']['instrumentSha256']:
        raise ValueError('verification cache source mismatch')
    old_source, new_source = source.read_text(), Path(engine.__file__).read_text()
    if source_node(old_source, 'polynomial_terms') != source_node(new_source, 'polynomial_terms'):
        raise ValueError('verification cache polynomial collection changed')
    old = next(n for n in ast.parse(old_source).body if getattr(n, 'name', None) == 'ClosedCurrentPairing')
    new = next(n for n in ast.parse(new_source).body if getattr(n, 'name', None) == 'ClosedCurrentPairing')
    methods = ('__init__', 'polarize', 'construct', 'rational_coefficient', 'carrier_expansion',
               'rational_reconstruction', 'wave_curve_reduction', 'split_balance_checks', 'split_derivative_checks')
    for name in methods:
        before = next(n for n in old.body if getattr(n, 'name', None) == name)
        after = next(n for n in new.body if getattr(n, 'name', None) == name)
        if ast.dump(before, include_attributes=False) != ast.dump(after, include_attributes=False):
            raise ValueError(('verification cache method changed', name))
    def diagonal_statements(source):
        return tuple(ast.dump(node, include_attributes=False) for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id in
                ('same_frequency', 'diagonal') for target in node.targets))
    if diagonal_statements(helper.read_text()) != diagonal_statements(Path(__file__).read_text()):
        raise ValueError('verification cache diagonal join changed')
    objects = base/'check_objects.pickle'
    if digest(objects) != finish['sha256']:
        raise ValueError('verification cache payload mismatch')
    checks = pickle.loads(objects.read_bytes())
    if derivatives != ('NORMAL_RADICAL_TRANSPORT' in checks):
        raise ValueError('verification cache derivative selection mismatch')
    provenance.update({'verificationCacheSha256':digest(objects), 'verificationCacheMethods':methods,
                       'verificationProducerEngineSha256':digest(source),
                       'verificationProducerInstrumentSha256':digest(helper)})
    return checks


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--current-manifest', type=Path)
    parser.add_argument('--cache-result', type=Path, required=True)
    parser.add_argument('--construction-cache', type=Path)
    parser.add_argument('--verification-cache', type=Path)
    parser.add_argument('--emit', action='store_true')
    parser.add_argument('--balances', action='store_true')
    parser.add_argument('--split-balances', action='store_true')
    parser.add_argument('--derivatives', action='store_true')
    parser.add_argument('--checks-json', type=Path)
    parser.add_argument('--require-zero-residuals', action='store_true')
    args = parser.parse_args()
    args.split_balances = args.split_balances or args.balances or args.derivatives or not args.emit
    args.end = 'REFERENCE'
    started = time.monotonic()
    def progress(record):
        record.setdefault('elapsedSeconds', time.monotonic()-started)
        line = json.dumps(record)
        if args.emit:
            with (args.cache_result.parent/'progress.jsonl').open('a') as stream:
                stream.write(line+'\n')
        else:
            print(line, flush=True)
    faulthandler.dump_traceback_later(300, repeat=True)
    builder, inputs, provenance = build(args)
    progress({'stage':'construction', 'provenance':provenance})
    if args.construction_cache:
        base = args.construction_cache
        records = [json.loads(line) for line in (base/'progress.jsonl').read_text().splitlines()
                   if line.startswith('{')]
        start = next(record for record in records if record.get('stage') == 'construction')
        finish = next(record for record in records if record.get('stage') == 'constructed')
        source = base/'source/scripts/S11c_d_mixing_scattering_sympy_audit.py'
        helper = base/'source/_measurements/S11c_d_modal_current_check.py'
        if digest(source) != start['provenance']['engineSha256'] or digest(helper) != start['provenance']['instrumentSha256']:
            raise ValueError('construction cache producer source mismatch')
        for key in ('producerManifestSha256', 'cacheSha256', 'currentCacheSha256', 'currentManifestSha256'):
            if start['provenance'][key] != provenance[key]:
                raise ValueError(('construction cache input mismatch', key))
        old = next(n for n in ast.parse(source.read_text()).body if getattr(n, 'name', None) == 'ClosedCurrentPairing')
        new = next(n for n in ast.parse(Path(engine.__file__).read_text()).body if getattr(n, 'name', None) == 'ClosedCurrentPairing')
        for name in ('__init__', 'polarize', 'construct'):
            old_method = next(n for n in old.body if getattr(n, 'name', None) == name)
            new_method = next(n for n in new.body if getattr(n, 'name', None) == name)
            if ast.dump(old_method, include_attributes=False) != ast.dump(new_method, include_attributes=False):
                raise ValueError(('construction cache method changed', name))
        objects = base/'objects.pickle'
        if digest(objects) != finish['cacheSha256']:
            raise ValueError('construction cache artifact mismatch')
        result, known = pickle.loads(objects.read_bytes())
        engine.PHYSICAL_METADATA.dimensions.known.update(known)
        builder.construct = lambda anchoring, endpoint:result
        provenance['constructionCacheSha256'] = digest(objects)
    else:
        result = builder.construct('LAB_HELD', None)
    args.cache_result.write_bytes(pickle.dumps((result, engine.PHYSICAL_METADATA.dimensions.known), protocol=5))
    progress({'stage':'constructed', 'cacheSha256':digest(args.cache_result),
              'wallSeconds':time.monotonic()-started,
              'peakRssKiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    modes = builder.modes
    algebraic, relation, _ = modes.analytic(builder.acoustic.strong)
    mapping = inputs.mapping(algebraic, relation, (modes.k, modes.q, modes.eta, modes.sigma, modes.r.omega))
    mapping.update({modes.eta:0, modes.sigma:0})
    cached_checks = verified_checks_cache(args.verification_cache, provenance, args.derivatives) if args.verification_cache else None
    if args.emit:
        engine.physical('MODAL_CURRENT_PREFLIGHT_PROVENANCE', provenance)
        builder.emit('LAB_HELD', None, 'REFERENCE_LAB_HELD_RHO4_CONSTANT')
        dims = engine.PHYSICAL_METADATA.dimensions
        binding_record = {str(key):value for key, value in mapping.items()}
        binding_units = {str(key):dims.measure(key) for key in mapping}
        engine.emit('MODAL_CURRENT_PREFLIGHT_BINDINGS', binding_record)
        engine.emit('METADATA_MODAL_CURRENT_PREFLIGHT_BINDINGS', modes.numeric_metadata(
            engine.cas(binding_record), lambda path:binding_units[path[0]]))
        engine.physical('MODAL_CURRENT_PREFLIGHT_DIMENSION_CONSTRAINTS', tuple(dims.constraints))
    if args.balances or args.split_balances or args.emit:
        if cached_checks is not None:
            checks = cached_checks
        elif args.split_balances:
            checks = builder.split_balance_checks(result, mapping, progress)
            if args.derivatives:
                checks.update(builder.split_derivative_checks(result, checks, mapping, progress))
        else:
            checks = {'SLAB_BALANCE_RESIDUAL':result['SLAB_BALANCE_RESIDUAL'].xreplace(mapping).applyfunc(
                builder.rational_coefficient)}
        same_frequency = dict.fromkeys(builder.frequencies, builder.r.omega)
        diagonal = (result['SLAB_CURRENT_MATRIX'].xreplace(same_frequency)-
                    builder.balance.construct('LAB_HELD', None)['SLAB_CURRENT_MATRIX']).applyfunc(sp.expand)
        checks['SLAB_CURRENT_DIAGONAL_JOIN_RESIDUAL'] = diagonal
        checks_cache = args.cache_result.parent/'check_objects.pickle'
        checks_cache.write_bytes(pickle.dumps(checks, protocol=5))
        progress({'stage':'checks_constructed', 'sha256':digest(checks_cache)})
        report = {'provenance':provenance, 'records':{}}
        for key, value in checks.items():
            fingerprint = engine.carrier_fingerprint(engine.cas(value))
            source_key = key.removesuffix('_DIVISION_RESIDUAL')
            units = (builder.check_output_units(key, value, checks) if args.split_balances else
                     builder.output_units(source_key, value))
            if key == 'SLAB_CURRENT_DIAGONAL_JOIN_RESIDUAL':
                units = builder.output_units('SLAB_CURRENT_MATRIX', value)
            if key.endswith('_DIVISION_RESIDUAL') and not args.split_balances:
                # The division reconstruction is a dimensionless polynomial
                # identity in the explicitly declared coefficient unit frame.
                units = {(i,):engine.PHYSICAL_METADATA.dimensions.zero for i in range(len(value))}
            if key in ('EQUAL_DEPTH_WAVE_ROWS', 'EQUAL_DEPTH_WAVE_ELIMINANT'):
                unit = tuple(2*v for v in engine.PHYSICAL_METADATA.dimensions.measure(modes.r.omega))
                units = {(i,):unit for i in range(2)} if key.endswith('_ROWS') else {():unit}
            if args.emit:
                if key.endswith('_RESIDUAL') and all(v == 0 for _, v in engine.leaves(engine.cas(value))):
                    engine.physical('MODAL_CURRENT_PREFLIGHT_'+key, value, zero_dimensions=units)
                elif key.startswith(('NORMAL_', 'FREQUENCY_')) and key.endswith((
                        'BALANCE_DERIVATIVE_OPERANDS', 'SOURCE_POWER_DERIVATIVE_MATRIX',
                        'SOURCE_PRODUCT_DERIVATIVE_TERMS', 'CLOSED_PENCIL_DERIVATIVE_MATRIX')):
                    engine.emit('MODAL_CURRENT_PREFLIGHT_'+key, fingerprint)
                    engine.emit('METADATA_MODAL_CURRENT_PREFLIGHT_'+key, modes.numeric_metadata(
                        engine.cas(value), lambda path:units[path]))
                else:
                    engine.fingerprinted('MODAL_CURRENT_PREFLIGHT_'+key, engine.cas(value), units)
            scalars = [v for _, v in engine.leaves(engine.cas(value))]
            report['records'][key] = {'scalars':len(scalars), 'nonzeroScalars':sum(v != 0 for v in scalars),
                                     'fingerprintSha256':hashlib.sha256(sp.srepr(fingerprint).encode()).hexdigest()}
        report['wallSeconds'] = time.monotonic()-started
        report['peakRssKiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        report['dimensionConstraints'] = [sp.srepr(v) for v in engine.PHYSICAL_METADATA.dimensions.constraints]
        if args.checks_json:
            args.checks_json.write_text(json.dumps(report, indent=2)+'\n')
        progress({'stage':'balance_checks', 'records':report['records']})
        if args.emit:
            resources = {key:report[key] for key in ('wallSeconds', 'peakRssKiB')}
            engine.emit('MODAL_CURRENT_PREFLIGHT_RESOURCES', resources)
            engine.emit('METADATA_MODAL_CURRENT_PREFLIGHT_RESOURCES', modes.numeric_metadata(
                engine.cas(resources), lambda path:(0, 1, 0) if path[0] == 'wallSeconds' else (0, 0, 0)))
            index = engine.emission_index(engine.EMISSION_LINES)
            engine.emit('MODAL_CURRENT_PREFLIGHT_EMISSION_LINES', index)
            engine.emit('METADATA_MODAL_CURRENT_PREFLIGHT_EMISSION_LINES', modes.numeric_metadata(
                engine.cas(index), lambda path:(0, 0, 0)))
        faulthandler.cancel_dump_traceback_later()
        if args.require_zero_residuals and any(record['nonzeroScalars'] for key, record in report['records'].items()
                                             if key.endswith('_RESIDUAL')):
            raise ValueError('current construction residual; see emitted operands')
        return


if __name__ == '__main__':
    run()
