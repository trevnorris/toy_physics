#!/usr/bin/env python3
"""Normalize saved Fourier reconstruction pairs without rebuilding their operands."""
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
import S11c_d_source_fourier_factorization_check as original
from S11c_d_source_fourier_factorization_check import ROOT, STORE, engine, digest, save, atomic_pickle
from S11c_d_output_codec import decoded_lines, restore_emission_index
from ledger_fold import _restore

PREFIX = original.PREFIX+'_EXACT_RECONSTRUCTION'
PLAN = ROOT/'_measurements/S11c_d_source_fourier_residual_plan.md'
SOURCES = (*original.SOURCES, Path(__file__).resolve(), PLAN,
           ROOT/'_measurements/S11c_d_source_fourier_residual_repair_checkpoint.json')


def load(previous):
    checks = json.loads((previous/'checks.json').read_text())
    for name, item in checks['artifacts'].items():
        if digest(previous/name) != item['sha256']:
            raise ValueError(('saved factorization artifact changed', name))
    for item in checks['savedRows']:
        if digest(previous/item['path']) != item['sha256']:
            raise ValueError(('saved factorization row changed', item['path']))
    engine_name = str(engine.HERE.relative_to(ROOT))
    for name, sha in checks['sourceFiles'].items():
        if digest(previous/'source'/name) != sha:
            raise ValueError(('frozen factorization source changed', name))
        if name != engine_name and digest(ROOT/name) != sha:
            raise ValueError(('consumed factorization helper changed', name))
    before = ast.parse((previous/'source'/engine_name).read_text())
    after = ast.parse(engine.HERE.read_text())
    selected = next(n for n in after.body if getattr(n, 'name', None) == 'BoundedSourceFourierAssembly')
    additions = [n for n in selected.body if getattr(n, 'name', None) == 'reconstruction_certificate']
    if len(additions) != 1:
        raise ValueError('expected exactly one reconstruction-certificate addition')
    selected.body.remove(additions[0])
    if ast.dump(before) != ast.dump(after):
        raise ValueError('engine changed beyond the new certificate helper')
    pencil, assembly, provenance = original.load()
    with (previous/'factorization.pickle').open('rb') as stream:
        packet = pickle.load(stream)
    if packet['provenance'] != provenance or checks['provenance'] != provenance:
        raise ValueError('saved factorization provenance does not join current accepted operands')
    engine.PHYSICAL_METADATA.dimensions.__dict__.update(packet['dimensionState'])
    result = packet['result']
    builder = engine.BoundedSourceFourierAssembly(pencil.r, assembly['NONLOCAL_INTEGRALS'])
    joins = []
    for i, row in enumerate(result['ROWS']):
        bounded, source_index, source_limit, remaining = builder.limit_layout(assembly['NONLOCAL_INTEGRALS'][i])
        joins.append(row['INDEX'] == i and row['ORIGINAL'] == assembly['NONLOCAL_INTEGRALS'][i]
            and row['BOUNDED'] == bounded and row['SOURCE_LIMIT_INDEX'] == source_index
            and row['SOURCE_LIMIT'] == source_limit and row['REMAINING_LIMITS'] == remaining
            and all(f['SOURCE_INTEGRAL'] == sp.Integral(f['SOURCE'], source_limit) for f in row['FACTORS']))
    if len(joins) != len(assembly['NONLOCAL_INTEGRALS']) or not all(joins):
        raise ValueError('full original integral/source-limit joins failed')
    return pencil, packet, checks, joins


def pairs(result):
    for row in result['ROWS']:
        yield row['INDEX'], 'INTEGRAND', row['BOUNDED'].function, sp.Add(*(
            f['COEFFICIENT']*f['SOURCE'] for f in row['FACTORS'])), row['RECONSTRUCTION_RESIDUAL']
        for index, f in enumerate(row['FACTORS']):
            yield row['INDEX'], 'AMPLITUDE_'+str(index), f['SOURCE'], f['CHARACTER']*f['AMPLITUDE'], f['AMPLITUDE_RECONSTRUCTION_RESIDUAL']


def emit_all(result, pencil, provenance, records, certificate_hash):
    # The original operand/residual values are retained. Fingerprint the large
    # raw residual DAGs even when their shared-node count is deceptively small.
    saved_physical = engine.physical
    def compact(name, value, **kwargs):
        if '_INTEGRAND_RECONSTRUCTION_RESIDUAL_' in name:
            return engine.fingerprinted(name, engine.cas(value), kwargs.get('zero_dimensions'))
        return saved_physical(name, value, **kwargs)
    engine.physical = compact
    try:
        original.emit_result(result, pencil, provenance)
    finally:
        engine.physical = saved_physical
    dimensions = engine.PHYSICAL_METADATA.dimensions
    engine.physical(PREFIX+'_CERTIFICATE_PACKET_SHA256', certificate_hash)
    for i, record in enumerate(records):
        certificate = record['certificate']; unit = dimensions.measure(certificate['LEFT'])
        suffix = str(i)
        engine.fingerprinted(PREFIX+'_OPERANDS_'+suffix,
            sp.Tuple(certificate['LEFT'], certificate['RIGHT']), {(0,): unit, (1,): unit})
        engine.fingerprinted(PREFIX+'_RAW_RESIDUAL_'+suffix, record['raw'], {(): unit})
        engine.physical(PREFIX+'_RATIONAL_RESIDUAL_'+suffix, certificate['RESIDUAL'], zero_dimensions={(): unit})
        engine.physical(PREFIX+'_SHARED_REPLAY_RESIDUALS_'+suffix, sp.Tuple(*certificate['REPLAY_RESIDUALS']),
            zero_dimensions={(j,): unit for j in range(len(certificate['REPLAY_RESIDUALS']))})
        phase_residuals = tuple(v[1] for v in certificate['PHASE_SPLITS'].values())
        root_residuals = tuple(v[2] for v in certificate['RADICAL_POWERS'].values())
        engine.physical(PREFIX+'_EXPONENT_RESIDUALS_'+suffix, sp.Tuple(*phase_residuals, *root_residuals),
            zero_dimensions={(j,): dimensions.zero for j in range(len(phase_residuals)+len(root_residuals))})
        engine.fingerprinted(PREFIX+'_ORIGINAL_DENOMINATOR_BASES_'+suffix,
            sp.Tuple(*certificate['ORIGINAL_DENOMINATOR_BASES']))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-directory', type=Path, required=True)
    parser.add_argument('--source-run-directory', type=Path, required=True)
    args = parser.parse_args()
    previous = args.source_run_directory.resolve(); previous.relative_to(STORE)
    base = args.run_directory.resolve(); base.relative_to(STORE); base.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); paths = tuple(dict.fromkeys(SOURCES))
    pins = {str(p.relative_to(ROOT)): digest(p) for p in paths}
    for p in paths:
        target = base/'source'/p.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(p, target)
    pencil, packet, old_checks, joins = load(previous)
    result, provenance = packet['result'], packet['provenance']
    original_packet = base/'source-factorization.pickle'; shutil.copyfile(previous/'factorization.pickle', original_packet)
    source_hash = digest(original_packet)
    save(base/'preflight.json', {'sourceFiles': pins, 'sourceRunDirectory': str(previous),
        'sourceChecksSha256': digest(previous/'checks.json'), 'sourcePacketSha256': source_hash,
        'fullConstructorAstJoin': True, 'originalIntegralLimitJoins': joins, 'provenance': provenance})
    records, residuals, inventory = [], [], []
    destination = base/'certificates'; destination.mkdir()
    for row, kind, left, right, raw in pairs(result):
        if raw == 0:
            residuals.append(raw)
            continue
        certificate = engine.BoundedSourceFourierAssembly.reconstruction_certificate(left, right,
            shared=kind == 'INTEGRAND')
        record = {'row': row, 'kind': kind, 'raw': raw, 'certificate': certificate}
        path = destination/(str(row).zfill(3)+'-'+kind.lower()+'.pickle')
        atomic_pickle(path, record)
        inventory.append({'path': str(path.relative_to(base)), 'bytes': path.stat().st_size, 'sha256': digest(path)})
        save(base/'certificate-inventory.json', inventory)
        records.append(record); residuals.append(certificate['RESIDUAL'])
        with (base/'progress.jsonl').open('a') as stream:
            stream.write(json.dumps({'row': row, 'kind': kind, 'certificates': len(records),
                'nonzeroNormalizedResidual': certificate['RESIDUAL'] != 0,
                'wallSeconds': time.monotonic()-started})+'\n')
    residuals.extend(f[k] for r in result['ROWS'] for f in r['FACTORS'] for k in
        ('CHARACTER_NORMALIZATION_RESIDUAL','CHARACTER_EQUATION_RESIDUAL'))
    residuals.extend(v[k] for _, v in result['PHASES'] for k in ('EXPONENT_RESIDUAL','SECOND_SOURCE_DERIVATIVE'))
    proof_residuals = [v for r in records for v in r['certificate']['REPLAY_RESIDUALS']]
    proof_residuals.extend(v[1] for r in records for v in r['certificate']['PHASE_SPLITS'].values())
    proof_residuals.extend(v[2] for r in records for v in r['certificate']['RADICAL_POWERS'].values())
    certificate_packet = base/'reconstruction-certificates.pickle'
    atomic_pickle(certificate_packet, {'records': records, 'sourcePacketSha256': source_hash,
        'provenance': provenance, 'dimensionState': dict(vars(engine.PHYSICAL_METADATA.dimensions))})
    before = digest(certificate_packet)
    with (base/'full.out').open('x') as stream, contextlib.redirect_stdout(stream):
        emit_all(result, pencil, provenance, records, before)
        keys = {tag: 's11cd'+''.join(w.title() for w in tag.removeprefix('PY_S11CD_').split('_'))
            for tag in engine.EMISSION_LINES if not tag.startswith('PY_S11CD_METADATA_')}
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        index = engine.emission_index(engine.EMISSION_LINES)
        index_units = {p: (0,0,0) for p,_ in engine.leaves(engine.cas(index))}
        engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=index_units)
    entries = {}
    for line in decoded_lines(base/'full.out'):
        tag, _, body = line.rstrip('\n').partition(': ')
        if tag in entries: raise ValueError('duplicate reconstruction-recovery emission')
        entries[tag] = _restore(body)
    seen = set(); old_emit = engine.emit
    def compare(name, value):
        tag = 'PY_S11CD_'+name
        if tag in seen or entries.get(tag) != engine.cas(value):
            raise ValueError(('reconstruction-recovery emission mismatch', tag))
        seen.add(tag)
    engine.emit = compare
    try:
        emit_all(result, pencil, provenance, records, before)
        engine.physical(PREFIX+'_WRITE_KEYS', keys)
        engine.physical(PREFIX+'_EMISSION_LINES', index, zero_dimensions=index_units)
    finally:
        engine.emit = old_emit
    if seen != set(entries) or len(keys) != len(set(keys.values())) or set(keys.values()) & set(engine.IMPORT_KEYS):
        raise ValueError('reconstruction-recovery emission/write-key census')
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    restore_emission_index({str(k):v for k,v in entries[final]}, list(entries)[:list(entries).index(final)])
    metadata_paths = 0
    for tag, body in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'): continue
        for path, descriptor in body:
            fields = {str(k):v for k,v in descriptor}; unit = fields['DIMENSION_L_T_M']
            if not isinstance(unit, sp.Tuple) or len(unit)!=3 or any(s.free_symbols for s in unit):
                raise ValueError(('unresolved reconstruction-recovery unit', tag, path))
            if any(k not in fields for k in ('MULTIGRADE','EPSILON_LAMBDA_SUPPORT')):
                raise ValueError(('missing reconstruction-recovery grades', tag))
            metadata_paths += 1
    atomic_pickle(base/'dimensions-after-emission.pickle', dict(vars(engine.PHYSICAL_METADATA.dimensions)))
    summary = {'runDirectory': str(base), 'sourceRunDirectory': str(previous), 'sourceFiles': pins,
        'sourceChecksSha256': digest(previous/'checks.json'), 'sourcePacketSha256': source_hash,
        'provenance': provenance, 'fullConstructorAstJoin': True, 'originalIntegralLimitJoins': joins,
        'factorizedIntegralCount': len(result['ROWS']), 'certificateCount': len(records),
        'rawResidualScalars': old_checks['residualScalars'], 'nonzeroRawResidualScalars': old_checks['nonzeroResidualScalars'],
        'normalizedResidualScalars': len(residuals), 'nonzeroNormalizedResidualScalars': sum(v!=0 for v in residuals),
        'proofResidualScalars': len(proof_residuals), 'nonzeroProofResidualScalars': sum(v!=0 for v in proof_residuals),
        'packetSha256BeforeEmission': before, 'packetSha256AfterEmission': digest(certificate_packet),
        'sourcePacketUnchanged': source_hash == digest(original_packet) == digest(previous/'factorization.pickle'),
        'tagCount': len(entries), 'metadataPaths': metadata_paths, 'writeKeyCount': len(keys),
        'certificateArtifacts': inventory, 'wallSeconds': time.monotonic()-started,
        'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'artifacts': {p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in base.iterdir() if p.suffix in ('.pickle','.out')},
        'scope': 'Exact reconstruction certificates for saved finite-domain source Fourier operands. Original raw residuals, all Piecewise operands and denominator domains are retained. No numerical action convergence, extension onto singular loci, unbounded integration interchange, Abel weak limit or scattering result.'}
    save(base/'checks.json', summary)
    if pins != {str(p.relative_to(ROOT)):digest(p) for p in paths} or before != digest(certificate_packet) or not summary['sourcePacketUnchanged']:
        raise ValueError('reconstruction-recovery source/packet changed')
    if any(digest(base/item['path'])!=item['sha256'] for item in inventory):
        raise ValueError('reconstruction certificate changed')
    if any(v!=0 for v in (*residuals,*proof_residuals)) or engine.PHYSICAL_METADATA.dimensions.constraints:
        raise ValueError('saved-pair reconstruction or dimension requires inspection')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
