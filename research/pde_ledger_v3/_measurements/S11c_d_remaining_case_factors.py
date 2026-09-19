#!/usr/bin/env python3
"""Factor only distinct new case integrals, retaining native rows and proofs."""
import argparse
import contextlib
import copy
import faulthandler
import json
import multiprocessing as mp
import os
from pathlib import Path
import resource
import shutil
import signal
import sys
import time
import sympy as sp
import S11c_d_remaining_case_sources as cases
import S11c_d_source_fourier_factorization_check as output
import S11c_d_source_fourier_residual_recover as proof_output

f, engine, grades = cases.f, cases.engine, cases.grades
PLAN = f.M/'S11c_d_remaining_case_factors_plan.md'
SOURCE = f.M/'S11c_d_remaining_case_sources_checkpoint.json'
FACTOR = f.M/'S11c_d_source_fourier_factorization_checkpoint.json'
PREFIX = 'REMAINING_CASE_FOURIER'


def artifact(path):
    return {'bytes': path.stat().st_size, 'sha256': f.digest(path)}


def context(source, actions, assembly):
    r, dimensions = cases.source.restore_context(source)
    dimensions.__dict__.update(assembly['dimensionState'])
    values = {key: engine.named(body, 'VALUE') for (key, _), body in source['payloads'].items()}
    pencil = engine.ReducedPencil(*(values[key] for key in engine.CLOSED_KEYS), r)
    f.require(pencil.strong == actions['strong'] and pencil.kernel == actions['kernel']
              and pencil.fields == actions['fields'] and pencil.probes == actions['probes'],
              'actual full case source/field/kernel identities')
    return r, dimensions, pencil


def load(base):
    cp = json.loads(SOURCE.read_text()); origin = Path(cp['runDirectory'])
    f.require(cp['status'] == 'PUBLISHED_ANNEX_VERIFIED', 'accepted all-case source checkpoint')
    f.require((f.ROOT/cp['publication']['path']).is_symlink() and
              f.digest(f.ROOT/cp['publication']['path']) == cp['publication']['sha256'], 'published source identity')
    f.require(f.digest(origin/'checks.json') == cp['checksSha256'], 'source final checks identity')
    for name, sha in cp['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(origin/'source'/name) == sha, ('current/frozen source', name))
    accepted, fc, fp, fcp = grades.packet(FACTOR.name, 'source-factorization.pickle')
    certificates_path = Path(fc['runDirectory'])/'reconstruction-certificates.pickle'
    f.require(f.digest(certificates_path) == fc['artifacts'][certificates_path.name]['sha256'], 'accepted factor certificates')
    certificate_packet = f.unpickle(certificates_path)
    f.require(f.digest(fp) == fc['sourcePacketSha256'] == certificate_packet['sourcePacketSha256'], 'factor/proof packet join')
    for record in certificate_packet['records']:
        cert = record['certificate']
        f.require(cert['RESIDUAL'] == 0 and all(v == 0 for v in proof_scalars(cert)), 'accepted exact factor proof')
    helper_join = cases.definition_joins(Path(fc['runDirectory'])/'source'/str(engine.HERE.relative_to(f.ROOT)),
                                         {'BoundedSourceFourierAssembly'})
    states, inputs = {}, {str(p): f.digest(p) for p in (fp, certificates_path, origin/'checks.json')}
    for case in cases.CASES:
        label = '__'.join(case); location = origin if case == cases.BASELINE else origin/'cases'/label
        prefix = 'accepted-' if case == cases.BASELINE else ''
        packets = {}
        for kind in ('reduced-action', 'actions', 'assembly'):
            path = location/(prefix+kind+'.pickle'); name = str(path.relative_to(origin))
            f.require(f.digest(path) == cp['artifacts'][name]['sha256'], ('accepted complete case packet', name))
            dest = base/'accepted-cases'/label/(kind+'.pickle'); dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, dest); f.require(f.digest(dest) == f.digest(path), 'byte-identical case copy')
            packets[kind] = f.unpickle(dest); inputs[str(path)] = f.digest(path)
        r, dimensions, pencil = context(packets['reduced-action'], packets['actions'], packets['assembly'])
        builder = engine.BoundedSourceFourierAssembly(r, packets['assembly']['result']['NONLOCAL_INTEGRALS'])
        states[label] = {'packets': packets, 'r': r, 'dimensions': dimensions, 'pencil': pencil, 'builder': builder}
    baseline = states['__'.join(cases.BASELINE)]; originals = baseline['packets']['assembly']['result']['NONLOCAL_INTEGRALS']
    f.require(len(originals) == len(accepted['result']['ROWS']), 'complete accepted factor/source census')
    exact = {}; unique = []
    for i, (original, row) in enumerate(zip(originals, accepted['result']['ROWS'])):
        f.require(row['ORIGINAL'] == original and row['INDEX'] == i, 'accepted original factor address')
        exact[original] = {'kind': 'accepted', 'index': i}
    catalogue = {}
    for label, state in states.items():
        r, dimensions, pencil = context(state['packets']['reduced-action'], state['packets']['actions'], state['packets']['assembly'])
        builder = state['builder']; engine.PHYSICAL_METADATA.dimensions.__dict__.update(vars(state['dimensions']))
        f.require(builder.cutoffs == accepted['result']['CUTOFFS'] and r.zp == baseline['r'].zp
                  and pencil.fields == baseline['pencil'].fields, 'inherited cutoff and field coordinates')
        addresses = []
        for i, original in enumerate(builder.integrals):
            bounded, si, sl, other = builder.limit_layout(original)
            if original not in exact:
                index = len(unique); exact[original] = {'kind': 'new', 'index': index}
                unique.append({'index': index, 'owner': label, 'caseIndex': i, 'original': original,
                               'bounded': bounded, 'sourceLimitIndex': si, 'sourceLimit': sl, 'remainingLimits': other})
            address = exact[original]
            if address['kind'] == 'accepted':
                row = accepted['result']['ROWS'][address['index']]
                f.require((row['BOUNDED'], row['SOURCE_LIMIT_INDEX'], row['SOURCE_LIMIT'], row['REMAINING_LIMITS'])
                          == (bounded, si, sl, other), 'exact baseline source and ordered-limit consumption')
            else:
                entry = unique[address['index']]
                f.require((entry['bounded'], entry['sourceLimitIndex'], entry['sourceLimit'], entry['remainingLimits'])
                          == (bounded, si, sl, other), 'exact inter-case source/limit consumption')
            addresses.append({'caseIndex': i, **address})
        catalogue[label] = addresses
        # The address lookup must reject an actual changed physical limit.
        original = builder.integrals[0]; variable, lower, upper = original.limits[0]
        wrong = sp.Integral(original.function, (variable, lower, sp.S.Zero), *original.limits[1:])
        f.require(upper != 0 and wrong not in exact, 'wrong-limit reuse rejection')
    for path, name in ((fp, 'accepted-factorization.pickle'), (certificates_path, 'accepted-factor-certificates.pickle')):
        shutil.copyfile(path, base/name); f.require(f.digest(base/name) == f.digest(path), 'byte-identical accepted factor/proof copy')
    paths = {Path(m.__file__).resolve() for m in tuple(sys.modules.values()) if getattr(m, '__file__', None)
             and Path(m.__file__).resolve().is_relative_to(f.ROOT) and Path(m.__file__).suffix == '.py'}
    paths.update((Path(__file__).resolve(), PLAN, SOURCE, FACTOR, engine.HERE, f.ACCEPTANCE,
                  f.M/'S11c_d_variable_profile_development_input.json', f.M/'S11c_d_focused_completion_plan.md',
                  f.ROOT/'directives/S11c_d_NONLINEAR_POLE_CONTRACT.md'))
    paths.update(f.ROOT/name for name in cp['sourceFiles'] if name.startswith('directives/') or name.endswith('_exports.py'))
    pins = {str(p.relative_to(f.ROOT)): f.digest(p) for p in paths}
    for name in pins:
        dest = base/'source'/name; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(f.ROOT/name, dest)
    manifest = {'sourceFiles': pins, 'inputPackets': inputs, 'nativeHelperJoins': helper_join,
                'acceptedFactorSourceSha256': f.digest(fp), 'newDistinctIntegrals': len(unique),
                'cases': {label: {'integrals': len(rows), 'baselineReused': sum(v['kind'] == 'accepted' for v in rows),
                                  'newFactorUses': sum(v['kind'] == 'new' for v in rows)} for label, rows in catalogue.items()},
                'scope': 'Only distinct new source factorizations; no numerical quadrature, source reduction, mode or response solve.'}
    f.atomic_pickle(base/'integral-catalogue.pickle', {'unique': unique, 'cases': catalogue})
    f.save(base/'inputs.json', manifest)
    return states, accepted, unique, catalogue, manifest


def proof_scalars(certificate):
    return (*certificate['REPLAY_RESIDUALS'], *(v[1] for v in certificate['PHASE_SPLITS'].values()),
            *(v[2] for v in certificate['RADICAL_POWERS'].values()))


def worker(directory, entry, state):
    with (directory/'stdout').open('x') as out, (directory/'stderr').open('x') as err:
        os.dup2(out.fileno(), 1); os.dup2(err.fileno(), 2)
        resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
        faulthandler.dump_traceback_later(150)
        started = time.monotonic()
        r, dimensions, pencil = context(state['packets']['reduced-action'], state['packets']['actions'], state['packets']['assembly'])
        f.atomic_pickle(directory/'operands.pickle', entry)
        builder = engine.BoundedSourceFourierAssembly(r, (entry['original'],))
        def checkpoint(row, phases):
            f.atomic_pickle(directory/'raw-factorization.pickle', {'row': row, 'phases': phases,
                           'dimensionState': dict(vars(dimensions))})
        result = builder.construct(checkpoint=checkpoint)
        row = result['ROWS'][0]
        f.require((row['BOUNDED'], row['SOURCE_LIMIT_INDEX'], row['SOURCE_LIMIT'], row['REMAINING_LIMITS']) ==
                  (entry['bounded'], entry['sourceLimitIndex'], entry['sourceLimit'], entry['remainingLimits']), 'native row/layout join')
        records = []; proofs = []
        for _, kind, left, right, raw in proof_output.pairs(result):
            pair = {'left': left, 'right': right, 'raw': raw, 'representations': tuple(sp.srepr(v) for v in (left, right))}
            pair['representationHashes'] = tuple(__import__('hashlib').sha256(v.encode()).hexdigest() for v in pair['representations'])
            f.atomic_pickle(directory/(kind.lower()+'-pair.pickle'), pair)
            certificate = engine.BoundedSourceFourierAssembly.reconstruction_certificate(left, right, shared=False)
            record = {'row': entry['index'], 'kind': kind, 'raw': raw, 'certificate': certificate}
            f.atomic_pickle(directory/(kind.lower()+'-certificate.pickle'), record)
            f.require(certificate['LEFT'] == left and certificate['RIGHT'] == right and certificate['RESIDUAL'] == 0,
                      'actual native uncompressed reconstruction certificate')
            checks = proof_scalars(certificate); f.require(all(v == 0 for v in checks), 'exact phase/root/replay proofs')
            records.append(record); proofs.extend(checks)
        first = row['FACTORS'][0]
        reconstructed = sp.Add(*(v['COEFFICIENT']*v['SOURCE'] for v in row['FACTORS']))
        mutation = engine.BoundedSourceFourierAssembly.reconstruction_certificate(
            row['BOUNDED'].function, reconstructed+first['COEFFICIENT']*first['SOURCE'], shared=False)
        f.atomic_pickle(directory/'coefficient-mutation.pickle', mutation)
        f.require(mutation['RESIDUAL'] != 0 and all(v == 0 for v in proof_scalars(mutation)), 'actual one-factor coefficient mutation')
        residuals = [v[k] for v in row['FACTORS'] for k in ('CHARACTER_NORMALIZATION_RESIDUAL', 'CHARACTER_EQUATION_RESIDUAL')]
        residuals += [v[k] for _, v in result['PHASES'] for k in ('EXPONENT_RESIDUAL', 'SECOND_SOURCE_DERIVATIVE')]
        f.require(all(v == 0 for v in residuals), 'actual source-character and phase residuals')
        f.atomic_pickle(directory/'factorization.pickle', {'result': result, 'certificates': records,
                        'dimensionState': dict(vars(dimensions)), 'entry': entry})
        checks = {'index': entry['index'], 'owner': entry['owner'], 'certificates': len(records),
                  'zeroProofScalars': len(proofs), 'zeroCharacterScalars': len(residuals), 'coefficientMutationNonzero': True,
                  'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  'artifacts': {p.name: artifact(p) for p in directory.glob('*.pickle')}}
        f.save(directory/'checks.json', checks); print(json.dumps(checks), flush=True)
        faulthandler.cancel_dump_traceback_later()


def run_row(base, entry, state):
    directory = base/'rows'/str(entry['index']).zfill(3); directory.mkdir(parents=True)
    child = mp.get_context('fork').Process(target=worker, args=(directory, entry, state))
    started = time.monotonic(); child.start()
    try:
        child.join(180)
        if child.is_alive():
            child.terminate(); child.join(5)
            if child.is_alive(): child.kill(); child.join()
            raise TimeoutError('bounded native factorization row; retain all completed operands')
    finally:
        if child.is_alive(): child.terminate(); child.join(5)
        f.save(directory/'outcome.json', {'exitCode': child.exitcode, 'wallSeconds': time.monotonic()-started,
                'stderrBytes': (directory/'stderr').stat().st_size if (directory/'stderr').exists() else None})
    f.require(child.exitcode == 0 and (directory/'stderr').read_bytes() == b'', 'clean completed factorization worker')
    checks = json.loads((directory/'checks.json').read_text())
    f.require(json.loads((directory/'stdout').read_text()) == checks, 'worker checks/stdout identity')
    return checks


def verify_row(directory, entry):
    checks = json.loads((directory/'checks.json').read_text()); outcome = json.loads((directory/'outcome.json').read_text())
    f.require(outcome['exitCode'] == 0 and (directory/'stderr').read_bytes() == b'' and
              json.loads((directory/'stdout').read_text()) == checks, 'completed saved worker')
    for name, item in checks['artifacts'].items():
        f.require(artifact(directory/name) == item, ('worker artifact', name))
    packet = f.unpickle(directory/'factorization.pickle')
    f.require(packet['entry'] == entry and packet['result']['ROWS'][0]['ORIGINAL'] == entry['original'], 'saved row full source address')
    for record in packet['certificates']:
        f.require(record['certificate']['RESIDUAL'] == 0 and all(v == 0 for v in proof_scalars(record['certificate'])), 'saved exact certificate')
    f.require(f.unpickle(directory/'coefficient-mutation.pickle')['RESIDUAL'] != 0, 'saved responding mutation')
    return packet, checks


def emit(base, states, unique, catalogue, packets, manifest):
    # The unchanged emitter receives each actual new native row in its owner context.
    # Accepted baseline factors remain joined to their verified published transcript.
    old_prefix = output.PREFIX; old_physical = engine.physical
    def compact(name, value, **kwargs):
        if '_INTEGRAND_RECONSTRUCTION_RESIDUAL_' in name:
            return engine.fingerprinted(name, engine.cas(value), kwargs.get('zero_dimensions'))
        return old_physical(name, value, **kwargs)
    try:
        engine.physical = compact
        for index, packet in sorted(packets.items()):
            entry = unique[index]; state = states[entry['owner']]
            r, dimensions, pencil = context(state['packets']['reduced-action'], state['packets']['actions'], state['packets']['assembly'])
            dimensions.__dict__.update(packet['dimensionState'])
            output.PREFIX = PREFIX+'_NEW_'+str(index)
            output.emit_result(packet['result'], pencil, {'source': manifest['acceptedFactorSourceSha256'], 'unionIndex': index})
            for j, record in enumerate(packet['certificates']):
                cert = record['certificate']; unit = dimensions.measure(cert['LEFT'])
                tag = output.PREFIX+'_CERTIFICATE_'+str(j)
                engine.fingerprinted(tag+'_PAIR', sp.Tuple(cert['LEFT'], cert['RIGHT']), {(0,): unit, (1,): unit})
                engine.physical(tag+'_NORMALIZED', cert['RESIDUAL'], zero_dimensions={(): unit})
                residuals = sp.Tuple(*proof_scalars(cert))
                engine.physical(tag+'_PROOF_RESIDUALS', residuals, zero_dimensions={
                    (i,): unit if i < len(cert['REPLAY_RESIDUALS']) else dimensions.zero for i in range(len(residuals))})
            mutation = f.unpickle(base/'rows'/str(index).zfill(3)/'coefficient-mutation.pickle')
            engine.fingerprinted(output.PREFIX+'_COEFFICIENT_MUTATION', mutation['RESIDUAL'],
                                 {(): dimensions.measure(packet['result']['ROWS'][0]['BOUNDED'].function)})
    finally:
        output.PREFIX = old_prefix; engine.physical = old_physical
    cases.boundary.structural_flags(PREFIX+'_CASE_ADDRESSES', catalogue)
    cases.boundary.structural_flags(PREFIX+'_INPUTS', manifest)


def output_replay(base, states, unique, catalogue, packets, manifest):
    structural = cases.boundary.structural_flags
    def emit_all(): emit(base, states, unique, catalogue, packets, manifest)
    engine.EMISSION_LINES.clear(); engine.PAYLOAD_ENCODER = grades.PayloadEncoder()
    with (base/'full.out').open('x') as stream, contextlib.redirect_stdout(stream):
        emit_all()
        keys = {tag: 's11cdRemainingCaseFactor'+str(i) for i, tag in enumerate(engine.EMISSION_LINES)
                if not tag.startswith('PY_S11CD_METADATA_')}
        structural(PREFIX+'_WRITE_KEYS', keys); index = engine.emission_index(engine.EMISSION_LINES)
        structural(PREFIX+'_EMISSION_LINES', index)
    entries = {}
    for line in grades.decoded_lines(base/'full.out'):
        tag, _, body = line.rstrip('\n').partition(': '); f.require(tag not in entries, 'unique emitted factor tag')
        entries[tag] = grades._restore(body)
    previous, seen = engine.emit, set()
    def compare(name, value):
        tag = 'PY_S11CD_'+name
        f.require(tag not in seen and entries.get(tag) == engine.cas(value), ('complete factor emission replay', tag)); seen.add(tag)
    engine.emit = compare
    try: emit_all(); structural(PREFIX+'_WRITE_KEYS', keys); structural(PREFIX+'_EMISSION_LINES', index)
    finally: engine.emit = previous
    f.require(seen == set(entries) and len(keys) == len(set(keys.values())) and not set(keys.values()) & set(engine.IMPORT_KEYS), 'full factor/key replay')
    final = 'PY_S11CD_'+PREFIX+'_EMISSION_LINES'
    grades.restore_emission_index({str(k): v for k, v in entries[final]}, list(entries)[:list(entries).index(final)])
    count = 0
    for tag, value in entries.items():
        if not tag.startswith('PY_S11CD_METADATA_'): continue
        for path, descriptor in value:
            d = {str(k): v for k, v in descriptor}; unit = d['DIMENSION_L_T_M']
            f.require(isinstance(unit, sp.Tuple) and len(unit) == 3 and all(not v.free_symbols for v in unit)
                      and 'MULTIGRADE' in d and 'EPSILON_LAMBDA_SUPPORT' in d, 'complete units and independent grades')
            count += 1
    return {'tagCount': len(entries), 'writeKeys': len(keys), 'metadataPaths': count}


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--mode', choices=('preflight', 'construct'), required=True)
    parser.add_argument('--run-directory', type=Path, required=True); parser.add_argument('--resume-from', type=Path)
    args = parser.parse_args(); base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); started = time.monotonic()
    def budget(*_): raise TimeoutError('remaining case factorization budget; keep every completed native row and proof')
    signal.signal(signal.SIGALRM, budget); signal.alarm(900)
    states, accepted, unique, catalogue, manifest = load(base)
    # One actual new source from each owning case, plus one two-momentum source.
    selected = list(range(len(unique)))
    if args.mode == 'preflight':
        selected = sorted({next(v['index'] for v in unique if v['owner'] == label) for label in {v['owner'] for v in unique}}
                          | {next(v['index'] for v in unique if len(v['remainingLimits']) == 2)})
    f.save(base/'preflight.json', {'mode': args.mode, 'selected': selected, **manifest})
    packets, records, reused = {}, {}, []
    if args.resume_from:
        old = args.resume_from.resolve(); old.relative_to(f.STORE)
        previous = json.loads((old/'checks.json').read_text()); old_inputs = json.loads((old/'inputs.json').read_text())
        f.require(previous['mode'] == 'preflight' and old_inputs == manifest, 'exact preflight source/helper/input reuse')
        for path in sorted((old/'rows').iterdir()):
            index = int(path.name); packet, check = verify_row(path, unique[index])
            destination = base/'rows'/path.name; destination.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(path, destination)
            packets[index], records[index] = packet, check; reused.append(index)
    for index in selected:
        if index not in packets:
            records[index] = run_row(base, unique[index], states[unique[index]['owner']])
            packets[index], _ = verify_row(base/'rows'/str(index).zfill(3), unique[index])
        f.save(base/'row-inventory.json', records)
    artifacts_before = {str(p.relative_to(base)): artifact(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    if args.mode == 'construct':
        for label, addresses in catalogue.items():
            state = states[label]; rows, phases, sources = [], {}, set()
            for address in addresses:
                result = accepted['result'] if address['kind'] == 'accepted' else packets[address['index']]['result']
                row = copy.copy(result['ROWS'][address['index'] if address['kind'] == 'accepted' else 0]); row['INDEX'] = address['caseIndex']
                f.require(row['ORIGINAL'] == state['builder'].integrals[address['caseIndex']], 'full assembled case factor address')
                rows.append(row); phases.update(dict(result['PHASES'])); sources.update(v['SOURCE_INTEGRAL'] for v in row['FACTORS'])
            case_result = {'ROWS': rows, 'CUTOFFS': state['builder'].cutoffs,
                           'SOURCE_INTEGRALS': tuple(sorted(sources, key=sp.default_sort_key)),
                           'PHASES': tuple(sorted(phases.items(), key=lambda v: sp.default_sort_key(v[0])))}
            target = base/'cases'/label; target.mkdir(parents=True)
            f.atomic_pickle(target/'factorization.pickle', {'result': case_result, 'dimensionState': dict(vars(state['dimensions'])),
                        'case': label, 'sourceFiles': manifest['sourceFiles'], 'inputPackets': manifest['inputPackets'], 'addresses': addresses})
    artifacts_before = {str(p.relative_to(base)): artifact(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    metadata = output_replay(base, states, unique, catalogue, packets, manifest)
    for name, record in artifacts_before.items(): f.require(artifact(base/name) == record, 'pre/post saved factor packet')
    for name, sha in manifest['sourceFiles'].items(): f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == sha, 'current/frozen helper/source identity')
    for name, sha in manifest['inputPackets'].items(): f.require(f.digest(Path(name)) == sha, 'unchanged input packet')
    checks = {'mode': args.mode, 'runDirectory': str(base), **manifest, **metadata, 'computedRows': sorted(packets), 'reusedPreflightRows': reused,
              'workers': records, 'artifacts': {str(p.relative_to(base)): artifact(p) for p in base.rglob('*')
                    if p.suffix in ('.pickle', '.out') and 'source' not in p.relative_to(base).parts},
              'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'scope': 'Preflight instrument evidence only.' if args.mode == 'preflight' else manifest['scope']}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__': main()
