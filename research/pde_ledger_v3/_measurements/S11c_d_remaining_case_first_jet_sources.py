#!/usr/bin/env python3
"""Actual case first-w-derivative controls, reusing completed source/grade work."""
import argparse
import ast
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import time

import sympy as sp
import S11c_d_remaining_case_bindings as native
import S11c_d_remaining_case_matrices as matrices
import S11c_d_remaining_case_modes as modes
import S11c_d_remaining_case_profile_bindings as provenance
import S11c_d_first_jet as first

f, engine, grades = native.f, native.engine, native.grades
BASELINE = native.BASELINE
PLAN = f.M/'S11c_d_remaining_case_first_jet_sources_plan.md'
BCP = f.M/'S11c_d_remaining_case_bindings_checkpoint.json'
JCP = f.M/'S11c_d_first_jet_checkpoint.json'
SCOPE = ('Three missing one-sided literal first-w-derivative source controls. '
         'The baseline control and exact expression/unit grade proofs are reused. '
         'No new consistent profile, isolated advection channel, numerical binding, '
         'quadrature, end/current/mode construction or response is computed here.')


def body(function):
    return hashlib.sha256(ast.dump(ast.parse(inspect.getsource(function))).encode()).hexdigest()


def require_record(item, address, value, unit):
    f.require(native.same(item['address'], address)
              and native.same(item['record']['ORIGINAL'], value)
              and tuple(item['record']['UNIT']) == tuple(unit),
              'actual source record, physical address and inherited unit')


def load(base, resume=None):
    origin, checks, _ = provenance.accepted(BCP, 'ACCEPTED_BINDINGS_AND_GRADES')
    jr, jc, publication = provenance.accepted(JCP, 'PUBLISHED_ANNEX_VERIFIED')
    published = f.ROOT/publication['publication']['path']
    f.require(published.is_symlink() and f.digest(published) == publication['publication']['sha256'],
              'accepted baseline derivative-control publication')
    pins = dict(checks['sourceFiles'])
    for name, digest in jc['sourceFiles'].items():
        f.require(name not in pins or pins[name] == digest, 'same consumed native sources')
        pins[name] = digest
    for path in (Path(__file__).resolve(), PLAN, BCP, JCP,
                 Path(matrices.__file__).resolve(), Path(modes.__file__).resolve(),
                 Path(provenance.__file__).resolve(), Path(first.__file__).resolve(),
                 Path(first.source.__file__).resolve(), Path(grades.__file__).resolve()):
        name = str(path.relative_to(f.ROOT))
        f.require(name not in pins or pins[name] == f.digest(path), 'unchanged helper join')
        pins[name] = f.digest(path)
    settings = f.unpickle(origin/'accepted-finite-system.pickle')['settings']
    manifest = {'runDirectory': str(base), 'sourceFiles': pins, 'inputPackets': {},
                'copiedInputs': {}, 'input': checks['input'], 'settings': settings,
                'nativeBodies': {name: body(fn) for name, fn in (
                    ('literalMutation', first.source.first_jet_mutation),
                    ('gradeSplit', grades.split), ('gradeCheck', grades.check),
                    ('gradeSpecs', grades.specs), ('termJoins', grades.term_joins),
                    ('recordJoin', require_record))}, 'scope': SCOPE}
    for root, result in ((origin, checks), (jr, jc)):
        manifest['inputPackets'][str(root/'checks.json')] = f.digest(root/'checks.json')
        for name, digest in result['inputPackets'].items():
            f.require(name not in manifest['inputPackets'] or manifest['inputPackets'][name] == digest,
                      'shared original source inputs')
            manifest['inputPackets'][name] = digest
    if resume is None:
        for label in checks['cases']:
            modes.retain(origin/'cases'/label/'case-binding.pickle',
                         base/'accepted-bindings'/label/'case-binding.pickle', manifest)
            for name in ('reduced-action', 'actions', 'assembly'):
                modes.retain(origin/'accepted-cases'/label/(name+'.pickle'),
                             base/'accepted-cases'/label/(name+'.pickle'), manifest)
            modes.retain(origin/'cases'/label/'factorization.pickle',
                         base/'accepted-cases'/label/'factorization.pickle', manifest)
        for name, item in jc['artifacts'].items():
            modes.retain(jr/name, base/'accepted-first-jet'/name, manifest, item['sha256'])
        modes.retain(jr/'inputs.json', base/'accepted-first-jet/inputs.json', manifest)
    else:
        old = json.loads((resume/'checks.json').read_text())
        f.require(old['status'] == 'VALIDATED_CASE_FIRST_JET_SOURCE_INPUTS', 'completed input preflight')
        f.require(old['sourceFiles'] == pins and old['nativeBodies'] == manifest['nativeBodies'],
                  'same complete source and validator bodies before input reuse')
        for name, digest in old['sourceFiles'].items():
            f.require(f.digest(resume/'source'/name) == digest, 'original focused source snapshot')
        for name, digest in old['inputPackets'].items():
            f.require(f.digest(Path(name)) == digest, 'original focused input')
            manifest['inputPackets'][name] = digest
        for name, item in old['artifacts'].items():
            modes.retain(resume/name, base/name, manifest, item['sha256'])
        modes.retain(resume/'inputs.json', base/'preflight-inputs.json', manifest)
        modes.retain(resume/'checks.json', base/'preflight-checks.json', manifest)
        manifest['completedInputReuse'] = {'runDirectory': str(resume),
            'checksSha256': f.digest(resume/'checks.json'),
            'validatorAstSha256': body(preflight), 'artifacts': old['artifacts']}
        f.require(old['preflightBodySha256'] == body(preflight), 'unchanged whole preflight validator')
    for name, digest in pins.items():
        dst = base/'source'/name; dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(f.ROOT/name, dst)
        f.require(f.digest(dst) == digest, 'frozen current source identity')
    f.save(base/'inputs.json', manifest)
    return manifest, tuple(checks['cases'])


def rejects(function):
    try:
        function()
    except (AssertionError, ValueError):
        return True
    return False


def preflight(base, labels, manifest):
    baseline = f.unpickle(base/'accepted-first-jet/first-jet-binding.pickle')
    summaries = {}; controls = []
    for label in labels:
        case = f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle')
        r, adapter, packets = matrices.context(base, label, case, manifest['input'])
        record = case['grades']; binding = case['binding']
        f.require(native.same(record['generators'], engine.PHYSICAL_METADATA.generators),
                  'actual independent material and shape generators')
        f.require(binding['settings'] == manifest['settings']
                  and len(binding['fieldUnits']) == len(binding['equationUnits']) == 5,
                  'approved finite settings and complete field/equation units')
        chosen = list(grades.specs(packets['assembly']['result'], packets['factorization']['result'],
                                  binding['fieldUnits'], binding['equationUnits'], r))
        pairs = [(key, address, value, unit, record['records'][key])
                 for key, value, unit, address in chosen]
        path = base/'input-joins'/label/'source-record-pairs.pickle'
        path.parent.mkdir(parents=True, exist_ok=True); f.atomic_pickle(path, pairs)
        f.require(len(chosen) == len(record['records']) and len({v[0] for v in chosen}) == len(chosen),
                  'complete actual native record census')
        for _, address, value, unit, item in pairs:
            require_record(item, address, value, unit)
        original = packets['factorization']['result']
        addresses = native.integral_addresses(packets['assembly']['result'], original)
        f.require(len(record['termJoins']) == sum(len(c['NONLOCAL']) for c in packets['assembly']['result']['ROWS']),
                  'all actual native nonlocal terms')
        by_index = {row['INDEX']: row for row in original['ROWS']}
        for term in record['termJoins']:
            row = by_index[term['integralIndex']]
            f.require(native.same((term['originalIntegral'], term['sourceLimit'], term['remainingLimits']),
                                  (row['ORIGINAL'], row['SOURCE_LIMIT'], row['REMAINING_LIMITS'])),
                      'complete integral and internally ordered source/momentum limits')
        if label == BASELINE:
            f.require(set(record['records']) == set(baseline['originalRecords']), 'all accepted baseline records')
            for key, item in record['records'].items():
                old = baseline['originalRecords'][key]
                require_record(old, item['address'], item['record']['ORIGINAL'], item['record']['UNIT'])
                f.require(native.same(baseline['records'][key]['address'], item['address'])
                          and baseline['records'][key]['record']['UNIT'] == item['record']['UNIT'],
                          'accepted baseline mutated record address and unit')
            f.require(native.same(record['termJoins'], baseline['termJoins'])
                      and baseline['settings'] == manifest['settings'], 'whole baseline native term and settings join')
        key, value, unit, address = next(v for v in chosen if v[1] != 0)
        item = record['records'][key]
        changed_address = (address[0], 999, *address[2:])
        changed_unit = (unit[0]+1, *unit[1:])
        tests = {'changedCoefficient': rejects(lambda: require_record(item, address, 2*value, unit)),
                 'changedPhysicalAddress': rejects(lambda: require_record(item, changed_address, value, unit)),
                 'changedUnit': rejects(lambda: require_record(item, address, value, changed_unit))}
        f.require(all(tests.values()), 'actual record coefficient/address/unit controls reject')
        source_index = next(i for i in range(len(original['SOURCE_INTEGRALS'])))
        source = binding['bound']['sources'][0, source_index]
        f.require(native.same(record['records'][f'sourceAmplitude{source_index}']['record']['ORIGINAL'],
                              source['symbolicAmplitude']), 'actual original field source before reversal')
        controls.append({'case': label, 'key': key, 'address': address, 'unit': unit, 'tests': tests})
        f.atomic_pickle(base/'input-joins'/label/'control-operands.pickle',
                        {'item': item, 'addressMutation': changed_address, 'unitMutation': changed_unit,
                         'coefficientMutation': 2*value, 'controls': tests})
        summaries[label] = {'records': len(chosen), 'rows': len(original['ROWS']),
                            'sources': len(original['SOURCE_INTEGRALS']), 'terms': len(record['termJoins']),
                            'wholeIntegralAddresses': len(addresses), 'controls': tests,
                            'fieldUnits': [list(map(str, v)) for v in binding['fieldUnits']],
                            'equationUnits': [list(map(str, v)) for v in binding['equationUnits']]}
    f.save(base/'preflight.json', {'cases': summaries, 'preflightBodySha256': body(preflight),
                                  'sourceFiles': manifest['sourceFiles'], 'scope': SCOPE})
    return summaries


def endpoint_operands(base, r, adapter):
    derivative = sp.diff(adapter.input.profiles['w'], r.xi)
    endpoints = []
    for direction in (-sp.oo, sp.oo):
        value = {'direction': direction, 'profile': adapter.input.profiles['w'],
                 'derivative': derivative, 'reversedDerivative': -derivative,
                 'profileLimit': sp.limit(adapter.input.profiles['w'], r.xi, direction),
                 'derivativeLimit': sp.limit(derivative, r.xi, direction),
                 'reversedDerivativeLimit': sp.limit(-derivative, r.xi, direction)}
        endpoints.append(value)
    f.atomic_pickle(base/'endpoint-operands.pickle', endpoints)
    f.require(all(v['derivativeLimit'] == v['reversedDerivativeLimit'] == 0 for v in endpoints),
              'computed unchanged first-derivative endpoints; prerequisite to future end-map reuse')
    return endpoints


def cached_grade(cache, value, unit):
    return next(((origin, record) for origin, record in cache.get((hash(value), tuple(unit)), ())
                 if native.same(record['ORIGINAL'], value) and tuple(record['UNIT']) == tuple(unit)), None)


def transformed_terms(base, packets, original, records):
    assembly, fourier = packets['assembly']['result'], packets['factorization']['result']
    addresses = native.integral_addresses(assembly, fourier)
    by_index = {v['INDEX']: v for v in fourier['ROWS']}
    by_address = {tuple(v['address']): v['record'] for v in records.values()}
    cells = [dict(cell, NONLOCAL=tuple((by_index[addresses[id(integral)]]['ORIGINAL'],
                 by_address['cell', cell['ROW'], cell['COLUMN'], ti]['ORIGINAL'])
                 for ti, (integral, _) in enumerate(cell['NONLOCAL']))) for cell in assembly['ROWS']]
    # Original whole integrals are immutable addresses of the closed operator.
    # Its actually reversed coefficient operands are stored separately and are
    # consumed by the original grade convolution, not claimed as new factorization.
    selected = dict(assembly, ROWS=cells)
    f.atomic_pickle(base/'term-inputs.pickle', {'originalCells': assembly['ROWS'], 'selectedCells': cells,
                    'originalFourier': fourier, 'addressRole': 'Original source identifiers; reversed factor/source operands are the selected grade records.'})
    terms = grades.term_joins(selected, fourier, records)
    f.atomic_pickle(base/'term-joins.pickle', terms)
    f.require(len(terms) == len(original), 'complete untruncated native term convolution')
    changed = 0
    for old, new in zip(original, terms):
        f.require(native.same({k: v for k, v in old.items() if k not in ('cellCoefficient', 'factors')},
                              {k: v for k, v in new.items() if k not in ('cellCoefficient', 'factors')}),
                  'all original physical term addresses and ordered limits retained')
        f.require(len(old['factors']) == len(new['factors']), 'all native factor branches retained')
        for a, b in zip(old['factors'], new['factors']):
            f.require(native.same({k: v for k, v in a.items() if k != 'combinations'},
                                  {k: v for k, v in b.items() if k != 'combinations'}),
                      'same actual source field, Fourier character and frequency')
        changed += int(not native.same(old, new))
    return terms, changed


def construct(base, labels, manifest, preflight_summary):
    baseline = f.unpickle(base/'accepted-first-jet/first-jet-binding.pickle')
    cases = {label: f.unpickle(base/'accepted-bindings'/label/'case-binding.pickle') for label in labels}
    generators = cases[BASELINE]['grades']['generators']; cache = {}
    for label, case in cases.items():
        f.require(native.same(case['grades']['generators'], generators), 'identical independent generators')
        for key, item in case['grades']['records'].items():
            record = item['record']; origin = {'packet': str(base/'accepted-bindings'/label/'case-binding.pickle'),
                                               'kind': 'accepted-original', 'case': label, 'key': key}
            cache.setdefault((hash(record['ORIGINAL']), tuple(record['UNIT'])), []).append((origin, record))
    for key, item in baseline['records'].items():
        record = item['record']; origin = {'packet': str(base/'accepted-first-jet/first-jet-binding.pickle'),
                                          'kind': 'accepted-first-jet', 'case': BASELINE, 'key': key}
        cache.setdefault((hash(record['ORIGINAL']), tuple(record['UNIT'])), []).append((origin, record))
    inventory = {BASELINE: dict(preflight_summary[BASELINE], reusedWholeControl=True,
                               changedRecords=baseline['changedRecordCounts'], newGradeProofs=0)}
    paths = {BASELINE: str(base/'accepted-first-jet/first-jet-binding.pickle')}
    for label in labels:
        if label == BASELINE:
            continue
        target = base/'cases'/label; target.mkdir(parents=True, exist_ok=False)
        case = cases[label]; r, adapter, packets = matrices.context(base, label, case, manifest['input'])
        f.require(native.same(engine.PHYSICAL_METADATA.generators, generators), 'current restored generator frame')
        records = {}; index = {}; counts = dict.fromkeys(('local', 'cell', 'factor', 'source'), 0)
        reused = {}; new_proofs = 0; coefficient_controls = []; mutation_controls = []
        for key, item in case['grades']['records'].items():
            original = item['record']; address = item['address']; unit = original['UNIT']
            mutation = first.source.first_jet_mutation(original['ORIGINAL'], r)
            raw = {'address': address, 'original': original['ORIGINAL'], 'unit': unit, 'mutation': mutation}
            raw_path = target/'raw-records'/(key+'.pickle'); raw_path.parent.mkdir(exist_ok=True)
            f.atomic_pickle(raw_path, raw)
            selected = mutation['mutated']; changed = not native.same(selected, original['ORIGINAL'])
            restored = first.source.first_jet_mutation(selected, r)['mutated'] if mutation['distinctAtomCount'] else selected
            proof = {'address': address, 'restored': restored, 'involutionResidual': restored-original['ORIGINAL'],
                     'omittedReversalResidual': selected-original['ORIGINAL'], 'changed': changed,
                     'occurrenceUnits': [(atom, engine.PHYSICAL_METADATA.dimensions.measure(atom),
                                          engine.PHYSICAL_METADATA.dimensions.measure(-atom)) for atom in mutation['occurrences']]}
            proof_path = target/'mutation-proofs'/(key+'.pickle'); proof_path.parent.mkdir(exist_ok=True)
            f.atomic_pickle(proof_path, proof)
            f.require(native.same(restored, original['ORIGINAL']) and proof['involutionResidual'] == 0,
                      'actual literal first-derivative reversal involution')
            for atom, before, after in proof['occurrenceUnits']:
                f.require(isinstance(atom, sp.Subs) and isinstance(atom.expr, sp.Derivative)
                          and atom.expr.expr.func == r.profiles['w']
                          and sum(n for _, n in atom.expr.variable_count) == 1 and before == after,
                          'only literal first w derivative changes, with unchanged unit')
            if changed:
                f.require(proof['omittedReversalResidual'] != 0, 'actual omitted reversal responds')
                mutation_controls.append((key, proof['omittedReversalResidual']))
            counts[address[0]] += int(changed)
            match = cached_grade(cache, selected, unit)
            if match:
                origin, record = match; reused[key] = origin
            else:
                record = grades.split(selected, generators, unit); new_proofs += 1
            out = {'address': address, 'record': record}; records[key] = out
            path = target/'records'/(key+'.pickle'); path.parent.mkdir(exist_ok=True)
            f.atomic_pickle(path, out)
            index[key] = {'address': address, 'changed': changed, 'path': str(path.relative_to(base)),
                          'sha256': f.digest(path), 'rawSha256': f.digest(raw_path),
                          'mutationProofSha256': f.digest(proof_path), 'reusedGradeFrom': reused.get(key)}
            f.save(target/'record-inventory.json', index)
            if not match:
                grades.check(record)
                cache.setdefault((hash(selected), tuple(unit)), []).append(
                    ({'packet': str(path), 'kind': 'completed-new', 'case': label, 'key': key}, record))
            f.require(native.same(record['ORIGINAL'], selected) and tuple(record['UNIT']) == tuple(unit),
                      'actual selected expression and unit before completed grade-proof reuse')
            f.require(not record['OUTSIDE_RECTANGLE'], 'actual selected independent retained rectangle')
            if record['COEFFICIENTS']:
                grade, coefficient = next(iter(record['COEFFICIENTS'].items()))
                control = coefficient*sp.prod(v**p for v, p in zip(generators, grade))
                f.require(control != 0, 'actual coefficient omission responds')
                coefficient_controls.append((key, grade, control))
            with (base/'progress.jsonl').open('a') as stream:
                stream.write(json.dumps({'case': label, 'record': key, 'completedRecords': len(records),
                                         'newGradeProofs': new_proofs})+'\n')
        f.atomic_pickle(target/'coefficient-controls.pickle', coefficient_controls)
        f.atomic_pickle(target/'omitted-reversal-controls.pickle', mutation_controls)
        terms, changed_terms = transformed_terms(target, packets, case['grades']['termJoins'], records)
        endpoints = endpoint_operands(target, r, adapter)
        result = {'records': records, 'originalRecords': case['grades']['records'], 'termJoins': terms,
                  'originalTermJoins': case['grades']['termJoins'], 'generators': generators,
                  'fieldUnits': case['binding']['fieldUnits'], 'equationUnits': case['binding']['equationUnits'],
                  'dimensionState': dict(vars(engine.PHYSICAL_METADATA.dimensions)),
                  'changedRecordCounts': counts, 'reusedGradeRecords': reused, 'newGradeProofs': new_proofs,
                  'endpointOperands': endpoints, 'coefficientControls': coefficient_controls,
                  'omittedReversalControls': mutation_controls, 'sourceFiles': manifest['sourceFiles'],
                  'inputPackets': manifest['inputPackets'], 'scope': SCOPE}
        f.atomic_pickle(target/'first-jet-sources.pickle', result); paths[label] = str(target/'first-jet-sources.pickle')
        inventory[label] = dict(preflight_summary[label], changedRecords=counts, changedTermOperands=changed_terms,
                                newGradeProofs=new_proofs, reusedGradeRecords=len(reused),
                                coefficientControls=len(coefficient_controls), omittedReversalControls=len(mutation_controls),
                                gradeSupport=sorted({str(g) for item in records.values() for g in item['record']['COEFFICIENTS']}))
        f.save(base/'case-inventory.json', inventory)
    f.atomic_pickle(base/'remaining-case-first-jet-sources.pickle',
                    {'cases': paths, 'inventory': inventory, 'sourceFiles': manifest['sourceFiles'],
                     'inputPackets': manifest['inputPackets'], 'scope': SCOPE})
    return inventory, paths


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--run-directory', type=Path, required=True)
    ap.add_argument('--preflight', action='store_true'); ap.add_argument('--resume-from', type=Path)
    args = ap.parse_args(); started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    def timeout(*_):
        raise TimeoutError('case first-derivative source budget; preserve completed operands')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(900)
    base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True, exist_ok=False)
    manifest, labels = load(base, args.resume_from)
    summary = (json.loads((base/'preflight.json').read_text())['cases'] if args.resume_from
               else preflight(base, labels, manifest))
    if args.preflight:
        cases, paths = summary, {}
    else:
        cases, paths = construct(base, labels, manifest, summary)
    for name, digest in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == digest, 'source/frozen pre/post identity')
    for name, digest in manifest['inputPackets'].items():
        f.require(f.digest(Path(name)) == digest, 'original input pre/post identity')
    for name, digest in manifest['copiedInputs'].items():
        f.require(f.digest(base/name) == digest, 'unchanged copied completed operands')
    artifacts = {str(p.relative_to(base)): {'sha256': f.digest(p), 'bytes': p.stat().st_size}
                 for p in base.rglob('*') if p.is_file() and 'source' not in p.relative_to(base).parts
                 and p not in (base/'inputs.json', base/'checks.json')}
    checks = {**manifest, 'status': ('VALIDATED_CASE_FIRST_JET_SOURCE_INPUTS' if args.preflight
                                   else 'COMPLETED_CASE_FIRST_JET_SOURCE_CONTROLS'),
              'cases': cases, 'casePackets': paths, 'preflightBodySha256': body(preflight),
              'newNumericalBindings': 0, 'newQuadratureNodes': 0, 'newModes': 0, 'newResponseSolves': 0,
              'artifacts': artifacts, 'wallSeconds': time.monotonic()-started,
              'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
