#!/usr/bin/env python3
"""Finish saved factor proofs with explicit Piecewise reachability records."""
import argparse
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import resource
import shutil
import signal
import textwrap
import time
import sympy as sp
import S11c_d_remaining_case_factors as c

f, engine, grades = c.f, c.engine, c.grades
ORIGIN = f.STORE/'s11c-remaining-case-factors-20260919/preflight/complete'
PLAN = f.M/'S11c_d_remaining_case_factor_recovery_plan.md'
BRANCH_RECORDS = []
FRACTION_RECORDS = []


def reachable_branches(branches):
    prior = sp.false; retained = []
    for value, condition in branches:
        effective = sp.And(condition, sp.Not(prior))
        simplified = sp.simplify_logic(effective, force=True)
        record = {'value': value, 'condition': condition, 'priorCoverage': prior,
                  'effectiveCondition': effective, 'simplifiedEffectiveCondition': simplified,
                  'shadowed': simplified == sp.false}
        BRANCH_RECORDS.append(record)
        if simplified != sp.false:
            retained.append((value, condition))
        prior = sp.Or(prior, condition)
    f.require(retained, 'nonempty reachable Piecewise branch set')
    return tuple(retained)


def rational_normalize(value):
    # Rational arithmetic on the same native encoded carriers; no multivariate
    # GCD is necessary to prove that an expanded numerator is zero.
    combined = sp.together(value)
    numerator, denominator = sp.fraction(combined)
    expanded = sp.expand(numerator)
    f.require(denominator != 0, 'nonzero formal denominator polynomial')
    result = sp.S.Zero if expanded == 0 else combined
    FRACTION_RECORDS.append({'value': value, 'combined': combined, 'numerator': numerator,
                            'denominator': denominator, 'expandedNumerator': expanded, 'result': result})
    return result


def certificate_constructor():
    original = ast.parse(textwrap.dedent(inspect.getsource(engine.BoundedSourceFourierAssembly.reconstruction_certificate))).body[0]
    changed = copy.deepcopy(original); counts = {'branches': 0, 'normalizer': 0}
    class Forward(ast.NodeTransformer):
        def visit_Assign(self, node):
            self.generic_visit(node)
            if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id == 'branches':
                counts['branches'] += 1
                node.value = ast.Call(func=ast.Name(id='reachable_branches', ctx=ast.Load()), args=[node.value], keywords=[])
            return node
        def visit_Call(self, node):
            self.generic_visit(node)
            if ast.unparse(node.func) == 'sp.cancel':
                counts['normalizer'] += 1; node.func = ast.Name(id='rational_normalize', ctx=ast.Load())
            return node
    Forward().visit(changed)
    class Reverse(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == 'reachable_branches': return node.args[0]
            if isinstance(node.func, ast.Name) and node.func.id == 'rational_normalize':
                node.func = ast.parse('sp.cancel', mode='eval').body
            return node
    f.require(counts == {'branches': 1, 'normalizer': 1} and ast.dump(Reverse().visit(copy.deepcopy(changed))) == ast.dump(original),
              'whole native certificate reverse AST join')
    joined = ast.dump(changed); changed.decorator_list = []
    namespace = dict(vars(engine), reachable_branches=reachable_branches, rational_normalize=rational_normalize)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[changed], type_ignores=[])), str(Path(__file__)), 'exec'), namespace)
    return namespace['reconstruction_certificate'], {'wholeCertificateReverseAstJoin': True, 'changes': counts,
        'originalAstSha256': hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'derivedAstSha256': hashlib.sha256(joined.encode()).hexdigest()}


def certificate(left, right, *, shared=False):
    f.require(shared is False, 'uncompressed actual certificate operands')
    BRANCH_RECORDS.clear(); FRACTION_RECORDS.clear()
    constructor, join = certificate_constructor(); result = constructor(left, right, shared=False)
    result['REACHABILITY_RECORDS'] = tuple(BRANCH_RECORDS)
    result['FRACTION_RECORDS'] = tuple(FRACTION_RECORDS)
    result['NORMALIZER_JOIN'] = join
    for item in result['REACHABILITY_RECORDS']:
        if item['shadowed']:
            f.require(sp.simplify_logic(sp.And(item['condition'], sp.Not(item['priorCoverage'])), force=True) == sp.false,
                      'every omitted branch has impossible effective condition')
    return result


def resumed_factorization(builder, checkpoint, directory, entry):
    path = directory/'retained-raw-factorization.pickle'
    if not path.exists(): return builder.construct(checkpoint=checkpoint)
    packet = f.unpickle(path); row = packet['row']
    f.require(row['ORIGINAL'] == entry['original'] and row['BOUNDED'] == entry['bounded'] and row['INDEX'] == 0,
              'actual saved native factor row and source address')
    builder.phases = packet['phases']
    engine.PHYSICAL_METADATA.dimensions.__dict__.update(packet['dimensionState'])
    # Native completed-row continuation traverses no separation or amplitude algebra.
    return builder.construct(completed=(row,))


def worker_adapter():
    original = ast.parse(inspect.getsource(c.worker)).body[0]; changed = copy.deepcopy(original)
    counts = {'construct': 0, 'certificate': 0}
    class Forward(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if ast.unparse(node.func) == 'builder.construct':
                counts['construct'] += 1
                return ast.parse('resumed_factorization(builder, checkpoint, directory, entry)', mode='eval').body
            if ast.unparse(node.func) == 'engine.BoundedSourceFourierAssembly.reconstruction_certificate':
                counts['certificate'] += 1; node.func = ast.Name(id='build_exact_certificate', ctx=ast.Load())
            return node
    Forward().visit(changed)
    class Reverse(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == 'resumed_factorization':
                return ast.parse('builder.construct(checkpoint=checkpoint)', mode='eval').body
            if isinstance(node.func, ast.Name) and node.func.id == 'build_exact_certificate':
                node.func = ast.parse('engine.BoundedSourceFourierAssembly.reconstruction_certificate', mode='eval').body
            return node
    f.require(counts == {'construct': 1, 'certificate': 2} and ast.dump(Reverse().visit(copy.deepcopy(changed))) == ast.dump(original),
              'whole worker reverse AST join')
    namespace = dict(vars(c), build_exact_certificate=certificate, resumed_factorization=resumed_factorization)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[changed], type_ignores=[])), str(Path(__file__)), 'exec'), namespace)
    return namespace['worker'], {'wholeWorkerReverseAstJoin': True, 'changes': counts,
        'originalAstSha256': hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'derivedAstSha256': hashlib.sha256(ast.dump(changed).encode()).hexdigest()}


def load_saved(base):
    manifest = json.loads((ORIGIN/'inputs.json').read_text()); outcome = json.loads((ORIGIN.parent/'factors_preflight.invocation.json').read_text())
    f.require(outcome == json.loads((ORIGIN.parent/'active.json').read_text()) and outcome['exitCode'] == 1, 'actual original preflight failure')
    for name, sha in manifest['sourceFiles'].items(): f.require(f.digest(f.ROOT/name) == f.digest(ORIGIN/'source'/name) == sha, 'unchanged original helper/source')
    for name, sha in manifest['inputPackets'].items(): f.require(f.digest(Path(name)) == sha, 'unchanged original input')
    frozen = {str(p.relative_to(ORIGIN)): f.digest(p) for p in ORIGIN.rglob('*') if p.is_file() and 'source' not in p.relative_to(ORIGIN).parts}
    for name in ('accepted-cases', 'source'):
        shutil.copytree(ORIGIN/name, base/name)
    for name in ('accepted-factorization.pickle', 'accepted-factor-certificates.pickle', 'integral-catalogue.pickle'):
        shutil.copyfile(ORIGIN/name, base/name)
    manifest = copy.deepcopy(manifest)
    manifest['inputPackets'].update({str(ORIGIN/name): sha for name, sha in frozen.items()})
    manifest['inputPackets'].update({str(p): f.digest(p) for p in (ORIGIN.parent/'factors_preflight.invocation.json', ORIGIN.parent/'factors_preflight.stderr')})
    for path in (Path(__file__).resolve(), PLAN):
        name = str(path.relative_to(f.ROOT)); manifest['sourceFiles'][name] = f.digest(path)
        dest = base/'source'/name; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(path, dest)
    manifest['originalOutcome'] = outcome; manifest['originalArtifacts'] = frozen
    catalogue = f.unpickle(base/'integral-catalogue.pickle'); states = {}
    for label in catalogue['cases']:
        packets = {kind: f.unpickle(base/'accepted-cases'/label/(kind+'.pickle')) for kind in ('reduced-action', 'actions', 'assembly')}
        r, dimensions, pencil = c.context(packets['reduced-action'], packets['actions'], packets['assembly'])
        builder = engine.BoundedSourceFourierAssembly(r, packets['assembly']['result']['NONLOCAL_INTEGRALS'])
        states[label] = {'packets': packets, 'r': r, 'dimensions': dimensions, 'pencil': pencil, 'builder': builder}
    # Validate saved addresses and retained full limits without re-factorization.
    accepted = f.unpickle(base/'accepted-factorization.pickle')
    for label, addresses in catalogue['cases'].items():
        for address in addresses:
            actual = states[label]['builder'].integrals[address['caseIndex']]
            expected = accepted['result']['ROWS'][address['index']]['ORIGINAL'] if address['kind'] == 'accepted' else catalogue['unique'][address['index']]['original']
            f.require(actual == expected and actual.limits == expected.limits, 'saved catalogue source/ordered-limit address')
    f.save(base/'inputs.json', manifest)
    return states, accepted, catalogue, manifest


def focused(base):
    _, join = certificate_constructor(); _, wjoin = worker_adapter(); results = []
    for index in (0, 3, 9, 28):
        folder = base/str(index).zfill(3); folder.mkdir()
        origin = ORIGIN/'rows'/str(index).zfill(3)
        p = f.unpickle(origin/'integrand-pair.pickle')
        f.require(tuple(hashlib.sha256(s.encode()).hexdigest() for s in p['representations']) == p['representationHashes'], 'original live representation hashes')
        before = f.digest(origin/'integrand-pair.pickle')
        result = certificate(p['left'], p['right']); f.atomic_pickle(folder/'certificate.pickle', result)
        f.require(result['RESIDUAL'] == 0 and all(v == 0 for v in c.proof_scalars(result)), 'actual saved reconstruction proof')
        raw = f.unpickle(origin/'raw-factorization.pickle')['row']; first = raw['FACTORS'][0]
        mutated = certificate(p['left'], p['right']+first['COEFFICIENT']*first['SOURCE'])
        f.atomic_pickle(folder/'mutation.pickle', mutated)
        f.require(mutated['RESIDUAL'] != 0 and all(v == 0 for v in c.proof_scalars(mutated)), 'actual changed coefficient response')
        f.require(f.digest(origin/'integrand-pair.pickle') == before, 'unchanged saved pair')
        results.append({'index': index, 'proofScalars': len(c.proof_scalars(result)),
                        'shadowedBranches': sum(v['shadowed'] for v in result['REACHABILITY_RECORDS']),
                        'mutationNonzero': mutated['RESIDUAL'] != 0, 'pairSha256': before})
        if index == 28:
            branches = result['REACHABILITY_RECORDS']; dropped = next(v for v in branches if v['shadowed'] and v['value'] != 0)
            # The previously shadowed nonzero branch becomes reachable when the
            # actual earlier guard is disabled. It must then remain and respond.
            original_branches = tuple((v['value'], v['condition']) for v in branches)
            changed_branches = ((original_branches[0][0], sp.false), *original_branches[1:])
            BRANCH_RECORDS.clear(); kept = reachable_branches(changed_branches)
            restored = sp.Piecewise(*((rational_normalize(v), cond) for v, cond in kept))
            f.atomic_pickle(folder/'guard-mutation.pickle', {'original': original_branches, 'changed': changed_branches,
                             'retained': kept, 'residual': restored, 'records': tuple(BRANCH_RECORDS)})
            f.require(any(value == dropped['value'] for value, _ in kept) and restored != 0, 'changed earlier guard exposes nonzero branch')
    return {'certificateJoin': join, 'workerJoin': wjoin, 'cases': results, 'guardMutationRejected': True}


def emit_domain_evidence(base, packets):
    for index, packet in sorted(packets.items()):
        for ordinal, record in enumerate(packet['certificates']):
            cert = record['certificate']
            records = [{key: (sp.srepr(value) if isinstance(value, sp.Basic) else value)
                        for key, value in branch.items() if key != 'value'}
                       for branch in cert.get('REACHABILITY_RECORDS', ())]
            # Full raw encoded branches/fractions and carrier definitions remain
            # in the hashed certificate packets. Original physical pairs and
            # normalized residuals use the unchanged native physical emitter.
            path = base/'rows'/str(index).zfill(3)/(record['kind'].lower()+'-certificate.pickle')
            c.cases.boundary.structural_flags(c.PREFIX+'_NEW_'+str(index)+'_DOMAIN_PROOF_'+str(ordinal),
                {'method': 'first-match conditions and rational numerator' if 'NORMALIZER_JOIN' in cert else 'accepted native certificate',
                 'records': records, 'certificateSha256': f.digest(path),
                 'fractionCount': len(cert.get('FRACTION_RECORDS', ())),
                 'originalDenominatorCount': len(cert['ORIGINAL_DENOMINATOR_BASES'])})


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--mode', choices=('focused', 'preflight', 'construct'), required=True)
    parser.add_argument('--run-directory', type=Path, required=True); parser.add_argument('--resume-from', type=Path)
    args = parser.parse_args(); base = args.run_directory.resolve(); base.relative_to(f.STORE); base.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3)); started = time.monotonic()
    def timeout(*_): raise TimeoutError('saved factor validation budget; preserve completed operands')
    signal.signal(signal.SIGALRM, timeout); signal.alarm(900)
    if args.mode == 'focused':
        result = {**focused(base), 'wallSeconds': time.monotonic()-started}; f.save(base/'checks.json', result); print(json.dumps(result, indent=2)); return
    states, accepted, catalogue, manifest = load_saved(base); unique = catalogue['unique']; addresses = catalogue['cases']
    worker, worker_join = worker_adapter(); _, certificate_join = certificate_constructor(); c.worker = worker
    packets, records, reused = {}, {}, []
    selected = json.loads((ORIGIN/'preflight.json').read_text())['selected'] if args.mode == 'preflight' else list(range(len(unique)))
    previous = args.resume_from.resolve() if args.resume_from else ORIGIN
    if args.resume_from:
        previous.relative_to(f.STORE); previous_checks = json.loads((previous/'checks.json').read_text())
        f.require(previous_checks['mode'] == 'preflight', 'complete accepted preflight reuse')
        f.require(json.loads((previous/'inputs.json').read_text()) == manifest, 'same recovery source/input manifest')
    for index in selected:
        source_dir = previous/'rows'/str(index).zfill(3); dest = base/'rows'/str(index).zfill(3)
        if (source_dir/'checks.json').exists():
            packet, check = c.verify_row(source_dir, unique[index]); dest.parent.mkdir(parents=True, exist_ok=True); shutil.copytree(source_dir, dest)
            packets[index], records[index] = packet, check; reused.append(index)
        else:
            if index == 28 and not args.resume_from:
                # c.run_row creates the row directory. The adapter uses the
                # preserved original raw row via a private copy made before launch.
                raw = source_dir/'raw-factorization.pickle'
                f.require(raw.exists(), 'completed native row exists')
                original_run = c.run_row
                def run_with_raw(b, entry, state):
                    directory = b/'rows'/str(entry['index']).zfill(3)
                    # Reuse the unchanged run_row body; only mkdir copies saved raw data.
                    original_resume = resumed_factorization
                    def from_origin(builder, checkpoint, directory, entry):
                        shutil.copyfile(raw, directory/'retained-raw-factorization.pickle')
                        f.require(f.digest(raw) == f.digest(directory/'retained-raw-factorization.pickle'), 'byte-identical completed raw row')
                        return original_resume(builder, checkpoint, directory, entry)
                    namespace = worker.__globals__; saved = namespace['resumed_factorization']; namespace['resumed_factorization'] = from_origin
                    try: return original_run(b, entry, state)
                    finally: namespace['resumed_factorization'] = saved
                records[index] = run_with_raw(base, unique[index], states[unique[index]['owner']])
            else: records[index] = c.run_row(base, unique[index], states[unique[index]['owner']])
            packets[index], _ = c.verify_row(dest, unique[index])
        f.save(base/'row-inventory.json', records)
    if args.mode == 'construct':
        for label, entries in addresses.items():
            state = states[label]; rows, phases, sources = [], {}, set()
            for address in entries:
                result = accepted['result'] if address['kind'] == 'accepted' else packets[address['index']]['result']
                row = copy.copy(result['ROWS'][address['index'] if address['kind'] == 'accepted' else 0]); row['INDEX'] = address['caseIndex']
                f.require(row['ORIGINAL'] == state['builder'].integrals[address['caseIndex']], 'actual case row source identity')
                rows.append(row); phases.update(dict(result['PHASES'])); sources.update(v['SOURCE_INTEGRAL'] for v in row['FACTORS'])
            target = base/'cases'/label; target.mkdir(parents=True)
            f.atomic_pickle(target/'factorization.pickle', {'result': {'ROWS': rows, 'CUTOFFS': state['builder'].cutoffs,
                'SOURCE_INTEGRALS': tuple(sorted(sources, key=sp.default_sort_key)), 'PHASES': tuple(sorted(phases.items(), key=lambda v: sp.default_sort_key(v[0])))},
                'dimensionState': dict(vars(state['dimensions'])), 'case': label, 'sourceFiles': manifest['sourceFiles'], 'inputPackets': manifest['inputPackets'], 'addresses': entries})
    before = {str(p.relative_to(base)): c.artifact(p) for p in base.rglob('*.pickle') if 'source' not in p.relative_to(base).parts}
    previous_emit = c.emit
    def emit_with_domains(*args):
        previous_emit(*args)
        emit_domain_evidence(base, packets)
    c.emit = emit_with_domains
    try: metadata = c.output_replay(base, states, unique, addresses, packets, manifest)
    finally: c.emit = previous_emit
    for name, record in before.items(): f.require(c.artifact(base/name) == record, 'unchanged pre/post artifact')
    for name, sha in manifest['sourceFiles'].items(): f.require(f.digest(f.ROOT/name) == f.digest(base/'source'/name) == sha, 'unchanged current/frozen source')
    for name, sha in manifest['inputPackets'].items(): f.require(f.digest(Path(name)) == sha, 'unchanged original inputs and artifacts')
    checks = {'mode': args.mode, 'runDirectory': str(base), **manifest, **metadata, 'workerJoin': worker_join, 'certificateJoin': certificate_join,
              'computedRows': sorted(packets), 'reusedRows': reused, 'workers': records,
              'artifacts': {str(p.relative_to(base)): c.artifact(p) for p in base.rglob('*') if p.suffix in ('.pickle', '.out') and 'source' not in p.relative_to(base).parts},
              'wallSeconds': time.monotonic()-started, 'peakRssKiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    f.save(base/'checks.json', checks); signal.alarm(0); print(json.dumps(checks, indent=2))


if __name__ == '__main__': main()
