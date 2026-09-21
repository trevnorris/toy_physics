#!/usr/bin/env python3
"""Continue saved character proofs through the native raw-zero branch."""
import ast
import copy
import hashlib
import inspect
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import S11c_d_remaining_case_frequency_characters_recover as r

h, f, native = r.h, r.f, r.native
PREVIOUS = r.PREVIOUS.parent/'characters-recovery-01'
COVERAGE = PREVIOUS.parent/'characters-certificate-coverage'
PLAN = f.M/'S11c_d_remaining_case_frequency_characters_zero_recovery_plan.md'
REPAIR = f.M/'S11c_d_remaining_case_frequency_characters_zero_repair.json'
original_checker = r.check_source_proofs
STATE = {'writes': {}, 'completedProofReuses': [], 'zeroControls': set()}
WRITE_ROUTES = {}
WRITES = []


def save_pickle(path, value):
    path = Path(path)
    name = str(path.relative_to(STATE['base']))
    if name in WRITE_ROUTES:
        item = WRITE_ROUTES[name]
        f.require(path.is_symlink() and f.digest(path) == item['sha256'], 'completed input reference bytes')
        f.require(native.same(f.unpickle(path), value), 'full typed completed/requested input identity')
        entry = {'path': name, 'original': item['original'], 'sha256': item['sha256'], 'typedIdentity': True}
        STATE['writes'][name] = entry
        WRITES.append(entry)
        f.save(STATE['base']/'zero-recovery-input-write-joins.json', WRITES)
        return
    f.require(not path.exists() and not path.is_symlink(), ('new packet creates once', name))
    f.atomic_pickle(path, value)


def save_json(path, value):
    path = Path(path)
    name = str(path.relative_to(STATE['base']))
    if name in WRITE_ROUTES:
        item = WRITE_ROUTES[name]
        f.require(path.is_symlink() and f.digest(path) == item['sha256'], 'completed JSON reference bytes')
        f.require(json.loads(path.read_text()) == json.loads(json.dumps(value)), 'complete requested/completed JSON identity')
        WRITES.append({'path': name, 'original': item['original'], 'sha256': item['sha256'], 'jsonIdentity': True})
        f.save(STATE['base']/'zero-recovery-input-write-joins.json', WRITES)
        return
    f.require(not path.is_symlink(), ('no write through a preserved reference', name))
    f.save(path, value)


def zero_join(factor, source, owner, ri, fi, evidence, certificates):
    f.require(owner == evidence['owner'] and ri == evidence['row'] and fi == evidence['factorIndex'],
              'actual raw-zero physical owner and factor address')
    f.require(owner['kind'] == 'accepted' and (owner['index'], 'AMPLITUDE_'+str(fi)) not in certificates,
              'native accepted raw-zero branch has no certificate record')
    f.require(evidence['certificatePresent'] is False and factor['AMPLITUDE_RECONSTRUCTION_RESIDUAL'] == 0
              and evidence['originalFactor']['AMPLITUDE_RECONSTRUCTION_RESIDUAL'] == 0,
              'actual saved original and requested literal raw zeros')
    # The diagnostic already joined the entire original owner factor. This
    # binds that saved evidence to the present case operand, without replaying
    # the original reconstruction or its certificate.
    f.require(native.same(factor, evidence['factor']), 'full requested factor equals completed diagnostic operand')
    f.require(native.same((factor['AMPLITUDE'], factor['FREQUENCY']),
                          (source['symbolicAmplitude'], source['symbolicFrequency'])),
              'actual raw-zero source amplitude and frequency')
    f.require(factor['CHARACTER_EQUATION_RESIDUAL'] == factor['CHARACTER_NORMALIZATION_RESIDUAL'] == 0,
              'unchanged native raw-zero character equations')


def raw_zero_route(base, label, si, ri, fi, factor, physical, owner):
    evidence = STATE['omitted'][label].get((ri, fi))
    if evidence is None:
        return None
    location = base/'factor-proofs/accepted-factor-certificates.pickle'
    if 'baseline' not in r.PROOFS:
        r.PROOFS['baseline'] = f.unpickle(location)
    addresses = {(v['row'], v['kind']) for v in r.PROOFS['baseline']['records']}
    joined = {'case': label, 'sourceIndex': si, 'row': ri, 'factorIndex': fi,
              'factor': factor, 'source': physical['source'], 'owner': owner,
              'nativeBranch': 'raw == 0; residuals.append(raw); continue',
              'certificatePresent': False, 'completedDiagnosticPair': evidence,
              'diagnosticPairPath': str(COVERAGE/(label+'-omitted-certificate-pairs.pickle')),
              'diagnosticPairSha256': STATE['coverageFiles'][label+'-omitted-certificate-pairs.pickle']['sha256'],
              'diagnosticChecksSha256': STATE['coverageChecksSha256'],
              'currentCertificatePacketPath': str(location), 'certificatePacketSha256': f.digest(location)}
    directory = base/'character-proof-joins'/label
    directory.mkdir(parents=True, exist_ok=True)
    save_pickle(directory/f'{si:03}-{ri:03}-{fi}-raw-zero.pickle', joined)
    if label == native.BASELINE and si == 5:
        name = 'character-cases/'+label+'/source-inputs/005.pickle'
        route = STATE['writes'][name]
        original = str(PREVIOUS/'complete'/name)
        f.require(route['typedIdentity'] and route['sha256'] == STATE['coverageChecks']['paths'][original]
                  and owner == evidence['owner'] and (ri, fi) in STATE['failedAddresses'],
                  'completed six-use diagnostic source proof reused through exact full input')
        source_route = 'completed-failed-source-diagnostic'
    else:
        zero_join(factor, physical['source'], owner, ri, fi, evidence, addresses)
        source_route = 'native-raw-zero-source-join'
    if label not in STATE['zeroControls']:
        wrong_owner = dict(owner, index=owner['index']+10000)
        wrong_raw = dict(factor, AMPLITUDE_RECONSTRUCTION_RESIDUAL=h.sp.Integer(1))
        wrong_coefficient = dict(factor, COEFFICIENT=factor['COEFFICIENT']+1)
        operands = {'factor': factor, 'source': physical['source'], 'owner': owner,
                    'changedOwner': wrong_owner, 'changedRaw': wrong_raw, 'changedCoefficient': wrong_coefficient,
                    'row': ri, 'factorIndex': fi, 'evidence': evidence}
        save_pickle(directory/'raw-zero-mutation-operands.pickle', operands)
        rejects = h.source.inputs.source.rejects
        controls = {
            'owner': rejects(lambda: zero_join(factor, physical['source'], wrong_owner, ri, fi, evidence, addresses)),
            'raw': rejects(lambda: zero_join(wrong_raw, physical['source'], owner, ri, fi, evidence, addresses)),
            'coefficient': rejects(lambda: zero_join(wrong_coefficient, physical['source'], owner, ri, fi, evidence, addresses))}
        save_json(directory/'raw-zero-mutation-controls.json', controls)
        f.require(all(controls.values()), 'actual raw-zero owner/raw/coefficient mutations reject')
        STATE['zeroControls'].add(label)
    return {'row': ri, 'factor': fi, 'owner': owner, 'route': source_route,
            'rawAmplitudeZero': True, 'certificatePresent': False,
            'nativeProofRoute': 'accepted literal raw-zero residual',
            'coverageChecksSha256': STATE['coverageChecksSha256']}


def checker_adapter(proxy):
    original = ast.parse(inspect.getsource(original_checker))
    changed = copy.deepcopy(original)
    loop = next(v for v in ast.walk(changed) if isinstance(v, ast.For)
                and ast.unparse(v.target) == '(ri, fi, factor)' and ast.unparse(v.iter) == "physical['uses']")
    position = next(i for i, v in enumerate(loop.body) if isinstance(v, ast.If)
                    and ast.unparse(v.test) == "owner['kind'] == 'accepted'")
    added = ast.parse('zero_route = raw_zero_route(base,label,si,ri,fi,factor,physical,owner)\n'
                     'if zero_route is not None:\n joined.append(zero_route)\n continue').body
    loop.body[position:position] = added
    reverse = copy.deepcopy(changed)
    restored_loop = next(v for v in ast.walk(reverse) if isinstance(v, ast.For)
                         and ast.unparse(v.target) == '(ri, fi, factor)' and ast.unparse(v.iter) == "physical['uses']")
    del restored_loop.body[position:position+len(added)]
    f.require(ast.dump(reverse) == ast.dump(original), 'whole prior certificate reader reverse AST: raw-zero dispatch only')
    namespace = dict(vars(r), raw_zero_route=raw_zero_route, safe_pickle=save_pickle, f=proxy)
    exec(compile(ast.fix_missing_locations(changed), str(Path(r.__file__))+':native-raw-zero-route', 'exec'), namespace)
    return namespace['check_source_proofs'], {
        'wholePriorCheckerReverseAST': True, 'rawZeroDispatches': 1,
        'originalCheckerAST': hashlib.sha256(ast.dump(original).encode()).hexdigest(),
        'originalCertificateReaderUnchanged': True, 'originalMutationTailUnchanged': True}


def check_source_proofs(base, label, si, physical, factors_packet):
    completed = STATE['completedProofs'].get((label, si))
    if completed is None:
        return STATE['checker'](base, label, si, physical, factors_packet)
    name = 'character-cases/'+label+'/source-inputs/'+f'{si:03}.pickle'
    write = STATE['writes'][name]
    f.require(write['typedIdentity'] and write['sha256'] == completed['sourceInputSha256'],
              'complete actual source input at original native proof call')
    f.require(f.digest(base/'accepted/accepted-cases'/label/'factorization.pickle') == r.STATE['factorCases'][label],
              'actual unchanged full factor producer for completed source proof')
    for name, value in completed['evidence'].items():
        f.require(f.digest(base/name) == value, 'completed proof pair/check/control unchanged')
    checks = json.loads((base/completed['checksPath']).read_text())
    f.require(checks['sourceIndex'] == si and len(checks['joins']) == len(physical['uses'])
              and all(v['certifiedResidualZero'] for v in checks['joins']), 'full completed source proof result')
    if si == 0:
        controls = json.loads((base/'character-proof-joins'/label/'mutation-controls.json').read_text())
        f.require(len(controls) == 4 and all(controls.values()), 'completed actual certificate mutations retained')
    STATE['completedProofReuses'].append({'case': label, 'sourceIndex': si,
        'input': write, 'evidence': completed['evidence'], 'originalCheckerAST': STATE['checkerJoin']['originalCheckerAST'],
        'proofOrMutationReexecuted': False})
    save_json(base/'completed-source-proof-reuse.json', STATE['completedProofReuses'])


def load(base):
    repair = json.loads(REPAIR.read_text())
    previous = PREVIOUS/'complete'
    f.require(f.digest(Path(__file__)) == repair['newHelperSha256'] and f.digest(PLAN) == repair['newPlanSha256'],
              'reviewed exact new continuation helper and plan')
    f.require(f.digest(previous/'inputs.json') == repair['previousInputsSha256'], 'unchanged previous recovery manifest')
    f.require(f.digest(Path(r.__file__)) == repair['previousHelperSha256'] == f.digest(PREVIOUS/'helper-source.py'),
              'whole original recovery helper/current/frozen identity')
    f.require(f.digest(Path(h.__file__)) == repair['originalHelperSha256'], 'whole original character helper unchanged')
    receipt = h.source.receipts.inspect_guard(COVERAGE, 'diagnose')
    f.require(f.digest(COVERAGE/'checks.json') == repair['coverageChecksSha256']
              and (COVERAGE/'checks.json').read_bytes() == (COVERAGE/'diagnose.stdout').read_bytes(),
              'clean completed saved certificate coverage diagnostic')
    old = json.loads((previous/'inputs.json').read_text())
    manifest = dict(old, runDirectory=str(base), sourceFiles=dict(old['sourceFiles']),
                    inputPackets=dict(old['inputPackets']), referencedInputs={})
    actual = {str(p.relative_to(previous)) for p in previous.rglob('*') if p.is_file()
              and 'source' not in p.relative_to(previous).parts and p != previous/'inputs.json'}
    f.require(actual == set(repair['completedFiles']), 'all and only previous completed files')
    STATE.update(base=base, coverageFiles=repair['coverageFiles'], coverageChecksSha256=repair['coverageChecksSha256'])
    for name, item in repair['completedFiles'].items():
        path = previous/name
        f.require(f.digest(path) == item['sha256'] and path.stat().st_size == item['bytes']
                  and (str(path.readlink()) if path.is_symlink() else None) == item['rawLink']
                  and str(path.resolve()) == item['resolved'], 'previous literal reference/regular file identity')
        h.source.reference(base, manifest, path, name, item['sha256'])
    for name in repair['reusedWritePaths']:
        WRITE_ROUTES[name] = manifest['referencedInputs'][name]
    h.source.reference(base, manifest, previous/'inputs.json', 'previous-character-recovery-inputs.json', repair['previousInputsSha256'])
    for root, prefix, inventory in ((PREVIOUS, 'previous-character-recovery-logs', repair['previousLogs']),
                                    (COVERAGE, 'accepted-certificate-coverage', repair['coverageFiles'])):
        for name, item in inventory.items():
            h.source.reference(base, manifest, root/name, prefix+'/'+name, item['sha256'])
    cp = json.loads(r.FACTOR_CP.read_text())
    origin = Path(cp['runDirectory'])
    f.require(f.digest(r.FACTOR_CP) == repair['factorProofCheckpointSha256']
              and f.digest(origin/'checks.json') == cp['checksSha256'], 'actual accepted factor proof producer')
    name = 'accepted-factorization.pickle'
    h.source.reference(base, manifest, origin/name, 'factor-proofs/'+name, cp['artifacts'][name]['sha256'])
    factors = {}
    for label in repair['cases']:
        factors[label] = cp['artifacts']['cases/'+label+'/factorization.pickle']['sha256']
        f.require(f.digest(base/'accepted/accepted-cases'/label/'factorization.pickle') == factors[label],
                  'whole case packet original accepted factor hash')
    for path in (Path(__file__).resolve(), PLAN, REPAIR):
        name = str(path.relative_to(f.ROOT))
        value = f.digest(path)
        f.require(name not in manifest['sourceFiles'], 'fresh recovery source namespace')
        manifest['sourceFiles'][name] = value
    for name, value in manifest['sourceFiles'].items():
        f.require(f.digest(f.ROOT/name) == value, 'current helper/source pin')
        target = base/'source'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        if name in old['sourceFiles']:
            frozen = previous/'source'/name
            f.require(f.digest(frozen) == value, 'previous frozen source unchanged')
            target.symlink_to(frozen)
            manifest['inputPackets'][str(frozen)] = value
        else:
            shutil.copyfile(f.ROOT/name, target)
        f.require(f.digest(target) == value, 'current frozen helper identity')
    for path, value in manifest['inputPackets'].items():
        f.require(f.digest(Path(path)) == value, 'all original/input/reference prehashes')
    coverage = json.loads((COVERAGE/'checks.json').read_text())
    for path, value in coverage['paths'].items():
        f.require(f.digest(Path(path)) == value, 'all diagnostic producer paths unchanged')
        manifest['inputPackets'][path] = value
    STATE['coverageChecks'] = coverage
    STATE['omitted'] = {}
    for label in repair['cases']:
        pairs = f.unpickle(base/'accepted-certificate-coverage'/(label+'-omitted-certificate-pairs.pickle'))
        STATE['omitted'][label] = {(v['row'], v['factorIndex']): v for v in pairs}
        f.require(len(STATE['omitted'][label]) == len(pairs), 'unique actual omitted factor addresses')
    failed = f.unpickle(base/'accepted-certificate-coverage/failed-source-omitted-certificate-pairs.pickle')
    STATE['failedAddresses'] = {(v['row'], v['factorIndex']) for v in failed['pairs']}
    f.require(len(STATE['failedAddresses']) == coverage['failedSourceUses'] == 6
              and coverage['rawZeroBranchExplainsFailure'], 'complete failed-source saved diagnostic')
    STATE['completedProofs'] = {(v['case'], v['sourceIndex']): v for v in repair['completedSourceProofs']}
    r.STATE.update(base=base, factorCases=factors,
        diagnosticChecks=json.loads((r.DIAGNOSTIC/'checks.json').read_text()),
        diagnosticPairs=f.unpickle(base/'accepted-proof-diagnostic/failed-guard-source-certificate-pairs.pickle'))
    manifest['rawZeroProofContinuation'] = {'previous': str(PREVIOUS), 'previousInputsSha256': repair['previousInputsSha256'],
        'coverage': str(COVERAGE), 'coverageChecksSha256': repair['coverageChecksSha256'], 'coverageGuard': receipt,
        'originalFiles': len(repair['completedFiles']), 'completedSourceProofsReused': len(repair['completedSourceProofs']),
        'sourcePrefixWrites': len(WRITE_ROUTES), 'newScientificOperations': 0,
        'originalPrepareJoin': r.STATE['adapterJoin'], 'certificateCheckerJoin': STATE['checkerJoin']}
    f.save(base/'inputs.json', manifest)
    f.save(base/'zero-recovery-reference-reuse.json', manifest['referencedInputs'])
    f.save(base/'raw-zero-reader-join.json', manifest['rawZeroProofContinuation'])
    return manifest, tuple(repair['cases'])


if __name__ == '__main__':
    proxy = SimpleNamespace(**vars(f))
    proxy.atomic_pickle = save_pickle
    proxy.save = save_json
    checker, checker_join = checker_adapter(proxy)
    STATE.update(checker=checker, checkerJoin=checker_join)
    r.check_source_proofs = check_source_proofs
    h.f = proxy
    prepare, prepare_join = r.adapter()
    repair = json.loads(REPAIR.read_text())
    f.require(prepare_join == repair['originalPrepareJoin'] and checker_join == repair['checkerJoin'],
              'reviewed whole prepare and certificate-reader reverse AST joins')
    r.STATE['adapterJoin'] = prepare_join
    h.load = load
    h.prepare = prepare
    r.original_main()
