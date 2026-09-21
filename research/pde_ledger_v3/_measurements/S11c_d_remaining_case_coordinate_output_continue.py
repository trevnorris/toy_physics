#!/usr/bin/env python3
"""Remaining saved-stream work and final audits, in separate bounded workers."""
import argparse
import ast
import copy
import importlib.util
import json
from pathlib import Path
import resource
import shutil
import signal
import sys
import time
from types import SimpleNamespace
import S11c_d_remaining_case_coordinate_output_audit_finish as audit

ROOT, REPO, M, STORE = audit.ROOT, audit.REPO, audit.M, audit.STORE
CP = M / 'S11c_d_remaining_case_coordinate_output_audit_checkpoint.json'
PLAN = M / 'S11c_d_remaining_case_coordinate_output_continue_plan.md'
PROOF = M / 'S11c_d_remaining_case_coordinate_output_continue_wiring.json'
COORDINATOR = M / 'S11c_d_remaining_case_coordinate_output_continue_staged.py'
STAGED = M / 'S11c_d_remaining_case_coordinate_output_staged.py'
CASE = 'MATERIAL_ADVECTED__RHOBR_CONSTANT'
sha, require, save = audit.digest, audit.require, audit.save


def staged_module():
    spec = importlib.util.spec_from_file_location('original_staged_receipt_reader', STAGED)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def partition():
    tree = ast.parse(audit.WORKER.read_text()); original = audit.node(tree, 'main')
    suffix_module, join = audit.final_audit_ast(); split = join['splitStatement']
    prefix = copy.deepcopy(original)
    prefix.name = 'native_stage_work'
    prefix.body = copy.deepcopy(original.body[:split]) + ast.parse('return locals()').body
    reverse = copy.deepcopy(prefix); reverse.name = original.name
    reverse.body = reverse.body[:-1] + copy.deepcopy(original.body[split:])
    require(ast.dump(reverse) == ast.dump(original), 'entire native main reverse split')
    module = ast.fix_missing_locations(ast.Module(body=[prefix], type_ignores=[]))
    join.update(nativeWorkPrefixUnchanged=True, soleWorkAddition='return locals()',
                nativeEmissionAndReplayBodiesUnchanged=True,
                contextRoute='existing physical source address, exact source and unit-proof hash joins')
    return module, suffix_module, join


def reviewed_sources():
    proof = json.loads(PROOF.read_text())
    require(proof['status'] == 'PASSED_SPLIT_SAVED_OUTPUT_CONTINUATION_WIRING' and
            proof['workerSha256'] == sha(__file__) and proof['coordinatorSha256'] == sha(COORDINATOR) and
            proof['planSha256'] == sha(PLAN), 'reviewed unchanged continuation sources')
    for path, digest in proof['originalSourcePins'].items(): require(sha(path) == digest, 'original helper/guard source join')
    _, _, join = partition(); require(join == proof['nativeJoin'], 'exact whole native prefix/suffix join')
    return proof, join


def prepare(args, started):
    cp = json.loads(CP.read_text()); origin = Path(cp['runDirectory']); outer = origin.parent
    require(cp['status'] == 'ACCEPTED_MATERIAL_OUTPUT_PHASE19_AUDIT_CONTINUATION' and
            not cp['wholeOutputAccepted'] and not cp['originalFailedPhaseAccepted'], 'accepted audit continuation only')
    reader = staged_module(); guard = reader.inspect_guard(outer, 'coordinate_output_audit')
    check_path = Path(cp['checksPath']); require(sha(check_path) == cp['checksSha256'] and
            reader.identical(check_path, outer / 'coordinate_output_audit.stdout'), 'actual accepted audit final checks')
    old = json.loads(check_path.read_text()); base = args.run_directory.resolve()
    require(not base.exists(), 'fresh reference workspace'); base.mkdir(parents=True)
    manifest = json.loads((origin / 'inputs.json').read_text())
    require(manifest['sourceFiles'] == old['sourceFiles'] and manifest['inputPackets'] == old['inputPackets'], 'accepted manifest identity')
    references = {}
    for p in sorted(origin.rglob('*')):
        if not p.is_file(): continue
        n = str(p.relative_to(origin))
        target = 'accepted-audit-inputs.json' if n == 'inputs.json' else 'accepted-audit-part-inventory.json' if n == 'part-inventory.json' else n
        digest = sha(p)
        if n in old['artifacts']: require(digest == old['artifacts'][n]['sha256'], 'all accepted artifact bytes')
        elif n.startswith('source/'): require(digest == old['sourceFiles'][n[len('source/'):]], 'all accepted frozen source bytes')
        else: require(n == 'inputs.json', 'complete accepted reference census')
        dest = base / target; dest.parent.mkdir(parents=True, exist_ok=True); dest.symlink_to(p)
        references[target] = {'source': str(p), 'resolvedSource': str(p.resolve()), 'sha256': digest, 'bytes': p.stat().st_size}
    require(set(old['artifacts']) <= {str(p.relative_to(origin)) for p in origin.rglob('*') if p.is_file()}, 'no missing accepted artifacts')
    manifest['runDirectory'] = str(base)
    for info in references.values(): manifest['inputPackets'][info['source']] = info['sha256']
    for p in (CP, check_path, outer / 'coordinate_output_audit.invocation.json', outer / 'resource-guard/outcome.json'):
        manifest['inputPackets'][str(p)] = sha(p)
    additions = {}
    for p in (Path(__file__).resolve(), COORDINATOR, PLAN, PROOF):
        n = str(p.relative_to(ROOT)); require(n not in manifest['sourceFiles'], 'new continuation sources only')
        additions[n] = sha(p); manifest['sourceFiles'][n] = additions[n]
        dest = base / 'source' / n; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(p, dest)
    manifest['savedOutputContinuation'] = {'acceptedAuditDirectory': str(origin), 'acceptedAuditChecksSha256': cp['checksSha256'],
        'physicalSourceBase': old['completedOutputAuditContinuation']['originalRunDirectory'],
        'referenceFiles': len(references), 'newSourceFiles': additions, 'workAndFinalAuditSeparate': True,
        'originalFailedPhaseStillFailed': True, 'noCompletedEmissionOrReplayRepeated': True}
    save(base / 'continuation-input-references.json', references)
    save(base / 'continuation-source-joins.json', manifest['savedOutputContinuation'])
    save(base / 'inputs.json', manifest); save(base / 'part-inventory.json', old['parts'])
    result = {'status': 'COMPLETED_MATERIAL_OUTPUT_REFERENCE_PREPARATION', 'sourceFiles': manifest['sourceFiles'],
        'manifestSha256': sha(base / 'inputs.json'), 'partInventorySha256': sha(base / 'part-inventory.json'),
        'acceptedParts': len(old['parts']), 'references': len(references), 'referenceBytes': sum(x['bytes'] for x in references.values()),
        'allAcceptedArtifactJoins': len(old['artifacts']), 'priorGuard': guard,
        'newScientificWork': 0, 'newEmission': 0, 'newPayloadReplay': 0, 'seconds': time.monotonic() - started}
    save(args.phase_directory / 'checks.json', result); signal.alarm(0); print(json.dumps(result, indent=2))


def reference_context(h, base, label, packet, folder):
    """Keep the producer's literal saved source address; no operand normalization."""
    manifest = json.loads((base / 'inputs.json').read_text())
    source_base = Path(manifest['savedOutputContinuation']['physicalSourceBase'])
    require(source_base == audit.ORIGIN / 'complete', 'exact original physical source base')
    rel = Path('numerical/boundary-inputs/bindings/sources/accepted-cases') / label / 'reduced-action.pickle'
    source, reference = source_base / rel, base / rel
    expected = manifest['copiedInputs'][str(rel)]
    require(reference.resolve() == source.resolve() and sha(source) == sha(reference) == expected, 'full actual source path/hash join')
    old_proof = folder / 'context-unit-join.pickle'
    require(old_proof.is_symlink() and old_proof.resolve() == (source_base / 'parts' / folder.name / old_proof.name).resolve(), 'original complete context/unit proof')
    joined = {'newReference': str(reference), 'originalSource': str(source), 'sourceSha256': expected,
              'rawProofReference': str(old_proof), 'rawProofSource': str(old_proof.resolve()),
              'rawProofSha256': sha(old_proof), 'nativeContextBodyUnchanged': True,
              'packetDimensionStateUnchanged': True, 'unitSupplementUnchanged': True}
    save(folder / 'continuation-context-address-join.json', joined)
    return h.context(source_base, label, packet, folder)


def work(args, started, join):
    require(args.stage in ('replay', 'aggregate'), 'only remaining saved output work')
    if args.stage == 'replay': require(args.case == CASE and args.kind in ('continuum', 'current', 'coordinate'), 'only three unfinished native replays')
    base = args.run_directory.resolve(); before_manifest = sha(base / 'inputs.json'); before_parts = sha(base / 'part-inventory.json')
    # Scientific libraries are imported only inside this guarded work process.
    import S11c_d_remaining_case_coordinate_output_recover as recovered
    h = recovered.h
    original_save, original_pickle = h.f.save, h.f.atomic_pickle
    def protected_save(path, value):
        require(not Path(path).is_symlink(), 'cannot rewrite preserved referenced JSON')
        return original_save(path, value)
    def protected_pickle(path, value, *a, **kw):
        require(not Path(path).is_symlink(), 'cannot rewrite preserved scientific packet')
        return original_pickle(path, value, *a, **kw)
    h.f.save = protected_save; h.f.atomic_pickle = protected_pickle
    namespace = dict(vars(h), context=lambda b,l,p,d: reference_context(h,b,l,p,d))
    prefix, _, actual_join = partition(); require(actual_join == join, 'unchanged native work source join')
    exec(compile(prefix, str(audit.WORKER) + ':saved-work-prefix', 'exec'), namespace)
    argv = [str(audit.WORKER), '--run-directory', str(base), '--stage', args.stage, '--phase-directory', str(args.phase_directory)]
    if args.stage == 'replay': argv += ['--case', args.case, '--kind', args.kind]
    original_argv = sys.argv; sys.argv = argv
    try: state = namespace['native_stage_work']()
    finally: sys.argv = original_argv
    require(state['manifest'] == json.loads((base / 'inputs.json').read_text()) and
            state['parts'] == json.loads((base / 'part-inventory.json').read_text()), 'all native prefix writes saved')
    result = {'status': 'COMPLETED_SAVED_OUTPUT_WORK_AWAITING_FINAL_AUDIT', 'stage': args.stage,
        'case': args.case, 'kind': args.kind, 'runDirectory': str(base),
        'manifestBeforeSha256': before_manifest, 'partInventoryBeforeSha256': before_parts,
        'manifestSha256': sha(base / 'inputs.json'), 'partInventorySha256': sha(base / 'part-inventory.json'),
        'nativeCases': state['cases'], 'nativeAggregate': state['combined'], 'nativeJoin': join,
        'newIndividualEmissions': 0, 'newScientificWork': 0, 'workSeconds': time.monotonic() - started,
        'finalAuditPending': True}
    if args.stage == 'replay':
        name = args.kind + '_' + args.case; part = state['parts'][name]; folder = base / part['directory']
        require(sha(folder / 'full.out') == part['sha256'], 'original stream preserved')
        result['part'] = part
        result['completedEvidence'] = {n: sha(folder / n) for n in ('full.out', 'emission-checks.json', 'replay-controls.json', 'validation-join.json')}
    else:
        result['completedEvidence'] = {n: sha(base / n) for n in ('full.out', 'original-combined.out', 'aggregation-checks.json')}
    save(args.phase_directory / 'work-state.json', result); save(args.phase_directory / 'checks.json', result)
    signal.alarm(0); print(json.dumps(result, indent=2))


def finish(args, started, join):
    base = args.run_directory.resolve(); previous = args.work_directory.resolve(); reader = staged_module()
    guard = reader.inspect_guard(previous, 'output_continue')
    require(reader.identical(previous / 'checks.json', previous / 'output_continue.stdout'), 'clean exact preceding work checks/stdout')
    state = json.loads((previous / 'work-state.json').read_text())
    require(state == json.loads((previous / 'checks.json').read_text()) and
            state['status'] == 'COMPLETED_SAVED_OUTPUT_WORK_AWAITING_FINAL_AUDIT' and state['finalAuditPending'], 'completed work, awaiting audit')
    require(state['stage'] == args.stage and state['case'] == args.case and state['kind'] == args.kind and state['runDirectory'] == str(base), 'exact preceding work address')
    require(state['nativeJoin'] == join and sha(base / 'inputs.json') == state['manifestSha256'] and
            sha(base / 'part-inventory.json') == state['partInventorySha256'], 'saved locals and original source identity')
    manifest = json.loads((base / 'inputs.json').read_text()); parts = json.loads((base / 'part-inventory.json').read_text())
    folder = base / 'parts' / (args.kind + '_' + args.case) if args.stage == 'replay' else base
    for n, digest in state['completedEvidence'].items(): require(sha(folder / n) == digest, 'unchanged completed work evidence')
    for p in sorted(previous.rglob('*')):
        if p.is_file(): manifest['inputPackets'][str(p)] = sha(p)
    manifest.setdefault('completedSplitOutputStages', []).append({'stage': args.stage, 'case': args.case, 'kind': args.kind,
        'workDirectory': str(previous), 'checksSha256': sha(previous / 'checks.json'), 'actualWorkGuard': guard,
        'workSeconds': state['workSeconds'], 'nativeMainSplitJoined': True})
    save(base / 'phase-receipts' / (previous.name + '-work.json'), state)
    save(base / 'inputs.json', manifest)
    _, suffix, actual_join = partition(); require(actual_join == join, 'unchanged native final audit source')
    namespace = dict(Path=Path, json=json, signal=signal, time=time,
                    f=SimpleNamespace(ROOT=ROOT, digest=sha, require=require, save=save))
    exec(compile(suffix, str(audit.WORKER) + ':split-final-audit', 'exec'), namespace)
    native_args = SimpleNamespace(stage=args.stage, phase_directory=args.phase_directory)
    namespace['finish_native_audit'](base, manifest, parts, state['nativeCases'], state['nativeAggregate'], native_args, started)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--mode', choices=('prepare', 'work', 'audit'), required=True)
    ap.add_argument('--run-directory', type=Path, required=True); ap.add_argument('--phase-directory', type=Path, required=True)
    ap.add_argument('--stage', choices=('replay', 'aggregate')); ap.add_argument('--case'); ap.add_argument('--kind', choices=('continuum', 'current', 'coordinate'))
    ap.add_argument('--work-directory', type=Path); args = ap.parse_args()
    args.run_directory.resolve().relative_to(STORE); args.phase_directory.resolve().relative_to(STORE)
    if args.mode != 'work': resource.setrlimit(resource.RLIMIT_AS, (512 * 1024**2, 512 * 1024**2))
    signal.alarm(900); started = time.monotonic(); _, join = reviewed_sources()
    if args.mode == 'prepare': prepare(args, started)
    elif args.mode == 'work': work(args, started, join)
    else: finish(args, started, join)


if __name__ == '__main__': main()
