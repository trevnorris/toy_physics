#!/usr/bin/env python3
"""Finish only the saved material-output main's remaining file audit (stdlib)."""
import argparse
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
STORE = REPO / '_scratch/s11c'
M = ROOT / '_measurements'
ORIGIN = STORE / 's11c-remaining-case-coordinate-20260921/response/output/production'
WORKER = M / 'S11c_d_remaining_case_coordinate_output.py'
PLAN = M / 'S11c_d_remaining_case_coordinate_output_audit_finish_plan.md'
PROOF = M / 'S11c_d_remaining_case_coordinate_output_audit_finish_wiring.json'
PART = 'coordinate_MATERIAL_ADVECTED__RHO4_CONSTANT'
PHASE = '19_replay_' + PART
CASE = 'MATERIAL_ADVECTED__RHO4_CONSTANT'


def digest(path):
    out = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            out.update(block)
    return out.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def save(path, value):
    path = Path(path)
    require(not path.is_symlink(), 'never write through a preserved-file reference')
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.new')
    require(not temp.exists(), 'no stale temporary file replacement')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


def node(tree, name):
    return next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)


def final_audit_ast():
    """Partition the literal whole main; alter no statement in its final suffix."""
    tree = ast.parse(WORKER.read_text())
    main = node(tree, 'main')
    split = next(i for i, n in enumerate(main.body)
                 if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
                 and isinstance(n.value.func, ast.Name) and n.value.func.id == 'hashes')
    prefix, suffix = copy.deepcopy(main.body[:split]), copy.deepcopy(main.body[split:])
    reverse = copy.deepcopy(main)
    reverse.body = prefix + suffix
    require(ast.dump(reverse) == ast.dump(main), 'whole main reverse partition')
    # These are the actual final durable prefix writes before the split.
    text = ast.unparse(main)
    require("f.save(base / 'part-inventory.json', parts)" in text and
            "parts[name] = recovered.finish_tail()" in text,
            'saved complete post-control result precedes final audit')
    args = ('base', 'manifest', 'parts', 'cases', 'combined', 'args', 'started')
    fn = ast.FunctionDef(name='finish_native_audit', args=ast.arguments(
        posonlyargs=[], args=[ast.arg(x) for x in args], vararg=None,
        kwonlyargs=[], kw_defaults=[], kwarg=None, defaults=[]),
        body=suffix, decorator_list=[])
    hashes = copy.deepcopy(node(tree, 'hashes'))
    module = ast.fix_missing_locations(ast.Module(body=[hashes, fn], type_ignores=[]))
    proof = {'wholeMainReversePartition': True, 'nativeFinalAuditUnchanged': True,
             'nativeHashesUnchanged': True, 'splitStatement': split,
             'mainAstSha256': hashlib.sha256(ast.dump(main).encode()).hexdigest(),
             'suffixAstSha256': hashlib.sha256(ast.dump(ast.Module(body=suffix, type_ignores=[])).encode()).hexdigest(),
             'originalOutputSha256': digest(WORKER),
             'noNativeEmissionReplayOrScientificImports': True}
    return module, proof


def validate_completed_part(part, done, check, replay, controls, join):
    require(part['case'] == CASE and part['kind'] == 'coordinate' and
            part['directory'] == 'parts/' + PART, 'exact completed physical part address')
    require(done['status'] == 'COMPLETED_EMISSION_AWAITING_REPLAY', 'actual saved emission')
    require(part['sha256'] == done['transcriptSha256'] == replay['transcriptSha256'], 'same saved transcript')
    require(done['packetSha256'] == check['packetSha256'] == replay['packetSha256'], 'same saved full packet')
    require(part['tags'] == done['tags'] == check['tags'] == replay['tags'] == replay['seen'] == 3500,
            'whole completed native payload replay')
    require(part['keys'] == done['keys'] == len(check['keys']) == replay['keys'] == 1748,
            'all original emitted keys')
    require(len(set(check['keys'].values())) == 1748, 'unique saved output keys')
    require(part['metadataPaths'] == check['metadataPaths'] == 159716, 'completed full metadata census')
    require(check['emitterJoin'] == join and check['replayControls'] == controls,
            'complete native post-control result and exact source join')
    require(set(controls) == {'payload', 'unit'} and
            controls['unit']['original'] != controls['unit']['changed'] and
            controls['payload']['tag'] in check['keys'] and
            controls['unit']['tag'] == controls['payload']['tag'].replace('PY_S11CD_', 'PY_S11CD_METADATA_', 1),
            'actual saved changed payload/unit controls')
    require(replay['nativeReplayPredicateUnchanged'] and
            replay['sourceHelperSha256'] == digest(WORKER), 'original completed replay source')


def completed_part_evidence():
    base = ORIGIN / 'complete'
    folder = base / 'parts' / PART
    parts = json.loads((base / 'part-inventory.json').read_text())
    require(len(parts) == 9 and PART in parts, 'eight earlier parts and one completed unaudited part')
    names = ('emission-complete.json', 'emission-checks.json', 'completed-payload-replay.json',
             'replay-controls.json', 'validation-join.json')
    records = [json.loads((folder / name).read_text()) for name in names]
    validate_completed_part(parts[PART], *records)
    return parts, records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-directory', type=Path, required=True)
    args0 = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS, (512 * 1024**2, 512 * 1024**2))
    signal.alarm(900)
    started = time.monotonic()
    base = args0.run_directory.resolve()
    base.relative_to(STORE)
    require(not base.exists(), 'fresh audit continuation only')
    reviewed = json.loads(PROOF.read_text())
    require(reviewed['status'] == 'PASSED_SAVED_PART_AUDIT_CONTINUATION_WIRING' and
            reviewed['helperSha256'] == digest(__file__) and reviewed['planSha256'] == digest(PLAN),
            'exact reviewed audit-only helper')
    module, join = final_audit_ast()
    require(join == reviewed['nativeJoin'], 'unchanged whole original final suffix')
    inv = json.loads((ORIGIN / 'active.json').read_text())
    require(inv == json.loads((ORIGIN / 'coordinate_output_construct.invocation.json').read_text()) and
            inv['exitCode'] == 1 and inv['status'] == 'failed' and
            not (ORIGIN / 'staged-completion.json').exists(), 'preserved failed sequence, never acceptance')
    outcomes = json.loads((ORIGIN / 'stage-outcomes.json').read_text())
    require(len(outcomes) == 20 and [x['index'] for x in outcomes] == list(range(20)) and
            all(x['accepted'] and x['actualGuardProcessExitCode'] == 0 for x in outcomes[:19]) and
            not outcomes[-1]['accepted'] and outcomes[-1]['actualGuardProcessExitCode'] == 124,
            'exact nineteen clean phases and final time-limit stop')
    phase = ORIGIN / 'phases' / PHASE
    guard = json.loads((phase / 'resource-guard/outcome.json').read_text())
    require(guard['exitCode'] == guard['childOutcome']['exitCode'] == 124 and
            guard['childOutcome']['guardReason'] == 'wall-time limit' and guard['stderrBytes'] == 0 and
            not (phase / 'checks.json').exists(), 'only unfinished final audit after timeout')
    for name in ('guard.stderr', 'output_phase.stderr', 'resource-guard/stderr'):
        require((phase / name).stat().st_size == 0, 'original strict child stderr empty')
    with (phase / 'resource-guard/resource-samples.jsonl').open() as samples:
        for line in samples:
            v = json.loads(line); events = dict(x.split() for x in v['memory.events'].splitlines())
            require(all(events.get(k, '0') == '0' for k in ('oom', 'oom_kill', 'oom_group_kill')) and
                    v['memory.swap.current'] == '0', 'no original OOM or job swap')
    parts, records = completed_part_evidence()
    oldbase = ORIGIN / 'complete'
    require(digest(oldbase / 'inputs.json') == reviewed['originalManifestSha256'] and
            digest(oldbase / 'part-inventory.json') == reviewed['originalPartInventorySha256'] and
            digest(oldbase / 'parts' / PART / 'emission-checks.json') == reviewed['originalCompletedChecksSha256'],
            'exact reviewed saved manifest, inventory and completed native checks')
    old = json.loads((oldbase / 'inputs.json').read_text())
    require(digest(oldbase / 'parts' / PART / 'full.out') == parts[PART]['sha256'] and
            digest(oldbase / 'parts' / PART / 'coordinate-output.pickle') == records[0]['packetSha256'],
            'actual completed stream and bundle content')
    # Preserve complete existing files by explicit references, without rewriting
    # them or allocating another multi-gigabyte copy. Mutable state gets new files.
    base.mkdir(parents=True)
    references = {}
    for p in sorted(oldbase.rglob('*')):
        if not p.is_file():
            continue
        require(not p.is_symlink(), 'original run regular-file provenance')
        n = str(p.relative_to(oldbase))
        target = 'original-production-inputs.json' if n == 'inputs.json' else 'original-production-part-inventory.json' if n == 'part-inventory.json' else n
        dest = base / target; dest.parent.mkdir(parents=True, exist_ok=True)
        dest.symlink_to(p)
        references[target] = {'original': str(p), 'sha256': digest(p), 'bytes': p.stat().st_size}
    manifest = copy.deepcopy(old)
    manifest['runDirectory'] = str(base)
    for n, info in references.items():
        manifest['inputPackets'][info['original']] = info['sha256']
    for p in sorted(phase.rglob('*')):
        if p.is_file(): manifest['inputPackets'][str(p)] = digest(p)
    for name in ('active.json', 'coordinate_output_construct.invocation.json', 'coordinate_output_construct.stderr',
                 'coordinate_output_construct.stdout', 'supervisor.stdout', 'supervisor.stderr',
                 'stage-outcomes.json', 'staged-schedule.json', 'launch.json'):
        p = ORIGIN / name; manifest['inputPackets'][str(p)] = digest(p)
    for p in (Path(__file__).resolve(), PLAN, PROOF):
        n = str(p.relative_to(ROOT)); require(n not in manifest['sourceFiles'], 'new continuation source only')
        manifest['sourceFiles'][n] = digest(p)
        dest = base / 'source' / n; dest.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(p, dest)
    continuation = {'originalRunDirectory': str(oldbase), 'originalPhase': 19,
                    'originalExitCode': 124, 'originalGuardReason': 'wall-time limit',
                    'completedNativeReplayAndPostControlsReused': True,
                    'originalFailedPhaseAccepted': False, 'newPhysicalEmission': 0,
                    'newPayloadReplay': 0, 'newScientificConstruction': 0,
                    'referenceFiles': len(references), 'referenceBytes': sum(v['bytes'] for v in references.values()),
                    'writesThroughReferencesForbidden': True, 'remainingOriginalPhases': [20, 21, 22, 23],
                    'originalChecksSha256': digest(oldbase / 'parts' / PART / 'emission-checks.json')}
    manifest['completedOutputAuditContinuation'] = continuation
    save(base / 'referenced-completed-artifacts.json', references)
    save(base / 'completed-part-audit-continuation.json', continuation)
    save(base / 'native-final-audit-join.json', join)
    save(base / 'inputs.json', manifest)
    save(base / 'part-inventory.json', parts)
    # The literal original final hash, artifact inventory and output/checks suffix.
    # Only saved locals and the fresh receipt destination are supplied.
    namespace = dict(Path=Path, json=json, signal=signal, time=time,
                     f=SimpleNamespace(ROOT=ROOT, digest=digest, require=require, save=save))
    exec(compile(module, str(WORKER) + ':saved-final-audit', 'exec'), namespace)
    native_args = SimpleNamespace(stage='replay', phase_directory=base.parent)
    namespace['finish_native_audit'](base, manifest, parts, {}, None, native_args, started)


if __name__ == '__main__':
    main()
