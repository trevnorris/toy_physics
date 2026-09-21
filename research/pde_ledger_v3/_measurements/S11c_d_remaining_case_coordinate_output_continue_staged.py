#!/usr/bin/env python3
"""One guarded worker at a time for the nine remaining split-output phases."""
import argparse
from datetime import datetime, timezone
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
M = ROOT / 'research/pde_ledger_v3/_measurements'
STORE = ROOT / '_scratch/s11c'
WORKER = M / 'S11c_d_remaining_case_coordinate_output_continue.py'
PLAN = M / 'S11c_d_remaining_case_coordinate_output_continue_plan.md'
PROOF = M / 'S11c_d_remaining_case_coordinate_output_continue_wiring.json'
CP = M / 'S11c_d_remaining_case_coordinate_output_audit_checkpoint.json'
GUARD = ROOT / 'scripts/s11c_guarded_run.py'
SUPERVISOR = M / 'S11c_d_end_normalization_run.py'
STAGED = M / 'S11c_d_remaining_case_coordinate_output_staged.py'
spec = importlib.util.spec_from_file_location('original_staged_protocol', STAGED)
prior = importlib.util.module_from_spec(spec); spec.loader.exec_module(prior)
sha, save, identical, inspect_guard = prior.sha, prior.save, prior.identical, prior.inspect_guard
THREADS = prior.THREADS
CASE = 'MATERIAL_ADVECTED__RHOBR_CONSTANT'


def schedule():
    stages = [{'stage': 'replay', 'case': CASE, 'kind': k} for k in ('continuum', 'current', 'coordinate')] + [{'stage': 'aggregate'}]
    return [{'mode': 'prepare'}] + [dict(v, mode=mode) for v in stages for mode in ('work', 'audit')]


def check_schedule(entries):
    assert entries == schedule() and len(entries) == 9, 'only remaining work followed by its audit'
    assert not any(x.get('stage') in ('emit', 'native') for x in entries)


def worker_command(base, directory, entry, work_directory=None):
    command = [sys.executable, '-u', str(WORKER), '--mode', entry['mode'], '--run-directory', str(base), '--phase-directory', str(directory)]
    if entry['mode'] != 'prepare':
        command += ['--stage', entry['stage']]
        if entry['stage'] == 'replay': command += ['--case', entry['case'], '--kind', entry['kind']]
    if entry['mode'] == 'audit':
        assert work_directory is not None
        command += ['--work-directory', str(work_directory)]
    return command


def coordinate(args):
    resource.setrlimit(resource.RLIMIT_AS, (512 * 1024**2, 512 * 1024**2)); os.nice(max(0, 15 - os.getpriority(os.PRIO_PROCESS, 0)))
    base = args.run_directory.resolve(); base.relative_to(STORE); outer = base.parent
    assert not base.exists() and not (STORE / 'PAUSED_HOST_FREEZE.json').exists()
    proof = json.loads(PROOF.read_text()); cp = json.loads(CP.read_text())
    assert proof['status'] == 'PASSED_SPLIT_SAVED_OUTPUT_CONTINUATION_WIRING'
    assert proof['coordinatorSha256'] == sha(__file__) and proof['workerSha256'] == sha(WORKER) and proof['planSha256'] == sha(PLAN)
    assert cp['status'] == 'ACCEPTED_MATERIAL_OUTPUT_PHASE19_AUDIT_CONTINUATION'
    origin = Path(cp['runDirectory']); check = Path(cp['checksPath'])
    assert sha(check) == cp['checksSha256'] and identical(check, origin.parent / 'coordinate_output_audit.stdout')
    inspect_guard(origin.parent, 'coordinate_output_audit')
    pins = dict(proof['originalSourcePins'])
    pins.update({str(p): sha(p) for p in (WORKER, Path(__file__).resolve(), PLAN, PROOF, CP)})
    entries = schedule(); check_schedule(entries)
    phases = outer / 'phases'; phases.mkdir(exist_ok=False); records = []; current = [None]; previous_work = None; started = time.monotonic()
    save(outer / 'continuation-schedule.json', {'entries': entries, 'sourcePins': pins, 'phaseSeconds': 900,
        'maximumConcurrentWorkers': 1, 'automaticRetries': 0, 'acceptedAuditChecksSha256': cp['checksSha256'],
        'originalFailedPhaseRemainsFailed': True})
    def stop(signum, frame):
        if current[0] is not None:
            current[0].send_signal(signal.SIGINT); current[0].wait()
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, stop); signal.signal(signal.SIGINT, stop)
    for index, entry in enumerate(entries):
        for p, digest in pins.items(): assert sha(p) == digest, ('changed source', p)
        directory = phases / (str(index).zfill(2) + '_' + entry['mode'] + '_' + entry.get('stage', 'inputs') + ('_' + entry['kind'] if 'kind' in entry else ''))
        directory.mkdir(); native = worker_command(base, directory, entry, previous_work)
        assert native == proof['workerCommands'][index], 'exact reviewed worker route'
        command = [sys.executable, str(GUARD), '--log-directory', str(directory / 'resource-guard'), '--seconds', '900', '--',
                   sys.executable, str(SUPERVISOR), '--run-root', str(directory), '--stage', 'output_continue', '--', *native]
        save(directory / 'command.json', {'entry': entry, 'command': command, 'workerCommand': native, 'sourcePins': pins})
        if (base / 'inputs.json').exists(): shutil.copyfile(base / 'inputs.json', directory / 'input-manifest.json')
        tick = time.monotonic()
        with (directory / 'guard.stdout').open('xb') as out, (directory / 'guard.stderr').open('xb') as err:
            process = subprocess.Popen(command, cwd=ROOT, env=dict(os.environ, **THREADS), stdin=subprocess.DEVNULL,
                                       stdout=out, stderr=err, start_new_session=True, close_fds=True)
            current[0] = process; code = process.wait(); current[0] = None
        record = {'index': index, 'entry': entry, 'directory': str(directory), 'command': command,
                  'actualGuardProcessExitCode': code, 'wallSeconds': time.monotonic() - tick, 'accepted': False}
        records.append(record); save(directory / 'outcome.json', record); save(outer / 'stage-outcomes.json', records)
        assert code == 0, ('failed continuation phase; preserve everything, no retry', record)
        record['resourceGuard'] = inspect_guard(directory, 'output_continue')
        final = entry['mode'] == 'audit' and entry.get('stage') == 'aggregate'
        check = base / 'checks.json' if final else directory / 'checks.json'
        assert identical(check, directory / 'output_continue.stdout'), 'actual checks/stdout identity'
        checks = json.loads(check.read_text())
        if entry['mode'] == 'prepare': assert checks['status'] == 'COMPLETED_MATERIAL_OUTPUT_REFERENCE_PREPARATION' and checks['acceptedParts'] == 9
        elif entry['mode'] == 'work':
            assert checks['status'] == 'COMPLETED_SAVED_OUTPUT_WORK_AWAITING_FINAL_AUDIT' and checks['finalAuditPending']
            previous_work = directory
        else:
            assert checks['status'] == ('COMPLETED_FOUR_CASE_MATERIAL_OUTPUT' if final else 'COMPLETED_MATERIAL_OUTPUT_PHASE')
            assert checks['stage'] == entry['stage']
            if final: assert len(checks['parts']) == 12 and checks['aggregate']['parts'] == 12
            else: assert entry['kind'] + '_' + entry['case'] in checks['parts']
        record.update(accepted=True, checksSha256=sha(check), stdoutSha256=sha(directory / 'output_continue.stdout'),
                      finalAuditPending=entry['mode'] == 'work')
        save(directory / 'outcome.json', record); save(outer / 'stage-outcomes.json', records)
        # Do not change work's saved manifest between its completion and audit.
        if entry['mode'] != 'work' and not final:
            manifest = json.loads((base / 'inputs.json').read_text())
            for p in directory.rglob('*'):
                if p.is_file(): manifest['inputPackets'][str(p)] = sha(p)
            manifest['inputPackets'][str(outer / 'continuation-schedule.json')] = sha(outer / 'continuation-schedule.json')
            save(base / 'phase-receipts' / (directory.name + '.json'), record)
            save(base / 'inputs.json', manifest)
        del checks
    assert len(records) == 9 and all(x['accepted'] for x in records)
    save(outer / 'staged-completion.json', {'status': 'COMPLETED_REMAINING_GUARDED_MATERIAL_OUTPUT_PHASES',
         'continuationPhases': 9, 'originalPassedPhases': 19, 'originalTimedOutPhase': 19,
         'acceptedPhase19AuditChecksSha256': cp['checksSha256'], 'allTwelvePartsCompleted': True,
         'originalFailedPhaseRemainsFailed': True, 'checksSha256': sha(base / 'checks.json'),
         'sourcePins': pins, 'wallSeconds': time.monotonic() - started})
    with (base / 'checks.json').open('rb') as source: shutil.copyfileobj(source, sys.stdout.buffer)


def supervise(args):
    outer = args.run_directory.resolve().parent; outer.relative_to(STORE)
    lock = (STORE / 'material-output-coordinator.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    assert not (outer / 'active.json').exists()
    command = [sys.executable, '-u', str(Path(__file__).resolve()), '--run-directory', str(args.run_directory.resolve())]
    stage = 'coordinate_output_continue'; stdout = outer / (stage + '.stdout'); stderr = outer / (stage + '.stderr')
    record = {'stage': stage, 'command': command, 'startedUtc': datetime.now(timezone.utc).isoformat(),
              'stdout': str(stdout), 'stderr': str(stderr), 'status': 'running', 'scientificWorkersGuardedIndividually': True}
    started = time.monotonic()
    with stdout.open('xb') as out, stderr.open('xb') as err:
        process = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=err, start_new_session=True, close_fds=True)
        record['childPid'] = process.pid; save(outer / 'active.json', record)
        def stop(signum, frame): process.send_signal(signum); process.wait(); raise SystemExit(128 + signum)
        signal.signal(signal.SIGTERM, stop); signal.signal(signal.SIGINT, stop); code = process.wait()
    record.update(status='completed' if code == 0 else 'failed', exitCode=code, stderrBytes=stderr.stat().st_size,
                  wallSeconds=time.monotonic() - started, finishedUtc=datetime.now(timezone.utc).isoformat())
    save(outer / 'active.json', record); save(outer / (stage + '.invocation.json'), record); print(json.dumps(record, indent=2))
    return code


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--run-directory', type=Path, required=True); parser.add_argument('--supervise', action='store_true')
    args = parser.parse_args()
    if args.supervise: sys.exit(supervise(args))
    coordinate(args)
