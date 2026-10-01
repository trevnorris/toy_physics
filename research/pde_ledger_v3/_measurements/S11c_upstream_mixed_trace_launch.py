#!/usr/bin/env python3
"""Launch the authorized selected response trace once, after its completion hook is armed."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import select
import shutil
import subprocess
import sys
import time
import traceback

ROOT = Path('/var/projects/toy_physics')
M = ROOT / 'research/pde_ledger_v3/_measurements'
PREFIX = 'S11c_upstream_mixed_trace'
GATE = M / (PREFIX + '_gate.json')
MANIFEST = M / (PREFIX + '_inputs.json')
RUN = Path('/var/projects/toy_physics/_scratch/s11c/s11c-mixed-trace-20261001/diagnostic-01')
HOOK = ROOT / 'scripts/codex_job_watch.py'
THREAD = '01a0e01b-ef84-7192-817f-584cda5d339b'
MESSAGE = 'The user-authorized selected mixed-term closed-response diagnostic has completed or reported an error.\nRoot: /var/projects/toy_physics/_scratch/s11c/s11c-mixed-trace-20261001/diagnostic-01.\nInspect coordinator/launch/stage-start, shared guard invocation/duration/enforced resources/outcomes/logs, supervisor upstream_mixed_trace invocation/active, strict scientific stderr, stdout/checks identity, failure/incomplete operation, all journal inputs/returns and exact source/closure/reference/domain/control JSON, snapshots and every source/prior posthash. Process exit is not science acceptance. Completion inspection is JSON/source/hash/opaque bytes only; no scientific restoration outside containment. Preserve failure and stop; no automatic retry or new job.\nBoth literal reports remain Claude NEEDS REVISION/Grok NEEDS REVISION, preserved with exact reviewed bytes at8305e075 and sibling build-review; packetd23734d7d3c61e058a12bf6820c45cef95b783aadb17d3c2af9a84d7651fb86c. Both supported the selected mathematics; the common blocker was leftmap/rightmap versus actual native left_map/right_map. Local namespace fix, evidence-before-guard reordering and per-unit-V labeling are authorized tooling corrections.11stdlib tests pass including actual fragment execution on stand-ins; all assignment ASTs and non-science top-level definitions match the reviewed version. No fresh independent CLEAR or rerun; gate records independentBuildClearance=false. Grok11969-byte stderr was reviewer CLI diagnostics, not scientific stderr.\nThe diagnostic uses the saved whole bare eta*sigma convolution once in the upper-face native closed-response comparison, then derives reference pressure and normal jet. Original first-shape composition remains once. Inspect actual native source joins, saved physical input/profile/units, both external factors, middle-independence/domain certificates, original closure/trace residuals, addressed omission/sign/one-sided controls and copied evidence. All reported factors/integrands are per unit prescribed normal velocity; chemical response and transverse excitation are not supplied. The native reference and jet are algebraically dependent outputs, not extra independent physics evidence. An already integrated convolution must not later pass unchanged through the native second-slot middle integral.\nSame strict-rest-bulk LAB_HELD/RHO4_CONSTANT,omega3,cs10,saved tangents,profile kin0/kout1/10. No boundary/profile/integral replay, slab-mode/b-row contraction, production repair, producer or defect sweep. Prior bare result and all history remain unchanged. Stop with selected response evidence and disposition; no leakage/loss/primitive calibration/drain/fullGreen/FORM/A11/A12 claim.\nNo deadline of any kind. Shared guard around existing normalization supervisor,4GiB native/cgroup within16GiB pool,zero swap,one assignedCPU/thread,32tasks,4GiBhost reserve,desktop priority,RuntimeMaxUSec=infinity,Restart=no. No scheduler change,fallback or retry. Keep Lean/S11_lean,shared guard,protected builder suffixf01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2 and all accepted/failed/incident/review history unchanged. Canonical records outside scratch;scratch never committed. Hook targets session01a0e01b-ef84-7192-817f-584cda5d339b;no model polling or recurring task.\n'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())


def verify():
    gate, manifest = read(GATE), read(MANIFEST)
    assert gate['status'] == 'READY_FOR_SELECTED_MIXED_TRACE_DIAGNOSTIC'
    assert gate['independentBuildClearance'] is False and gate['pooledExecution'] is True
    assert gate['launcherSha256'] == sha(__file__)
    assert gate['workerSha256'] == sha(M / (PREFIX + '.py'))
    assert gate['manifestSha256'] == sha(MANIFEST)
    assert gate['completionMessageSha256'] == hashlib.sha256(MESSAGE.encode()).hexdigest()
    assert gate['outputDirectory'] == str(RUN / 'complete')
    assert gate['sourcePins'] == manifest['sourcePins']
    for path, digest in gate['sourcePins'].items():
        assert sha(path) == digest, path
    review = read(M / (PREFIX + '_review_record.json'))
    repair = read(M / (PREFIX + '_repair_record.json'))
    assert review['allMetadataChecksPassed'] is True and review['noScienceRun'] is True
    assert review['reports']['claude']['literalVerdict'] == 'NEEDS REVISION'
    assert review['reports']['grok']['literalVerdict'] == 'NEEDS REVISION'
    assert repair['status'] == 'LOCAL_TOOLING_REPAIRS_TESTED_NO_FRESH_INDEPENDENT_CLEAR'
    assert repair['workerSha256'] == gate['workerSha256'] and repair['tests']['exitCode'] == 0
    assert repair['independentBuildClearance'] is False
    assert gate['reviewRecordSha256'] == sha(M / (PREFIX + '_review_record.json'))
    assert gate['repairRecordSha256'] == sha(M / (PREFIX + '_repair_record.json'))
    assert manifest['resources']['memoryBytes'] == 4*1024**3
    assert manifest['resources']['durationLimits'] is None
    assert gate['command'] == [sys.executable, str(ROOT / 'scripts/s11c_guarded_run.py'),
        '--pool', 's11c-near-unity', '--memory-gib', '4', '--tasks-max', '32',
        '--log-directory', str(RUN / 'resource-guard'), '--', sys.executable,
        str(M / 'S11c_d_end_normalization_run.py'), '--parallel-prerequisite-read',
        '--run-root', str(RUN), '--stage', 'upstream_mixed_trace', '--', sys.executable,
        '-u', str(M / (PREFIX + '.py')), '--out', str(RUN / 'complete'),
        '--inputs', str(MANIFEST)]
    return gate, manifest


def coordinate(descriptor):
    try:
        # Startup handshake only; this is not a computation deadline.
        ready, _, _ = select.select([descriptor], [], [], 30)
        if not ready or os.read(descriptor, 1) != b'1':
            raise RuntimeError('Completion hook not armed; refusing scientific launch')
        os.close(descriptor)
        gate, _ = verify()
        start = time.monotonic()
        save(RUN / 'stage-start.json', {'utc': datetime.now(timezone.utc).isoformat(),
             'gateSha256': sha(GATE), 'command': gate['command']})
        with (RUN / 'guard-launch.stdout').open('x') as out, (RUN / 'guard-launch.stderr').open('x') as err:
            result = subprocess.run(gate['command'], cwd=ROOT, stdin=subprocess.DEVNULL,
                                    stdout=out, stderr=err)
        save(RUN / 'coordinator-outcome.json', {'exitCode': result.returncode,
             'wallSeconds': time.monotonic() - start, 'scientificAcceptance': 'PENDING_ACTUAL_INSPECTION'})
        return result.returncode
    except BaseException:
        save(RUN / 'coordinator-failure.json', {'traceback': traceback.format_exc(), 'automaticRetry': False})
        return 1


def launch():
    gate, _ = verify()
    assert shutil.which('codex'), 'completion queue unavailable'
    RUN.mkdir(exist_ok=False)
    copies = {}
    for source, digest in {**gate['sourcePins'], str(GATE): sha(GATE), str(MANIFEST): sha(MANIFEST)}.items():
        path = Path(source)
        target = RUN / 'source' / path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        assert sha(target) == digest
        copies[source] = {'path': str(target), 'sha256': digest, 'bytes': target.stat().st_size}
    save(RUN / 'source-snapshots.json', copies)
    (RUN / 'completion-message.txt').write_text(MESSAGE)
    reader, writer = os.pipe()
    with (RUN / 'coordinator.stdout').open('x') as out, (RUN / 'coordinator.stderr').open('x') as err:
        job = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '--coordinator', str(reader)],
              cwd=ROOT, stdin=subprocess.DEVNULL, stdout=out, stderr=err, pass_fds=(reader,), start_new_session=True)
    os.close(reader)
    try:
        with (RUN / 'watcher.stdout').open('x') as out, (RUN / 'watcher.stderr').open('x') as err:
            watcher = subprocess.Popen([sys.executable, str(HOOK), '--pid', str(job.pid),
                '--thread', THREAD, '--directory', str(RUN / 'completion-watcher'),
                '--message-file', str(RUN / 'completion-message.txt')], cwd=ROOT,
                stdin=subprocess.DEVNULL, stdout=out, stderr=err, start_new_session=True)
        save(RUN / 'launch.json', {'jobPid': job.pid, 'watcherPid': watcher.pid,
             'gateSha256': sha(GATE), 'command': gate['command'], 'thread': THREAD, 'automaticRetry': False})
        for _ in range(50):
            path = RUN / 'completion-watcher/state.json'
            if path.exists() and read(path)['status'] == 'waiting':
                os.write(writer, b'1')
                print(json.dumps({'runDirectory': str(RUN), 'jobPid': job.pid,
                      'watcherPid': watcher.pid, 'hookStatus': 'waiting'}))
                return
            if job.poll() is not None or watcher.poll() is not None:
                raise RuntimeError('Coordinator/watcher exited before arming')
            time.sleep(.1)
        raise RuntimeError('Hook did not arm')
    finally:
        os.close(writer)


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--coordinator':
        sys.exit(coordinate(int(sys.argv[2])))
    if len(sys.argv) != 1:
        raise ValueError('Unexpected launch arguments')
    launch()
