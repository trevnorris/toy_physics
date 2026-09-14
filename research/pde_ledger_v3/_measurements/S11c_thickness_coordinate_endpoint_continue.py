#!/usr/bin/env python3
"""Execute the prepared endpoint plan serially, with validated local checkpoints.

This controller constructs no physical objects. Native instruments emit them;
their validators check and publish them. A failed command stops the queue with
its output intact, including a pairing emitter's retained-discrepancy failure.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
PLAN = ROOT/'_measurements/S11c_thickness_coordinate_endpoint_plan.json'
STATUS = ROOT/'_measurements/S11c_thickness_coordinate_endpoint_status.md'
CONTRACT = b'## Retained user-approved solver/export contract'
CONTRACT_SHA = 'f01024512c9a028e7e524e5a95b5c455be42e125d10c36f94e20493109b013c2'


def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def write_json(path, value):
    temporary = path.with_name(path.name+'.new')
    with temporary.open('w') as stream:
        stream.write(json.dumps(value, indent=2)+'\n')
        stream.flush(); os.fsync(stream.fileno())
    temporary.replace(path)


def checked_file(path, record):
    if path.stat().st_size != record['bytes'] or sha(path) != record['sha256']:
        raise ValueError(('artifact mismatch', str(path)))


def stage_result(row):
    documents = {name: json.loads((ROOT/name).read_text()) for name in row['ordinary']}
    primary = documents[row['ordinary'][-1]]
    required = (('nonzeroResidualCount', 'failures') if row['name'] == 'reference_source'
                else ('retainedNonzeroScalars', 'nonzeroCancellationIdentities') if row['name'].endswith('_source')
                else ('retainedNonzeroScalars', 'retainedResidualScalars') if row['name'].endswith('_pairing')
                else ('summary', 'checks', 'normalReality'))
    if any(name not in primary for name in required):
        raise ValueError(('incomplete validator checkpoint', row['name']))
    for name in ('retainedNonzeroScalars', 'nonzeroCancellationIdentities', 'nonzeroResidualCount'):
        if primary.get(name, 0):
            raise ValueError(('recorded discrepancy', row['name'], name))
    if primary.get('failures'):
        raise ValueError(('recorded failures', row['name']))
    outputs = {}
    for name in row['annex']:
        expected = primary.get('publication')
        if not isinstance(expected, dict):
            expected = primary['artifacts']['full.out']
        checked_file(ROOT/name, expected)
        outputs[name] = {'bytes': (ROOT/name).stat().st_size, 'sha256': sha(ROOT/name)}
    keys = ('tagCount', 'objectCount', 'checkedSourceObjects', 'checkedMetadataPaths',
            'metadataPaths', 'literalSourceResidualScalars', 'residualScalarCount',
            'residualScalars', 'retainedResidualScalars', 'retainedNonzeroScalars',
            'cancellationIdentityCount', 'nonzeroCancellationIdentities',
            'nonzeroResidualCount', 'wallSeconds', 'peakRssKiB')
    return {'stage': row['name'], 'measurements': {k: primary[k] for k in keys if k in primary},
            'evidence': {name: sha(ROOT/name) for name in row['ordinary']}, 'outputs': outputs}


def committed_paths(row):
    paths = [str((ROOT/name).relative_to(REPO)) for name in row['ordinary']+row['annex']]
    for path in paths:
        subprocess.run(['git', 'ls-files', '--error-unmatch', '--', path], cwd=REPO,
                       stdout=subprocess.DEVNULL, check=True)
    if subprocess.check_output(['git', 'diff', 'HEAD', '--name-only', '--', *paths], cwd=REPO, text=True).strip():
        raise ValueError(('uncommitted predecessor', row['name']))
    for name in row['annex']:
        if not (ROOT/name).is_symlink():
            raise ValueError(('predecessor is not annexed', name))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--start-at', default='right_source')
    parser.add_argument('--state-directory', type=Path, required=True)
    parser.add_argument('--adopt-construction', action='store_true',
                        help='Wait for the selected stage\'s existing invocation; never repeat it.')
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    plan = json.loads(PLAN.read_text()); base = Path(plan['runRoot'])
    rows = plan['stages']; index = next(i for i, row in enumerate(rows) if row['name'] == args.start_at)
    state = args.state_directory.resolve(); state.relative_to(base)
    for row in rows[index:]:
        for step in row['steps']:
            if step['command'][0] != 'python' or not (ROOT/step['command'][1]).is_file():
                raise ValueError(('invalid instrument', step))
        for name in row['ordinary']+row['annex']:
            (ROOT/name).resolve().relative_to(REPO)
    if args.plan_only:
        print(json.dumps({'stages': [r['name'] for r in rows[index:]],
                          'adoptConstruction': args.adopt_construction, 'planSha256': sha(PLAN)}))
        return
    lock = (base/'endpoint_controller.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    state.mkdir(exist_ok=False)
    manifest_path = Path(plan['nativeProducerManifest'])
    manifest = json.loads(manifest_path.read_text())
    publication_path = ROOT/'_measurements/S11c_thickness_coordinate_d_publication.json'
    publication = json.loads(publication_path.read_text())
    checks_path = ROOT/publication['checksPath']; checks = json.loads(checks_path.read_text())
    if (manifest['exit_code'] != 0 or manifest['source_hashes_before'] != manifest['source_hashes_after']
            or checks['failures'] or sha(manifest_path) != publication['producerManifestSha256']
            or sha(checks_path) != publication['checksSha256']):
        raise ValueError('native producer/publication prerequisite')
    checked_file(ROOT/publication['path'], publication)
    if not (ROOT/publication['path']).is_symlink():
        raise ValueError('native transcript is not annexed')
    native_paths = [str(p.relative_to(REPO)) for p in (publication_path, checks_path, ROOT/publication['path'])]
    if subprocess.check_output(['git', 'diff', 'HEAD', '--name-only', '--', *native_paths], cwd=REPO, text=True).strip():
        raise ValueError('native publication is uncommitted')
    if sha(ROOT/plan['input']) != plan['inputSha256']:
        raise ValueError('supplied input changed')
    pins = {ROOT/name: pin for name, pin in manifest['source_hashes_after'].items()}
    files = [Path(__file__), PLAN, manifest_path, publication_path, checks_path, ROOT/plan['input']]
    files += [ROOT/name for name in plan['instrumentPathsChecked']]
    files += [ROOT/'_measurements'/name for name in (
        'S11c_d_end_pairing_check.py', 'S11c_d_end_pairing_workers.py',
        'S11c_d_current_runtime_metadata.py', 'S11c_d_joint_sheet_check.py',
        'S11c_d_modal_current_check.py', 'S11c_d_nonlocal_current_check.py')]
    for path in files:
        pins[path] = sha(path)
    def stable():
        if any(sha(path) != pin for path, pin in pins.items()):
            raise ValueError('endpoint plan/source changed')
        report = (ROOT/'_measurements/S11c_d_sympy_builder_report.md').read_bytes()
        if hashlib.sha256(report[report.index(CONTRACT):]).hexdigest() != CONTRACT_SHA:
            raise ValueError('retained solver/export contract changed')
    stable()
    write_json(state/'source_pins.json', {str(p): v for p, v in pins.items()})
    write_json(state/'plan.json', rows[index:])
    completed = []
    for row in rows[:index]:
        committed_paths(row); completed.append(stage_result(row))
    events = []
    pointer = {'controllerPid': os.getpid(), 'stateDirectory': str(state),
               'planSha256': sha(PLAN), 'controllerSha256': sha(Path(__file__))}
    def progress(event):
        event = {**event, 'utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
        events.append(event); write_json(state/'progress.json', events)
        write_json(base/'active_endpoint_controller.json', {**pointer, 'latest': event})
        print(json.dumps(event), flush=True)
    resource.setrlimit(resource.RLIMIT_STACK, (256*1024*1024, resource.getrlimit(resource.RLIMIT_STACK)[1]))
    progress({'event': 'started', 'stage': args.start_at})
    for row in rows[index:]:
        stable()
        if any(name != 'validated_native_regeneration' and name not in {v['stage'] for v in completed}
               for name in row['requires']):
            raise ValueError(('unfinished prerequisite', row['name']))
        run = Path(row.get('createRunDirectoryBefore', row.get('runnerCreatesFreshDirectory')))
        for name in row['ordinary']+row['annex']:
            if (ROOT/name).exists() or (ROOT/name).is_symlink():
                raise ValueError(('stage publication already exists', name))
        steps = row['steps']
        if row['name'] == args.start_at and args.adopt_construction:
            invocation = run/'invocation.json'
            progress({'event': 'waiting_for_existing_construction', 'stage': row['name'], 'invocation': str(invocation)})
            while True:
                prior = json.loads(invocation.read_text())
                if prior['planSha256'] != sha(PLAN) or prior['command'] != steps[0]['command'] or prior['stage'] != row['name']:
                    raise ValueError('existing construction invocation differs')
                if prior['status'] != 'running':
                    break
                os.kill(prior['supervisorPid'], 0)
                time.sleep(15)
            if prior['status'] != 'construction_complete' or prior['exitCode'] or prior['stderrBytes']:
                raise ValueError('existing construction failed; inspect preserved output')
            write_json(state/(row['name']+'_adopted.json'), prior)
            steps = steps[1:]
        elif 'createRunDirectoryBefore' in row:
            run.mkdir(exist_ok=False)
        for step in steps:
            stable()
            stem = row['name']+'_'+step['name']
            out = Path(step.get('stdout', str(state/(stem+'.stdout'))))
            err = Path(step.get('stderr', str(state/(stem+'.stderr'))))
            progress({'event': 'operation_started', 'stage': row['name'], 'operation': step['name'],
                      'command': step['command'], 'stdout': str(out), 'stderr': str(err), 'runDirectory': str(run)})
            with out.open('xb') as stdout, err.open('xb') as stderr:
                child = subprocess.Popen(step['command'], cwd=ROOT, stdout=stdout, stderr=stderr)
                progress({'event': 'child_running', 'stage': row['name'], 'operation': step['name'], 'childPid': child.pid})
                child.wait()
            progress({'event': 'operation_exited', 'stage': row['name'], 'operation': step['name'],
                      'exitCode': child.returncode, 'stdoutSha256': sha(out), 'stderrBytes': err.stat().st_size,
                      'stderrSha256': sha(err)})
            if child.returncode or err.stat().st_size:
                raise ValueError(('endpoint operation failed; inspect emitted operands', stem))
        stable(); result = stage_result(row)
        completed.append(result)
        write_json(state/(row['name']+'_validated.json'), result)
        report = '# S11c thickness-coordinate endpoint checkpoints\n\n'
        report += 'Fresh calculations against the regenerated native producer. Each stage uses its own validator; numerical residuals retain their recorded domains.\n\n'
        for item in completed:
            report += '## '+item['stage']+'\n\n```json\n'+json.dumps(item['measurements'], indent=2)+'\n```\n\n'
        report += 'Remaining planned stages: '+', '.join(r['name'] for r in rows if r['name'] not in {v['stage'] for v in completed})+'.\n\n'
        report += 'Full endpoint current/adjoint maps and outward orientations precede variable-profile matching. No global exceptional coverage, complete scattering or profile-frequency bound-pole claim.\n'
        STATUS.write_text(report)
        ordinary = [ROOT/name for name in row['ordinary']]+[STATUS]
        annex = [ROOT/name for name in row['annex']]
        if subprocess.check_output(['git', 'diff', '--cached', '--name-only'], cwd=REPO, text=True).strip():
            raise ValueError('preexisting staged changes require inspection')
        message = state/(row['name']+'_commit.txt')
        message.write_text('S11c thickness coordinate: validate and publish fresh '+row['name']+'\n\n'
            'COMPUTED CHECKPOINT\n'+json.dumps(result, indent=2)+'\n\n'
            'Execute the prepared source/frequency/pairing stage against the regenerated d producer. Preserve emitted operands, raw/retained/remainder distinctions, numerical residuals and complete computed subspaces. Validator and command exits succeeded before publication.\n\n'
            'STORAGE AND BOUNDARIES\nOrdinary evidence/status are in Git; .out files use DataLad/git-annex with exact post-save content checks. Full endpoint current/adjoint normalization and variable-profile matching remain subsequent work. This is not global spectral coverage, complete scattering or a section 3b frequency-bound-pole result. No S10/Lean/authority changes, review/comparator/Wolfram, downstream run or push.\n')
        subprocess.run(['git', 'add', '--', *[str(p.relative_to(REPO)) for p in ordinary]], cwd=REPO, check=True)
        subprocess.run(['datalad', 'save', '-F', str(message), '--',
                        *[str(p.relative_to(REPO)) for p in ordinary+annex]], cwd=REPO, check=True)
        for path in annex:
            if not path.is_symlink():
                raise ValueError('published output is not an annex pointer')
            checked_file(path, result['outputs'][str(path.relative_to(ROOT))])
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip()
        progress({'event': 'committed', 'stage': row['name'], 'commit': commit, 'result': result})
    progress({'event': 'endpoint_plan_complete', 'remaining': plan['remainingAfterPlan']})


if __name__ == '__main__':
    main()
