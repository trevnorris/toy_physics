"""Sequential fixed-packet, read-only reviews; explicit submission gate."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import traceback
import uuid

ROOT = Path('/var/projects/toy_physics')
RUN = ROOT / '_scratch/s11c/s11c-d-localized-thickness-20260928/build-review'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix+'.new')
    with temporary.open('w') as f:
        json.dump(value, f, indent=2)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())
    temporary.replace(path)


def now():
    return datetime.now(timezone.utc).isoformat()


def verify(state):
    packet = Path(state['packetDirectory'])
    assert sorted(str(p.relative_to(packet)) for p in packet.rglob('*') if p.is_file()) == sorted(state['fileHashes'])
    assert all(sha(packet/n) == h for n, h in state['fileHashes'].items())
    assert hashlib.sha256(json.dumps(state['fileHashes'], sort_keys=True).encode()).hexdigest() == state['packetSha256']
    assert sha(RUN/'packet.tar.gz') == state['archiveSha256']
    assert all(sha(ROOT/r['path']) == r['sourceSha256'] for r in state['sourceRecords'])
    approval = json.loads((RUN/'explicit-packet-approval.json').read_text())
    assert approval['status'] == 'EXPLICITLY_APPROVED_FIXED_PACKET'
    assert approval['packetSha256'] == state['packetSha256']
    assert approval['archiveSha256'] == state['archiveSha256']
    assert approval['reviewers'] == ['claude', 'grok']


def main():
    state = json.loads((RUN/'state.json').read_text())
    assert state['status'] == 'PREPARED_NOT_LAUNCHED'
    verify(state)
    deadline = time.monotonic()+30
    while not (RUN/'start-authorized').exists():
        if time.monotonic() >= deadline:
            raise RuntimeError('No armed completion hook; no review launched')
        time.sleep(.1)
    packet = Path(state['packetDirectory'])
    prompt_path = packet/'review-prompt.md'
    state.update(status='RUNNING', startedUtc=now(), supervisorPid=os.getpid(),
                 concurrentReviewWorkers=1, scienceCalls=0)
    save(RUN/'state.json', state)
    for engine in state['reviewers']:
        verify(state)
        executable = shutil.which(engine)
        if not executable:
            raise RuntimeError('Missing reviewer executable: '+engine)
        session = str(uuid.uuid4())
        if engine == 'claude':
            command = [executable, '--print', '--output-format', 'json', '--safe-mode',
                '--restricted', '--tools', 'Read,Glob,Grep', '--allowedTools', 'Read,Glob,Grep',
                '--permission-mode', 'dontAsk', '--strict-mcp-config', '--mcp-config',
                '{"mcpServers":{}}', '--setting-sources', '', '--disable-slash-commands',
                '--session-id', session]
        else:
            command = [executable, '--cwd', str(packet), '--prompt-file', str(prompt_path),
                '--output-format', 'json', '--permission-mode', 'dontAsk', '--no-plan',
                '--tools', 'Read', '--no-subagents', '--disable-web-search', '--verbatim',
                '--session-id', session]
        paths = [RUN/(engine+suffix) for suffix in ('.json', '.stderr', '-run.json')]
        assert not any(p.exists() for p in paths)
        out_path, error_path, record_path = paths
        record = {'engine': engine, 'sessionId': session, 'command': command,
            'cwd': str(packet), 'startedUtc': now(), 'packetSha256': state['packetSha256'],
            'promptSha256': sha(prompt_path), 'scientificExecutionEnabled': False,
            'reviewVerdict': 'NOT_ADJUDICATED'}
        save(record_path, record)
        started = time.monotonic()
        with out_path.open('x') as out, error_path.open('x') as err:
            child = subprocess.Popen(command, cwd=packet,
                stdin=subprocess.PIPE if engine == 'claude' else subprocess.DEVNULL,
                stdout=out, stderr=err, text=True, start_new_session=True,
                env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
            record['pid'] = child.pid
            save(record_path, record)
            # Plain local process wait: no model polling and no automatic retry.
            child.communicate(input=prompt_path.read_text() if engine == 'claude' else None)
        record.update(exitStatus=child.returncode, finishedUtc=now(),
            wallSeconds=time.monotonic()-started, outputSha256=sha(out_path),
            stderrSha256=sha(error_path), stderrBytes=error_path.stat().st_size)
        save(record_path, record)
        verify(state)
        if child.returncode:
            with (RUN/'issues.log').open('a') as f:
                f.write(engine+' exited nonzero; inspect the actual record.\n')
    state.update(status='AWAITING_ADJUDICATION', finishedUtc=now(),
                 productionAuthorized=False, reviewVerdict='NOT_ADJUDICATED')
    save(RUN/'state.json', state)
    save(RUN/'coordinator-outcome.json', {'status':'BOTH_ATTEMPTS_FINISHED_AWAITING_LITERAL_ADJUDICATION', 'finishedUtc':now(), 'scienceCalls':0})


if __name__ == '__main__':
    try:
        main()
    except BaseException:
        save(RUN/'coordinator-outcome.json', {'status':'COORDINATOR_ERROR_NO_AUTOMATIC_RETRY', 'finishedUtc':now(), 'traceback':traceback.format_exc()})
        with (RUN/'issues.log').open('a') as f:
            f.write(traceback.format_exc())
        raise
