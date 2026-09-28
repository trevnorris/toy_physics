"""Launch only the explicitly approved review packet, with a local event hook."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path('/var/projects/toy_physics')
RUN = ROOT / '_scratch/s11c/s11c-d-localized-thickness-20260928/build-review'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    state = json.loads((RUN/'state.json').read_text())
    approval = json.loads((RUN/'explicit-packet-approval.json').read_text())
    static = json.loads((RUN/'static-review.json').read_text())
    assert sha(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_localized_thickness_review_run.py') == static['runnerSha256']
    assert sha(Path(__file__)) == static['launcherSha256']
    assert sha(RUN/'completion_message.txt') == static['completionMessageSha256']
    assert state['status'] == 'PREPARED_NOT_LAUNCHED'
    assert approval['status'] == 'EXPLICITLY_APPROVED_FIXED_PACKET'
    assert approval['packetSha256'] == state['packetSha256']
    assert approval['archiveSha256'] == state['archiveSha256']
    assert approval['reviewers'] == ['claude', 'grok']
    assert not any((RUN/n).exists() for n in ('launch.json', 'start-authorized', 'coordinator.stdout'))
    assert state['currentThread'] == '01a0e01b-ef84-7192-817f-584cda5d339b'
    for name in ('claude', 'grok', 'codex'):
        assert shutil.which(name), 'Missing required executable: '+name
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'}
    with (RUN/'coordinator.stdout').open('x') as out, (RUN/'coordinator.stderr').open('x') as err:
        job = subprocess.Popen([sys.executable, str(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_localized_thickness_review_run.py')], cwd=ROOT,
            stdin=subprocess.DEVNULL, stdout=out, stderr=err, env=env, start_new_session=True)
    try:
        with (RUN/'watcher.stdout').open('x') as out, (RUN/'watcher.stderr').open('x') as err:
            watcher = subprocess.Popen([sys.executable, str(ROOT/'scripts/codex_job_watch.py'),
                '--pid', str(job.pid), '--thread', state['currentThread'],
                '--directory', str(RUN/'completion-watcher'),
                '--message-file', str(RUN/'completion_message.txt'),
                '--error-log', str(RUN/'issues.log')], cwd=ROOT,
                stdin=subprocess.DEVNULL, stdout=out, stderr=err, env=env,
                start_new_session=True)
        record = {'jobPid': job.pid, 'watcherPid': watcher.pid,
            'packetSha256': state['packetSha256'], 'reviewers': ['Claude', 'Grok'],
            'concurrentReviewWorkers': 1, 'scienceCalls': 0,
            'runnerSha256': sha(ROOT/'research/pde_ledger_v3/_measurements/S11c_d_localized_thickness_review_run.py'), 'launcherSha256': sha(Path(__file__))}
        with (RUN/'launch.json').open('x') as f:
            json.dump(record, f, indent=2)
            f.write('\n')
        for _ in range(30):
            path = RUN/'completion-watcher/state.json'
            if path.exists() and json.loads(path.read_text())['status'] == 'waiting':
                (RUN/'start-authorized').write_text('Fixed packet explicitly approved; completion hook armed.\n')
                print(json.dumps({'launch': record, 'hookStatus': 'waiting'}))
                return
            if watcher.poll() is not None or job.poll() is not None:
                raise RuntimeError('Review coordinator or hook exited before arming')
            time.sleep(.1)
        raise RuntimeError('Completion hook did not arm')
    except BaseException:
        if job.poll() is None:
            os.killpg(job.pid, signal.SIGTERM)
        raise


if __name__ == '__main__':
    main()
