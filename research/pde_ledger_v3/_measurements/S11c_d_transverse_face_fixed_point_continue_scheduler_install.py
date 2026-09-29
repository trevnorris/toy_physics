#!/usr/bin/python3
"""Install the explicitly approved one-job scheduler exception; no science work.

Run once with administrator authentication. Cleanup is a transient root service
whose complete code is passed by value, never loaded later from this workspace.
The existing worker remains stopped until a separate priority verification.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
from datetime import datetime, timezone

SOURCE = Path('/var/projects/toy_physics/research/pde_ledger_v3/_measurements/S11c_d_transverse_face_fixed_point_continue_scheduler_exception.kdl')
SHA256 = 'ec3a8b3f2514520c3a4f153bf7785777cc2ac995528ba529a411e06616ac7520'
GROUP = '/user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-37f533cfe2e7.service'
WORKER = 2497607
DIRECTORY = Path('/etc/system76-scheduler/process-scheduler')
DESTINATION = DIRECTORY / 's11c-37f533cfe2e7.kdl'
CLEANUP_UNIT = 's11c-priority-exception-cleanup-37f533cfe2e7'


def main():
    if os.geteuid() != 0:
        raise SystemExit('Administrator authentication is required.')
    payload = SOURCE.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == SHA256
    proc = Path('/proc') / str(WORKER)
    assert proc.stat().st_uid == 1000
    assert proc.joinpath('cgroup').read_text().strip() == '0::' + GROUP
    assert '\nState:\tT ' in proc.joinpath('status').read_text()
    cgroup = Path('/sys/fs/cgroup' + GROUP)
    inode = cgroup.stat().st_ino
    assert not DESTINATION.exists() and not DESTINATION.is_symlink()
    created = []
    for parent in (DIRECTORY.parent, DIRECTORY):
        if not parent.exists():
            parent.mkdir(mode=0o755)
            created.append(str(parent))
        assert not parent.is_symlink() and parent.stat().st_uid == 0
    cleanup = f'''import hashlib, pathlib, subprocess, time
p = pathlib.Path({str(DESTINATION)!r})
c = pathlib.Path({str(cgroup)!r})
while True:
    try:
        if c.stat().st_ino != {inode}:
            break
        if 'populated 0' in (c / 'cgroup.events').read_text():
            break
    except FileNotFoundError:
        break
    time.sleep(10)
if p.exists():
    if p.is_symlink() or hashlib.sha256(p.read_bytes()).hexdigest() != {SHA256!r}:
        raise SystemExit('Cleanup refused: exception file was changed.')
    p.unlink()
    subprocess.run(['/usr/bin/systemctl', 'reload', 'com.system76.Scheduler.service'], check=True)
for name in {list(reversed(created))!r}:
    try:
        pathlib.Path(name).rmdir()
    except OSError:
        pass
print('Exact-job scheduler exception removed.', flush=True)
'''
    fd = os.open(DESTINATION, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o644)
    with os.fdopen(fd, 'wb') as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    try:
        subprocess.run(['/usr/bin/systemd-run', '--quiet', '--collect',
                        '--unit=' + CLEANUP_UNIT, '--property=Restart=no',
                        '--property=MemoryMax=64M', '--property=TasksMax=8',
                        '/usr/bin/python3', '-I', '-c', cleanup], check=True)
        subprocess.run(['/usr/bin/systemctl', 'reload', 'com.system76.Scheduler.service'], check=True)
    except BaseException:
        if DESTINATION.exists() and hashlib.sha256(DESTINATION.read_bytes()).hexdigest() == SHA256:
            DESTINATION.unlink()
            subprocess.run(['/usr/bin/systemctl', 'reload', 'com.system76.Scheduler.service'], check=False)
        raise
    print(json.dumps({'utc': datetime.now(timezone.utc).isoformat(),
                      'installed': str(DESTINATION), 'sha256': SHA256,
                      'exactCgroup': GROUP, 'cleanupUnit': CLEANUP_UNIT,
                      'schedulerReloaded': True, 'workerResumed': False}, indent=2), flush=True)


if __name__ == '__main__':
    main()
