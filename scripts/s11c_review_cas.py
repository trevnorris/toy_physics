#!/usr/bin/env python3
"""Run a reviewer's own Python probe inside the S11c resource pool.

This dispatcher imports only the standard library. Each invocation uses a new
durable directory, preserves the script and output, and never retries a probe.
It is the sole computational command exposed to the external repair reviewers.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT / '_scratch' / 's11c'
GUARD = ROOT / 'scripts' / 's11c_guarded_run.py'
SUPERVISOR = ROOT / 'research/pde_ledger_v3/_measurements/S11c_d_end_normalization_run.py'


def validate_script(workspace, script):
    workspace = workspace.resolve(strict=True)
    workspace.relative_to(STORE)
    script = script.resolve(strict=True)
    script.relative_to(workspace)
    if not workspace.is_dir() or not script.is_file() or script.suffix != '.py':
        raise ValueError('a Python probe inside the assigned review workspace is required')
    return workspace, script


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--script', type=Path, required=True)
    args = parser.parse_args()
    workspace, script = validate_script(args.workspace, args.script)
    base = workspace / 'cas-runs' / uuid.uuid4().hex
    base.mkdir(parents=True, exist_ok=False)
    source = base / 'probe.py'
    source.write_bytes(script.read_bytes())
    # Execute the immutable copy. The submitted probe must use explicit paths
    # to any staged source operands; hidden imports from the repository are not
    # added by this dispatcher.
    execute = base / 'execute_probe.py'
    execute.write_text('import os, runpy, sys\n'
                       + 'os.chdir(' + repr(str(workspace)) + ')\n'
                       + 'sys.path[0] = ' + repr(str(workspace)) + '\n'
                       + 'runpy.run_path(' + repr(str(source)) + ', run_name="__main__")\n')
    command = [sys.executable, str(GUARD), '--pool', 's11c-review',
               '--memory-gib', '4', '--log-directory', str(base / 'resource-guard'), '--',
               sys.executable, str(SUPERVISOR), '--parallel-prerequisite-read',
               '--run-root', str(base / 'supervisor'), '--stage', 'review_probe', '--',
               sys.executable, str(execute)]
    record = {'startedUtc': datetime.now(timezone.utc).isoformat(),
              'script': str(script), 'scriptSha256': digest(script),
              'sourceCopy': str(source), 'executionWrapperSha256': digest(execute), 'command': command,
              'guardSha256': digest(GUARD), 'supervisorSha256': digest(SUPERVISOR),
              'wallDeadlineSeconds': None, 'retry': False}
    (base / 'invocation.json').write_text(json.dumps(record, indent=2) + '\n')
    with (base / 'launcher.stdout').open('xb') as out, (base / 'launcher.stderr').open('xb') as err:
        result = subprocess.run(command, cwd=workspace, stdin=subprocess.DEVNULL,
                                stdout=out, stderr=err)
    record.update(exitCode=result.returncode,
                  finishedUtc=datetime.now(timezone.utc).isoformat(),
                  originalScriptUnchanged=digest(script) == record['scriptSha256'],
                  sourceCopyUnchanged=digest(source) == record['scriptSha256'])
    (base / 'outcome.json').write_text(json.dumps(record, indent=2) + '\n')
    for label, name in [('stdout', 'review_probe.stdout'), ('stderr', 'review_probe.stderr')]:
        path = base / 'supervisor' / name
        print(label + ': ' + str(path))
        if path.exists():
            # Complete literal output stays on disk; avoid flooding the leg's
            # context with a large symbolic dump.
            raw = path.read_bytes()
            print(raw[:24000].decode('utf-8', errors='replace'))
            if len(raw) > 24000:
                print('[display truncated; complete output retained at the path above]')
    print(json.dumps({'runDirectory': str(base), **record}, indent=2))
    return result.returncode or (0 if record['sourceCopyUnchanged'] and record['originalScriptUnchanged'] else 1)


if __name__ == '__main__':
    sys.exit(main())
