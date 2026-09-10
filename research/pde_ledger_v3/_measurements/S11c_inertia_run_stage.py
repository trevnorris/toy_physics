"""Run one authorized rebuild stage, measure it, and install completed stdout."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
STAGES = {
    'c1': 'S11c_c1_bulk_closure_sympy_audit',
    'c2': 'S11c_c2_selfenergy_fold_sympy_audit',
    'd': 'S11c_d_mixing_scattering_sympy_audit',
    'diagnostic': 'S11c_d_transverse_sign_diagnostic',
}
parser = argparse.ArgumentParser()
parser.add_argument('stage', choices=STAGES)
args = parser.parse_args()
name = STAGES[args.stage]
scratch = Path('/tmp/s11c-inertia-rebuild')
scratch.mkdir(exist_ok=True)
out, err, metrics = (scratch / (args.stage + suffix) for suffix in ('.out', '.err', '.time'))
command = ['/usr/bin/time', '-v', '-o', str(metrics), sys.executable,
           '-u', str(ROOT / 'scripts' / (name + '.py'))]
started = time.monotonic()
with out.open('wb') as stdout, err.open('wb') as stderr:
    child = subprocess.Popen(command, cwd=ROOT, stdout=stdout, stderr=stderr)
    while child.poll() is None:
        try:
            child.wait(timeout=45)
        except subprocess.TimeoutExpired:
            print(json.dumps({'stage': args.stage, 'elapsed_seconds': round(time.monotonic()-started),
                              'stdout_bytes': out.stat().st_size, 'stderr_bytes': err.stat().st_size}), flush=True)
result = {'stage': args.stage, 'command': command, 'exit_code': child.returncode,
          'elapsed_seconds': time.monotonic()-started, 'stdout_bytes': out.stat().st_size,
          'stderr_bytes': err.stat().st_size, 'time_file': str(metrics)}
if child.returncode == 0 and out.stat().st_size:
    destination = ROOT / 'scripts/out' / (name + '.out')
    staged = destination.with_name(destination.name + '.inertia-staged')
    shutil.copyfile(out, staged)
    os.replace(staged, destination)
    result['stdout_path'] = str(destination)
(ROOT / '_measurements' / ('S11c_inertia_run_' + args.stage + '.json')).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result), flush=True)
raise SystemExit(child.returncode or (0 if 'stdout_path' in result else 1))
