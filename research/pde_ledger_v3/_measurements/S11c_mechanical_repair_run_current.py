#!/usr/bin/env python3
"""Pin and measure the repaired reference-current check without publishing it."""
import hashlib
import json
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/tmp/s11c-mechanical-repair-20260912')


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def run():
    destination = BASE/'d_current'
    destination.mkdir(exist_ok=False)
    reference = BASE/'d_reference/manifest.json'
    paths = [Path(__file__).relative_to(ROOT),
             Path('_measurements/S11c_d_nonlocal_current_check.py'),
             Path('_measurements/S11c_d_joint_sheet_check.py'),
             Path('_measurements/S11c_d_channel_preflight_input.json'),
             Path('scripts/S11c_d_mixing_scattering_sympy_audit.py'),
             Path('scripts/S11c_d_output_codec.py'), Path('scripts/ledger_fold.py'),
             Path('directives/S11c_d_SHARED_PHYSICS.md'),
             Path('directives/S11b_SHARED_PHYSICS.md'),
             *(Path('scripts')/('S11c_'+stage+'_exports.py') for stage in ('b','c1','c2'))]
    pins = {str(path): digest(ROOT/path) for path in paths}
    for path in paths:
        target = destination/'source'/path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT/path, target)
    command = [sys.executable, '-u', '_measurements/S11c_d_nonlocal_current_check.py',
               '--manifest', str(reference),
               '--input', '_measurements/S11c_d_channel_preflight_input.json',
               '--acoustic', '--require-zero-residuals',
               '--cache-result', str(destination/'objects.pickle')]
    manifest = {'run_directory': str(destination), 'command': command,
                'scope': 'reduced LAB_HELD/RHO4_CONSTANT reference current and acoustic load',
                'source_hashes_before': pins, 'reference_manifest_sha256_before': digest(reference)}
    manifest_path = destination/'manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    started = time.monotonic()
    with (destination/'full.out').open('wb') as out, (destination/'stderr.txt').open('wb') as err:
        child = subprocess.Popen(command, cwd=ROOT, stdout=out, stderr=err)
        while child.poll() is None:
            try:
                child.wait(timeout=45)
            except subprocess.TimeoutExpired:
                print(json.dumps({'stage': 'd_current', 'wallSeconds': time.monotonic()-started,
                                  'stdoutBytes': out.tell(), 'stderrBytes': err.tell()}), flush=True)
    manifest.update(exit_code=child.returncode, wall_seconds=time.monotonic()-started,
                    peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                    source_hashes_after={str(path): digest(ROOT/path) for path in paths},
                    reference_manifest_sha256_after=digest(reference))
    manifest['artifacts'] = {path.name: {'bytes': path.stat().st_size, 'sha256': digest(path)}
                             for path in destination.iterdir() if path.is_file() and path != manifest_path}
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({key: value for key, value in manifest.items()
                      if key not in ('source_hashes_before','source_hashes_after','artifacts')}), flush=True)
    stable = (pins == manifest['source_hashes_after'] and
              manifest['reference_manifest_sha256_before'] == manifest['reference_manifest_sha256_after'])
    raise SystemExit(child.returncode or (0 if stable else 1))


if __name__ == '__main__':
    run()
