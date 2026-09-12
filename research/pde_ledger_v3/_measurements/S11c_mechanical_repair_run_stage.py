#!/usr/bin/env python3
"""Source-pinned serial producer run; transcript publication is a separate step.

Native upstream producers write their own exports. This runner preserves the
completed native export and stdout in a unique run directory without replacing
an annex-managed transcript. Exit status records execution, not physics checks.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/tmp/s11c-mechanical-repair-20260912')
STAGES = {
    'b': 'S11c_b_brane_operator_sympy_audit',
    'c1': 'S11c_c1_bulk_closure_sympy_audit',
    'c2': 'S11c_c2_selfenergy_fold_sympy_audit',
    'd': 'S11c_d_mixing_scattering_sympy_audit',
}
INPUTS = {'b': ('a',), 'c1': ('b',), 'c2': ('b', 'c1'), 'd': ('b', 'c1', 'c2')}


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=STAGES)
    parser.add_argument('--run-name')
    parser.add_argument('producer_arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    name = STAGES[args.stage]
    destination = BASE / (args.run_name or args.stage + '_full')
    destination.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__).relative_to(ROOT), Path('scripts') / (name + '.py'),
               Path('scripts/ledger_fold.py'),
               Path('directives') / ('S11c_' + args.stage + '_SHARED_PHYSICS.md')]
    sources += [Path('scripts') / ('S11c_' + stage + '_exports.py') for stage in INPUTS[args.stage]]
    sources += [path.relative_to(ROOT) for path in (ROOT / 'directives').glob('S11c_' + args.stage + '_sympy_build*.md')]
    if args.stage == 'd':
        sources += [Path('scripts/S11c_d_output_codec.py'),
                    Path('_measurements/S11c_d_channel_preflight_input.json')]
    pins = {str(path): digest(ROOT / path) for path in sources}
    for path in sources:
        frozen = destination / 'source' / path
        frozen.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / path, frozen)
    command = [sys.executable, '-u', str(Path('scripts') / (name + '.py'))]
    command += args.producer_arguments[1:] if args.producer_arguments[:1] == ['--'] else args.producer_arguments
    environment = os.environ.copy()
    controls = {}
    if args.stage == 'b':
        controls = {'S11CB_PRIMARIES_ONLY': '1', 'S11CB_PROJECTION_WORKERS': '1'}
        environment.update(controls)
    started = time.monotonic()
    manifest = {'stage': args.stage, 'cwd': str(ROOT), 'run_directory': str(destination),
                'command': command, 'environment_overrides': controls, 'source_hashes_before': pins,
                'transcript_published': False}
    manifest_path = destination / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    with (destination / 'full.out').open('wb') as stdout, (destination / 'stderr.txt').open('wb') as stderr:
        child = subprocess.Popen(command, cwd=ROOT, env=environment, stdout=stdout, stderr=stderr)
        manifest['child_pid'] = child.pid
        manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
        while child.poll() is None:
            try:
                child.wait(timeout=45)
            except subprocess.TimeoutExpired:
                print(json.dumps({'stage': args.stage, 'wall_seconds': time.monotonic()-started,
                                  'stdout_bytes': stdout.tell(), 'stderr_bytes': stderr.tell()}), flush=True)
    manifest.update(exit_code=child.returncode, wall_seconds=time.monotonic()-started,
                    peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
                    source_hashes_after={str(path): digest(ROOT / path) for path in sources})
    export = ROOT / 'scripts' / ('S11c_' + args.stage + '_exports.py')
    if args.stage != 'd' and export.exists():
        shutil.copyfile(export, destination / export.name)
    manifest['artifacts'] = {str(path.relative_to(destination)): {'bytes': path.stat().st_size, 'sha256': digest(path)}
                             for path in destination.rglob('*') if path.is_file()
                             and 'source' not in path.relative_to(destination).parts and path != manifest_path}
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({key: value for key, value in manifest.items()
                      if key not in ('source_hashes_before', 'source_hashes_after', 'artifacts')}), flush=True)
    raise SystemExit(child.returncode or (0 if pins == manifest['source_hashes_after'] else 1))


if __name__ == '__main__':
    run()
