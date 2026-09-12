#!/usr/bin/env python3
"""Run the prepared, read-only d inventories on a completed native producer.

This records execution and immutable inputs. Physical coverage and residuals
remain the outputs of the existing inventory instruments, inspected separately.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / '_measurements/S11c_mechanical_repair_d_recheck_plan.json'


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def run():
    plan = json.loads(PLAN.read_text())
    producer_path = Path(plan['nativeProducerManifest'])
    producer = json.loads(producer_path.read_text())
    pins = producer['source_hashes_before']
    current = {path: digest(ROOT / path) for path in pins}
    if (producer.get('exit_code') != 0 or pins != producer.get('source_hashes_after')
            or pins != current):
        raise RuntimeError('completed source-stable native producer required')
    base = Path(producer['run_directory'])
    for name, record in producer['artifacts'].items():
        if digest(base / name) != record['sha256']:
            raise RuntimeError(('producer artifact changed', name))
    if (base / 'stderr.txt').stat().st_size:
        raise RuntimeError('inspect native stderr before rechecking')
    destination = base / 'recheck'
    destination.mkdir(exist_ok=False)
    instruments = {Path(stage['command'][1]) for stage in plan['stages']}
    instruments.update(path.relative_to(ROOT) for path in
                       (ROOT / '_measurements').glob('S11c_d_*inventory.py'))
    instruments.update((Path(__file__).relative_to(ROOT), PLAN.relative_to(ROOT)))
    instrument_pins = {str(path): digest(ROOT / path) for path in sorted(instruments)}
    for path in instruments:
        frozen = destination / 'source' / path
        frozen.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / path, frozen)
    report = {'producerManifest': str(producer_path),
              'producerManifestSha256': digest(producer_path),
              'producerSourceHashes': pins, 'instrumentHashesBefore': instrument_pins,
              'stages': [], 'scope': plan['scope']}
    report_path = ROOT / '_measurements/S11c_mechanical_repair_d_rechecks.json'
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    for stage in plan['stages']:
        name = stage['name']
        output = destination / (name + '.stdout')
        errors = destination / (name + '.stderr')
        result_path = ROOT / stage['result']
        if result_path.exists():
            raise RuntimeError(('inventory result already exists', str(result_path)))
        started = time.monotonic()
        with output.open('xb') as stdout, errors.open('xb') as stderr:
            process = subprocess.Popen(stage['command'], cwd=ROOT, stdout=stdout, stderr=stderr)
            while process.poll() is None:
                try:
                    process.wait(timeout=45)
                except subprocess.TimeoutExpired:
                    print(json.dumps({'inventory': name, 'wallSeconds': time.monotonic()-started,
                                      'stdoutBytes': stdout.tell(), 'stderrBytes': stderr.tell()}), flush=True)
        if process.returncode == 0 and stage.get('captureStdoutAsResult'):
            value = json.loads(output.read_text())
            with result_path.open('x') as stream:
                stream.write(json.dumps(value, indent=2) + '\n')
        record = {'name': name, 'command': stage['command'], 'exitCode': process.returncode,
                  'wallSeconds': time.monotonic()-started, 'stderrBytes': errors.stat().st_size,
                  'stdout': str(output), 'stdoutSha256': digest(output),
                  'stderr': str(errors), 'stderrSha256': digest(errors)}
        if result_path.exists():
            record['result'] = {'path': str(result_path), 'bytes': result_path.stat().st_size,
                                'sha256': digest(result_path)}
        report['stages'].append(record)
        report_path.write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps(record), flush=True)
        if process.returncode or errors.stat().st_size:
            raise SystemExit(process.returncode or 1)
    report['instrumentHashesAfter'] = {path: digest(ROOT / path) for path in instrument_pins}
    report['producerSourceHashesAfter'] = {path: digest(ROOT / path) for path in pins}
    report['producerManifestSha256After'] = digest(producer_path)
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    if (report['instrumentHashesAfter'] != instrument_pins or report['producerSourceHashesAfter'] != pins
            or report['producerManifestSha256After'] != report['producerManifestSha256']):
        raise RuntimeError('recheck source pin changed')


if __name__ == '__main__':
    run()
