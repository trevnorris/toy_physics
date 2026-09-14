#!/usr/bin/env python3
"""Inventory native producer/export artifacts and optionally install stdout.

Publication here checks execution and provenance. Physical comparison records
remain separate and must be inspected before choosing --publish.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shutil

from S11c_inertia_artifact_audit import export_data, locate
from S11c_thickness_coordinate_run_stage import STAGES

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('b', 'c1', 'c2'))
    parser.add_argument('--run-directory', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--physical-check', type=Path)
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    run = args.run_directory
    manifest_path = run / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    transcript = run / 'full.out'
    export = ROOT / 'scripts' / ('S11c_'+args.stage+'_exports.py')
    records = {}
    tags = []
    for line in transcript.open():
        if not line.startswith('PY_'):
            continue
        tag, separator, payload = line.partition(': ')
        if not separator:
            continue
        tags.append(tag)
        if tag.endswith(('_OPERATIONAL_EXCEPTIONS', '_RUN_TASKS', '_SKIPPED_TASKS', '_TASK_TIMING_SECONDS')):
            records[tag] = payload.rstrip()
    values, pins, imports = export_data(export)
    before, _, _ = export_data(args.baseline / 'scripts' / export.name)
    result = {'stage': args.stage, 'producerExitCode': manifest['exit_code'],
              'producerSourcePinsStable': manifest['source_hashes_before'] == manifest['source_hashes_after'],
              'currentSourcePinMismatches': [name for name, expected in manifest['source_hashes_after'].items()
                                           if digest(ROOT/name) != expected],
              'transcriptSha256': digest(transcript), 'transcriptBytes': transcript.stat().st_size,
              'exportSha256': digest(export), 'exportBytes': export.stat().st_size,
              'stderrBytes': (run/'stderr.txt').stat().st_size,
              'tagCount': len(tags), 'duplicateTags': [tag for tag, count in Counter(tags).items() if count > 1],
              'nativeExecutionRecords': records, 'exportRows': len(values), 'directImportKeys': list(imports),
              'changedValueSerializations': sorted(key for key in values.keys() & before.keys() if values[key] != before[key]),
              'addedKeys': sorted(values.keys()-before.keys()), 'removedKeys': sorted(before.keys()-values.keys()),
              'exportPinMismatches': [name for name, expected in pins.items() if digest(locate(name)) != expected],
              'published': False}
    result['transcriptMatchesProducerArtifact'] = result['transcriptSha256'] == manifest['artifacts']['full.out']['sha256']
    result['exportMatchesProducerArtifact'] = result['exportSha256'] == manifest['artifacts'][export.name]['sha256']
    if args.stage == 'b':
        _, _, c1_imports = export_data(args.baseline/'scripts/S11c_c1_exports.py')
        result['changedC1DirectInputs'] = sorted(set(c1_imports) & set(result['changedValueSerializations']))
    if args.stage == 'c2':
        evidence_path = run/'native_guard_evidence.json'
        evidence = json.loads(evidence_path.read_text())
        result['nativeGuardEvidenceSha256'] = digest(evidence_path)
        result['publicationSemanticResults'] = dict(Counter(line.rsplit('=', 1)[-1].strip()
                                                for line in evidence['export']['publication_semantic']))
        result['publicationSemanticFailures'] = [line for line in evidence['export']['publication_semantic']
            if line.rsplit('=',1)[-1].strip() != ('0' if ' expanded_difference = ' in line else 'True')]
        result['nativeDirectLookupsEqualImports'] = sorted(evidence['lookups']) == sorted(evidence['import_keys'])
    failures = [key for key in ('producerExitCode', 'currentSourcePinMismatches', 'stderrBytes', 'duplicateTags',
                               'exportPinMismatches') if result[key]]
    failures += [key for key in ('producerSourcePinsStable', 'transcriptMatchesProducerArtifact',
                                 'exportMatchesProducerArtifact') if not result[key]]
    failures += [tag for tag, value in records.items() if tag.endswith('_OPERATIONAL_EXCEPTIONS') and value != 'Tuple()']
    if args.stage == 'c2':
        if not result['nativeDirectLookupsEqualImports'] or not result['publicationSemanticResults'] or result['publicationSemanticFailures']:
            failures.append('nativeSemanticChecks')
    if args.stage == 'c1' and (result['changedValueSerializations'] or result['addedKeys'] or result['removedKeys']):
        failures.append('unexpectedC1ValueChanges')
    if args.physical_check:
        proof=json.loads(args.physical_check.read_text())
        result['physicalCheckSha256']=digest(args.physical_check)
        if proof['stage']!=args.stage or proof['newExportSha256']!=result['exportSha256'] or proof['nonzeroResidualScalars']:
            failures.append('physicalCheckMismatch')
    elif args.publish and args.stage in ('b','c2'):
        failures.append('physicalCheckRequired')
    result['operationalOrProvenanceFailures'] = failures
    if args.publish and not failures:
        destination = ROOT/'scripts/out'/(STAGES[args.stage]+'.out')
        import tempfile
        with tempfile.NamedTemporaryFile(dir=destination.parent,prefix='.coordinate-producer-',delete=False) as f:
            staging=Path(f.name)
            with transcript.open('rb') as source:shutil.copyfileobj(source,f)
            f.flush();os.fsync(f.fileno())
        if digest(staging)!=result['transcriptSha256']:raise ValueError('publication digest')
        os.replace(staging,destination)
        result['published'] = True
        result['publishedPath'] = str(destination.relative_to(ROOT))
    result['producerManifestSha256'] = digest(manifest_path)
    destination = ROOT/'_measurements'/('S11c_thickness_coordinate_'+args.stage+'_stage_inventory.json')
    destination.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    raise SystemExit(1 if failures else 0)


if __name__ == '__main__':
    run()
