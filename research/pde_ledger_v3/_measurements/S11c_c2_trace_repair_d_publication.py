#!/usr/bin/env python3
"""Inventory the completed native run and read-only checks; publish atomically."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MEASUREMENTS = ROOT / '_measurements'
PREFIX = 'S11c_c2_trace_repair_d_'


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--publish', action='store_true')
    args = parser.parse_args()
    plan_path = MEASUREMENTS / (PREFIX + 'recheck_plan.json')
    plan = json.loads(plan_path.read_text())
    manifest_path = Path(plan['nativeProducerManifest'])
    manifest = json.loads(manifest_path.read_text())
    base = Path(manifest['run_directory'])
    checks_path = MEASUREMENTS / (PREFIX + 'rechecks.json')
    checks = json.loads(checks_path.read_text())
    failures = []

    def require(condition, name):
        if not condition:
            failures.append(name)

    pins = manifest['source_hashes_before']
    current_pins = {name: digest(ROOT / name) for name in pins}
    require(manifest.get('exit_code') == 0, 'native exit')
    require(pins == manifest.get('source_hashes_after') == current_pins, 'native sources')
    require((base / 'stderr.txt').stat().st_size == 0, 'native stderr')
    require(checks['producerManifestSha256'] == checks.get('producerManifestSha256After')
            == digest(manifest_path), 'inventory producer pin')
    require(checks['instrumentHashesBefore'] == checks.get('instrumentHashesAfter'), 'inventory source stability')
    require(checks['instrumentHashesBefore'] == {
        name: digest(ROOT / name) for name in checks['instrumentHashesBefore']}, 'inventory current sources')
    require([v['name'] for v in checks['stages']] == [v['name'] for v in plan['stages']], 'inventory stage census')
    for name, artifact in manifest['artifacts'].items():
        require(digest(base / name) == artifact['sha256'], 'native artifact: ' + name)
    inventories = {}
    for stage in checks['stages']:
        require(stage['exitCode'] == 0 and stage['stderrBytes'] == 0, 'inventory execution: ' + stage['name'])
        for key in ('stdout', 'stderr'):
            require(digest(Path(stage[key])) == stage[key + 'Sha256'], 'inventory log: ' + stage['name'] + '/' + key)
        result = stage['result']
        require(digest(Path(result['path'])) == result['sha256'], 'inventory result: ' + stage['name'])
        if stage['name'] != 'codec_expand':
            value = json.loads(Path(result['path']).read_text())
            inventories[stage['name']] = value['results'][0] if 'results' in value else value
    for name, value in inventories.items():
        for key in ('duplicates', 'duplicate_tags', 'missingMetadata', 'metadataGaps',
                    'unmatchedResidualMetadata', 'unmatchedObjects', 'nonfiniteObjects',
                    'payloadMismatches', 'rectangular_coverage_gaps', 'changed_source_pins', 'zero_map_metadata'):
            if key in value:
                require(not value[key], name + '/' + key)
        for key in ('dimensionConstraints', 'dimensionalConstraints', 'dimension_records'):
            for tag, record in value.get(key, {}).items():
                # INPUT_DIMENSION_CONSTRAINTS is the unsolved equation set;
                # only the later residual/unresolved records test closure.
                if tag != 'INPUT_DIMENSION_CONSTRAINTS':
                    require(record in ([], (), 'Tuple()'), name + '/' + tag)
        if 'completionMarkers' in value:
            require(value['completionMarkers'] == 1, name + '/completion')
    spectrum = inventories['spectrum']
    require(not spectrum['inputDigestMismatch'], 'spectrum input digest')
    modes = [record for packet in spectrum['packets'].values() for record in packet['modes']]
    require(spectrum['packets'].keys() == spectrum['coverage'].keys(), 'spectrum packet/coverage census')
    for name, record in spectrum['coverage'].items():
        require(record['finitePolynomialRootCoverage'], 'finite root coverage: ' + name)
        require(len(spectrum['packets'][name]['modes']) == 2 * record['distinctRadicalRoots'],
                'isolated root/lift census: ' + name)
    coverage = Counter(tuple(v[k] for k in ('polynomialDegree', 'distinctRadicalRoots',
        'countWithMultiplicity', 'degreeCountResidual', 'allDisksDisjoint', 'allDisksIsolated'))
        for v in spectrum['coverage'].values())
    inverse = inventories['inverse']
    require(not inverse['coverage_gaps'], 'inverse carrier coverage')
    for tag, value in inverse['literal_residuals'].items():
        require(not value['nonzero'], 'inverse literal residual: ' + tag)
    jets = inventories['mode_jets']
    for tag, value in jets['coverage_gaps'].items():
        require(not value, 'mode jets/' + tag)
    codec = inventories['codec']
    require(codec['unencodedPayloadSha256'] == codec['decodedPayloadSha256'], 'lossless codec payload')
    require(len(codec['sourceIndexChecks']) == 1 and
            all(v['differentAssignments'] == 0 for v in codec['sourceIndexChecks']), 'lossless codec source index')
    require(len(codec['malformedReferenceRejections']) == 2 and
            all(codec['malformedReferenceRejections']), 'malformed codec rejection')
    require(codec['encodedBytes'] == (base / 'full.out').stat().st_size, 'codec native size')
    exceptional = inventories['exceptional']
    report = {
        'nativeManifest': str(manifest_path), 'nativeManifestSha256': digest(manifest_path),
        'instrumentSha256': digest(Path(__file__)), 'sourcePinsStable': pins == current_pins,
        'nativeWallSeconds': manifest['wall_seconds'], 'nativePeakRssKiB': manifest['peak_rss_kib'],
        'nativeArtifact': manifest['artifacts']['full.out'],
        'rechecks': {'path': str(checks_path.relative_to(ROOT)), 'sha256': digest(checks_path),
                     'stages': len(checks['stages'])},
        'spectrum': {'packets': len(spectrum['packets']), 'candidates': len(modes),
                     'nullities': dict(Counter(v['nullity'] for v in modes)),
                     'coverage': {str(k): v for k, v in coverage.items()},
                     'maximumResidualByDimension': spectrum['maximumResidualByDimension']},
        'jointSheet': {'packets': len(inventories['joint']['packets']),
                       'pathStatuses': inventories['joint']['pathStatuses']},
        'exceptional': {'packets': exceptional['packetCount'],
                        'inputDimensionEquationCount': len(exceptional['dimensionConstraints'].get('INPUT_DIMENSION_CONSTRAINTS', [])),
                        'stageLocalStatuses': dict(Counter(str(v['statuses']) for v in exceptional['packets'].values())),
                        'sourceIndexChecks': exceptional['sourceIndexChecks']},
        'resolvent': {key: inventories['resolvent'][key] for key in
                      ('totals', 'poleStatuses', 'contourCount', 'residualMaximaByDimension')},
        'thresholdModes': {key: inventories['threshold_modes'][key] for key in
                           ('counts', 'residualMaximaByDimension', 'sourceIndexChecks')},
        'inverseFourier': {key: inverse[key] for key in
                           ('carrier_counts_by_case', 'literal_residuals', 'residual_projections', 'dimension_records')},
        'modeJets': {key: jets[key] for key in
                     ('rectangular_jet_states', 'nullspace_dimensions', 'fingerprint_residual_projections')},
        'codec': codec, 'outstandingConstructions': jets['outstanding_constructions'],
        'failures': failures,
        'scope': ('Computed point, slice, bank, path and contour domains only. Preserve exceptional-domain records '
                  'and nonzero numerical residuals in the inventories. No global parameter/sheet coverage, '
                  'physical flux normalization, complete scattering or profile-frequency bound-pole claim.')}
    report_path = MEASUREMENTS / (PREFIX + 'full_checks.json')
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ('nativeWallSeconds', 'nativePeakRssKiB', 'nativeArtifact', 'failures')}, indent=2))
    if failures:
        raise SystemExit(1)
    if args.publish:
        target = ROOT / 'scripts/out/S11c_d_mixing_scattering_sympy_audit.out'
        old_payload = target.resolve() if target.is_symlink() else None
        old_hash = digest(target) if target.exists() else None
        with tempfile.NamedTemporaryFile(dir=target.parent, prefix='.s11cd-trace-repair-', delete=False) as stream:
            temporary = Path(stream.name)
            with (base / 'full.out').open('rb') as source:
                shutil.copyfileobj(source, stream)
            stream.flush()
            os.fsync(stream.fileno())
        if digest(temporary) != manifest['artifacts']['full.out']['sha256']:
            raise ValueError('publication content hash')
        os.replace(temporary, target)
        publication = {'path': str(target.relative_to(ROOT)), 'bytes': target.stat().st_size,
                       'sha256': digest(target), 'previousSha256': old_hash,
                       'previousAnnexPayload': str(old_payload) if old_payload else None,
                       'previousAnnexPayloadUnchanged': digest(old_payload) == old_hash if old_payload else None,
                       'producerManifestSha256': digest(manifest_path),
                       'checksPath': str(report_path.relative_to(ROOT)), 'checksSha256': digest(report_path)}
        (MEASUREMENTS / (PREFIX + 'publication.json')).write_text(json.dumps(publication, indent=2) + '\n')


if __name__ == '__main__':
    run()
