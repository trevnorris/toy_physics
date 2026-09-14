#!/usr/bin/env python3
"""Check recorded reference-current evidence and publish a fresh transcript."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--checks', type=Path, required=True)
    parser.add_argument('--publication-suffix', default='trace_repair')
    parser.add_argument('--checkpoint', type=Path,
                        default=ROOT/'_measurements/S11c_c2_trace_repair_d_current_checkpoint.json')
    args = parser.parse_args()
    if not re.fullmatch('[a-z][a-z0-9_]*',args.publication_suffix):
        raise ValueError('publication suffix must be a lowercase identifier')
    manifest = json.loads(args.manifest.read_text())
    checks = json.loads(args.checks.read_text())
    base = Path(manifest['run_directory'])
    failures = []

    def require(condition, name):
        if not condition:
            failures.append(name)

    pins = manifest['source_hashes_before']
    require(manifest['exit_code'] == 0, 'producer exit')
    require(pins == manifest['source_hashes_after'] == {
        name: digest(ROOT / name) for name in pins}, 'producer sources')
    command = manifest['command']
    full_manifest = Path(command[command.index('--manifest') + 1])
    require(manifest['reference_manifest_sha256_before'] ==
            manifest['reference_manifest_sha256_after'] == digest(full_manifest),
            'full producer manifest')
    for name, artifact in manifest['artifacts'].items():
        require(digest(base / name) == artifact['sha256'] and
                (base / name).stat().st_size == artifact['bytes'], 'artifact: ' + name)
    require((base / 'stderr.txt').stat().st_size == 0, 'producer stderr')
    require(checks['producerManifestSha256'] == digest(args.manifest), 'inventory producer pin')
    require(checks['instrumentSha256'] == digest(
        ROOT / '_measurements/S11c_d_nonlocal_current_inventory.py'), 'inventory source pin')
    for key in ('metadataGaps', 'nonfiniteTags', 'zeroMapMetadataTags', 'unknownDimensionSymbols'):
        require(not checks[key], key)
    require(checks['dimensionConstraints'] == 'Tuple()', 'dimension closure')
    for name, residual in checks['residuals'].items():
        require(residual['nonzeroCount'] == 0, 'residual: ' + name)
    require(checks['mechanicalOrientationAnchor']['residual'] == 'Integer(0)', 'row orientation')
    report = {
        'manifest': str(args.manifest), 'manifestSha256': digest(args.manifest),
        'instrumentSha256': digest(Path(__file__)),
        'checksPath': str(args.checks), 'checksSha256': digest(args.checks),
        'wallSeconds': manifest['wall_seconds'], 'peakRssKiB': manifest['peak_rss_kib'],
        'tagCount': checks['tagCount'],
        'residualScalarCount': sum(v['scalarCount'] for v in checks['residuals'].values()),
        'nonzeroResidualCount': sum(v['nonzeroCount'] for v in checks['residuals'].values()),
        'mechanicalOrientationAnchor': checks['mechanicalOrientationAnchor'],
        'failures': failures,
        'scope': 'LAB_HELD/RHO4_CONSTANT reference current and acoustic source; no endpoint flux normalization.'}
    print(json.dumps(report, indent=2), flush=True)
    if failures:
        raise SystemExit(1)
    target = ROOT / 'scripts/out' / ('S11c_d_nonlocal_current_reference_'+args.publication_suffix+'.out')
    for path in (target,args.checkpoint):
        if path.exists() or path.is_symlink():raise FileExistsError(path)
    with tempfile.NamedTemporaryFile(dir=target.parent, prefix='.s11cd-current-trace-', delete=False) as stream:
        temporary = Path(stream.name)
        with (base / 'full.out').open('rb') as source:
            shutil.copyfileobj(source, stream)
        stream.flush()
        os.fsync(stream.fileno())
    if digest(temporary) != manifest['artifacts']['full.out']['sha256']:
        raise ValueError('publication content hash')
    os.replace(temporary, target)
    report['publication'] = {'path': str(target.relative_to(ROOT)),
                             'bytes': target.stat().st_size, 'sha256': digest(target)}
    args.checkpoint.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    run()
