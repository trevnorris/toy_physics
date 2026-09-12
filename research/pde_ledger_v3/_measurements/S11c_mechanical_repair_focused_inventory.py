#!/usr/bin/env python3
"""Inventory completed native-case diagnostic records without expanding them."""
import ast
import hashlib
import json
import os
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/tmp/s11c-mechanical-repair-20260912')


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def samples(record):
    # The producer's scalar carrier fingerprint contains (path, SHA, PIT).
    node = ast.parse(record['value'], mode='eval').body
    row = node.args[0].args[1].args[0]
    return [ast.unparse(value) for value in row.args[2].args]


def run():
    checks = {}
    for stage in ('b_case', 'c2_case'):
        run = BASE / stage
        manifest = json.loads((run / 'manifest.json').read_text())
        records = [json.loads(line) for line in (run / 'full.out').open()]
        keys = [record['writeKey'] for record in records]
        check = {'records': len(keys), 'duplicateKeys': len(keys)-len(set(keys)),
                 'sourcePinsStable': manifest['sourcePins'] == manifest['sourcePinsAfter'],
                 'stderrBytes': (run / 'stderr.txt').stat().st_size,
                 'nonfiniteOrUnknownMetadata': sum(any(word in record['dimensionLTM'] for word in ('nan', 'zoo', 'oo'))
                                                    for record in records)}
        if stage == 'c2_case':
            check['canonicalPowerResidualZero'] = manifest['canonicalResidualZero']
            check['canonicalPowerResidualPit'] = samples(next(record for record in records
                                                if record['writeKey'].endswith('CanonicalPowerResidual')))
            check['kineticNormalizationResidual'] = next(record['value'] for record in records
                                                  if record['writeKey'].endswith('KineticNormalizationResidual'))
            check['controlPitSamples'] = {record['writeKey']: samples(record) for record in records
                                         if record['writeKey'].endswith(('TractionSignResidual', 'FaceRoutingResidual'))}
        else:
            check['literalResiduals'] = {record['writeKey']: record['value'] for record in records
                                        if record['writeKey'].endswith('Residual')}
        for name, expected in manifest['sourcePins'].items():
            path = Path(name)
            if not path.is_absolute():
                path = ROOT / path
            frozen = run / 'source' / path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
            source = frozen if frozen.exists() else path
            if digest(source) != expected:
                raise ValueError(('source pin', stage, name))
            if path.is_relative_to(ROOT) and not frozen.exists():
                frozen.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, frozen)
        if not check['sourcePinsStable'] or check['duplicateKeys'] or check['nonfiniteOrUnknownMetadata'] or check['stderrBytes']:
            raise ValueError(('diagnostic publication metadata', stage, check))
        destination = ROOT / 'scripts/out' / ('S11c_mechanical_repair_' + stage + '.out')
        staging = destination.with_suffix('.out.mechanical-staged')
        shutil.copyfile(run / 'full.out', staging)
        os.replace(staging, destination)
        check.update(publishedPath=str(destination.relative_to(ROOT)), outputBytes=destination.stat().st_size,
                     outputSha256=digest(destination), producerManifestSha256=digest(run / 'manifest.json'))
        checks[stage] = check
        (run / 'checked_inventory.json').write_text(json.dumps(check, indent=2) + '\n')
    (ROOT / '_measurements/S11c_mechanical_repair_focused_checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    run()
