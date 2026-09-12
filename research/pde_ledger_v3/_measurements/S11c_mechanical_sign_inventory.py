#!/usr/bin/env python3
"""Verify the sign-audit transcript's metadata, literal residuals, and pins."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.core.symbol import Str

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from ledger_fold import _restore


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def leaves(value, path=()):
    if isinstance(value, (sp.Tuple, sp.MatrixBase)):
        association = all(isinstance(v, sp.Tuple) and len(v) == 2 and isinstance(v[0], Str) for v in value)
        for i, item in enumerate(value):
            yield from leaves(item[1] if association else item, path+(str(item[0]) if association else i,))
    elif not isinstance(value, Str):
        yield path, value


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    base = Path(manifest['run_directory'])
    native = json.loads((base/'inventory.json').read_text())
    records = [json.loads(line) for line in (base/'full.out').open()]
    native_records = {row['writeKey']:row for row in native['records']}
    issues, residuals, anchors, residual_pit = [], {}, {}, {}
    if manifest['exit_code'] != 0:
        issues.append('producer exit')
    if manifest['source_hashes_before'] != manifest['source_hashes_after']:
        issues.append('producer source stability')
    for name, expected in manifest['source_hashes_after'].items():
        if digest(ROOT/name) != expected or digest(base/'source'/name) != expected:
            issues.append('source pin '+name)
    for name, expected in manifest['artifacts'].items():
        if digest(base/name) != expected['sha256']:
            issues.append('artifact pin '+name)
    if len({row['writeKey'] for row in records}) != len(records):
        issues.append('duplicate write-key')
    if {row['writeKey'] for row in records} != set(native_records):
        issues.append('native/transcript record set')
    for row in records:
        key = row['writeKey']
        value = _restore(row['object'])
        metadata = {tuple(item['path']):item for item in row['metadata']}
        if len(metadata) != len(row['metadata']):
            issues.append('duplicate metadata path '+key)
        if row['representation'] == 'carrierSha256AndNumericPit':
            entries = dict(value)[Str('OBJECT_SHA_AND_NUMERIC_PIT')]
            paths = {tuple(str(p) if isinstance(p,Str) else int(p) for p in path) for path,_,_ in entries}
            numbers = [number for _,_,points in entries for number in points]
            if native_records[key]['residual']:
                residual_pit[key] = [{'path':list(map(str,path)),
                                     'nonzeroPoints':sum(number != 0 for number in points),
                                     'pointCount':len(points)} for path,_,points in entries]
        else:
            actual = list(leaves(value))
            paths = {path for path,_ in actual}
            numbers = [value for _,value in actual]
            if native_records[key]['residual']:
                count = sum(value != 0 for _,value in actual)
                if count != native_records[key]['nonzeroLeaves']:
                    issues.append('literal residual census '+key)
        if paths != set(metadata):
            issues.append('metadata coverage '+key)
        if any(value.has(sp.nan,sp.zoo,sp.oo,-sp.oo) for value in numbers):
            issues.append('nonfinite '+key)
        if any(len(item['dimensionLTM']) != 3 or 'nan' in item['dimensionLTM'] or
               not item['epsEtaSigmaOrder'] for item in metadata.values()):
            issues.append('dimension/order '+key)
        if native_records[key]['residual']:
            residuals[key] = {name:native_records[key][name] for name in ('scalarLeaves','zeroLeaves','nonzeroLeaves')}
        if key.endswith('ActionToRowScale') or key.endswith(('S11bStiffnessCoefficient','S11bInertiaCoefficient')):
            anchors[key] = str(value)
    result = {'instrumentSha256':digest(Path(__file__)), 'producerManifestSha256':digest(args.manifest),
              'transcriptSha256':digest(base/'full.out'), 'recordCount':len(records), 'issues':issues,
              'residuals':residuals, 'residualPit':residual_pit, 'orientationAnchors':anchors,
              'nonfiniteLeaves':sum(row['nonfiniteLeaves'] for row in native['records']),
              'unknownDimensions':sum(row['unknownDimensions'] for row in native['records'])}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items()
                     if key not in ('residuals','residualPit','orientationAnchors')},indent=2))
    if issues:
        raise ValueError('sign-audit artifact inventory; see emitted issues')


if __name__ == '__main__':
    run()
