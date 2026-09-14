#!/usr/bin/env python3
"""Reconstruct trace-repair transcript payloads, units, grades and source pins."""
import argparse
import json
from pathlib import Path
import pickle
import sys
import tempfile
import sympy as sp
from S11c_c2_trace_repair_check import ROOT,c,d,digest,TraceMetadata,canonical_power
from ledger_fold import _restore


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-directory',type=Path,required=True)
    parser.add_argument('--publish-name');args=parser.parse_args();base=args.run_directory
    checks=json.loads((base/'checks.json').read_text())
    if checks['sourcePins']!=checks['sourcePinsAfter']:raise ValueError('unstable trace producer')
    for name,pin in checks['sourcePins'].items():
        if digest(ROOT/name)!=pin or digest(base/'source'/name)!=pin:raise ValueError(('trace source pin',name))
    if digest(base/'objects.pickle')!=checks['objectsSha256']:raise ValueError('trace objects pin')
    payload=pickle.loads((base/'objects.pickle').read_bytes())
    c.NEW_DIMENSIONS.update(payload['dimensions']);c.DIMENSION_SCHEMA.update(payload['schema']);c.dimension.cache_clear()
    metadata=TraceMetadata.__new__(TraceMetadata);metadata.generators=payload['metadataGenerators']
    lines=(base/'full.out').read_text().splitlines()
    if len(lines)!=len(payload['records']):raise ValueError('trace record census')
    keys=set();unit_paths=0;residual_scalars=0;nonzero=0
    for line,record in zip(lines,payload['records']):
        restored=_restore(line);expected=record['payload'];value=record['value'];units=record['units']
        fingerprint=d.carrier_fingerprint(value) if d.dag_size(value)>1200 else value
        if restored!=c.cas(expected) or expected['value']!=fingerprint:raise ValueError('trace transcript reconstruction')
        key=expected['writeKey']
        if key in keys:raise ValueError('trace duplicate key')
        keys.add(key)
        if units.has(sp.nan,sp.zoo):raise ValueError('trace unknown units')
        unitmap=d.payload_units(value,units);grades=set();support=set()
        for path,leaf in d.leaves(value):
            if isinstance(leaf,c.Str):continue
            unit_paths+=1
            if leaf!=0 and tuple(c.dimension(leaf))!=tuple(unitmap[path]):raise ValueError(('trace units',key,path))
            coefficients=metadata.coefficients(leaf.xreplace(payload['profiles']));grades.update(coefficients)
            homotopy={}
            for (e,a,b),coefficient in coefficients.items():
                homotopy[e,a+b]=homotopy.get((e,a+b),sp.S.Zero)+coefficient*payload['homotopyRatio']**b
            support.update(g for g,v in homotopy.items() if v!=0)
        if expected['multigrade']!=sp.Tuple(*(sp.Tuple(*g) for g in sorted(grades))):raise ValueError('trace multigrade')
        if expected['epsilonLambdaSupport']!=sp.Tuple(*(sp.Tuple(*g) for g in sorted(support))):raise ValueError('trace lambda support')
        if record['name'].endswith(('_RESIDUAL','CanonicalResidual')) or record['name'].startswith('SourceTraceResidual'):
            leaves=[v for _,v in d.leaves(value) if not isinstance(v,c.Str)]
            residual_scalars+=len(leaves);nonzero+=sum(v!=0 for v in leaves)
    if (residual_scalars,nonzero)!=(checks['residualScalars'],checks['nonzeroResidualScalars']):raise ValueError('trace residual inventory')
    for case,faces in payload['results'].items():
        if 'power' in faces:
            if canonical_power(faces['power']['raw'])!=faces['power']['canonical']:raise ValueError('power canonical replay')
            faces={k:v for k,v in faces.items() if k!='power'}
        if set(faces)!=set(c.FACES):raise ValueError('trace face coverage')
        for face,result in faces.items():
            for direct,factorized,residual in result['sourceJoins']:
                eta,sigma=metadata.generators[1:]
                rebuilt=sum(v*eta**a*sigma**b for (a,b),v in c.shape_coefficients((direct-factorized).xreplace(payload['profiles']),eta,sigma).items() if a<=1 and b<=1)
                if sp.cancel(rebuilt)!=residual:raise ValueError('trace source coefficient reconstruction')
    checks.update(checkedObjects=len(lines),checkedMetadataPaths=unit_paths,checkedKeys=len(keys),
                  validatorSha256=digest(Path(__file__)),cases=['/'.join(case) for case in payload['results']],
                  transcript={'bytes':(base/'full.out').stat().st_size,'sha256':digest(base/'full.out')},runDirectory=str(base))
    if args.publish_name:
        if not args.publish_name.startswith('S11c_c2_trace_repair_') or '/' in args.publish_name:raise ValueError('trace publication name')
        target=ROOT/'scripts/out'/(args.publish_name+'.out')
        if target.exists() or target.is_symlink():raise FileExistsError(target)
        with tempfile.NamedTemporaryFile(dir=target.parent,delete=False) as f:
            f.write((base/'full.out').read_bytes());stage=Path(f.name)
        stage.replace(target);checks['publication']=str(target)
        (ROOT/'_measurements'/(args.publish_name+'_checkpoint.json')).write_text(json.dumps(checks,indent=2)+'\n')
    (base/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps({k:checks[k] for k in ('cases','residualScalars','nonzeroResidualScalars','checkedObjects','checkedMetadataPaths','checkedKeys','wallSeconds','peakRssKiB','transcript')},indent=2))


if __name__=='__main__':main()
